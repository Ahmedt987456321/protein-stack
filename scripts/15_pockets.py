"""Step 15 (v0.4) - pocket detection / druggability triage over the dark set.

Runs fpocket (inside WSL) on every dark-protein AlphaFold model. Structures
are copied to WSL-native disk first - fpocket over the /mnt/c bridge is an
order of magnitude slower.

Outputs:
  data/dark/pockets.tsv        (accession, n_pockets, best_druggability, best_score)
  results/druggability_top.md  top candidates joined with their function hypotheses
"""
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pis.common import data_dir, load_config, results_dir, to_wsl_path

BASH_TEMPLATE = r"""
set -e
SRC={src}
OUT={out}
rm -rf /root/work/pk && mkdir -p /root/work/pk
: > "$OUT"
# Per-structure worker: run fpocket, extract the summary line, delete output
# immediately so peak disk stays at a few structures regardless of set size.
process_one() {{
  local pdb="$1" src="$2" out="$3"
  local acc="${{pdb%.pdb}}"
  local work="/root/work/pk/$acc"
  rm -rf "$work" && mkdir -p "$work"
  cp "$src/$pdb" "$work/"
  ( cd "$work" && fpocket -f "$pdb" > /dev/null 2>&1 || true )
  local info="$work/${{acc}}_out/${{acc}}_info.txt"
  if [ -f "$info" ]; then
    awk -v acc="$acc" '
      /^Pocket/            {{ n++ }}
      /Druggability Score/ {{ if ($4+0 > d) d = $4+0 }}
      /^\tScore/           {{ if ($3+0 > s) s = $3+0 }}
      END {{ printf "%s\t%d\t%.3f\t%.3f\n", acc, n+0, d+0, s+0 }}' "$info" >> "$out"
  fi
  rm -rf "$work"
}}
export -f process_one
ls "$SRC" | grep '\.pdb$' | \
  xargs -P {workers} -I XX bash -c 'process_one "XX" "$0" "$1"' "$SRC" "$OUT"
sort -o "$OUT" "$OUT"
rm -rf /root/work/pk
"""


def main():
    cfg = load_config()
    d = data_dir(cfg)
    dd = d / "dark"
    r = results_dir(cfg)

    script = BASH_TEMPLATE.format(
        src=to_wsl_path(dd / "structures"),
        out=to_wsl_path(dd / "pockets.tsv"),
        workers=cfg["tools"]["threads"],
    )
    sh_path = dd / "run_pockets.sh"
    with open(sh_path, "w", encoding="utf-8", newline="\n") as f:
        f.write(script)

    pk_path = dd / "pockets.tsv"
    if pk_path.exists() and pk_path.stat().st_size > 0:
        print("pockets.tsv present; skipping fpocket run (delete to recompute).")
    else:
        print("Running fpocket over the dark set inside WSL ...")
        subprocess.run(["wsl", "-u", "root", "bash", to_wsl_path(sh_path)], check=True)

    seen, rows = set(), []
    with open(pk_path, encoding="utf-8") as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 4:
                continue  # skip blanks/partial lines from parallel appends
            acc, n, drug, score = cols
            if acc in seen:
                continue
            seen.add(acc)
            try:
                rows.append((acc, int(n), float(drug), float(score)))
            except ValueError:
                continue
    rows.sort(key=lambda x: -x[2])
    print("Pockets computed for {} proteins; {} with druggability >= 0.5".format(
        len(rows), sum(1 for _a, _n, dg, _s in rows if dg >= 0.5)))

    # join top candidates with their best MF hypothesis for context
    top_mf = {}
    hyp_path = r / "dark_hypotheses.tsv"
    if hyp_path.exists():
        with open(hyp_path, encoding="utf-8") as f:
            next(f)
            for line in f:
                acc, branch, term, score, tier, streams, _so = line.rstrip("\n").split("\t")
                if branch == "molecular_function" and (
                        acc not in top_mf or float(score) > top_mf[acc][1]):
                    top_mf[acc] = (term, float(score), tier)

    with open(r / "druggability_top.md", "w", encoding="utf-8", newline="\n") as f:
        f.write("# Druggability triage - dark proteins, fpocket over AlphaFold models\n\n")
        f.write("fpocket druggability score ranges 0-1; >= 0.5 is conventionally "
                "considered druggable. All values are COMPUTATIONAL evidence on "
                "predicted structures - triage input, not conclusions.\n\n")
        f.write("| accession | druggability | pockets | best MF hypothesis | conf | tier |\n")
        f.write("|---|---:|---:|---|---:|---|\n")
        for acc, n, drug, _score in rows[:25]:
            hyp = top_mf.get(acc)
            f.write("| {} | {:.3f} | {} | {} | {} | {} |\n".format(
                acc, drug, n,
                hyp[0] if hyp else "-",
                "{:.2f}".format(hyp[1]) if hyp else "-",
                hyp[2] if hyp else "-"))
    print("Wrote data/dark/pockets.tsv and results/druggability_top.md")


if __name__ == "__main__":
    main()
