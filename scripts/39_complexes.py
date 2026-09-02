"""Exploration - structural template validation of the interaction layer.

For each high-confidence STRING pair whose two partners both have AlphaFold
models, we ask whether those two models each match a DIFFERENT chain of the
same experimental PDB complex (both TM >= 0.5). Such a shared-complex
template is structural evidence that the two proteins can form a complex,
and it names a concrete interface (the PDB entry) - the GPU-free analogue of
predicting the dimer.

The result only means something against a control: random protein pairs
drawn from the same proteins. If STRING pairs are templated far more often
than random pairs, structure independently corroborates the interaction
network. Dark-protein interactions with template support are novel,
structurally backed interaction hypotheses.

Outputs: results/explore/complexes.md, templated_interactions.tsv
"""
import random
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pis.common import data_dir, load_config, norm_hit_id

TM_MIN = 0.5
# Foldseek PDB target -> (pdb_id, chain). Handles '1abc_A', 'pdb1abc.ent.gz_A',
# '1abc-assembly1_A', etc.: take trailing _<chain>, pull a 4-char PDB code.
PDB_RE = re.compile(r"([0-9][a-z0-9]{3})", re.IGNORECASE)


def parse_target(t):
    if "_" not in t:
        return None
    base, chain = t.rsplit("_", 1)
    m = PDB_RE.search(base.replace("pdb", "", 1))
    if not m or not chain:
        return None
    return m.group(1).lower(), chain


def main():
    cfg = load_config()
    d = data_dir(cfg)
    out_dir = Path("results/explore")
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(cfg["seed"] + 39)

    # protein -> pdb_id -> {chain: best_tm}
    hits = defaultdict(lambda: defaultdict(dict))
    n_lines = 0
    with open(d / "explore" / "involved_vs_pdb.tsv", encoding="utf-8") as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 3:
                continue
            q = norm_hit_id(c[0])
            pt = parse_target(c[1])
            if not pt:
                continue
            try:
                tm = float(c[2])
            except ValueError:
                continue
            if tm < TM_MIN:
                continue
            pdb_id, chain = pt
            prev = hits[q][pdb_id].get(chain, 0.0)
            if tm > prev:
                hits[q][pdb_id][chain] = tm
            n_lines += 1
    print("proteins with >=1 PDB match (TM>={}): {} (from {} qualifying hit rows)".format(
        TM_MIN, len(hits), n_lines))

    def template(a, b):
        """Return (pdb_id, min_tm) of the best shared complex where a and b
        map to different chains, or None."""
        ha, hb = hits.get(a), hits.get(b)
        if not ha or not hb:
            return None
        best = None
        for pdb_id in set(ha) & set(hb):
            for ca, ta in ha[pdb_id].items():
                for cb, tb in hb[pdb_id].items():
                    if ca != cb:
                        m = min(ta, tb)
                        if best is None or m > best[1]:
                            best = (pdb_id, m)
        return best

    # STRING pairs with both partners structured
    structured = set(hits)
    string_pairs = []
    dark = set((d / "dark" / "accessions.txt").read_text(encoding="utf-8").split())
    with open(d / "string_edges.tsv", encoding="utf-8") as f:
        next(f)
        for line in f:
            a, b, s = line.rstrip("\n").split("\t")
            string_pairs.append((a, b))
    # restrict to pairs where both have any PDB match (fair denominator)
    eligible = [(a, b) for a, b in string_pairs if a in structured and b in structured]
    print("STRING pairs eligible (both partners have a PDB match): {}".format(len(eligible)))

    templated = []
    for a, b in eligible:
        t = template(a, b)
        if t:
            templated.append((a, b, t[0], t[1], (a in dark or b in dark)))
    frac_string = len(templated) / len(eligible) if eligible else 0.0

    # control: random pairs from the same eligible protein pool, same count
    pool = sorted({p for pair in eligible for p in pair})
    n_ctrl = min(len(eligible), 20000)
    ctrl_string = {frozenset(p) for p in eligible}
    ctrl_hit = 0
    tries = 0
    seen = set()
    while len(seen) < n_ctrl and tries < n_ctrl * 20:
        tries += 1
        a, b = rng.sample(pool, 2)
        key = frozenset((a, b))
        if key in ctrl_string or key in seen:
            continue
        seen.add(key)
        if template(a, b):
            ctrl_hit += 1
    frac_ctrl = ctrl_hit / len(seen) if seen else 0.0
    enrichment = (frac_string / frac_ctrl) if frac_ctrl > 0 else float("inf")

    dark_templated = [t for t in templated if t[4]]
    templated.sort(key=lambda x: -x[3])
    with open(out_dir / "templated_interactions.tsv", "w", encoding="utf-8", newline="\n") as f:
        f.write("protein_a\tprotein_b\ttemplate_pdb\tmin_TM\tinvolves_dark\n")
        for a, b, pdb, tm, isd in templated:
            f.write("{}\t{}\t{}\t{:.3f}\t{}\n".format(a, b, pdb, tm, int(isd)))

    with open(out_dir / "complexes.md", "w", encoding="utf-8", newline="\n") as f:
        f.write("# Structural template validation of interactions\n\n")
        f.write("For each STRING pair whose partners both have an AlphaFold "
                "model, we test whether the two models match different chains "
                "of the same experimental PDB complex (TM >= {}). This is a "
                "GPU-free proxy for predicting the dimer: a shared-complex "
                "template both corroborates the interaction and names an "
                "interface.\n\n".format(TM_MIN))
        f.write("| set | pairs | with structural template | fraction |\n")
        f.write("|---|---|---|---|\n")
        f.write("| STRING (high-confidence) | {} | {} | {:.1%} |\n".format(
            len(eligible), len(templated), frac_string))
        f.write("| random control (same proteins) | {} | {} | {:.1%} |\n\n".format(
            len(seen), ctrl_hit, frac_ctrl))
        f.write("**Enrichment of real interactions over random: {:.1f}x.** "
                "If large, structure independently supports the STRING "
                "network and supplies candidate interfaces.\n\n".format(enrichment))
        f.write("Interactions involving a dark (unannotated) protein that "
                "have structural template support: **{}** - novel, "
                "structurally backed interaction hypotheses. Top 15 by "
                "confidence:\n\n".format(len(dark_templated)))
        f.write("| protein A | protein B | template PDB | min TM |\n|---|---|---|---|\n")
        for a, b, pdb, tm, _ in sorted(dark_templated, key=lambda x: -x[3])[:15]:
            f.write("| {} | {} | {} | {:.2f} |\n".format(a, b, pdb, tm))
        f.write("\nCaveat: template-based, not de novo AlphaFold-Multimer; a "
                "shared fold-pair template is evidence of feasibility and a "
                "likely interface, not proof of a specific cellular complex.\n")

    print("STRING templated: {:.1%} | control: {:.1%} | enrichment {:.1f}x".format(
        frac_string, frac_ctrl, enrichment))
    print("dark templated interactions: {}".format(len(dark_templated)))
    print("Wrote results/explore/complexes.md")


if __name__ == "__main__":
    main()
