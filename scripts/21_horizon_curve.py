"""Step 21 (analysis) - temporal validation at multiple horizons.

Repeats the leakage-controlled time-machine protocol (sequence + structure
streams only, fitted fusion gate) for several past GOA releases, giving a
top-1 precision versus prediction-horizon curve. The 2023 horizon reproduces
the Section 3.5 protocol and serves as a consistency check.

Outputs:
  data/timesplit_h<year>/...      (cached per-horizon corpora and searches)
  results/horizon_curve.csv, results/horizon_curve.md
"""
import csv
import datetime
import json
import random
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pis.common import (
    data_dir,
    download,
    http_session,
    link_or_copy,
    load_config,
    read_fasta,
    results_dir,
    run_tool,
    write_fasta,
)
from pis.go import GoDag, parse_gaf
from pis.streams import (
    BRANCHES,
    fuse,
    parse_raw_seq,
    parse_raw_struct,
    specific_terms,
    split_by_branch,
    transfer_stream,
)

OLD_BASE = "https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/{DIR}/"
HORIZONS = ["2019-07-01", "2021-07-01", "2023-07-12"]
N_PERM = 1000


def pick_old_release(session, species, target):
    url = OLD_BASE.format(DIR=species.upper())
    html = session.get(url, timeout=120).text
    pat = re.compile(
        r'href="(goa_{}\.gaf\.\d+\.gz)".*?(\d{{4}}-\d{{2}}-\d{{2}})'.format(species))
    best = None
    for fname, datestr in pat.findall(html):
        d = datetime.date.fromisoformat(datestr)
        gap = abs((d - target).days)
        if best is None or gap < best[0]:
            best = (gap, fname, d)
    if best is None:
        raise RuntimeError("no archived releases found for " + species)
    return url + best[1], best[1], best[2]


def run_horizon(cfg, d, dag, session, params, universe, kb, current_exp,
                all_seqs, date_str):
    target = datetime.date.fromisoformat(date_str)
    ts = d / "timesplit_h{}".format(target.year)
    ts.mkdir(exist_ok=True)

    past_exp = defaultdict(set)
    releases = []
    for sp in cfg["species"]:
        try:
            url, fname, date = pick_old_release(session, sp["name"], target)
        except RuntimeError:
            # species without archived per-species GOA releases (e.g. E. coli)
            # are excluded from the temporal cohort, as in step 19.
            releases.append("{} (no archive, excluded)".format(sp["name"]))
            continue
        releases.append("{} {} ({})".format(sp["name"], fname, date))
        gaf = download(session, url, ts / fname)
        for acc, term, _aspect in parse_gaf(gaf, cfg["evidence_codes"]):
            past_exp[acc].add(term)

    cohort = sorted(a for a in universe if current_exp.get(a) and not past_exp.get(a))
    kb_past = sorted(a for a in kb - set(cohort) if past_exp.get(a))
    print("{}: cohort {} | past KB {} | {}".format(
        date_str, len(cohort), len(kb_past), "; ".join(releases)))

    if not (ts / "cohort.fasta").exists():
        write_fasta(ts / "cohort.fasta", {a: all_seqs[a] for a in cohort})
        write_fasta(ts / "kbpast.fasta", {a: all_seqs[a] for a in kb_past})
    for dirname, members in (("struct_cohort", cohort), ("struct_kbpast", kb_past)):
        out = ts / dirname
        if not out.exists():
            out.mkdir()
            for acc in members:
                for src_dir in (d / "structures", d / "dark" / "structures"):
                    src = src_dir / (acc + ".pdb")
                    if src.exists():
                        link_or_copy(src, out / (acc + ".pdb"))
                        break
    th = str(cfg["tools"]["threads"])
    if not (ts / "seq_hits_raw.tsv").exists():
        run_tool(cfg, "mmseqs", [
            "easy-search", ts / "cohort.fasta", ts / "kbpast.fasta",
            ts / "seq_hits_raw.tsv", ts / "tmp_seq",
            "--format-output", "query,target,pident,evalue",
            "-e", "10", "-s", "7.5",
            "--max-seqs", str(cfg["search"]["max_seqs"]), "--threads", th,
        ])
        shutil.rmtree(ts / "tmp_seq", ignore_errors=True)
    if not (ts / "struct_hits_raw.tsv").exists():
        run_tool(cfg, "foldseek", [
            "easy-search", ts / "struct_cohort", ts / "struct_kbpast",
            ts / "struct_hits_raw.tsv", ts / "tmp_struct",
            "--format-output", "query,target,alntmscore,bits,evalue",
            "-e", str(cfg["search"]["evalue"]),
            "--max-seqs", str(cfg["search"]["max_seqs"]), "--threads", th,
        ])
        shutil.rmtree(ts / "tmp_struct", ignore_errors=True)

    kb_past_set = set(kb_past)
    past_terms = {a: dag.propagate(past_exp[a]) for a in kb_past}
    s_seq, seq_transfer = parse_raw_seq(ts / "seq_hits_raw.tsv", kb_past_set,
                                        cfg["search"]["evalue"])
    struct_hits = parse_raw_struct(ts / "struct_hits_raw.tsv", kb_past_set)
    seq_stream = split_by_branch(transfer_stream(seq_transfer, past_terms), dag)
    str_stream = split_by_branch(transfer_stream(struct_hits, past_terms), dag)
    empty = {b: {} for b in BRANCHES}
    fused = {b: fuse(cohort, seq_stream[b], str_stream[b], empty[b], s_seq,
                     params["tau"], params["alpha"], params["beta"],
                     params["gamma"], params["tau_d"], params["delta"])
             for b in BRANCHES}

    truth = {b: {} for b in BRANCHES}
    for acc in cohort:
        prop = dag.propagate(current_exp[acc])
        for b in BRANCHES:
            sub = {t for t in prop if dag.branch(t) == b}
            if sub:
                truth[b][acc] = sub

    rng = random.Random(cfg["seed"] + 8)
    rows = []
    min_conf = cfg["timesplit"]["min_confidence"]
    for b in BRANCHES:
        tops = {}
        for acc in cohort:
            picks = specific_terms(fused[b].get(acc, {}), dag, min_conf, 1)
            if picks and acc in truth[b]:
                tops[acc] = picks[0][0]
        accs = sorted(tops)
        if not accs:
            rows.append([date_str, b, 0, 0.0, 1.0, len(cohort)])
            continue
        hits = sum(1 for a in accs if tops[a] in truth[b][a])
        obs = hits / len(accs)
        terms_list = [tops[a] for a in accs]
        ge = 0
        for _ in range(N_PERM):
            perm = terms_list[:]
            rng.shuffle(perm)
            if sum(1 for a, t in zip(accs, perm) if t in truth[b][a]) / len(accs) >= obs:
                ge += 1
        rows.append([date_str, b, len(accs), round(obs, 3),
                     round((ge + 1) / (N_PERM + 1), 4), len(cohort)])
    return rows


def main():
    cfg = load_config()
    d = data_dir(cfg)
    r = results_dir(cfg)
    session = http_session()
    dag = GoDag(d / "go-basic.obo")
    params = json.load(open(r / "fusion_params.json", encoding="utf-8"))

    kb = set((d / "final_accessions.txt").read_text(encoding="utf-8").split())
    dark = set((d / "dark" / "accessions.txt").read_text(encoding="utf-8").split())
    universe = kb | dark

    current_exp = defaultdict(set)
    for sp in cfg["species"]:
        gaf = d / "gaf" / "goa_{}.gaf.gz".format(sp["name"])
        for acc, term, _aspect in parse_gaf(gaf, cfg["evidence_codes"]):
            current_exp[acc].add(term)

    all_seqs = read_fasta(d / "final.fasta")
    all_seqs.update(read_fasta(d / "dark" / "dark.fasta"))

    all_rows = []
    for date_str in HORIZONS:
        all_rows += run_horizon(cfg, d, dag, session, params, universe, kb,
                                current_exp, all_seqs, date_str)

    with open(r / "horizon_curve.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["horizon", "branch", "n_scored", "top1_precision",
                    "perm_p", "cohort_size"])
        w.writerows(all_rows)
    with open(r / "horizon_curve.md", "w", encoding="utf-8", newline="\n") as f:
        f.write("# Top-1 precision vs prediction horizon (arm SS)\n\n")
        f.write("| horizon | branch | n scored | top-1 precision | perm p | cohort |\n")
        f.write("|---|---|---|---|---|---|\n")
        for row in all_rows:
            f.write("| " + " | ".join(str(x) for x in row) + " |\n")
    for row in all_rows:
        print(row)
    print("Wrote results/horizon_curve.csv|md")


if __name__ == "__main__":
    main()
