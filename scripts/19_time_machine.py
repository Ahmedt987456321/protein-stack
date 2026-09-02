"""Step 19 (v1.0) - the experiment loop, closed against time.

Rebuild the annotation world as of an archived GOA release (~3 years back),
generate function hypotheses for proteins that were dark THEN, and score them
against the experimental annotations that have accumulated SINCE. Outcomes -
supported / contradicted / unconfirmed - are written back for the graph.

Leakage control: the gated arm (SS) transfers from sequence + structure
neighbours only. InterPro/interpro2go are current-day databases, so the
domain-augmented arm (FULL) is reported for information with that caveat.

Outcome semantics (open world): `supported` = a new experimental annotation
implies the predicted term; `contradicted` = an experimental NOT-qualifier
annotation contradicts the exact term; otherwise `unconfirmed` - absence of
annotation is not absence of function.

Outputs:
  data/timesplit/...                    (past GAFs, cohort inputs, searches)
  results/hypothesis_outcomes.tsv     (accession, term, branch, confidence,
                                       arm, outcome, source)
  results/report_v10.md
"""
import datetime
import gzip
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
    domain_stats,
    domain_stream,
    fuse,
    load_domains,
    load_interpro2go,
    parse_raw_seq,
    parse_raw_struct,
    specific_terms,
    split_by_branch,
    transfer_stream,
)

OLD_BASE = "https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/{DIR}/"
N_PERM = 1000
TIER_HIGH, TIER_MED = 0.7, 0.45


def tier_of(score):
    return "HIGH" if score >= TIER_HIGH else ("MEDIUM" if score >= TIER_MED else "LOW")


def pick_old_release(session, species, target):
    """Return (filename, date) of the archived GAF nearest the target date."""
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


def parse_not_experimental(gaf_gz, evidence_codes):
    """Experimental NOT-qualifier rows: acc -> {term} (genuine negatives)."""
    from pis.go import _open_gaf
    codes = set(evidence_codes)
    neg = defaultdict(set)
    with _open_gaf(gaf_gz) as f:
        for line in f:
            if line.startswith("!"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 15 or cols[0] != "UniProtKB":
                continue
            if "NOT" not in cols[3].split("|"):
                continue
            if cols[6] in codes and cols[8] in ("F", "P", "C"):
                neg[cols[1].split("-")[0]].add(cols[4])
    return neg


def main():
    cfg = load_config()
    d = data_dir(cfg)
    ts = d / "timesplit"
    ts.mkdir(exist_ok=True)
    r = results_dir(cfg)
    session = http_session()
    target = datetime.date.fromisoformat(cfg["timesplit"]["target_date"])
    dag = GoDag(d / "go-basic.obo")

    # ---- past vs current annotation worlds --------------------------------
    past_exp = defaultdict(set)
    releases = []
    for sp in cfg["species"]:
        try:
            url, fname, date = pick_old_release(session, sp["name"], target)
        except Exception as e:
            print("{}: no archived release ({}); excluded from time-machine".format(
                sp["name"], e))
            continue
        print("{}: past release {} ({})".format(sp["name"], fname, date))
        releases.append("{} {} ({})".format(sp["name"], fname, date))
        gaf = download(session, url, ts / fname)
        for acc, term, _aspect in parse_gaf(gaf, cfg["evidence_codes"]):
            past_exp[acc].add(term)

    current_exp = defaultdict(set)
    current_not = defaultdict(set)
    for sp in cfg["species"]:
        gaf = d / "gaf" / "goa_{}.gaf.gz".format(sp["name"])
        for acc, term, _aspect in parse_gaf(gaf, cfg["evidence_codes"]):
            current_exp[acc].add(term)
        for acc, terms in parse_not_experimental(gaf, cfg["evidence_codes"]).items():
            current_not[acc] |= terms
    print("Past experimental proteins: {}; current: {}".format(
        len(past_exp), len(current_exp)))

    # ---- cohort: dark then, annotated now, inside our universe ------------
    kb = set((d / "final_accessions.txt").read_text(encoding="utf-8").split())
    dark = set((d / "dark" / "accessions.txt").read_text(encoding="utf-8").split())
    universe = kb | dark
    cohort = sorted(a for a in universe
                    if current_exp.get(a) and not past_exp.get(a))
    kb_past = sorted(a for a in kb - set(cohort) if past_exp.get(a))
    print("Cohort (dark in {}, experimentally annotated now): {}".format(
        target.year, len(cohort)))
    print("Past knowledge base: {} proteins".format(len(kb_past)))
    if len(cohort) < 30:
        print("Cohort too small for a meaningful gate; aborting.")
        sys.exit(1)

    # ---- searches: cohort vs past KB --------------------------------------
    all_seqs = read_fasta(d / "final.fasta")
    all_seqs.update(read_fasta(d / "dark" / "dark.fasta"))
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

    # ---- streams from PAST knowledge only ---------------------------------
    import json
    params = json.load(open(r / "fusion_params.json", encoding="utf-8"))
    kb_past_set = set(kb_past)
    past_terms = {a: dag.propagate(past_exp[a]) for a in kb_past}
    s_seq, seq_transfer = parse_raw_seq(ts / "seq_hits_raw.tsv", kb_past_set,
                                        cfg["search"]["evalue"])
    struct_hits = parse_raw_struct(ts / "struct_hits_raw.tsv", kb_past_set)
    seq_stream = split_by_branch(transfer_stream(seq_transfer, past_terms), dag)
    str_stream = split_by_branch(transfer_stream(struct_hits, past_terms), dag)

    # FULL arm only: current-day domains (leakage-caveated)
    domains = load_domains(d / "domains.tsv")
    for acc, doms in load_domains(d / "dark" / "domains.tsv").items():
        domains.setdefault(acc, set()).update(doms)
    ip2go = load_interpro2go(d / "interpro2go.txt", dag)
    kb_probs = domain_stats(kb_past, past_terms, domains)
    dom_stream = split_by_branch(domain_stream(cohort, domains, kb_probs, ip2go), dag)

    empty = {b: {} for b in BRANCHES}
    arms = {
        "SS": {b: fuse(cohort, seq_stream[b], str_stream[b], empty[b], s_seq,
                       params["tau"], params["alpha"], params["beta"],
                       params["gamma"], params["tau_d"], params["delta"])
               for b in BRANCHES},
        "FULL": {b: fuse(cohort, seq_stream[b], str_stream[b], dom_stream[b], s_seq,
                         params["tau"], params["alpha"], params["beta"],
                         params["gamma"], params["tau_d"], params["delta"])
                 for b in BRANCHES},
    }

    truth = {b: {} for b in BRANCHES}
    for acc in cohort:
        prop = dag.propagate(current_exp[acc])
        for b in BRANCHES:
            sub = {t for t in prop if dag.branch(t) == b}
            if sub:
                truth[b][acc] = sub

    # ---- top-1 precision + permutation null -------------------------------
    min_conf = cfg["timesplit"]["min_confidence"]
    rng = random.Random(cfg["seed"] + 7)
    prec = {}
    pvals = {}
    for arm in arms:
        for b in BRANCHES:
            tops = {}
            for acc in cohort:
                picks = specific_terms(arms[arm][b].get(acc, {}), dag, min_conf, 1)
                if picks and acc in truth[b]:
                    tops[acc] = picks[0][0]
            if not tops:
                prec[(arm, b)] = (0.0, 0)
                pvals[(arm, b)] = 1.0
                continue
            accs = sorted(tops)
            hits = sum(1 for a in accs if tops[a] in truth[b][a])
            obs = hits / len(accs)
            prec[(arm, b)] = (obs, len(accs))
            terms_list = [tops[a] for a in accs]
            ge = 0
            for _ in range(N_PERM):
                perm = terms_list[:]
                rng.shuffle(perm)
                ph = sum(1 for a, t in zip(accs, perm) if t in truth[b][a])
                if ph / len(accs) >= obs:
                    ge += 1
            pvals[(arm, b)] = (ge + 1) / (N_PERM + 1)

    # ---- outcome writeback (arm SS) ---------------------------------------
    counts = defaultdict(int)
    tier_stats = defaultdict(lambda: [0, 0])  # tier -> [supported, total]
    source = "GOA-timesplit-{}".format(target.isoformat())
    with open(r / "hypothesis_outcomes.tsv", "w", encoding="utf-8", newline="\n") as f:
        f.write("accession\tgo_term\tbranch\tconfidence\tarm\toutcome\tsource\n")
        for b in BRANCHES:
            for acc in cohort:
                picks = specific_terms(arms["SS"][b].get(acc, {}), dag, min_conf,
                                       cfg["timesplit"]["top_per_branch"])
                tr = truth[b].get(acc, set())
                for term, score in picks:
                    if term in tr:
                        outcome = "supported"
                    elif term in current_not.get(acc, set()):
                        outcome = "contradicted"
                    else:
                        outcome = "unconfirmed"
                    counts[outcome] += 1
                    tstat = tier_stats[tier_of(score)]
                    tstat[1] += 1
                    if outcome == "supported":
                        tstat[0] += 1
                    f.write("{}\t{}\t{}\t{:.4f}\t{}\t{}\t{}\n".format(
                        acc, term, b, score, "SS", outcome, source))

    g_a = pvals[("SS", "molecular_function")] < 0.01 and prec[("SS", "molecular_function")][0] > 0
    verdict = "PASSED" if g_a else "NOT PASSED"

    # ---- report ------------------------------------------------------------
    with open(r / "report_v10.md", "w", encoding="utf-8", newline="\n") as f:
        f.write("# v1.0 - the experiment loop (time-machine study) - report\n\n")
        f.write("Past world: {}. Cohort: **{} proteins** that had zero "
                "experimental GO then and have it now; past KB: {} proteins "
                "(annotations truncated to the past release).\n\n".format(
                    "; ".join(releases), len(cohort), len(kb_past)))
        f.write("Arm SS = sequence+structure transfer only (leakage-controlled, "
                "gated). Arm FULL adds current-day InterPro domains - reported "
                "with that caveat, not gated.\n\n")
        f.write("## Top-1 precision vs permutation null (n = scored proteins)\n\n")
        f.write("| arm | branch | n | precision@1 | perm. p |\n|---|---|---:|---:|---:|\n")
        for arm in ("SS", "FULL"):
            for b in BRANCHES:
                p, n = prec[(arm, b)]
                f.write("| {} | {} | {} | {:.3f} | {:.4f} |\n".format(
                    arm, b, n, p, pvals[(arm, b)]))
        f.write("\n## Outcomes written back (arm SS, top {} per branch)\n\n".format(
            cfg["timesplit"]["top_per_branch"]))
        f.write("supported {} | contradicted {} | unconfirmed {} - all stored as "
                "typed edges in kg.db (rebuild step 13). `unconfirmed` is weak "
                "negative evidence only: GO is open-world.\n\n".format(
                    counts["supported"], counts["contradicted"], counts["unconfirmed"]))
        f.write("## Confidence calibration (arm SS, supported rate by tier)\n\n")
        f.write("| tier | supported | total | rate |\n|---|---:|---:|---:|\n")
        for t in ("HIGH", "MEDIUM", "LOW"):
            s, n = tier_stats[t]
            f.write("| {} | {} | {} | {} |\n".format(
                t, s, n, "{:.1%}".format(s / n) if n else "-"))
        rates = [tier_stats[t][0] / tier_stats[t][1] if tier_stats[t][1] else 0.0
                 for t in ("HIGH", "MEDIUM", "LOW")]
        f.write("\nCalibration monotone (HIGH >= MEDIUM >= LOW): **{}**\n\n".format(
            "yes" if rates[0] >= rates[1] >= rates[2] else "NO - flagged"))
        f.write("**Gate 1.0** (MF precision@1 of arm SS beats permutation null "
                "at p < 0.01, and outcomes stored): **{}**\n".format(verdict))

    for arm in ("SS", "FULL"):
        for b in BRANCHES:
            p, n = prec[(arm, b)]
            print("{:4s} {:22s} n={:4d}  p@1={:.3f}  perm-p={:.4f}".format(
                arm, b, n, p, pvals[(arm, b)]))
    print("Outcomes: " + ", ".join("{}={}".format(k, counts[k])
                                   for k in ("supported", "contradicted", "unconfirmed")))
    print("PHASE 1.0 GATE: " + verdict)
    print("Wrote results/hypothesis_outcomes.tsv and results/report_v10.md")


if __name__ == "__main__":
    main()
