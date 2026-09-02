"""Step 11 (v0.3) - build the dark-protein set.

"Dark" here = a protein appearing in the GOA GAFs with NO experimentally
supported GO annotation (only electronic/other non-experimental evidence).
Their electronic (IEA-class) annotations are kept aside - never used for
prediction - as an independent corroboration signal for step 12.

Fetches sequences, AlphaFold models (same pLDDT gate), and InterPro entries
for a sampled subset.

Outputs (all under data/dark/):
  accessions.txt, dark.fasta, structures/<ACC>.pdb, domains.tsv
  iea_annotations.tsv   (accession, go_term, aspect) non-experimental only
"""
import random
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tqdm import tqdm

from pis.common import data_dir, http_session, load_config, write_fasta
from pis.fetch import UNIPROT_BATCH, fetch_alphafold, fetch_interpro, fetch_uniprot_batch
from pis.go import parse_gaf_full


def main():
    cfg = load_config()
    d = data_dir(cfg)
    dd = d / "dark"
    dd.mkdir(exist_ok=True)
    (dd / "structures").mkdir(exist_ok=True)
    session = http_session()
    exp_codes = set(cfg["evidence_codes"])

    # ---- find IEA-only proteins -------------------------------------------
    has_exp = set()
    nonexp = defaultdict(set)  # acc -> {(term, aspect)} non-experimental
    for sp in cfg["species"]:
        gaf = d / "gaf" / "goa_{}.gaf.gz".format(sp["name"])
        for acc, term, aspect, evidence in parse_gaf_full(gaf):
            if evidence in exp_codes:
                has_exp.add(acc)
            else:
                nonexp[acc].add((term, aspect))
    dark = sorted(set(nonexp) - has_exp)
    print("Proteins with GO annotations but zero experimental evidence: {}".format(len(dark)))

    rng = random.Random(cfg["seed"] + 3)
    rng.shuffle(dark)
    sample = dark[: cfg["dark"]["oversample"]]

    # ---- sequences ---------------------------------------------------------
    print("Fetching sequences for {} dark proteins ...".format(len(sample)))
    seqs = {}
    for i in tqdm(range(0, len(sample), UNIPROT_BATCH), unit="batch"):
        seqs.update(fetch_uniprot_batch(session, sample[i : i + UNIPROT_BATCH]))
    lo, hi = cfg["sampling"]["min_len"], cfg["sampling"]["max_len"]
    seqs = {a: s for a, s in seqs.items() if lo <= len(s) <= hi}
    print("Kept {} after length filter".format(len(seqs)))

    # ---- structures --------------------------------------------------------
    min_plddt = cfg["structures"]["min_plddt"]
    accs = sorted(seqs)
    kept = []
    with ThreadPoolExecutor(max_workers=cfg["structures"]["workers"]) as pool:
        futures = [pool.submit(fetch_alphafold, session, a, dd / "structures", min_plddt)
                   for a in accs]
        for fut in tqdm(as_completed(futures), total=len(futures), unit="protein"):
            acc, _plddt, ok = fut.result()
            if ok:
                kept.append(acc)
    print("Models kept after pLDDT gate: {} / {}".format(len(kept), len(accs)))

    cap = cfg["dark"]["max_proteins"]
    if len(kept) > cap:
        rng.shuffle(kept)
        kept = kept[:cap]
    kept = sorted(kept)

    # ---- domains -----------------------------------------------------------
    print("Fetching InterPro entries for {} dark proteins ...".format(len(kept)))
    rows = []
    with ThreadPoolExecutor(max_workers=cfg["structures"]["workers"]) as pool:
        futures = {pool.submit(fetch_interpro, session, a): a for a in kept}
        for fut in tqdm(as_completed(futures), total=len(futures), unit="protein"):
            acc = futures[fut]
            try:
                entries = fut.result()
            except Exception as e:
                print("  warn: {}: {}".format(acc, e), file=sys.stderr)
                entries = []
            for ipr, etype in entries:
                rows.append([acc, ipr, etype])
            if not entries:
                rows.append([acc, "-", "-"])

    # ---- write -------------------------------------------------------------
    with open(dd / "accessions.txt", "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(kept) + "\n")
    write_fasta(dd / "dark.fasta", {a: seqs[a] for a in kept})
    rows.sort()
    with open(dd / "domains.tsv", "w", encoding="utf-8", newline="\n") as f:
        f.write("accession\tinterpro_id\tentry_type\n")
        for cols in rows:
            f.write("\t".join(cols) + "\n")
    with open(dd / "iea_annotations.tsv", "w", encoding="utf-8", newline="\n") as f:
        f.write("accession\tgo_term\taspect\n")
        for acc in kept:
            for term, aspect in sorted(nonexp.get(acc, ())):
                f.write("{}\t{}\t{}\n".format(acc, term, aspect))

    print("Dark set: {} proteins -> data/dark/".format(len(kept)))


if __name__ == "__main__":
    main()
