"""Exploration 3 - druggability-to-disease target triage.

Human dark proteins with a druggable pocket (fpocket best score >= 0.5) are
mapped UniProt -> Ensembl gene (UniProt ID-mapping API) and crossed with
Open Targets disease associations. Output is a ranked triage table joining:
no experimental annotation + druggable pocket + disease association +
predicted function.

If the Open Targets API is unreachable the table degrades to
druggability + hypothesis only, and the degradation is recorded.

Outputs: results/explore/target_triage.md, target_triage.tsv
"""
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tqdm import tqdm

from pis.common import data_dir, http_session, load_config
from pis.go import GoDag, parse_gaf_full

OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"
OT_QUERY = """query($id: String!) {
  target(ensemblId: $id) {
    approvedSymbol
    associatedDiseases(page: {index: 0, size: 3}) {
      count
      rows { disease { name } score }
    }
  }
}"""


def uniprot_to_ensembl(session, accs):
    """Map via mygene.info in batches (reliable, no job queue)."""
    mapping = {}
    for i in range(0, len(accs), 500):
        chunk = accs[i:i + 500]
        for attempt in range(3):
            try:
                r = session.post("https://mygene.info/v3/query",
                                 data={"q": ",".join(chunk),
                                       "scopes": "uniprot",
                                       "fields": "ensembl.gene,symbol",
                                       "species": "human"},
                                 timeout=120)
                r.raise_for_status()
                for row in r.json():
                    ens = row.get("ensembl")
                    if isinstance(ens, list):
                        ens = ens[0]
                    gene = (ens or {}).get("gene")
                    if gene and row.get("query") not in mapping:
                        mapping[row["query"]] = gene
                break
            except Exception:
                time.sleep(5 * (attempt + 1))
    return mapping


def main():
    cfg = load_config()
    d = data_dir(cfg)
    out_dir = Path("results/explore")
    out_dir.mkdir(parents=True, exist_ok=True)
    session = http_session()
    dag = GoDag(d / "go-basic.obo")

    # human dark accessions
    dark = set((d / "dark" / "accessions.txt").read_text(encoding="utf-8").split())
    human = set()
    for acc, _t, _a, _e in parse_gaf_full(d / "gaf" / "goa_human.gaf.gz"):
        if acc in dark:
            human.add(acc)

    pockets = {}
    with open(d / "dark" / "pockets.tsv", encoding="utf-8") as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            pockets[c[0]] = float(c[2])
    cands = sorted(a for a in human if pockets.get(a, 0) >= 0.5)
    print("human dark proteins with druggable pocket: {}".format(len(cands)))

    top_mf = {}
    with open(Path("results") / "dark_hypotheses.tsv", encoding="utf-8") as f:
        next(f)
        for line in f:
            acc, br, term, score, tier, _st, _so = line.rstrip("\n").split("\t")
            if br == "molecular_function" and acc in human:
                if acc not in top_mf or float(score) > top_mf[acc][1]:
                    top_mf[acc] = (term, float(score), tier)

    print("mapping {} accessions to Ensembl ...".format(len(cands)))
    mapping = uniprot_to_ensembl(session, cands)
    print("mapped: {}".format(len(mapping)))

    rows = []
    ot_failures = 0
    for acc in tqdm(cands, unit="target"):
        ens = mapping.get(acc)
        symbol, dis_name, dis_score, dis_count = "-", "-", 0.0, 0
        if ens:
            for attempt in range(3):
                try:
                    r = session.post(OT_URL, json={"query": OT_QUERY,
                                                   "variables": {"id": ens}},
                                     timeout=30)
                    js = r.json().get("data", {}).get("target")
                    if js:
                        symbol = js.get("approvedSymbol") or "-"
                        ad = js.get("associatedDiseases") or {}
                        dis_count = ad.get("count", 0)
                        rws = ad.get("rows") or []
                        if rws:
                            dis_name = rws[0]["disease"]["name"]
                            dis_score = round(rws[0]["score"], 3)
                    break
                except Exception:
                    if attempt == 2:
                        ot_failures += 1
                    else:
                        time.sleep(3 * (attempt + 1))
        hyp = top_mf.get(acc)
        rows.append([acc, symbol, ens or "-", "{:.2f}".format(pockets[acc]),
                     str(dis_count), "{:.3f}".format(dis_score), dis_name,
                     hyp[0] if hyp else "-",
                     dag.name(hyp[0]) if hyp else "-",
                     hyp[2] if hyp else "-"])

    rows.sort(key=lambda r: (-float(r[5]), -float(r[3])))
    with open(out_dir / "target_triage.tsv", "w", encoding="utf-8", newline="\n") as f:
        f.write("accession\tsymbol\tensembl\tdruggability\tdisease_assoc_count\t"
                "top_disease_score\ttop_disease\ttop_MF_term\tterm_name\ttier\n")
        for r in rows:
            f.write("\t".join(r) + "\n")

    n_dis = sum(1 for r in rows if float(r[5]) > 0)
    with open(out_dir / "target_triage.md", "w", encoding="utf-8", newline="\n") as f:
        f.write("# Druggability-to-disease triage (human dark proteins)\n\n")
        f.write("{} human dark proteins carry a druggable pocket (>= 0.5); "
                "{} mapped to Ensembl; {} have at least one Open Targets "
                "disease association.{}\n\n".format(
                    len(cands), len(mapping), n_dis,
                    "" if ot_failures == 0 else
                    " NOTE: {} targets could not be queried after retries.".format(ot_failures)))
        f.write("Top 20 by disease association score:\n\n")
        f.write("| accession | symbol | druggability | top disease | score "
                "| predicted MF | tier |\n|---|---|---|---|---|---|---|\n")
        for r in rows[:20]:
            f.write("| {} | {} | {} | {} | {} | {} | {} |\n".format(
                r[0], r[1], r[3], r[6][:40], r[5], r[8][:40], r[9]))
    print("with disease association: {} | query failures: {}".format(n_dis, ot_failures))
    print("Wrote results/explore/target_triage.md")


if __name__ == "__main__":
    main()
