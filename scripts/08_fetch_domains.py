"""Step 08 (v0.2) - fetch InterPro entries for the final protein set, plus the
curated InterPro->GO mapping.

Uses the InterPro API per protein (paginated). Resumable: accessions already
present in data/domains.tsv are skipped on re-run.

Outputs:
  data/domains.tsv       (accession, interpro_id, entry_type)
  data/interpro2go.txt   curated mapping (raw file)
"""
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tqdm import tqdm

from pis.common import data_dir, download, http_session, load_config
from pis.fetch import fetch_interpro as fetch_domains

IP2GO_URL = "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/interpro2go"


def main():
    cfg = load_config()
    d = data_dir(cfg)
    session = http_session()

    print("Downloading interpro2go mapping ...")
    download(session, IP2GO_URL, d / "interpro2go.txt")

    accs = (d / "final_accessions.txt").read_text(encoding="utf-8").split()

    out_path = d / "domains.tsv"
    done = set()
    rows = []
    if out_path.exists():
        with open(out_path, encoding="utf-8") as f:
            next(f, None)
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if cols:
                    done.add(cols[0])
                    rows.append(cols)
    todo = [a for a in accs if a not in done]
    print("InterPro entries: {} cached, fetching {} ...".format(len(done), len(todo)))

    failures = 0
    with ThreadPoolExecutor(max_workers=cfg["structures"]["workers"]) as pool:
        futures = {pool.submit(fetch_domains, session, a): a for a in todo}
        for fut in tqdm(as_completed(futures), total=len(futures), unit="protein"):
            acc = futures[fut]
            try:
                entries = fut.result()
            except Exception as e:
                print("  warn: {}: {}".format(acc, e), file=sys.stderr)
                failures += 1
                continue
            if entries:
                for ipr, etype in entries:
                    rows.append([acc, ipr, etype])
            else:
                rows.append([acc, "-", "-"])  # mark as fetched-but-empty

    rows.sort()
    with open(out_path, "w", encoding="utf-8", newline="\n") as f:
        f.write("accession\tinterpro_id\tentry_type\n")
        for cols in rows:
            f.write("\t".join(cols) + "\n")

    n_with = len({c[0] for c in rows if c[1] != "-"})
    print("Proteins with >=1 InterPro entry: {} / {} ({} fetch failures)".format(
        n_with, len(accs), failures))
    if failures:
        print("Re-run this script to retry failed accessions.")


if __name__ == "__main__":
    main()
