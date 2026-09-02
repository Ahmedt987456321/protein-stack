"""Exploration B2 prep - fetch predicted aligned error (PAE) matrices for
the 1,832 test-set proteins. Resumable; files land in data/explore/pae/.
"""
import json
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tqdm import tqdm

from pis.common import data_dir, http_session, load_config

API = "https://alphafold.ebi.ac.uk/api/prediction/{}"


def fetch_one(session, acc, out_dir):
    dest = out_dir / (acc + ".json")
    if dest.exists() and dest.stat().st_size > 0:
        return True
    try:
        r = session.get(API.format(acc), timeout=60)
        if r.status_code != 200:
            return False
        entries = r.json()
        url = entries[0].get("paeDocUrl") if entries else None
        if not url:
            return False
        pae = session.get(url, timeout=300)
        pae.raise_for_status()
        tmp = dest.with_suffix(".json.part")
        tmp.write_bytes(pae.content)
        tmp.replace(dest)
        return True
    except Exception:
        return False


def main():
    cfg = load_config()
    d = data_dir(cfg)
    out_dir = d / "explore" / "pae"
    out_dir.mkdir(parents=True, exist_ok=True)
    session = http_session()

    test = []
    with open(d / "split.csv", encoding="utf-8") as f:
        next(f)
        for line in f:
            acc, split, _b = line.rstrip("\n").split(",")
            if split == "test":
                test.append(acc)
    print("fetching PAE for {} test proteins ...".format(len(test)))
    ok = 0
    with ThreadPoolExecutor(max_workers=cfg["structures"]["workers"]) as pool:
        futs = [pool.submit(fetch_one, session, a, out_dir) for a in test]
        for fut in tqdm(as_completed(futs), total=len(futs), unit="pae"):
            ok += bool(fut.result())
    print("fetched {}/{}".format(ok, len(test)))


if __name__ == "__main__":
    main()
