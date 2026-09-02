"""Shared fetchers: UniProt sequences, AlphaFold models, InterPro entries."""
import sys
from pathlib import Path

from pis.common import fasta_accession

UNIPROT_BATCH = 100  # UniProt accessions endpoint limit
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/accessions"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{}"
INTERPRO_API = "https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/{}?page_size=200"


def fetch_uniprot_batch(session, accs):
    """Fetch canonical FASTA for up to UNIPROT_BATCH accessions -> {acc: seq}."""
    r = session.get(
        UNIPROT_URL,
        params={"accessions": ",".join(accs), "format": "fasta"},
        timeout=180,
    )
    r.raise_for_status()
    records = {}
    acc, chunks = None, []
    for line in r.text.splitlines():
        if line.startswith(">"):
            if acc is not None:
                records[acc] = "".join(chunks)
            acc, chunks = fasta_accession(line[1:]), []
        elif line:
            chunks.append(line.strip())
    if acc is not None:
        records[acc] = "".join(chunks)
    return records


def mean_ca_plddt(pdb_path: Path):
    total, n = 0.0, 0
    with open(pdb_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    total += float(line[60:66])
                    n += 1
                except ValueError:
                    pass
    return (total / n) if n else 0.0


def fetch_alphafold(session, acc, out_dir: Path, min_plddt: float):
    """Download the AlphaFold model for acc into out_dir/<acc>.pdb.

    Returns (acc, mean_plddt, kept). mean_plddt < 0 means no model available.
    """
    dest = out_dir / (acc + ".pdb")
    try:
        if not dest.exists():
            r = session.get(ALPHAFOLD_API.format(acc), timeout=60)
            if r.status_code == 404:
                return acc, -1.0, False
            r.raise_for_status()
            entries = r.json()
            if not entries or "pdbUrl" not in entries[0]:
                return acc, -1.0, False
            pdb = session.get(entries[0]["pdbUrl"], timeout=180)
            pdb.raise_for_status()
            tmp = dest.with_suffix(".pdb.part")
            tmp.write_bytes(pdb.content)
            tmp.replace(dest)
        plddt = mean_ca_plddt(dest)
        if plddt < min_plddt:
            dest.unlink(missing_ok=True)
            return acc, plddt, False
        return acc, plddt, True
    except Exception as e:
        print("  warn: {}: {}".format(acc, e), file=sys.stderr)
        dest.unlink(missing_ok=True)
        return acc, -1.0, False


def fetch_interpro(session, acc):
    """Return list of (interpro_id, entry_type); empty list if none."""
    out = []
    url = INTERPRO_API.format(acc)
    while url:
        r = session.get(url, timeout=60)
        if r.status_code in (204, 404):
            return out
        r.raise_for_status()
        js = r.json()
        for res in js.get("results", []):
            meta = res.get("metadata", {})
            if meta.get("accession"):
                out.append((meta["accession"], meta.get("type", "")))
        url = js.get("next")
    return out
