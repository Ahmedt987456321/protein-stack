"""Shared plumbing: config, HTTP with retries, FASTA IO, external tool runner."""
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import requests
import yaml
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

ROOT = Path(__file__).resolve().parent.parent


def load_config():
    with open(ROOT / "config.yaml", encoding="utf-8") as f:
        return yaml.safe_load(f)


def data_dir(cfg) -> Path:
    p = ROOT / cfg["dirs"]["data"]
    p.mkdir(parents=True, exist_ok=True)
    return p


def results_dir(cfg) -> Path:
    p = ROOT / cfg["dirs"]["results"]
    p.mkdir(parents=True, exist_ok=True)
    return p


def http_session() -> requests.Session:
    s = requests.Session()
    retry = Retry(
        total=5,
        backoff_factor=1.0,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
    )
    adapter = HTTPAdapter(max_retries=retry, pool_connections=16, pool_maxsize=16)
    s.mount("https://", adapter)
    s.mount("http://", adapter)
    s.headers["User-Agent"] = "protein-stack/0.1 (research pipeline)"
    return s


def download(session: requests.Session, url: str, dest: Path) -> Path:
    """Download url to dest, resuming nothing but skipping files already present."""
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with session.get(url, stream=True, timeout=180) as r:
        r.raise_for_status()
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(1 << 20):
                f.write(chunk)
    tmp.replace(dest)
    return dest


# ---------------------------------------------------------------- FASTA ----

_ACC_HEADER = re.compile(r"^(?:sp|tr)\|([^|]+)\|")


def fasta_accession(header: str) -> str:
    """'sp|P12345|NAME ...' -> 'P12345'; otherwise first whitespace token."""
    m = _ACC_HEADER.match(header)
    if m:
        return m.group(1).split("-")[0]
    return header.split()[0].split("-")[0]


def read_fasta(path: Path):
    """Return {accession: sequence}."""
    records = {}
    acc, chunks = None, []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if acc is not None:
                    records[acc] = "".join(chunks)
                acc, chunks = fasta_accession(line[1:]), []
            elif line:
                chunks.append(line.strip())
    if acc is not None:
        records[acc] = "".join(chunks)
    return records


def write_fasta(path: Path, records) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        for acc, seq in records.items():
            f.write(">" + acc + "\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i : i + 60] + "\n")


def norm_hit_id(s: str) -> str:
    """Normalize an mmseqs/foldseek hit name to a bare accession.

    Handles 'P12345', 'P12345.pdb', 'P12345.pdb_A', and stray path prefixes.
    UniProt accessions never contain '.' or '_', so splitting on either is safe.
    """
    s = s.replace("\\", "/").split("/")[-1]
    return re.split(r"[._]", s)[0]


# ------------------------------------------------------------- tooling ----


def to_wsl_path(p: Path) -> str:
    """C:\\dev\\x -> /mnt/c/dev/x"""
    p = Path(p).resolve()
    drive = p.drive.rstrip(":").lower()
    rest = "/".join(p.parts[1:])
    return "/mnt/{}/{}".format(drive, rest)


def tool_uses_wsl(cfg, tool_key: str) -> bool:
    """Per-tool override (e.g. foldseek_wsl: true) falls back to use_wsl."""
    per_tool = cfg["tools"].get(tool_key + "_wsl")
    if per_tool is not None:
        return bool(per_tool)
    return bool(cfg["tools"].get("use_wsl"))


def run_tool(cfg, tool_key: str, args) -> None:
    """Run mmseqs/foldseek with the configured binary, optionally through WSL.

    Path arguments are converted to /mnt/... form automatically when the tool
    runs inside WSL, so callers just pass Path objects.
    """
    wsl = tool_uses_wsl(cfg, tool_key)
    binary = cfg["tools"][tool_key]
    conv = to_wsl_path if wsl else str
    cmd = [binary] + [conv(a) if isinstance(a, Path) else str(a) for a in args]
    if wsl:
        cmd = ["wsl"] + cmd
    print("+ " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def tool_available(cfg, tool_key: str) -> bool:
    binary = cfg["tools"][tool_key]
    cmd = (["wsl"] if tool_uses_wsl(cfg, tool_key) else []) + [binary, "version"]
    try:
        subprocess.run(cmd, capture_output=True, check=True, timeout=60)
        return True
    except Exception:
        return False


def link_or_copy(src: Path, dst: Path) -> None:
    """Hardlink when possible (same NTFS volume, no admin needed), else copy."""
    if dst.exists():
        return
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)


def fail(msg: str) -> None:
    print("ERROR: " + msg, file=sys.stderr)
    sys.exit(1)
