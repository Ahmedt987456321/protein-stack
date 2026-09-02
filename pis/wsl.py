"""Safe WSL invocation for this project.

Consolidates the hard-won lessons from the GPU/complexes run so they cannot
recur (see docs/WSL_GPU_NOTES.md):

  1. Path mangling. Git Bash rewrites POSIX paths like /root and /mnt when
     they are passed to `wsl`. Everything here goes through subprocess.run
     (Python, PowerShell-equivalent), which never mangles. NEVER call `wsl`
     from a Git Bash command line with a /root or /mnt argument.

  2. tmpfs /tmp. On this machine WSL mounts /tmp as a ~3.8 GB RAM disk.
     Staging thousands of files there exhausts RAM. Use WORK (=/root/work),
     which is on the 900+ GB ext4 root.

  3. RAM contention. 7.6 GB total system RAM cannot hold a large Foldseek
     search AND AlphaFold-Multimer at once; the OOM-killer terminates both.
     heavy() serializes such jobs with an flock so a second waits for the
     first, even across separately launched processes.

  4. WSL instability. run() can recover once from a crashed WSL service by
     `wsl --shutdown` and retrying.
"""
import subprocess
from pathlib import Path

WORK = "/root/work"                 # disk-backed staging (never /tmp)
HEAVY_LOCK = "/root/.protein_heavy.lock"


def to_wsl_path(winpath) -> str:
    """C:\\dev\\x -> /mnt/c/dev/x. Accepts str or Path."""
    p = Path(winpath).resolve()
    drive = p.drive.rstrip(":").lower()
    return "/mnt/{}/{}".format(drive, "/".join(p.parts[1:]))


def run(argv, timeout=None, check=False, recover=True):
    """Run a command inside WSL as root. argv is a list of Linux-side tokens.

    Returns subprocess.CompletedProcess. Never routed through a shell that
    mangles paths. On a WSL service crash, optionally shuts WSL down and
    retries once.
    """
    # errors="replace": ColabFold stdout carries non-cp1252 bytes (progress
    # glyphs); without this the reader thread throws UnicodeDecodeError on the
    # Windows default codec. utf-8/replace decodes it cleanly.
    cmd = ["wsl", "-u", "root"] + [str(a) for a in argv]
    dec = dict(encoding="utf-8", errors="replace")
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, **dec)
    if recover and r.returncode != 0 and "E_UNEXPECTED" in (r.stderr or ""):
        subprocess.run(["wsl", "--shutdown"], capture_output=True, text=True, **dec)
        subprocess.run(["wsl", "-u", "root", "true"], capture_output=True,
                       timeout=60, **dec)  # re-init
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, **dec)
    if check and r.returncode != 0:
        raise RuntimeError("WSL command failed ({}): {}".format(
            r.returncode, (r.stderr or r.stdout or "")[-400:]))
    return r


def run_script(bash_text, timeout=None, check=True, heavy=False,
               require_mb=0, out_path=None):
    """Write bash to a LF file on the Windows side, run it in WSL.

    Avoids inline-quoting bugs entirely. If heavy=True the whole script is
    wrapped in an flock on HEAVY_LOCK so it serializes with other heavy
    jobs. If require_mb>0, refuses to start unless that much RAM is free.
    """
    if require_mb and available_ram_mb() < require_mb:
        raise RuntimeError("insufficient free RAM: need {} MB, have {} MB "
                           "(close other heavy jobs)".format(
                               require_mb, available_ram_mb()))
    sh = out_path or (Path(__file__).resolve().parent.parent
                      / "data" / "explore" / "_wsl_job.sh")
    sh.parent.mkdir(parents=True, exist_ok=True)
    body = "set -e\nmkdir -p {}\n".format(WORK) + bash_text
    with open(sh, "w", encoding="utf-8", newline="\n") as f:
        f.write(body)
    wsl_sh = to_wsl_path(sh)
    if heavy:
        argv = ["bash", "-c",
                "exec 9>{lock}; flock 9; bash {sh}".format(lock=HEAVY_LOCK, sh=wsl_sh)]
    else:
        argv = ["bash", wsl_sh]
    return run(argv, timeout=timeout, check=check)


def available_ram_mb() -> int:
    """Free RAM in MB. Parses `free -m` in Python (WSL mangles inline awk)."""
    r = run(["free", "-m"])
    for line in (r.stdout or "").splitlines():
        parts = line.split()
        if parts and parts[0].rstrip(":") == "Mem":
            # total used free shared buff/cache available
            return int(parts[6]) if len(parts) >= 7 else int(parts[3])
    return 0


def gpu_ok() -> bool:
    """True if a CUDA device is visible to the ColabFold JAX."""
    py = "/root/localcolabfold/colabfold-conda/bin/python"
    r = run([py, "-c",
             "import jax;print(any('cuda' in str(d).lower() for d in jax.devices()))"])
    return "True" in (r.stdout or "")
