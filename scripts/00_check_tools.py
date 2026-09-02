"""Preflight: verify Python deps and external tools before running the pipeline."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

ok = True

for mod in ("requests", "yaml", "tqdm"):
    try:
        __import__(mod)
        print("[ok]   python module: " + mod)
    except ImportError:
        print("[MISS] python module: {}  ->  pip install -r requirements.txt".format(mod))
        ok = False

if ok:
    from pis.common import load_config, tool_available

    cfg = load_config()
    for key in ("mmseqs", "foldseek"):
        if tool_available(cfg, key):
            print("[ok]   tool: " + key)
        else:
            where = "inside WSL" if cfg["tools"].get("use_wsl") else "on PATH"
            print("[MISS] tool: {} not found {} (needed from step 04 onward)".format(key, where))
            print("       install: conda install -c bioconda mmseqs2 foldseek")
            print("       on Windows: install in WSL and set tools.use_wsl: true in config.yaml")
            ok = False

    # WSL/GPU preflight - surfaces the operational hazards (see
    # docs/WSL_GPU_NOTES.md) before a heavy run rather than mid-crash.
    if cfg["tools"].get("foldseek_wsl") or cfg["tools"].get("use_wsl"):
        try:
            from pis import wsl
            r = wsl.run(["bash", "-c", "findmnt -n -o FSTYPE /tmp || true"])
            if "tmpfs" in (r.stdout or ""):
                print("[note] WSL /tmp is tmpfs (RAM disk) - staging is on "
                      "/root/work instead (handled in code)")
            free = wsl.available_ram_mb()
            print(("[ok]   " if free >= 4000 else "[warn] ")
                  + "WSL free RAM: {} MB".format(free)
                  + ("" if free >= 4000 else " - heavy jobs may OOM; serialised via flock"))
            disk = wsl.run(["bash", "-c", "df -BG --output=avail /root | tail -1"])
            print("[ok]   WSL /root free:" + (disk.stdout or "").strip().replace("\n", " "))
            if Path("/root").exists() or True:  # GPU is optional
                print(("[ok]   " if wsl.gpu_ok() else "[note] ")
                      + ("GPU (CUDA JAX) available" if wsl.gpu_ok()
                         else "GPU not set up (run scripts/setup_gpu.sh for AF2-Multimer)"))
        except Exception as e:
            print("[note] WSL preflight skipped: {}".format(str(e)[:80]))

print()
print("Ready." if ok else "Fix the [MISS] items above. Steps 01-03 need only the Python modules.")
sys.exit(0 if ok else 1)
