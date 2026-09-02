# WSL / GPU operational notes

Four failures hit this project's GPU and heavy-search runs on a Windows +
WSL2 + laptop-GPU machine. Each is now prevented in code (`pis/wsl.py`) or
documented here so it does not recur.

## 1. Git Bash mangles POSIX paths passed to `wsl`

Calling `wsl ... /root/x` or `wsl ... /mnt/c/x` from a Git Bash command line
rewrites the path to `C:/Program Files/Git/root/x`. Symptom: "No such file
or directory" for a path that plainly exists.

**Fix:** never invoke `wsl` with a `/root` or `/mnt` argument from Git Bash.
Use `pis.wsl.run(argv)` / `run_script(...)` (Python subprocess) or PowerShell.
Both pass arguments to `wsl.exe` verbatim.

## 2. `/tmp` is a RAM disk (tmpfs)

This WSL mounts `/tmp` as a ~3.8 GB tmpfs. Staging thousands of structures
there exhausts RAM and fails with "No space left on device" (and can crash
WSL). The ext4 root (`/`, `/root`) has 900+ GB.

**Fix:** stage on `pis.wsl.WORK` (`/root/work`), never `/tmp`. All heavy
copy-then-process steps use it.

## 3. Concurrent heavy jobs trigger the OOM-killer

7.6 GB system RAM cannot hold a large Foldseek search AND an
AlphaFold-Multimer run at once. `dmesg` shows `Out of memory: Killed process
... foldseek`; both jobs die.

**Fix:** `pis.wsl.run_script(..., heavy=True)` wraps the job in an `flock` on
`/root/.protein_heavy.lock`, so heavy jobs serialize even when launched from
separate processes. `require_mb=N` additionally refuses to start below a RAM
floor.

## 4. localColabFold ships CPU-only JAX

The installer leaves a CPU jaxlib; the GPU is visible to WSL but unused
("a CUDA-enabled jaxlib is not installed. Falling back to cpu").

**Fix:** `scripts/setup_gpu.sh` installs `jax[cuda12]` matched to the
installed jax version and verifies `jax.devices()` shows a `CudaDevice`.
`pis.wsl.gpu_ok()` checks this at runtime.

## Card-specific limits (RTX 3050 Ti, 4 GB)

AlphaFold-Multimer memory scales with the square of combined sequence
length. 4 GB fits complexes up to roughly 450-600 combined residues; beyond
that it OOMs on the GPU. Set `TF_FORCE_UNIFIED_MEMORY=1` and
`XLA_PYTHON_CLIENT_MEM_FRACTION=3.0` to allow spillover into system RAM at a
speed cost. Keep candidate complexes small and fold sequentially.

## Recovery

`pis.wsl.run(...)` retries once through `wsl --shutdown` if a call fails with
`E_UNEXPECTED` (a crashed WSL service). If WSL is wedged, from Windows:
`wsl --shutdown` then re-run.
