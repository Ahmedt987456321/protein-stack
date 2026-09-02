"""Run the full v0.1 pipeline in order. Stops at the first failing step.

Steps are idempotent where it matters: downloads are cached, so re-running
after a failure resumes cheaply.
"""
import subprocess
import sys
from pathlib import Path

STEPS = [
    "scripts/00_check_tools.py",
    "scripts/01_fetch_annotations.py",
    "scripts/02_fetch_sequences.py",
    "scripts/03_fetch_structures.py",
    "scripts/04_split.py",
    "scripts/05_search.py",
    "scripts/06_predict.py",
    "scripts/07_evaluate.py",
    "scripts/08_fetch_domains.py",
    "scripts/09_predict_v2.py",
    "scripts/10_evaluate_v2.py",
    "scripts/11_dark_candidates.py",
    "scripts/12_dark_sweep.py",
    "scripts/14_fetch_interactions.py",
    "scripts/15_pockets.py",
    "scripts/16_interaction_experiment.py",
    "scripts/19_time_machine.py",
    "scripts/13_build_graph.py",  # graph assembly last: folds in every layer present
    "scripts/18_agent_eval.py",   # agent grounding gate (defers without LLM credentials)
]

root = Path(__file__).resolve().parent
for step in STEPS:
    print("\n" + "=" * 70)
    print("RUNNING " + step)
    print("=" * 70, flush=True)
    result = subprocess.run([sys.executable, str(root / step)], cwd=root)
    if result.returncode != 0:
        print("\nStep failed: {} (exit {})".format(step, result.returncode))
        sys.exit(result.returncode)

print("\nPipeline complete. See results/report.md")
