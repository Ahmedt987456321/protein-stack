"""LLM reasoning pass for the v0.5 agent - provider-neutral.

Set the LLM_CLI environment variable to any command-line LLM client that
reads a prompt on stdin and prints the completion to stdout (most vendor
CLIs have such a mode; include any flags in the variable, e.g.
LLM_CLI="somecli --plain"). No provider is assumed or required.

llm_available() reports whether a working CLI is configured;
llm_answer() returns the model's text or an "[unavailable: ...]" marker
string (never raises).
"""
import os
import shutil
import subprocess

SYSTEM_PROMPT = """You are a protein-science assistant reasoning over an \
evidence dossier assembled from a provenance-typed knowledge graph.

Hard rules:
- Cite only GO terms, domains, proteins, and scores that appear in the dossier.
- Respect evidence tiers: EXPERIMENTAL findings may be stated as fact; \
CURATED as strong support; COMPUTATIONAL as supporting signals; HYPOTHESIS \
must always be presented as a hypothesis needing experimental validation, \
never as established function.
- For a DARK protein, function statements must be framed as hypotheses.
- If evidence is thin or conflicting, say so plainly.
- End with a short "Suggested validation" line proposing one concrete \
experiment when hypotheses are involved.
Answer concisely (under 250 words)."""

UNAVAILABLE = "[unavailable: "


def _cli():
    """Resolve LLM_CLI into an argv list, or None if unset/not found."""
    raw = os.environ.get("LLM_CLI", "").strip()
    if not raw:
        return None
    parts = raw.split()
    exe = shutil.which(parts[0])
    if exe is None:
        return None
    return [exe] + parts[1:]


def llm_available():
    """Return ('cli'|None, detail)."""
    cmd = _cli()
    if cmd is None:
        return None, ("no LLM configured - set the LLM_CLI environment "
                      "variable to a CLI that reads a prompt on stdin and "
                      "prints the completion")
    try:
        probe = subprocess.run(
            cmd, input="reply with exactly: PONG",
            capture_output=True, text=True, encoding="utf-8", timeout=300,
        )
    except Exception as e:
        return None, "LLM_CLI probe failed: {}".format(e)
    if probe.returncode == 0 and "PONG" in probe.stdout:
        return "cli", "LLM CLI: " + " ".join(cmd)
    return None, "LLM_CLI is set but the probe call failed (exit {})".format(probe.returncode)


def llm_answer(dossier_md: str, question: str) -> str:
    cmd = _cli()
    if cmd is None:
        return UNAVAILABLE + "no LLM_CLI configured]"
    prompt = "{}\n\n---\n\nQuestion: {}".format(dossier_md, question)
    try:
        result = subprocess.run(
            cmd, input=SYSTEM_PROMPT + "\n\n" + prompt,
            capture_output=True, text=True, encoding="utf-8", timeout=600,
        )
    except Exception as e:
        return UNAVAILABLE + "LLM CLI error: {}]".format(e)
    if result.returncode != 0:
        detail = (result.stderr.strip() or result.stdout.strip())[:300]
        return UNAVAILABLE + "LLM CLI: " + detail + "]"
    return result.stdout.strip()
