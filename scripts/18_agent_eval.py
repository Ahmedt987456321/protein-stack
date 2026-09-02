"""Step 18 (v0.5) - the agent grounding gate.

Samples 5 annotated + 5 dark proteins, asks the agent the standard question,
and mechanically checks each answer:

  G-a (grounding): every GO id the answer cites appears in the dossier the
      agent was given.
  G-b (tier discipline): for dark proteins the answer frames function as a
      hypothesis (mentions "hypothes-"/"predicted"/"putative"); it never
      claims experimental confirmation the dossier does not contain.

Gate: >= 9 of 10 answers pass both checks.

If no LLM is configured (LLM_CLI environment variable unset or not
responding), the gate is recorded as DEFERRED - re-run this script after
configuring LLM_CLI (see pis/agent_llm.py).

Outputs: results/report_v05.md (+ per-answer transcripts in results/agent_eval/)
"""
import random
import re
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pis.agent_llm import UNAVAILABLE, llm_answer, llm_available
from pis.common import data_dir, load_config, results_dir
from pis.dossier import build_dossier, render_markdown
from pis.go import GoDag

QUESTION = "What does this protein likely do, and how confident should I be?"
GO_RE = re.compile(r"GO:\d{7}")
HYP_RE = re.compile(r"hypothes|predict|putativ|suggest|candidate", re.IGNORECASE)
OVERCLAIM_RE = re.compile(r"experimentally (?:confirmed|verified|shown|demonstrated)",
                          re.IGNORECASE)


def dossier_go_ids(d):
    ids = set()
    for key in ("experimental_functions", "electronic_annotations", "hypotheses"):
        ids |= {x["term"] for x in d[key]}
    ids |= {x["term"] for x in d["curated_domain_implications"]}
    for key in ("sequence_neighbours", "structure_neighbours", "interaction_partners"):
        for n in d[key]:
            ids |= {f["term"] for f in n["partner_functions"]}
    return ids


def main():
    cfg = load_config()
    d = data_dir(cfg)
    r = results_dir(cfg)
    out_dir = r / "agent_eval"
    out_dir.mkdir(exist_ok=True)

    mode, detail = llm_available()
    if mode is None:
        with open(r / "report_v05.md", "w", encoding="utf-8", newline="\n") as f:
            f.write("# v0.5 - agent - report\n\n")
            f.write("Agent infrastructure complete (dossier layer + LLM "
                    "reasoning pass). **Gate: DEFERRED** - {}.\n\nRe-run "
                    "`python scripts/18_agent_eval.py` after configuring an "
                    "LLM; the gate then completes in ~5 minutes.\n".format(detail))
        print("GATE DEFERRED: " + detail)
        return

    print("LLM path: {} ({})".format(mode, detail))
    dag = GoDag(d / "go-basic.obo")
    con = sqlite3.connect("file:{}?mode=ro".format((d / "kg.db").as_posix()), uri=True)

    annotated = [row[0] for row in con.execute(
        "SELECT DISTINCT subject FROM edges WHERE predicate='has_function'")]
    dark = [row[0] for row in con.execute(
        "SELECT DISTINCT subject FROM edges WHERE predicate='predicted_function'")]
    rng = random.Random(cfg["seed"] + 6)
    sample = sorted(rng.sample(sorted(annotated), 5)) + sorted(rng.sample(sorted(dark), 5))

    rows = []
    for acc in sample:
        dossier = build_dossier(con, acc, dag)
        md = render_markdown(dossier)
        answer = llm_answer(md, QUESTION)
        (out_dir / (acc + ".md")).write_text(
            md + "\n\n---\n# Answer\n\n" + answer, encoding="utf-8")

        if answer.startswith(UNAVAILABLE):
            rows.append((acc, dossier["is_dark"], False, False, "LLM error"))
            continue
        cited = set(GO_RE.findall(answer))
        allowed = dossier_go_ids(dossier)
        grounded = cited <= allowed
        if dossier["is_dark"]:
            tiered = bool(HYP_RE.search(answer)) and not OVERCLAIM_RE.search(answer)
        else:
            tiered = not (OVERCLAIM_RE.search(answer) and not dossier["experimental_functions"])
        note = "" if grounded else "cites unknown terms: {}".format(
            sorted(cited - allowed)[:5])
        rows.append((acc, dossier["is_dark"], grounded, tiered, note))
        print("{}  dark={}  grounded={}  tier-ok={}  {}".format(
            acc, dossier["is_dark"], grounded, tiered, note))

    passed = sum(1 for _a, _d, g, t, _n in rows if g and t)
    verdict = "PASSED" if passed >= 9 else "NOT PASSED"

    with open(r / "report_v05.md", "w", encoding="utf-8", newline="\n") as f:
        f.write("# v0.5 - agent - report\n\n")
        f.write("LLM path: {} ({}). Question: \"{}\"\n\n".format(mode, detail, QUESTION))
        f.write("| accession | dark | grounded | tier discipline | note |\n")
        f.write("|---|---|---|---|---|\n")
        for acc, is_dark, g, t, note in rows:
            f.write("| {} | {} | {} | {} | {} |\n".format(
                acc, "yes" if is_dark else "no",
                "PASS" if g else "FAIL", "PASS" if t else "FAIL", note))
        f.write("\n**Gate ({}/10 passed, threshold 9): {}**\n\n".format(passed, verdict))
        f.write("Grounding = every cited GO id exists in the dossier given to "
                "the model. Tier discipline = dark-protein function framed as "
                "hypothesis, no fabricated experimental confirmation. Answer "
                "transcripts: results/agent_eval/.\n")
    print("\nGate: {}/10 -> {}".format(passed, verdict))
    print("Wrote results/report_v05.md")


if __name__ == "__main__":
    main()
