"""Step 17 (v0.5) - the protein agent.

Deterministic engines gather an evidence dossier from the knowledge graph;
an LLM reasons over it. Usage:

  python scripts/17_agent.py P69905                 # print the dossier only
  python scripts/17_agent.py P69905 --ask "what does this protein likely do?"

LLM access is provider-neutral: set the LLM_CLI environment variable to any
command-line LLM client that reads a prompt on stdin and prints the
completion (see pis/agent_llm.py). Without it, the dossier still prints.
"""
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pis.agent_llm import llm_answer
from pis.common import data_dir, load_config
from pis.dossier import build_dossier, render_markdown
from pis.go import GoDag


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if not args:
        print("usage: 17_agent.py ACCESSION [--ask \"question\"]")
        sys.exit(2)
    acc = args[0]
    question = None
    if "--ask" in sys.argv:
        question = sys.argv[sys.argv.index("--ask") + 1]

    cfg = load_config()
    d = data_dir(cfg)
    dag = GoDag(d / "go-basic.obo")
    con = sqlite3.connect("file:{}?mode=ro".format((d / "kg.db").as_posix()), uri=True)

    dossier = build_dossier(con, acc, dag)
    if dossier is None:
        print("No evidence in the knowledge graph for {}".format(acc))
        sys.exit(1)
    md = render_markdown(dossier)
    print(md)

    if question:
        print("\n---\n# Agent answer\n")
        print(llm_answer(md, question))


if __name__ == "__main__":
    main()
