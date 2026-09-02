"""Step 13 (v0.3) - build the provenance-first protein knowledge graph.

Every edge carries its evidence type - the architecture's non-negotiable:

  has_function            EXPERIMENTAL   GOA (experimental codes only)
  electronic_function     COMPUTATIONAL  GOA (IEA-class, dark proteins)
  has_domain              COMPUTATIONAL  InterPro
  domain_implies_function CURATED        interpro2go
  sequence_similar_to     COMPUTATIONAL  MMseqs2 (score = identity)
  structure_similar_to    COMPUTATIONAL  Foldseek (score = TM-score)
  predicted_function      HYPOTHESIS  fusion v0.2 (score = fused confidence)

Outputs:
  data/kg.db            SQLite; table edges(subject, predicate, object,
                        evidence, source, score)
  results/kg_stats.md
"""
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pis.common import data_dir, load_config, results_dir
from pis.streams import load_interpro2go_raw


def iter_tsv(path, ncols):
    with open(path, encoding="utf-8") as f:
        next(f)
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= ncols:
                yield cols


def main():
    cfg = load_config()
    d = data_dir(cfg)
    dd = d / "dark"
    r = results_dir(cfg)

    db_path = d / "kg.db"
    if db_path.exists():
        db_path.unlink()
    con = sqlite3.connect(str(db_path))
    con.execute("CREATE TABLE edges (subject TEXT, predicate TEXT, object TEXT, "
                "evidence TEXT, source TEXT, score REAL)")

    def add(rows):
        con.executemany("INSERT INTO edges VALUES (?,?,?,?,?,?)", rows)

    # experimental function annotations (the 14K sampled annotated proteins)
    add((acc, "has_function", term, "EXPERIMENTAL", "GOA", 1.0)
        for acc, term, _a in iter_tsv(d / "annotations.tsv", 3))

    # electronic annotations of dark proteins (held-aside corroboration set)
    add((acc, "electronic_function", term, "COMPUTATIONAL", "GOA-IEA", 1.0)
        for acc, term, _a in iter_tsv(dd / "iea_annotations.tsv", 3))

    # domains, annotated + dark
    for path in (d / "domains.tsv", dd / "domains.tsv"):
        add((acc, "has_domain", ipr, "COMPUTATIONAL", "InterPro", 1.0)
            for acc, ipr, _t in iter_tsv(path, 3) if ipr != "-")

    # curated domain->GO rules
    ip2go = load_interpro2go_raw(d / "interpro2go.txt")
    add((ipr, "domain_implies_function", term, "CURATED", "InterPro2GO", 1.0)
        for ipr, terms in ip2go.items() for term in terms)

    # similarity edges (test-vs-train from v0.1, dark-vs-KB from v0.3)
    add((q, "sequence_similar_to", t, "COMPUTATIONAL", "MMseqs2", float(s))
        for q, t, s in iter_tsv(d / "seq_hits.tsv", 3))
    add((q, "structure_similar_to", t, "COMPUTATIONAL", "Foldseek", float(s))
        for q, t, s in iter_tsv(d / "struct_hits.tsv", 3))
    for raw, pred in ((dd / "seq_hits_raw.tsv", "sequence_similar_to"),
                      (dd / "struct_hits_raw.tsv", "structure_similar_to")):
        # raw files have no header and raw ids; normalise like the streams do
        from pis.common import norm_hit_id
        rows = []
        with open(raw, encoding="utf-8") as f:
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 3:
                    continue
                sim = float(cols[2])
                if pred == "sequence_similar_to" and sim > 1.0:
                    sim /= 100.0
                rows.append((norm_hit_id(cols[0]), pred, norm_hit_id(cols[1]),
                             "COMPUTATIONAL", "MMseqs2" if "seq" in raw.name else "Foldseek",
                             sim))
        add(rows)

    # AI hypotheses from the dark sweep
    add((acc, "predicted_function", term, "HYPOTHESIS", "fusion-v0.2", float(score))
        for acc, _b, term, score, _tier, _streams, _so
        in iter_tsv(r / "dark_hypotheses.tsv", 7))

    # v0.4 layers, when present: STRING interactions and fpocket pockets
    string_path = d / "string_edges.tsv"
    if string_path.exists():
        add((a, "interacts_with", b, "COMPUTATIONAL", "STRING", int(s) / 1000.0)
            for a, b, s in iter_tsv(string_path, 3))
    # v1.0: experiment-loop outcomes (supported = validated by experimental
    # annotation; contradicted = experimental NOT-qualifier; unconfirmed = open world)
    outcomes_path = r / "hypothesis_outcomes.tsv"
    if outcomes_path.exists():
        ev_of = {"supported": "EXPERIMENTAL", "contradicted": "EXPERIMENTAL",
                 "unconfirmed": "COMPUTATIONAL"}
        add((acc, "hypothesis_" + outcome, term, ev_of[outcome], source, float(conf))
            for acc, term, _b, conf, _arm, outcome, source
            in iter_tsv(outcomes_path, 7) if outcome in ev_of)

    pockets_path = dd / "pockets.tsv"
    if pockets_path.exists():
        rows = []
        with open(pockets_path, encoding="utf-8") as f:
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if len(cols) >= 3 and float(cols[2]) > 0:
                    rows.append((cols[0], "has_druggable_pocket", "pocket",
                                 "COMPUTATIONAL", "fpocket", float(cols[2])))
        add(rows)

    con.execute("CREATE INDEX idx_subj ON edges(subject)")
    con.execute("CREATE INDEX idx_obj ON edges(object)")
    con.commit()

    stats = con.execute(
        "SELECT predicate, evidence, COUNT(*) FROM edges "
        "GROUP BY predicate, evidence ORDER BY 3 DESC").fetchall()
    n_edges = sum(s[2] for s in stats)
    n_nodes = con.execute(
        "SELECT COUNT(*) FROM (SELECT subject FROM edges UNION SELECT object FROM edges)"
    ).fetchone()[0]

    # demo query: best structure-only hypothesis and its 1-hop neighbourhood
    demo = con.execute(
        "SELECT subject, object, score FROM edges "
        "WHERE predicate='predicted_function' AND evidence='HYPOTHESIS' "
        "AND subject IN (SELECT subject FROM edges WHERE predicate='structure_similar_to') "
        "AND subject NOT IN (SELECT subject FROM edges WHERE predicate='sequence_similar_to') "
        "AND subject NOT IN (SELECT subject FROM edges WHERE predicate='has_domain') "
        "ORDER BY score DESC LIMIT 1").fetchone()

    with open(r / "kg_stats.md", "w", encoding="utf-8", newline="\n") as f:
        f.write("# Knowledge graph - data/kg.db\n\n")
        f.write("{} edges, {} nodes.\n\n".format(n_edges, n_nodes))
        f.write("| predicate | evidence | edges |\n|---|---|---:|\n")
        for pred, ev, n in stats:
            f.write("| {} | {} | {} |\n".format(pred, ev, n))
        if demo:
            acc, term, score = demo
            f.write("\n## Demo: top structure-only hypothesis\n\n")
            f.write("`{}` predicted_function `{}` (score {:.2f}) - supported only by "
                    "structural neighbours:\n\n".format(acc, term, score))
            for s, p, o, ev, src, sc in con.execute(
                    "SELECT * FROM edges WHERE subject=? ORDER BY predicate, score DESC LIMIT 15",
                    (acc,)):
                f.write("- {} **{}** {}  ({}, {}, {:.2f})\n".format(s, p, o, ev, src, sc))

    for pred, ev, n in stats:
        print("{:24s} {:14s} {:>9d}".format(pred, ev, n))
    print("Graph: {} edges, {} nodes -> data/kg.db (stats in results/kg_stats.md)".format(
        n_edges, n_nodes))
    con.close()


if __name__ == "__main__":
    main()
