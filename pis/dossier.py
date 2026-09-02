"""Evidence dossier: everything the knowledge graph knows about one protein,
grouped by evidence tier, rendered as markdown for a human or an LLM.

This is the deterministic half of the v0.5 agent: engines and the graph
gather, the LLM reasons over the result. Nothing here predicts anything.
"""
import sqlite3
from collections import defaultdict

TIER_ORDER = ["EXPERIMENTAL", "CURATED", "COMPUTATIONAL", "HYPOTHESIS"]


def _rows(con, sql, args=()):
    return con.execute(sql, args).fetchall()


def build_dossier(con: sqlite3.Connection, acc: str, dag, max_per_section: int = 8):
    """Return a dict dossier for one accession, or None if the graph has
    nothing at all about it."""
    edges = _rows(con, "SELECT predicate, object, evidence, source, score "
                       "FROM edges WHERE subject=? ORDER BY score DESC", (acc,))
    if not edges:
        return None

    by_pred = defaultdict(list)
    for pred, obj, ev, src, score in edges:
        by_pred[pred].append((obj, ev, src, score))

    # similarity and interaction edges are stored once per pair; pick up the
    # ones where this protein is the object so KB-side proteins get neighbours
    reverse = _rows(con, "SELECT predicate, subject, evidence, source, score "
                         "FROM edges WHERE object=? AND predicate IN "
                         "('interacts_with','sequence_similar_to','structure_similar_to') "
                         "ORDER BY score DESC", (acc,))
    for pred, subj, ev, src, score in reverse:
        by_pred[pred].append((subj, ev, src, score))
    for pred in ("interacts_with", "sequence_similar_to", "structure_similar_to"):
        if pred in by_pred:
            seen = set()
            merged = []
            for item in sorted(by_pred[pred], key=lambda x: -x[3]):
                if item[0] not in seen:
                    seen.add(item[0])
                    merged.append(item)
            by_pred[pred] = merged

    def go_list(pred):
        return [{"term": t, "name": dag.name(t), "branch": dag.branch(t),
                 "evidence": ev, "source": src, "score": sc}
                for t, ev, src, sc in by_pred.get(pred, [])][:max_per_section * 3]

    def neighbour_list(pred):
        out = []
        for target, _ev, _src, score in by_pred.get(pred, [])[:max_per_section]:
            funcs = _rows(con,
                          "SELECT object FROM edges WHERE subject=? AND "
                          "predicate='has_function' ORDER BY score DESC LIMIT 4",
                          (target,))
            out.append({"partner": target, "score": score,
                        "partner_functions": [
                            {"term": f[0], "name": dag.name(f[0])} for f in funcs]})
        return out

    pocket = by_pred.get("has_druggable_pocket")
    dossier = {
        "accession": acc,
        "is_dark": not by_pred.get("has_function"),
        "experimental_functions": go_list("has_function"),
        "electronic_annotations": go_list("electronic_function"),
        "domains": [{"interpro": t, "evidence": ev, "source": src}
                    for t, ev, src, _ in by_pred.get("has_domain", [])],
        "sequence_neighbours": neighbour_list("sequence_similar_to"),
        "structure_neighbours": neighbour_list("structure_similar_to"),
        "interaction_partners": neighbour_list("interacts_with"),
        "hypotheses": go_list("predicted_function"),
        "druggability": pocket[0][3] if pocket else None,
        "loop_outcomes": [
            {"outcome": pred.replace("hypothesis_", ""), "term": t,
             "name": dag.name(t), "evidence": ev, "source": src, "score": sc}
            for pred in ("hypothesis_contradicted", "hypothesis_supported",
                         "hypothesis_unconfirmed")
            for t, ev, src, sc in by_pred.get(pred, [])
        ],
    }

    # curated implications of this protein's domains
    implied = []
    for dom in dossier["domains"][:max_per_section]:
        for t, _ev, _src, _sc in by_pred_lookup(con, dom["interpro"],
                                                "domain_implies_function")[:4]:
            implied.append({"interpro": dom["interpro"], "term": t,
                            "name": dag.name(t)})
    dossier["curated_domain_implications"] = implied
    return dossier


def by_pred_lookup(con, subject, predicate):
    return _rows(con, "SELECT object, evidence, source, score FROM edges "
                      "WHERE subject=? AND predicate=?", (subject, predicate))


def render_markdown(d, max_per_section: int = 8) -> str:
    """Render a dossier dict as markdown."""
    lines = ["# Evidence dossier: {}".format(d["accession"]), ""]
    lines.append("Status: **{}**".format(
        "DARK PROTEIN -- no experimentally verified function"
        if d["is_dark"] else "annotated -- has experimental GO evidence"))
    lines.append("")

    def go_section(title, items, tier_note):
        if not items:
            return
        lines.append("## {} [{}]".format(title, tier_note))
        for it in items[:max_per_section]:
            lines.append("- {} -- {} ({}, score {:.2f})".format(
                it["term"], it["name"], (it["branch"] or "?")[:4], it["score"]))
        if len(items) > max_per_section:
            lines.append("- ... and {} more".format(len(items) - max_per_section))
        lines.append("")

    go_section("Experimental GO annotations", d["experimental_functions"],
               "EXPERIMENTAL -- trusted ground truth")

    if d["loop_outcomes"]:
        lines.append("## Experiment-loop outcomes [validation record]")
        for it in d["loop_outcomes"][:max_per_section]:
            lines.append("- {}: {} -- {} (hypothesis conf {:.2f}; {}; {})".format(
                it["outcome"].upper(), it["term"], it["name"], it["score"],
                it["evidence"], it["source"]))
        n_more = len(d["loop_outcomes"]) - max_per_section
        if n_more > 0:
            lines.append("- ... and {} more".format(n_more))
        lines.append("")
    go_section("Electronic annotations (never used for prediction)",
               d["electronic_annotations"], "COMPUTATIONAL -- unverified, IEA-class")
    go_section("Function hypotheses (fusion v0.2)", d["hypotheses"],
               "HYPOTHESIS -- requires validation")

    if d["domains"]:
        lines.append("## InterPro entries [COMPUTATIONAL]")
        lines.append("- " + ", ".join(x["interpro"] for x in d["domains"][:12]))
        lines.append("")
    if d["curated_domain_implications"]:
        lines.append("## Curated domain->function implications [CURATED, interpro2go]")
        for it in d["curated_domain_implications"][:max_per_section]:
            lines.append("- {} implies {} -- {}".format(
                it["interpro"], it["term"], it["name"]))
        lines.append("")

    def neigh_section(title, items, unit):
        if not items:
            return
        lines.append("## {} [COMPUTATIONAL]".format(title))
        for n in items[:max_per_section]:
            funcs = "; ".join("{} ({})".format(f["name"], f["term"])
                              for f in n["partner_functions"]) or "no experimental functions on record"
            lines.append("- {} ({} {:.2f}): {}".format(
                n["partner"], unit, n["score"], funcs))
        lines.append("")

    neigh_section("Structural neighbours (Foldseek over AlphaFold)",
                  d["structure_neighbours"], "TM")
    neigh_section("Sequence neighbours (MMseqs2)", d["sequence_neighbours"], "id")
    neigh_section("Interaction partners (STRING)", d["interaction_partners"], "conf")

    if d["druggability"] is not None:
        lines.append("## Druggability [COMPUTATIONAL, fpocket on predicted structure]")
        lines.append("- best pocket druggability score: {:.2f} (>= 0.5 is "
                     "conventionally druggable)".format(d["druggability"]))
        lines.append("")
    return "\n".join(lines)
