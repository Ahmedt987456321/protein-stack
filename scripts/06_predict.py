"""Step 06 - the two prediction arms.

Both arms transfer GO terms from training neighbours to test proteins with
score(term) = max similarity among neighbours annotated with the term
(annotations propagated to ancestors first).

  Arm A (baseline):   sequence neighbours only        (sim = identity)
  Arm B (treatment):  sequence + structure neighbours (sim = identity or TM)

Outputs:
  results/pred_armA.tsv, results/pred_armB.tsv   (accession, go_term, score)
"""
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pis.common import data_dir, load_config, results_dir
from pis.go import GoDag


def load_hits(path: Path):
    hits = defaultdict(list)  # query -> [(target, sim)]
    with open(path, encoding="utf-8") as f:
        next(f)  # header
        for line in f:
            q, t, sim = line.rstrip("\n").split("\t")
            hits[q].append((t, float(sim)))
    return hits


def load_split(path: Path):
    train, test = set(), set()
    with open(path, encoding="utf-8") as f:
        next(f)
        for line in f:
            acc, split, _bin = line.rstrip("\n").split(",")
            (train if split == "train" else test).add(acc)
    return train, test


def load_annotations(path: Path):
    ann = defaultdict(set)
    with open(path, encoding="utf-8") as f:
        next(f)
        for line in f:
            acc, term, _aspect = line.rstrip("\n").split("\t")
            ann[acc].add(term)
    return ann


def transfer(hit_sources, train_terms):
    """hit_sources: list of {query: [(target, sim)]}. Returns query -> {term: score}."""
    preds = defaultdict(dict)
    for hits in hit_sources:
        for q, neighbours in hits.items():
            scores = preds[q]
            for target, sim in neighbours:
                for term in train_terms.get(target, ()):
                    if sim > scores.get(term, 0.0):
                        scores[term] = sim
    return preds


def write_preds(path: Path, preds):
    n = 0
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("accession\tgo_term\tscore\n")
        for acc in sorted(preds):
            for term, score in sorted(preds[acc].items()):
                f.write("{}\t{}\t{:.4f}\n".format(acc, term, score))
                n += 1
    return n


def main():
    cfg = load_config()
    d = data_dir(cfg)
    r = results_dir(cfg)

    dag = GoDag(d / "go-basic.obo")
    train, test = load_split(d / "split.csv")
    raw_ann = load_annotations(d / "annotations.tsv")

    # propagated ground-truth terms for training proteins only
    train_terms = {acc: dag.propagate(raw_ann[acc]) for acc in train if acc in raw_ann}
    print("Training proteins with propagated annotations: {}".format(len(train_terms)))

    seq_hits = load_hits(d / "seq_hits.tsv")
    struct_hits = load_hits(d / "struct_hits.tsv")
    # keep only test-set queries and train-set targets (defence in depth)
    for hits in (seq_hits, struct_hits):
        for q in list(hits):
            if q not in test:
                del hits[q]
            else:
                hits[q] = [(t, s) for t, s in hits[q] if t in train]

    pred_a = transfer([seq_hits], train_terms)
    pred_b = transfer([seq_hits, struct_hits], train_terms)

    na = write_preds(r / "pred_armA.tsv", pred_a)
    nb = write_preds(r / "pred_armB.tsv", pred_b)
    covered_a = len(pred_a)
    covered_b = len(pred_b)
    print("Arm A (sequence only):      {} proteins covered, {} scored terms".format(covered_a, na))
    print("Arm B (+ structure):        {} proteins covered, {} scored terms".format(covered_b, nb))
    print("Structure-only rescues (proteins with hits in B but not A): {}".format(
        len(set(pred_b) - set(pred_a))))


if __name__ == "__main__":
    main()
