"""Unit tests for the pure-Python core: OBO parsing, propagation, hit
normalization, and the Fmax computation. Run: python tests/test_core.py
"""
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pis.common import norm_hit_id
from pis.go import GoDag, ROOTS

OBO = """format-version: 1.2

[Term]
id: GO:0008150
name: biological_process
namespace: biological_process

[Term]
id: GO:0003674
name: molecular_function
namespace: molecular_function

[Term]
id: GO:0016787
name: hydrolase activity
namespace: molecular_function
is_a: GO:0003674 ! molecular_function

[Term]
id: GO:0004252
name: serine-type endopeptidase activity
namespace: molecular_function
alt_id: GO:9999999
is_a: GO:0016787 ! hydrolase activity

[Term]
id: GO:0000001
name: obsolete thing
namespace: molecular_function
is_obsolete: true

[Term]
id: GO:0009987
name: cellular process
namespace: biological_process
is_a: GO:0008150 ! biological_process
relationship: part_of GO:0008150 ! biological_process

[Typedef]
id: part_of
name: part of
"""


def check(name, cond):
    status = "ok  " if cond else "FAIL"
    print("[{}] {}".format(status, name))
    return cond


def main():
    passed = True

    with tempfile.TemporaryDirectory() as td:
        obo = Path(td) / "test.obo"
        obo.write_text(OBO, encoding="utf-8")
        dag = GoDag(obo)

        anc = dag.ancestors("GO:0004252")
        passed &= check("ancestors climb is_a chain",
                        anc == frozenset({"GO:0004252", "GO:0016787", "GO:0003674"}))
        passed &= check("propagate drops roots",
                        dag.propagate(["GO:0004252"]) == {"GO:0004252", "GO:0016787"})
        passed &= check("alt_id resolves", dag.canonical("GO:9999999") == "GO:0004252")
        passed &= check("obsolete rejected", dag.canonical("GO:0000001") is None)
        passed &= check("part_of followed",
                        "GO:0008150" in dag.ancestors("GO:0009987"))
        passed &= check("branch lookup", dag.branch("GO:0016787") == "molecular_function")

    passed &= check("norm plain accession", norm_hit_id("P12345") == "P12345")
    passed &= check("norm .pdb suffix", norm_hit_id("P12345.pdb") == "P12345")
    passed &= check("norm foldseek chain suffix", norm_hit_id("P12345.pdb_A") == "P12345")
    passed &= check("norm path prefix", norm_hit_id("data/struct_test/A0A024R161.pdb") == "A0A024R161")

    # ---- Fmax on a hand-computable case -------------------------------
    import pis.eval as ev

    truth = {"p1": {"t1", "t2"}, "p2": {"t1"}}
    # p1: t1 scored 0.9 (right), t3 scored 0.9 (wrong); p2: no predictions.
    preds = {"p1": {"t1": 0.9, "t3": 0.9}}
    f, t = ev.fmax(["p1", "p2"], preds, truth)
    # at any t <= 0.9: precision = 0.5 (over p1 only), recall = (0.5 + 0)/2 = 0.25
    expected = 2 * 0.5 * 0.25 / 0.75
    passed &= check("fmax hand-computed value", abs(f - expected) < 1e-9)

    # perfect predictor sanity
    preds2 = {"p1": {"t1": 1.0, "t2": 1.0}, "p2": {"t1": 1.0}}
    f2, _ = ev.fmax(["p1", "p2"], preds2, truth)
    passed &= check("fmax perfect predictor = 1", abs(f2 - 1.0) < 1e-9)

    # empty predictions -> 0
    f3, _ = ev.fmax(["p1", "p2"], {}, truth)
    passed &= check("fmax no predictions = 0", f3 == 0.0)

    # ---- interface enrichment (Fisher) --------------------------------
    import math
    from pis.stats import interface_enrichment, kabsch_rmsf

    # regression against the published flagship numbers (80 complexes scored:
    # pathogenic 63/274 at interface vs benign 18/138)
    odds, p = interface_enrichment(63, 274, 18, 138)
    passed &= check("enrichment odds ratio ~ 1.99", abs(odds - 1.99) < 0.05)
    passed &= check("enrichment Fisher p ~ 0.0103", abs(p - 0.0103) < 5e-4)
    # undefined when a class is empty
    o2, p2 = interface_enrichment(5, 10, 0, 0)
    passed &= check("enrichment nan on empty class", math.isnan(o2) and math.isnan(p2))
    # stronger separation -> smaller p (monotone sanity)
    _, p_weak = interface_enrichment(30, 100, 10, 100)
    _, p_strong = interface_enrichment(50, 100, 10, 100)
    passed &= check("enrichment p shrinks with separation", p_strong < p_weak)
    # input validation: hits cannot exceed totals
    try:
        interface_enrichment(11, 10, 1, 5)
        raised = False
    except ValueError:
        raised = True
    passed &= check("enrichment rejects hit > total", raised)

    # ---- Kabsch RMSF --------------------------------------------------
    import numpy as np
    A = [(1, 0.0, 0.0, 0.0), (2, 1.0, 0.0, 0.0), (3, 1.0, 1.0, 0.0),
         (4, 0.0, 1.0, 0.5), (5, 2.0, 0.5, 1.0)]
    # a pure rotation + translation of A must yield ~zero RMSF: this is the
    # property that was broken before (centroid-only alignment reported huge
    # flexibility for models that only differed in orientation).
    Rz = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    shift = np.array([5.0, 3.0, -2.0])
    B = [(r, *(Rz @ np.array([x, y, z]) + shift)) for r, x, y, z in A]
    rmsf_rigid = kabsch_rmsf([A, B])
    passed &= check("kabsch: rigid motion gives ~0 RMSF",
                    max(rmsf_rigid.values()) < 1e-6)
    # a genuine single-residue displacement must localize to that residue
    C = [(r, x, y, z + (2.0 if r == 3 else 0.0)) for r, x, y, z in A]
    rmsf_local = kabsch_rmsf([A, C])
    passed &= check("kabsch: displaced residue has max RMSF",
                    rmsf_local[3] == max(rmsf_local.values()) and rmsf_local[3] > 0.5)
    # guards
    passed &= check("kabsch: <2 traces -> empty", kabsch_rmsf([A]) == {})
    passed &= check("kabsch: <3 shared residues -> empty",
                    kabsch_rmsf([A[:2], A[:2]]) == {})

    # ---- provenance: paper number reproduces from committed raw data ---
    # Re-derive the flagship enrichment straight from the per-complex tsv and
    # assert it still matches the published statistic. This fails loudly if the
    # report and the raw data ever drift apart (e.g. a hand-edited number, or a
    # truncated result file). Skips cleanly if the data file is absent.
    tsv = ROOT / "results" / "explore" / "variant_interface.tsv"
    if tsv.exists():
        ph = pt = bh = bt = 0
        for i, line in enumerate(tsv.read_text(encoding="utf-8").splitlines()):
            if i == 0:
                continue
            c = line.split("\t")
            if len(c) < 6 or c[4] == "-":
                continue
            a, tot = c[4].split("/"); ph += int(a); pt += int(tot)
            a, tot = c[5].split("/"); bh += int(a); bt += int(tot)
        odds_d, p_d = interface_enrichment(ph, pt, bh, bt)
        passed &= check("provenance: tsv reproduces published OR ~ 1.99",
                        abs(odds_d - 1.99) < 0.05)
        passed &= check("provenance: tsv reproduces published p ~ 0.0103",
                        abs(p_d - 0.0103) < 5e-4)
    else:
        print("[skip] provenance check (variant_interface.tsv not present)")

    print()
    if passed:
        print("ALL TESTS PASSED")
        return 0
    print("TESTS FAILED")
    return 1


if __name__ == "__main__":
    sys.exit(main())
