"""Causal layer (build2, D) - genetic x physical interaction crossing.

DepMap CRISPR co-dependency is behind a verification wall and not scriptable
here (logged). BioGRID genetic interactions are the scriptable substitute
and arguably more direct: a genetic interaction (synthetic lethality,
negative/positive genetic, dosage rescue) is a measured functional
dependency between two genes' perturbations. Crossing BioGRID GENETIC
interactions with our PHYSICAL interaction evidence (STRING high-confidence
edges + AF-template-supported pairs) yields pairs that both physically
associate AND functionally co-depend - a stronger, more causal statement
than either alone, and a functional-genomics evidence stream for the dark
proteins.

Outputs: results/explore/causal_biogrid.md, causal_pairs.tsv
"""
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pis.common import data_dir, load_config

# BioGRID tab3 1-indexed columns
C_SYS = 11        # Experimental System (specific assay)
C_TYPE = 12       # Experimental System Type: physical | genetic
C_SPA = 23        # SWISS-PROT Accessions Interactor A
C_SPB = 26        # SWISS-PROT Accessions Interactor B


def accs(cell):
    out = set()
    for tok in cell.replace(";", "|").split("|"):
        tok = tok.strip().split("-")[0]
        if tok and tok != "-":
            out.add(tok)
    return out


def main():
    cfg = load_config()
    d = data_dir(cfg)
    out_dir = Path("results/explore")
    out_dir.mkdir(parents=True, exist_ok=True)

    universe = set((d / "final_accessions.txt").read_text().split())
    universe |= set((d / "dark" / "accessions.txt").read_text().split())
    dark = set((d / "dark" / "accessions.txt").read_text().split())

    # physical interactions we already have: STRING (>=700) + AF-template
    physical = set()
    with open(d / "string_edges.tsv") as f:
        next(f)
        for line in f:
            a, b, s = line.rstrip("\n").split("\t")
            physical.add(frozenset((a, b)))
    templ = out_dir / "templated_interactions.tsv"
    if templ.exists():
        with open(templ) as f:
            next(f)
            for line in f:
                a, b, pdb, tm, isd = line.rstrip("\n").split("\t")
                physical.add(frozenset((a, b)))
    print("physical interaction pairs (STRING + template):", len(physical))

    # stream BioGRID; collect genetic interactions among our proteins
    bg = sorted(d.glob("explore/BIOGRID-ALL-*.tab3.txt"))
    if not bg:
        print("BioGRID file not found"); return
    genetic = defaultdict(set)  # frozenset(pair) -> set of assay types
    n_gen = 0
    with open(bg[-1], encoding="utf-8", errors="replace") as f:
        next(f)  # header
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) <= C_SPB or c[C_TYPE] != "genetic":
                continue
            aset = accs(c[C_SPA]) & universe
            bset = accs(c[C_SPB]) & universe
            if not aset or not bset:
                continue
            for a in aset:
                for b in bset:
                    if a != b:
                        genetic[frozenset((a, b))].add(c[C_SYS])
                        n_gen += 1
    print("BioGRID genetic-interaction rows mapped to our proteins:", n_gen)
    print("distinct genetic pairs among our proteins:", len(genetic))

    # cross: genetic AND physical
    both = [pair for pair in genetic if pair in physical]
    dark_both = [p for p in both if any(x in dark for x in p)]
    print("pairs both genetic AND physical:", len(both))
    print("of those involving a dark protein:", len(dark_both))

    with open(out_dir / "causal_pairs.tsv", "w", encoding="utf-8", newline="\n") as f:
        f.write("protein_a\tprotein_b\tinvolves_dark\tgenetic_assays\n")
        for pair in sorted(both, key=lambda p: (any(x in dark for x in p), sorted(p)),
                           reverse=True):
            a, b = sorted(pair)
            f.write("{}\t{}\t{}\t{}\n".format(
                a, b, int(any(x in dark for x in pair)),
                ";".join(sorted(genetic[pair]))))

    # enrichment: are physical pairs more often genetic than random?
    # (physical pairs that are genetic) vs (genetic rate among all our pairs)
    phys_genetic = len(both)
    phys_total = len(physical)
    with open(out_dir / "causal_biogrid.md", "w", encoding="utf-8", newline="\n") as f:
        f.write("# Genetic x physical interaction crossing (causal layer)\n\n")
        f.write("DepMap CRISPR co-dependency is not scriptable from this "
                "environment (verification wall); BioGRID genetic interactions "
                "are used instead - a measured functional dependency between "
                "two genes' perturbations.\n\n")
        f.write("- Physical interaction pairs (STRING >= 700 + AF-template): "
                "{}\n".format(phys_total))
        f.write("- Distinct genetic-interaction pairs among our proteins "
                "(BioGRID): {}\n".format(len(genetic)))
        f.write("- Pairs that are BOTH physical and genetic: **{}** "
                "({:.1%} of physical pairs also genetically interact)\n".format(
                    phys_genetic, phys_genetic / phys_total if phys_total else 0))
        f.write("- Of those, involving a dark (unannotated) protein: **{}** - "
                "these carry both physical and functional-genetic evidence "
                "for a protein with no experimental annotation.\n\n".format(
                    len(dark_both)))
        f.write("A pair with both evidence types is a stronger, more causal "
                "functional link than either alone: they physically associate "
                "AND perturbing one modifies the other's phenotype. Full list: "
                "causal_pairs.tsv.\n\n")
        f.write("## Dark-protein pairs with both physical and genetic evidence "
                "(top 20)\n\n")
        f.write("| dark-involving pair | genetic assay(s) |\n|---|---|\n")
        for pair in dark_both[:20]:
            a, b = sorted(pair)
            f.write("| {} - {} | {} |\n".format(
                a, b, ";".join(sorted(genetic[pair]))[:60]))
    print("Wrote results/explore/causal_biogrid.md")


if __name__ == "__main__":
    main()
