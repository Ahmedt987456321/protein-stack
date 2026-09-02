"""Exploration 4 - fold-space map analysis.

Consumes the Foldseek structural clustering of all 11,726 models and
cross-tabulates clusters with InterPro entries and experimental GO:
cluster count and size distribution, InterPro purity of large clusters,
dark-enriched clusters, and fold-function disagreements (clusters whose
annotated members have disjoint molecular-function profiles).

Outputs: results/explore/foldspace.md, foldspace_clusters.tsv
"""
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pis.common import data_dir, load_config, norm_hit_id
from pis.go import GoDag
from pis.streams import load_domains


def main():
    cfg = load_config()
    d = data_dir(cfg)
    out_dir = Path("results/explore")
    out_dir.mkdir(parents=True, exist_ok=True)
    dag = GoDag(d / "go-basic.obo")

    clusters = defaultdict(list)
    with open(d / "explore" / "foldclust_cluster.tsv", encoding="utf-8") as f:
        for line in f:
            rep, mem = line.rstrip("\n").split("\t")
            clusters[norm_hit_id(rep)].append(norm_hit_id(mem))

    domains = load_domains(d / "domains.tsv")
    for acc, ds in load_domains(d / "dark" / "domains.tsv").items():
        domains.setdefault(acc, set()).update(ds)
    dark = set((d / "dark" / "accessions.txt").read_text(encoding="utf-8").split())

    mf_terms = defaultdict(set)
    with open(d / "annotations.tsv", encoding="utf-8") as f:
        next(f)
        for line in f:
            acc, term, aspect = line.rstrip("\n").split("\t")
            if aspect == "F":
                mf_terms[acc].add(term)
    mf_prop = {a: dag.propagate(ts) for a, ts in mf_terms.items()}

    sizes = sorted((len(v) for v in clusters.values()), reverse=True)
    singletons = sum(1 for s in sizes if s == 1)

    rows = []
    disagreements = []
    dark_enriched = []
    for rep, members in clusters.items():
        n = len(members)
        if n < 2:
            continue
        iprs = Counter(ipr for m in members for ipr in domains.get(m, ()))
        top_ipr, top_n = (iprs.most_common(1)[0] if iprs else ("-", 0))
        with_dom = sum(1 for m in members if domains.get(m))
        purity = top_n / with_dom if with_dom else 0.0
        n_dark = sum(1 for m in members if m in dark)
        annotated = [m for m in members if m in mf_prop]
        # disjoint MF profiles among annotated members?
        disjoint_pairs = 0
        for i in range(min(len(annotated), 12)):
            for j in range(i + 1, min(len(annotated), 12)):
                a, b = mf_prop[annotated[i]], mf_prop[annotated[j]]
                if a and b and not (a & b):
                    disjoint_pairs += 1
        rows.append([rep, n, n_dark, top_ipr, "{:.2f}".format(purity),
                     len(annotated), disjoint_pairs])
        if disjoint_pairs > 0 and n >= 3:
            disagreements.append((rep, n, disjoint_pairs, top_ipr))
        if n >= 5 and n_dark / n >= 0.8:
            dark_enriched.append((rep, n, n_dark, top_ipr))

    rows.sort(key=lambda r: -r[1])
    with open(out_dir / "foldspace_clusters.tsv", "w", encoding="utf-8", newline="\n") as f:
        f.write("representative\tsize\tn_dark\ttop_interpro\tinterpro_purity\t"
                "n_annotated\tdisjoint_MF_pairs\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")

    big = [r for r in rows if r[1] >= 5]
    mean_purity = sum(float(r[4]) for r in big) / len(big) if big else 0.0
    with open(out_dir / "foldspace.md", "w", encoding="utf-8", newline="\n") as f:
        f.write("# Fold-space map of 11,726 AlphaFold models\n\n")
        f.write("Structural clustering (Foldseek, coverage 0.8): {} clusters; "
                "largest {}; median size {}; {} singletons ({:.0f}%).\n\n".format(
                    len(clusters), sizes[0], sizes[len(sizes) // 2], singletons,
                    100 * singletons / len(clusters)))
        f.write("Clusters with >= 5 members: {}; their mean InterPro purity "
                "(members sharing the most common entry) is {:.2f}, i.e. "
                "structural clusters largely recover curated families.\n\n".format(
                    len(big), mean_purity))
        f.write("## Ten largest clusters\n\n")
        f.write("| representative | size | dark | top InterPro | purity |\n")
        f.write("|---|---|---|---|---|\n")
        for r in rows[:10]:
            f.write("| {} | {} | {} | {} | {} |\n".format(r[0], r[1], r[2], r[3], r[4]))
        f.write("\n## Dark-enriched fold clusters (>= 80% dark, >= 5 members)\n\n")
        f.write("{} clusters. Top 10 by size:\n\n".format(len(dark_enriched)))
        f.write("| representative | size | dark | top InterPro |\n|---|---|---|---|\n")
        for rep, n, nd, ipr in sorted(dark_enriched, key=lambda x: -x[1])[:10]:
            f.write("| {} | {} | {} | {} |\n".format(rep, n, nd, ipr))
        f.write("\n## Fold-function disagreements\n\n")
        f.write("{} clusters (>= 3 members) contain annotated member pairs "
                "with fully disjoint experimental molecular-function "
                "profiles. Top 15 by disjoint pairs:\n\n".format(len(disagreements)))
        f.write("| representative | size | disjoint MF pairs | top InterPro |\n")
        f.write("|---|---|---|---|\n")
        for rep, n, dp, ipr in sorted(disagreements, key=lambda x: -x[2])[:15]:
            f.write("| {} | {} | {} | {} |\n".format(rep, n, dp, ipr))
    print("clusters {} | singletons {} | big {} | mean purity {:.2f} | "
          "disagreement clusters {} | dark-enriched {}".format(
              len(clusters), singletons, len(big), mean_purity,
              len(disagreements), len(dark_enriched)))
    print("Wrote results/explore/foldspace.md")


if __name__ == "__main__":
    main()
