"""Exploration B2 - PAE-based domain-level structural transfer.

Test proteins are segmented into structural domains by clustering the
predicted aligned error matrix (10-residue windows; two windows join when
their mean inter-window PAE is below 10 A; union-find components of at
least 40 residues become domains). Multi-domain proteins have each domain
searched separately against the training structures; per-domain transfers
are unioned. Single-domain proteins keep their whole-chain transfer, so
the comparison isolates the multi-domain treatment.

Arm C-dom (gated fusion with the domain-level structure stream) is
evaluated against arm C on the standard test cells, overall and restricted
to the treated (multi-domain) subgroup.

Outputs: results/explore/pae_domains.md
"""
import json
import shutil
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tqdm import tqdm

from pis.common import data_dir, load_config, norm_hit_id, results_dir, run_tool
from pis.eval import fmax
from pis.go import GoDag
from pis.streams import (
    BIN_LABELS,
    BRANCHES,
    domain_stats,
    domain_stream,
    fuse,
    load_domains,
    load_interpro2go,
    parse_raw_seq,
    parse_raw_struct,
    split_by_branch,
    transfer_stream,
)

WIN = 10
PAE_CUT = 10.0
MIN_DOM = 40


def load_pae(path):
    js = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(js, list):
        js = js[0]
    return js.get("predicted_aligned_error") or js.get("pae")


class UF:
    def __init__(self, n):
        self.p = list(range(n))
    def find(self, x):
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x
    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[rb] = ra


def segment(pae):
    n = len(pae)
    nw = (n + WIN - 1) // WIN
    uf = UF(nw)
    for i in range(nw):
        ri = range(i * WIN, min(n, (i + 1) * WIN))
        for j in range(i + 1, nw):
            rj = range(j * WIN, min(n, (j + 1) * WIN))
            tot = cnt = 0
            for a in ri:
                row = pae[a]
                for b in rj:
                    tot += row[b]
                    cnt += 1
            for a in rj:
                row = pae[a]
                for b in ri:
                    tot += row[b]
                    cnt += 1
            if cnt and tot / cnt < PAE_CUT:
                uf.union(i, j)
    comps = defaultdict(set)
    for w in range(nw):
        comps[uf.find(w)].update(range(w * WIN, min(n, (w + 1) * WIN)))
    doms = [sorted(c) for c in comps.values() if len(c) >= MIN_DOM]
    return doms


def write_domain_pdb(src, residues, dest):
    keep = set(r + 1 for r in residues)  # PDB is 1-based
    with open(src, encoding="utf-8", errors="replace") as f, \
         open(dest, "w", encoding="utf-8", newline="\n") as g:
        for line in f:
            if line.startswith(("ATOM", "TER", "END")):
                if line.startswith("ATOM"):
                    try:
                        if int(line[22:26]) not in keep:
                            continue
                    except ValueError:
                        continue
                g.write(line)


def main():
    cfg = load_config()
    d = data_dir(cfg)
    r = results_dir(cfg)
    out_dir = Path("results/explore")
    out_dir.mkdir(parents=True, exist_ok=True)
    dag = GoDag(d / "go-basic.obo")
    params = json.load(open(r / "fusion_params.json", encoding="utf-8"))

    train, test, bins = [], [], {}
    with open(d / "split.csv", encoding="utf-8") as f:
        next(f)
        for line in f:
            acc, split, b = line.rstrip("\n").split(",")
            (train if split == "train" else test).append(acc)
            if split == "test":
                bins[acc] = b

    # ---- segmentation -----------------------------------------------------
    dom_dir = d / "explore" / "domains_test"
    seg_cache = d / "explore" / "pae_segments.json"
    if seg_cache.exists():
        segments = {k: v for k, v in json.loads(seg_cache.read_text()).items()}
    else:
        dom_dir.mkdir(parents=True, exist_ok=True)
        segments = {}
        for acc in tqdm(test, unit="protein"):
            pae_path = d / "explore" / "pae" / (acc + ".json")
            if not pae_path.exists():
                continue
            try:
                doms = segment(load_pae(pae_path))
            except Exception:
                continue
            segments[acc] = [[dm[0], dm[-1], len(dm)] for dm in doms]
            if len(doms) >= 2:
                src = d / "structures" / (acc + ".pdb")
                if src.exists():
                    for i, dm in enumerate(doms):
                        write_domain_pdb(src, dm, dom_dir / "{}_d{}.pdb".format(acc, i))
        seg_cache.write_text(json.dumps(segments))
    multi = sorted(a for a, doms in segments.items() if len(doms) >= 2)
    print("test proteins segmented: {} | multi-domain: {}".format(
        len(segments), len(multi)))

    # ---- domain-level search ----------------------------------------------
    hits_path = d / "explore" / "domain_struct_hits.tsv"
    if not hits_path.exists():
        run_tool(cfg, "foldseek", [
            "easy-search", dom_dir, d / "struct_train",
            hits_path, d / "explore" / "tmp_domsearch",
            "--format-output", "query,target,alntmscore,bits,evalue",
            "-e", str(cfg["search"]["evalue"]),
            "--max-seqs", str(cfg["search"]["max_seqs"]),
            "--threads", str(cfg["tools"]["threads"]),
        ])
        shutil.rmtree(d / "explore" / "tmp_domsearch", ignore_errors=True)

    train_set = set(train)
    dom_hits = defaultdict(list)  # acc -> [(target, tm)]
    with open(hits_path, encoding="utf-8") as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 3:
                continue
            q = norm_hit_id(c[0])          # strips _dN as well (underscore split)
            t = norm_hit_id(c[1])
            if t in train_set and q != t:
                dom_hits[q].append((t, float(c[2])))

    # ---- streams: reuse manuscript arms, swap structure stream ------------
    raw_ann = defaultdict(set)
    with open(d / "annotations.tsv", encoding="utf-8") as f:
        next(f)
        for line in f:
            acc, term, _a = line.rstrip("\n").split("\t")
            raw_ann[acc].add(term)
    prop = {a: dag.propagate(ts) for a, ts in raw_ann.items()}
    train_terms = {a: prop[a] for a in train if a in prop}

    def load_clean(path):
        hits = defaultdict(list)
        with open(path, encoding="utf-8") as f:
            next(f)
            for line in f:
                q, t, sim = line.rstrip("\n").split("\t")
                if t in train_set and q != t:
                    hits[q].append((t, float(sim)))
        return hits
    seq_hits = load_clean(d / "seq_hits.tsv")
    chain_struct = load_clean(d / "struct_hits.tsv")
    s_seq, _ = parse_raw_seq(d / "bin_hits.tsv", train_set, cfg["search"]["evalue"])

    struct_dom = dict(chain_struct)
    for acc in multi:
        if acc in dom_hits:
            struct_dom[acc] = dom_hits[acc]

    domains = load_domains(d / "domains.tsv")
    ip2go = load_interpro2go(d / "interpro2go.txt", dag)
    train_probs = domain_stats(train, prop, domains)

    seq_stream = split_by_branch(transfer_stream(seq_hits, train_terms), dag)
    dstr = split_by_branch(transfer_stream(struct_dom, train_terms), dag)
    cstr = split_by_branch(transfer_stream(chain_struct, train_terms), dag)
    dm_stream = split_by_branch(domain_stream(test, domains, train_probs, ip2go), dag)

    arms = {}
    for name, stream in (("C", cstr), ("Cdom", dstr)):
        arms[name] = {b: fuse(test, seq_stream[b], stream[b], dm_stream[b],
                              s_seq, params["tau"], params["alpha"],
                              params["beta"], params["gamma"],
                              params["tau_d"], params["delta"])
                      for b in BRANCHES}

    truth = {b: {} for b in BRANCHES}
    for acc in test:
        for b in BRANCHES:
            sub = {t for t in prop.get(acc, ()) if dag.branch(t) == b}
            if sub:
                truth[b][acc] = sub

    lines = ["# PAE domain-level structural transfer", "",
             "{} of {} test proteins are multi-domain by PAE segmentation "
             "(windows {} residues, cut {} A, min domain {} residues); their "
             "structure stream is replaced by unioned per-domain transfer. "
             "Single-domain proteins are untouched, so all differences below "
             "come from the treated subgroup.".format(
                 len(multi), len(segments), WIN, PAE_CUT, MIN_DOM), ""]
    for scope, prots_of in (("all test proteins", lambda b: sorted(truth[b])),
                            ("multi-domain subgroup",
                             lambda b: sorted(a for a in truth[b] if a in set(multi)))):
        lines += ["## " + scope, "",
                  "| branch | n | C | C-dom | delta |", "|---|---|---|---|---|"]
        for b in BRANCHES:
            prots = prots_of(b)
            fc, _ = fmax(prots, arms["C"][b], truth[b])
            fd, _ = fmax(prots, arms["Cdom"][b], truth[b])
            lines.append("| {} | {} | {:.4f} | {:.4f} | {:+.4f} |".format(
                b, len(prots), fc, fd, fd - fc))
            print(scope, b, len(prots), "C={:.4f} Cdom={:.4f}".format(fc, fd))
        lines.append("")
    (out_dir / "pae_domains.md").write_text("\n".join(lines) + "\n",
                                            encoding="utf-8")
    print("Wrote results/explore/pae_domains.md")


if __name__ == "__main__":
    main()
