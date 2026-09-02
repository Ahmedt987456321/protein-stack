"""Evidence streams and fusion - shared by the v0.2 experiment (script 09)
and the v0.3 dark-proteome sweep (script 12).

A "stream" is {accession: {go_term: score}} with scores in [0, 1].
"""
from collections import Counter, defaultdict

BRANCHES = ["molecular_function", "biological_process", "cellular_component"]
BIN_LABELS = ["lt30", "30to50", "50to80", "ge80"]
MIN_SUPPORT = 3     # proteins a domain needs before P(term|domain) is trusted
MIN_PROB = 0.02     # drop weaker domain->term associations (below Fmax range)

from pis.common import norm_hit_id  # noqa: E402


def bin_of(identity, edges):
    """identity as fraction 0..1; edges e.g. [0.3, 0.5, 0.8]."""
    for i, edge in enumerate(edges):
        if identity < edge:
            return BIN_LABELS[i]
    return BIN_LABELS[len(edges)]


def parse_raw_seq(path, allowed_targets, transfer_evalue):
    """Raw mmseqs output (query,target,pident,evalue) at permissive e-value.

    Returns (s_seq: {query: max identity}, transfer: {query: [(target, ident)]}).
    """
    s_seq = defaultdict(float)
    transfer = defaultdict(list)
    with open(path, encoding="utf-8") as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4:
                continue
            q, t = norm_hit_id(cols[0]), norm_hit_id(cols[1])
            if t not in allowed_targets or q == t:
                continue
            ident = float(cols[2])
            if ident > 1.0:
                ident /= 100.0
            s_seq[q] = max(s_seq[q], ident)
            if float(cols[3]) <= transfer_evalue:
                transfer[q].append((t, ident))
    return s_seq, transfer


def parse_raw_struct(path, allowed_targets):
    """Raw foldseek output (query,target,alntmscore,...) -> {query: [(target, tm)]}."""
    hits = defaultdict(list)
    with open(path, encoding="utf-8") as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            q, t = norm_hit_id(cols[0]), norm_hit_id(cols[1])
            if t in allowed_targets and q != t:
                hits[q].append((t, float(cols[2])))
    return hits


def load_domains(path):
    domains = defaultdict(set)
    with open(path, encoding="utf-8") as f:
        next(f)
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 2 and cols[1] != "-":
                domains[cols[0]].add(cols[1])
    return domains


def load_interpro2go(path, dag):
    """interpro2go lines: 'InterPro:IPR000971 Globin > GO:... ; GO:0015671'.

    Returns {interpro_id: propagated GO term set}.
    """
    mapping = defaultdict(set)
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith("InterPro:IPR"):
                continue
            _head, _, tail = line.partition(";")
            ipr = line.split()[0].split(":")[1]
            go = tail.strip()
            if go.startswith("GO:"):
                mapping[ipr].add(go)
    return {ipr: dag.propagate(terms) for ipr, terms in mapping.items()}


def load_interpro2go_raw(path):
    """Same file, but the unpropagated {interpro_id: {go_term}} pairs."""
    mapping = defaultdict(set)
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith("InterPro:IPR"):
                continue
            _head, _, tail = line.partition(";")
            ipr = line.split()[0].split(":")[1]
            go = tail.strip()
            if go.startswith("GO:"):
                mapping[ipr].add(go)
    return mapping


def transfer_stream(hits, source_terms):
    """{query: {term: max sim among annotated neighbours}}"""
    out = {}
    for q, neighbours in hits.items():
        scores = {}
        for target, sim in neighbours:
            for term in source_terms.get(target, ()):
                if sim > scores.get(term, 0.0):
                    scores[term] = sim
        out[q] = scores
    return out


def domain_stats(accs, prop_terms, domains):
    """P(term | domain) over the given annotated proteins."""
    nd = Counter()
    dt = defaultdict(Counter)
    for acc in accs:
        terms = prop_terms.get(acc)
        doms = domains.get(acc)
        if not terms or not doms:
            continue
        for dm in doms:
            nd[dm] += 1
            dt[dm].update(terms)
    probs = {}
    for dm, n in nd.items():
        if n < MIN_SUPPORT:
            continue
        p = {t: c / n for t, c in dt[dm].items() if c / n >= MIN_PROB}
        if p:
            probs[dm] = p
    return probs


def domain_stream(accs, domains, probs, ip2go):
    out = {}
    for acc in accs:
        scores = {}
        for dm in domains.get(acc, ()):
            for t, p in probs.get(dm, {}).items():
                if p > scores.get(t, 0.0):
                    scores[t] = p
            for t in ip2go.get(dm, ()):
                scores[t] = 1.0
        if scores:
            out[acc] = scores
    return out


def split_by_branch(stream, dag):
    """{acc: {term: s}} -> {branch: {acc: {term: s}}}"""
    out = {b: {} for b in BRANCHES}
    for acc, scores in stream.items():
        per = {b: {} for b in BRANCHES}
        for t, s in scores.items():
            b = dag.branch(t)
            if b in per:
                per[b][t] = s
        for b in BRANCHES:
            if per[b]:
                out[b][acc] = per[b]
    return out


def specific_terms(scores, dag, min_score, top_n):
    """Keep the most specific high-scoring terms: drop a term if a descendant
    of it is predicted at (nearly) the same score. Returns [(term, score)]."""
    kept = []
    items = [(t, s) for t, s in scores.items() if s >= min_score]
    items.sort(key=lambda x: (-x[1], x[0]))  # deterministic tie-break
    for t, s in items:
        dominated = False
        for t2, s2 in items:
            if t2 != t and s2 >= s - 0.01 and t in dag.ancestors(t2):
                dominated = True
                break
        if not dominated:
            kept.append((t, s))
        if len(kept) >= top_n:
            break
    return kept


def fuse(accs, seq_b, str_b, dom_b, s_seq, tau, alpha, beta, gamma, tau_d, delta):
    """Noisy-OR fusion for one branch. Returns {acc: {term: score}}."""
    out = {}
    for acc in accs:
        strength = s_seq.get(acc, 0.0)
        w_str = alpha if strength < tau else beta
        w_dom = gamma if strength < tau_d else delta
        seq_s = seq_b.get(acc, {})
        str_s = str_b.get(acc, {})
        dom_s = dom_b.get(acc, {})
        terms = set(seq_s)
        if w_str > 0:
            terms |= set(str_s)
        if w_dom > 0:
            terms |= set(dom_s)
        scores = {}
        for t in sorted(terms):  # deterministic iteration order
            keep = (1.0 - seq_s.get(t, 0.0))
            keep *= (1.0 - w_str * str_s.get(t, 0.0))
            keep *= (1.0 - w_dom * dom_s.get(t, 0.0))
            s = 1.0 - keep
            if s > 0.0:
                scores[t] = s
        out[acc] = scores
    return out
