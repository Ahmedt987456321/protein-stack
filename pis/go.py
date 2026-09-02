"""Gene Ontology: OBO parsing, GAF parsing, ancestor propagation.

Self-contained (no obonet/networkx). Propagation follows is_a and part_of,
which is the CAFA convention.
"""
import gzip
from pathlib import Path

ROOTS = {"GO:0008150", "GO:0003674", "GO:0005575"}

NAMESPACE_OF_ASPECT = {
    "F": "molecular_function",
    "P": "biological_process",
    "C": "cellular_component",
}


class GoDag:
    def __init__(self, obo_path: Path):
        self.parents = {}     # term -> set(parent terms) via is_a / part_of
        self.namespace = {}   # term -> namespace string
        self.names = {}       # term -> human-readable name
        self.alt = {}         # alt_id -> canonical id
        self.obsolete = set()
        self._anc_cache = {}
        self._parse(obo_path)

    def _parse(self, obo_path: Path) -> None:
        term_id = None
        in_term = False
        with open(obo_path, encoding="utf-8") as f:
            for raw in f:
                line = raw.strip()
                if line == "[Term]":
                    in_term, term_id = True, None
                    continue
                if line.startswith("["):  # [Typedef] etc.
                    in_term = False
                    continue
                if not in_term or not line:
                    continue
                if line.startswith("id: "):
                    term_id = line[4:]
                    self.parents.setdefault(term_id, set())
                elif term_id is None:
                    continue
                elif line.startswith("name: "):
                    self.names[term_id] = line[6:]
                elif line.startswith("alt_id: "):
                    self.alt[line[8:]] = term_id
                elif line.startswith("namespace: "):
                    self.namespace[term_id] = line[11:]
                elif line.startswith("is_a: "):
                    self.parents[term_id].add(line[6:].split(" ")[0])
                elif line.startswith("relationship: part_of "):
                    self.parents[term_id].add(line[22:].split(" ")[0])
                elif line == "is_obsolete: true":
                    self.obsolete.add(term_id)

    def canonical(self, term: str):
        term = self.alt.get(term, term)
        if term not in self.parents or term in self.obsolete:
            return None
        return term

    def ancestors(self, term: str):
        """All ancestors of term including itself, excluding nothing."""
        cached = self._anc_cache.get(term)
        if cached is not None:
            return cached
        out = {term}
        stack = list(self.parents.get(term, ()))
        while stack:
            t = stack.pop()
            if t in out:
                continue
            out.add(t)
            stack.extend(self.parents.get(t, ()))
        out = frozenset(out)
        self._anc_cache[term] = out
        return out

    def propagate(self, terms):
        """Union of ancestor closures, minus the three root terms."""
        out = set()
        for t in terms:
            c = self.canonical(t)
            if c is not None:
                out |= self.ancestors(c)
        return out - ROOTS

    def branch(self, term: str):
        return self.namespace.get(term)

    def name(self, term: str):
        return self.names.get(term, term)


def _open_gaf(path: Path):
    """Open a GAF whether gzip-compressed or plain, by magic bytes."""
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def parse_gaf_full(gaf_gz_path: Path):
    """Yield (accession, go_term, aspect, evidence_code) for UniProtKB rows."""
    with _open_gaf(gaf_gz_path) as f:
        for line in f:
            if line.startswith("!"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 15:
                continue
            if cols[0] != "UniProtKB":
                continue
            if "NOT" in cols[3].split("|"):
                continue
            if cols[8] not in ("F", "P", "C"):
                continue
            acc = cols[1].split("-")[0]  # canonical accession, no isoform
            yield acc, cols[4], cols[8], cols[6]


def parse_gaf(gaf_gz_path: Path, evidence_codes):
    """Yield (accession, go_term, aspect) for rows matching evidence_codes."""
    codes = set(evidence_codes)
    for acc, term, aspect, evidence in parse_gaf_full(gaf_gz_path):
        if evidence in codes:
            yield acc, term, aspect
