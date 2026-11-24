"""
CS1005 Discrete Structures - Fall 2026 - Course Project
Finding Network Motifs in an Email Network (SNAP email-Eu-core)

BENCHMARK / REFERENCE IMPLEMENTATION (instructor use)

Usage:
    python motifs.py --digits 1234 4567 0890        # last-4-digit blocks of each roll number
    python motifs.py --R 6691                       # or give R directly
    python motifs.py --R 6691 --networks 5 --attempts-factor 1   # quick smoke test

Data files (downloaded automatically if missing):
    email-Eu-core.txt.gz
    email-Eu-core-department-labels.txt.gz

Only the standard library (+ matplotlib for the chart) is used for the core tasks.
networkx is imported ONLY inside verify_with_networkx(), as a cross-check of T4.
"""

import argparse
import gzip
import itertools
import os
import random
import statistics
import shutil
import ssl
import sys
import tempfile
import urllib.error
import urllib.request
from collections import defaultdict

EDGES_URL = "https://snap.stanford.edu/data/email-Eu-core.txt.gz"
DEPTS_URL = "https://snap.stanford.edu/data/email-Eu-core-department-labels.txt.gz"
# Data files and the default chart live next to this script, so the program behaves the
# same on Windows, macOS and Linux no matter which folder it is launched from.
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
EDGES_FILE = os.path.join(BASE_DIR, "email-Eu-core.txt.gz")
DEPTS_FILE = os.path.join(BASE_DIR, "email-Eu-core-department-labels.txt.gz")

# ---------------------------------------------------------------------------
# Triad classes (Section 2 of the handout)
# ---------------------------------------------------------------------------
# Bit k of the code is 1 when the k-th edge below is present.
#   bit:   0     1     2     3     4     5
#   edge: 0->1  0->2  1->0  1->2  2->0  2->1
EDGE_BIT = {(0, 1): 0, (0, 2): 1, (1, 0): 2, (1, 2): 3, (2, 0): 4, (2, 1): 5}

CONNECTED_CODES = [3, 6, 7, 10, 11, 15, 21, 23, 25, 27, 30, 31, 63]
CODE_TO_CLASS = {code: k + 1 for k, code in enumerate(CONNECTED_CODES)}

PERMS = list(itertools.permutations(range(3)))  # the 6 orderings


def code_of_ordering(present, order):
    """6-bit code of the triad when the nodes are listed in `order`.

    present : set of directed pairs (i, j) over positions {0,1,2}
              describing which edges exist among the three nodes,
              indexed by the *original* position.
    order   : a permutation; position p of the new ordering holds original node order[p].
    """
    code = 0
    for (p, q), bit in EDGE_BIT.items():
        if (order[p], order[q]) in present:
            code |= 1 << bit
    return code


def canonical_code(present):
    """Minimum 6-bit code over all 6 orderings (the orbit representative)."""
    return min(code_of_ordering(present, order) for order in PERMS)


def triad_class(a, b, c, out):
    """Return the class 1..13 of the triad on nodes a, b, c, or None if disconnected.

    `out` maps node -> set of out-neighbours (the adjacency of G).
    Steps: 64 labelled configurations -> canonical code (16 orbits) -> 13 connected classes.
    """
    nodes = (a, b, c)
    present = set()
    for i in range(3):
        for j in range(3):
            if i != j and nodes[j] in out.get(nodes[i], ()):
                present.add((i, j))
    code = canonical_code(present)
    return CODE_TO_CLASS.get(code)  # None -> disconnected orbit


# ---------------------------------------------------------------------------
# T1: load / clean / parameters
# ---------------------------------------------------------------------------
def _ssl_context():
    """TLS context that also works on python.org macOS builds (no system certs)."""
    try:
        import certifi  # optional; installed alongside many packages
        return ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        return ssl.create_default_context()


def download_if_missing(url, path):
    """Download `url` to `path` unless it already exists. Exits with help on failure."""
    if os.path.exists(path):
        return
    print(f"Downloading {url} ...")
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(path) or ".", suffix=".part")
    try:
        with os.fdopen(fd, "wb") as out, \
                urllib.request.urlopen(url, context=_ssl_context(), timeout=60) as resp:
            shutil.copyfileobj(resp, out)
        os.replace(tmp, path)  # atomic on every platform; no half-written file is left behind
    except (urllib.error.URLError, OSError, ssl.SSLError) as e:
        if os.path.exists(tmp):
            os.remove(tmp)
        sys.exit(f"Could not download {url}\n  ({e})\n"
                 f"Download it manually in a browser and save it as:\n  {path}")


def load_graph(edges_path=EDGES_FILE, depts_path=DEPTS_FILE):
    """Read SNAP files; drop self-loops and duplicate edges."""
    raw = 0
    edge_set = set()
    with gzip.open(edges_path, "rt", encoding="utf-8") as f:
        for line in f:
            parts = line.split()
            if len(parts) != 2:
                continue
            u, v = int(parts[0]), int(parts[1])
            raw += 1
            if u != v:
                edge_set.add((u, v))
    dept_of = {}
    with gzip.open(depts_path, "rt", encoding="utf-8") as f:
        for line in f:
            parts = line.split()
            if len(parts) == 2:
                dept_of[int(parts[0])] = int(parts[1])
    return edge_set, dept_of, raw


def top20_departments(edge_set, dept_of):
    """Rank departments by internal edges (largest first, ties -> smaller dept id)."""
    internal = defaultdict(int)
    for u, v in edge_set:
        if u in dept_of and v in dept_of and dept_of[u] == dept_of[v]:
            internal[dept_of[u]] += 1
    # include departments with zero internal edges too (they can never be in top 20 here)
    for d in set(dept_of.values()):
        internal.setdefault(d, 0)
    ranked = sorted(internal.items(), key=lambda kv: (-kv[1], kv[0]))
    return [d for d, _ in ranked[:20]], dict(internal)


def derive_parameters(R):
    pairs = [(i, j) for i in range(20) for j in range(i + 1, 20)]  # lexicographic, 190 pairs
    assert len(pairs) == 190
    pair_no = R % 190
    return {
        "R": R,
        "pair_number": pair_no,
        "pair": pairs[pair_no],
        "num_networks": 20 + (R % 31),
        "attempts_factor": 5 + (R % 6),
        "seed": R,
    }


# ---------------------------------------------------------------------------
# T2: induced subgraph
# ---------------------------------------------------------------------------
def induced_subgraph(edge_set, dept_of, dept):
    members = {n for n, d in dept_of.items() if d == dept}
    edges = {(u, v) for (u, v) in edge_set if u in members and v in members}
    return members, edges


# ---------------------------------------------------------------------------
# Adjacency helpers
# ---------------------------------------------------------------------------
def build_adjacency(edges):
    out = defaultdict(set)
    inn = defaultdict(set)
    for u, v in edges:
        out[u].add(v)
        inn[v].add(u)
    return out, inn


# ---------------------------------------------------------------------------
# T4: triad census (no C(|V|,3) enumeration)
# ---------------------------------------------------------------------------
def triad_census(edges):
    """Return [N1..N13]. Each connected 3-subset is counted exactly once.

    Every connected triple has a node adjacent to the other two. For each node v we
    enumerate pairs of its neighbours (in or out) and collect the sorted triple in a set,
    which removes duplicates produced by the other 'centres' of the same triple.
    """
    out, inn = build_adjacency(edges)
    nodes = set(out) | set(inn)
    triples = set()
    for v in nodes:
        nbrs = sorted(out.get(v, set()) | inn.get(v, set()))
        for a, b in itertools.combinations(nbrs, 2):
            triples.add(tuple(sorted((v, a, b))))
    counts = [0] * 13
    for a, b, c in triples:
        k = triad_class(a, b, c, out)
        assert k is not None, "a triple built from neighbour pairs must be connected"
        counts[k - 1] += 1
    return counts


# ---------------------------------------------------------------------------
# T5: double-edge-swap randomisation
# ---------------------------------------------------------------------------
def degree_sequences(edges):
    outd, ind = defaultdict(int), defaultdict(int)
    for u, v in edges:
        outd[u] += 1
        ind[v] += 1
    return dict(outd), dict(ind)


def randomize(edges, attempts, seed):
    """Return a randomised copy of `edges` after `attempts` swap attempts.

    One attempt: pick edges (a,b),(c,d) uniformly from the current edge list and
    replace them with (a,d),(c,b). Reject (no change) if a==d, c==b, (a,d) in E or
    (c,b) in E. Rejected attempts still count. `seed` may be an int or a random.Random.
    The two edges are drawn as two distinct list positions.
    """
    rng = seed if isinstance(seed, random.Random) else random.Random(seed)
    edge_list = list(edges)
    edge_set = set(edge_list)
    m = len(edge_list)
    orig_out, orig_in = degree_sequences(edge_list)

    if m >= 2:
        for _ in range(attempts):
            i = rng.randrange(m)
            j = rng.randrange(m - 1)
            if j >= i:
                j += 1  # uniform over positions != i
            a, b = edge_list[i]
            c, d = edge_list[j]
            if a == d or c == b:
                continue
            if (a, d) in edge_set or (c, b) in edge_set:
                continue
            edge_set.discard((a, b))
            edge_set.discard((c, d))
            edge_set.add((a, d))
            edge_set.add((c, b))
            edge_list[i] = (a, d)
            edge_list[j] = (c, b)

    new_out, new_in = degree_sequences(edge_list)
    assert len(edge_set) == len(edge_list) == len(edges), "graph must stay simple, same size"
    assert new_out == orig_out, "out-degrees changed!"
    assert new_in == orig_in, "in-degrees changed!"
    assert all(u != v for u, v in edge_list), "self-loop introduced!"
    return edge_list


# ---------------------------------------------------------------------------
# T6: Z-scores
# ---------------------------------------------------------------------------
def z_scores(real_counts, random_counts):
    """random_counts: list of 13-vectors (one per random graph). Population std."""
    z, mu, sd = [], [], []
    for k in range(13):
        sample = [rc[k] for rc in random_counts]
        m = statistics.fmean(sample)
        s = statistics.pstdev(sample)
        mu.append(m)
        sd.append(s)
        z.append(0.0 if s == 0 else (real_counts[k] - m) / s)
    return z, mu, sd


def run_department(edge_set, dept_of, dept, params, rng, label):
    members, edges = induced_subgraph(edge_set, dept_of, dept)
    print(f"\n[{label}] department {dept}: {len(members)} nodes, {len(edges)} edges")
    real = triad_census(edges)
    attempts = params["attempts_factor"] * len(edges)
    rand_counts = []
    for r in range(params["num_networks"]):
        g = randomize(edges, attempts, rng)  # starts from the ORIGINAL edges each time
        rand_counts.append(triad_census(g))
        if (r + 1) % 10 == 0 or r + 1 == params["num_networks"]:
            print(f"   random networks done: {r + 1}/{params['num_networks']}")
    z, mu, sd = z_scores(real, rand_counts)
    return {"dept": dept, "nodes": len(members), "edges": len(edges),
            "real": real, "mu": mu, "sd": sd, "z": z}


# ---------------------------------------------------------------------------
# T7: chart
# ---------------------------------------------------------------------------
def plot_z(resA, resB, R, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = list(range(1, 14))
    w = 0.38
    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.bar([i - w / 2 for i in x], resA["z"], w, label=f"Dept A = {resA['dept']}")
    ax.bar([i + w / 2 for i in x], resB["z"], w, label=f"Dept B = {resB['dept']}")
    ax.axhline(2, linestyle="--", color="black", linewidth=1)
    ax.axhline(-2, linestyle="--", color="black", linewidth=1)
    ax.axhline(0, color="grey", linewidth=0.6)
    ax.set_xticks(x)
    ax.set_xlabel("Triad class")
    ax.set_ylabel("Z-score")
    ax.set_title(f"Triad-class Z-scores, Dept {resA['dept']} vs Dept {resB['dept']} (R = {R})")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    print(f"Saved chart -> {path}")


# ---------------------------------------------------------------------------
# Cross-check against networkx (allowed ONLY to verify T4; class order differs)
# ---------------------------------------------------------------------------
def verify_with_networkx(edges, counts):
    try:
        import networkx as nx
    except ImportError:
        print("networkx not installed; skipping cross-check.")
        return None
    # Map each of our 13 classes to networkx's triad name via a representative graph.
    inv_bit = {bit: e for e, bit in EDGE_BIT.items()}
    class_to_name = {}
    for k, code in enumerate(CONNECTED_CODES, start=1):
        g = nx.DiGraph()
        g.add_nodes_from(range(3))
        for bit in range(6):
            if code >> bit & 1:
                g.add_edge(*inv_bit[bit])
        class_to_name[k] = nx.triad_type(g)
    G = nx.DiGraph()
    G.add_edges_from(edges)
    census = nx.triadic_census(G)
    ok = all(census[class_to_name[k]] == counts[k - 1] for k in range(1, 14))
    print("nx.triadic_census cross-check:", "PASS" if ok else "FAIL")
    return ok


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    # Windows consoles may use a legacy code page; never crash on printing.
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--R", type=int, help="group number R")
    g.add_argument("--digits", nargs="+", help="last-4-digit part of each member's roll number")
    ap.add_argument("--networks", type=int, help="override #random networks (smoke test only)")
    ap.add_argument("--attempts-factor", type=int, help="override swap factor (smoke test only)")
    ap.add_argument("--out", default=os.path.join(BASE_DIR, "zscores.png"),
                    help="chart file (default: zscores.png next to this script)")
    ap.add_argument("--no-verify", action="store_true", help="skip networkx cross-check")
    args = ap.parse_args()

    R = args.R if args.R is not None else sum(int(d) for d in args.digits)

    # --- T1 ---
    download_if_missing(EDGES_URL, EDGES_FILE)
    download_if_missing(DEPTS_URL, DEPTS_FILE)
    edge_set, dept_of, raw = load_graph()
    nodes = {n for e in edge_set for n in e}
    print(f"Raw edge lines: {raw}; after removing self-loops/duplicates: {len(edge_set)} edges, "
          f"{len(nodes)} non-isolated nodes, {len(set(dept_of.values()))} departments")

    top20, internal = top20_departments(edge_set, dept_of)
    params = derive_parameters(R)
    if args.networks:
        params["num_networks"] = args.networks
    if args.attempts_factor:
        params["attempts_factor"] = args.attempts_factor
    i, j = params["pair"]
    deptA, deptB = top20[i], top20[j]

    print("\n=== PARAMETERS ===")
    print(f"R                      = {R}")
    print(f"Pair number (R mod 190)= {params['pair_number']} -> (i, j) = ({i}, {j})")
    print(f"Department A           = top20[{i}] = {deptA}  ({internal[deptA]} internal edges)")
    print(f"Department B           = top20[{j}] = {deptB}  ({internal[deptB]} internal edges)")
    print(f"Random networks        = {params['num_networks']}")
    print(f"Swap attempts / network= {params['attempts_factor']} x (#edges)")
    print(f"Seed                   = {params['seed']}")
    print("Top-20 (dept: internal edges):", [(d, internal[d]) for d in top20])

    # --- T3 hand example (printed so it can be shown in the video) ---
    demo_out = {0: {1, 2}, 1: set(), 2: set()}
    print("\nT3 example: edges 0->1, 0->2 -> class", triad_class(0, 1, 2, demo_out),
          "(code 3 -> class 1)")

    # --- T2, T4-T6 ---
    rng = random.Random(params["seed"])  # one generator seeded with R for the whole run
    resA = run_department(edge_set, dept_of, deptA, params, rng, "A")
    resB = run_department(edge_set, dept_of, deptB, params, rng, "B")

    if not args.no_verify:
        for res in (resA, resB):
            _, edges = induced_subgraph(edge_set, dept_of, res["dept"])
            verify_with_networkx(edges, res["real"])

    print("\n=== T4: real counts ===")
    print(f"{'Class':>5} {'Dept '+str(deptA):>10} {'Dept '+str(deptB):>10}")
    for k in range(13):
        print(f"{k+1:>5} {resA['real'][k]:>10} {resB['real'][k]:>10}")

    print("\n=== T6: Z-scores (motif: Z>2, anti-motif: Z<-2) ===")
    print(f"{'Class':>5} {'Z(A)':>9} {'Z(B)':>9}  A-flag  B-flag")
    flag = lambda z: "motif" if z > 2 else ("anti" if z < -2 else "-")
    for k in range(13):
        print(f"{k+1:>5} {resA['z'][k]:>9.2f} {resB['z'][k]:>9.2f}  {flag(resA['z'][k]):>6}  {flag(resB['z'][k]):>6}")

    plot_z(resA, resB, R, args.out)


if __name__ == "__main__":
    main()
