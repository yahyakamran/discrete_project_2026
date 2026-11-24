"""
CS1005 Discrete Structures - Fall 2026 - Course Project
Finding Network Motifs in an Email Network (SNAP email-Eu-core)

Step 1 (T1, T2): load and clean the data, derive parameters, build subgraphs.

Usage:
    python motifs.py --digits 1234 4567 0890
    python motifs.py --R 6691
"""

import argparse
import gzip
import os
import shutil
import ssl
import sys
import tempfile
import urllib.error
import urllib.request
from collections import defaultdict

EDGES_URL = "https://snap.stanford.edu/data/email-Eu-core.txt.gz"
DEPTS_URL = "https://snap.stanford.edu/data/email-Eu-core-department-labels.txt.gz"

# Data files live next to this script, so the program behaves the same on Windows,
# macOS and Linux no matter which folder it is launched from.
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
EDGES_FILE = os.path.join(BASE_DIR, "email-Eu-core.txt.gz")
DEPTS_FILE = os.path.join(BASE_DIR, "email-Eu-core-department-labels.txt.gz")


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

    # --- T2 ---
    for label, dept in (("A", deptA), ("B", deptB)):
        members, edges = induced_subgraph(edge_set, dept_of, dept)
        print(f"[{label}] department {dept}: {len(members)} nodes, {len(edges)} edges")


if __name__ == "__main__":
    main()
