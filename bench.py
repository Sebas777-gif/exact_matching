#!/usr/bin/env python3
"""Regression tests and benchmarks for play, which builds the A0 arrangement.

  ./bench.py small [--fast]    exact A0 counts for small pairs; exit 1 on a mismatch
  ./bench.py medium [--aida]   larger pairs with known distances, Hera as reference
  ./bench.py finite            the 18 random pairs with a finite matching distance
  ./bench.py random            every pair of tests/random_bifiltrations/*.bif
  ./bench.py kerber            feasible pairs of data_benchmarks/small_instances

finite, random, and kerber append one row per pair to a results file and
resume where they stopped. ./bench.py SET --help lists the options of a set.
Ctrl+N skips the running pair. Ctrl+C stops and keeps all completed rows.
"""

import argparse
import atexit
import os
import re
import shutil
import signal
import subprocess
import sys
import termios
import threading
import time
import tty
from collections import defaultdict
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BIF_DIR = ROOT / "tests" / "random_bifiltrations"
KERBER_DIR = ROOT / "data_benchmarks" / "small_instances"


# --- pairs ------------------------------------------------------------------

def v10(sample):
    return f"tests/random_bifiltrations/r_bif_v_10_dim_2_n_top_simplices_50_sample_{sample}.bif"


# small: name, expected "lines vertices edges faces", play arguments. Always
# run with --no-aida. The comments give the matching distance d each pair
# was built for; it is not checked, since play does not evaluate F yet.
SMALL = [
    # Rank-1 intervals (d = 1). Grades (0,0), (1,1), (2,2): all three
    # duals and s = 1 meet at (s, t) = (1, 0).
    ("interval scc", "5 4 15 12", ["tests/interval_long.scc", "tests/interval_short.scc"]),
    ("interval bif", "5 4 15 12", ["tests/interval_long.bif", "tests/interval_short.bif"]),
    # Two origin vertices merging at (51,0) vs a single infinite bar
    # (d = 25.5). The .scc version adds unused grades at x, y in
    # {10, 20, 30} that only cut A0; the .bif version does not.
    ("far_merge vs A0 anchors", "10 24 62 39", ["tests/far_merge.scc", "tests/a0_anchors.scc"]),
    ("far_merge bif", "5 4 15 12", ["tests/far_merge.bif", "tests/a0_anchors.bif"]),
    # Horizontal merge at (50,0) vs vertical merge at (0,50) (d = 25).
    ("horizontal vs vertical merge", "5 6 18 13",
     ["tests/merge_horizontal.scc", "tests/merge_vertical.scc"]),
    # Different numbers of infinite bars (d = infinity).
    ("unmatched infinite bars", "4 3 11 9", ["tests/free_one.scc", "tests/free_two.scc"]),
    # Raw H_2 presentations of two prism 2-skeleta (d = infinity).
    ("prism H_2", "7 12 34 23", ["tests/prism_1.bif", "tests/prism_2.bif", "--dim", "2"]),
    # L of two intervals vs one spanning interval (d = 1.5). The maximum of
    # F sits at (s, t) = (1, 1/2), on an A0 edge, not an A0 vertex.
    ("ell vs span", "7 9 29 21", ["tests/ell.scc", "tests/span.scc"]),
]

# small, skipped with --fast.
SMALL_SLOW = [
    # far_merge with 140 unused relation grades (d = 25.5).
    ("far_merge vs 140 A0 anchors", "144 5183 10581 5399",
     ["tests/far_merge.scc", "tests/a0_anchors_medium.scc"]),
    # Merge pair with ~140 distinct grades (d = 25).
    ("horizontal vs vertical merge medium", "145 5493 11222 5730",
     ["tests/merge_horizontal_medium.scc", "tests/merge_vertical_medium.scc"]),
    # Real 10-vertex 1-skeleta, H_0 (d = 3 and 2.5).
    ("v10 bif sample 1 vs 2", "74 1966 4221 2256", [v10(1), v10(2), "--dim", "0"]),
    ("v10 bif sample 1 vs 5", "71 1806 3883 2078", [v10(1), v10(5), "--dim", "0"]),
]

# medium: name, known matching distance, the two inputs. A bare name is a
# synthetic pair: tests/<name>.scc (empty dummy relation columns that only
# add grades) without AIDA, and tests/<name>.bif with AIDA, because AIDA
# removes the dummy columns but keeps the extra short H_0 bars of the .bif
# 1-skeleta. Hera always reads the .bif version.
MEDIUM = [
    ("far_merge vs 140 A0 anchors", "25.5", "far_merge", "a0_anchors_medium"),
    ("far_merge_41 vs 140 A0 anchors", "20.5", "far_merge_41", "a0_anchors_medium"),
    ("far_merge vs 10x10 A0 grid", "25.5", "far_merge", "a0_anchors_grid"),
    ("horizontal vs vertical merge L=50", "25", "merge_horizontal_medium", "merge_vertical_medium"),
    ("horizontal vs vertical merge L=40", "20", "merge_h40_medium", "merge_v40_medium"),
    ("horizontal vs vertical merge L=60", "30", "merge_h60_medium", "merge_v60_medium"),
] + [
    (f"v10 n50 sample {i} vs {j}", d, v10(i), v10(j))
    for i, j, d in [(1, 2, "3"), (1, 5, "2.5"), (2, 8, "3"), (3, 7, "4.5"),
                    (4, 9, "6"), (7, 10, "3"), (8, 10, "3")]
]


def rbif(v, dim, n, sample):
    return f"r_bif_v_{v}_dim_{dim}_n_top_simplices_{n}_sample_{sample}.bif"


# finite: the 18 pairs of tests/results.txt (dim 0) whose distance was
# finite; all other pairs among the first 1577 were infinity. Third entry:
# Hera's distance, for reference. Ordered by the old running time.
FINITE = [
    (rbif(100, 2, 100, 3), rbif(100, 5, 100, 9), "34.4853"),
    (rbif(100, 2, 100, 3), rbif(100, 5, 100, 5), "39.2855"),
    (rbif(100, 2, 100, 10), rbif(100, 5, 50, 9), "44.4961"),
    (rbif(100, 2, 100, 3), rbif(50, 3, 50, 3), "34.4945"),
    (rbif(100, 2, 100, 3), rbif(50, 2, 50, 10), "20.9998"),
    (rbif(100, 2, 100, 3), rbif(100, 3, 100, 2), "44.4919"),
    (rbif(100, 2, 100, 1), rbif(50, 2, 50, 5), "34.987"),
    (rbif(100, 2, 100, 3), rbif(50, 2, 50, 4), "28.9896"),
    (rbif(100, 2, 100, 2), rbif(50, 2, 50, 5), "25.4931"),
    (rbif(100, 2, 100, 3), rbif(50, 2, 50, 6), "25.9974"),
    (rbif(100, 2, 100, 3), rbif(100, 3, 100, 1), "34.4977"),
    (rbif(100, 2, 100, 3), rbif(100, 3, 100, 8), "30.4884"),
    (rbif(100, 2, 100, 1), rbif(100, 3, 100, 3), "41.4983"),
    (rbif(100, 2, 100, 2), rbif(100, 3, 100, 3), "37.9846"),
    (rbif(100, 2, 100, 3), rbif(100, 3, 100, 9), "23.4991"),
    (rbif(100, 2, 100, 3), rbif(100, 3, 100, 4), "30.9924"),
    (rbif(100, 2, 100, 3), rbif(100, 2, 150, 10), "24.9894"),
    (rbif(100, 2, 100, 1), rbif(100, 2, 100, 2), "45.9965"),
]


# --- running commands ---------------------------------------------------------

class Run:
    """Outcome of one command. status is ok, fail, timeout, or skip (Ctrl+N)."""

    def __init__(self, status, code, stdout, stderr, seconds):
        self.status, self.code, self.stdout, self.stderr = status, code, stdout, stderr
        self.seconds = f"{seconds:.3f}"


current = None                     # the running process, for Ctrl+N
skip_requested = threading.Event()


def kill(proc):
    """Kill proc and everything it started (play starts mpfree and AIDA)."""
    try:
        os.killpg(proc.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass


def run(cmd, timeout=None):
    global current
    skip_requested.clear()
    start = time.monotonic()
    # A new session is a new process group: Ctrl+C reaches only this script,
    # which then kills the group, and kill() also gets mpfree and AIDA.
    proc = subprocess.Popen([str(c) for c in cmd], stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True, start_new_session=True)
    current = proc
    try:
        out, err = proc.communicate(timeout=timeout)
        status = "ok" if proc.returncode == 0 else "fail"
    except subprocess.TimeoutExpired:
        kill(proc)
        out, err = proc.communicate()
        status = "timeout"
    except BaseException:          # Ctrl+C
        kill(proc)
        proc.wait()
        raise
    finally:
        current = None
    if skip_requested.is_set():
        status = "skip"
    return Run(status, proc.returncode, out, err, time.monotonic() - start)


def listen_for_ctrl_n():
    """Read single keys from the terminal; Ctrl+N kills the running command."""
    try:
        fd = os.open("/dev/tty", os.O_RDONLY)
        if os.tcgetpgrp(fd) != os.getpgrp():   # in the background: leave the terminal alone
            return
        old = termios.tcgetattr(fd)
    except OSError:
        return
    tty.setcbreak(fd)              # keys arrive one by one; Ctrl+C still works
    atexit.register(termios.tcsetattr, fd, termios.TCSADRAIN, old)

    def loop():
        while key := os.read(fd, 1):   # b"" when the terminal is gone
            if key == b"\x0e":
                skip_requested.set()
                if current:
                    kill(current)

    threading.Thread(target=loop, daemon=True).start()


def run_play(args, files, dim=None):
    cmd = [args.play, *files]
    if dim is not None:
        cmd += ["--dim", dim]
    if not args.aida:
        cmd.append("--no-aida")
    return run(cmd, args.timeout)


def a0_stats(r):
    """(lines, vertices, edges, faces) from play's output, or None."""
    lines = re.search(r"^A0 lines: (\d+)", r.stdout, re.M)
    arr = re.search(r"^Arrangement has (\d+) vertices, (\d+) edges, and (\d+) faces",
                    r.stdout, re.M)
    if r.status != "ok" or not (lines and arr):
        return None
    return (lines[1], *arr.groups())


def a0_fields(r, stats):
    """The a0_lines, vertices, edges, faces columns of a results row."""
    if stats:
        return list(stats)
    if r.status == "skip":
        return ["???"] * 4
    return ["timeout" if r.status == "timeout" else "error", "", "", ""]


def outcome(r, stats=None, lines=8):
    """One-line summary of a play run; on failure, the last lines of its output."""
    if stats:
        return "{} lines, {} vertices, {} edges, {} faces".format(*stats) + f"  ({r.seconds}s)"
    if r.status == "ok":
        return f"FAIL no A0 counts in the output  ({r.seconds}s)"
    if r.status != "fail":
        return f"{r.status.upper()}  ({r.seconds}s)"
    out = "\n".join((r.stderr + r.stdout).rstrip().splitlines()[-lines:])
    return f"FAIL exit {r.code}  ({r.seconds}s)\n{out}"


def find_hera(args):
    if not args.use_hera:
        return None
    if args.hera:
        return args.hera
    for c in [Path.home() / "hera_aida_matching/extern/hera/build/matching/matching_dist",
              ROOT.parent / "hera_aida_matching/extern/hera/build/matching/matching_dist",
              ROOT.parent / "hera/build/matching/matching_dist"]:
        if os.access(c, os.X_OK):
            return str(c)
    found = shutil.which("matching_dist")
    if not found:
        print("warning: Hera's matching_dist not found; reference skipped. Pass --hera PATH.",
              file=sys.stderr)
    return found


def run_hera(args, a, b, dim):
    """Hera's distance and running time, as strings."""
    r = run([args.hera, "-d", dim, a, b], args.hera_timeout)
    if r.status != "ok":
        return {"timeout": "timeout", "skip": "???"}.get(r.status, "error"), r.seconds
    last = (r.stdout.strip().splitlines() or [""])[-1].strip()
    if last.lower().lstrip("+") in ("inf", "infinity"):
        last = "infinity"
    return last or "error", r.seconds


# --- results files --------------------------------------------------------------

class Results:
    """Tab-separated file: '# key: value' header lines, a column line, one row
    per pair. Rows are appended right away, so a crash loses nothing."""

    def __init__(self, path, title, settings, info, columns, fresh):
        self.path = path
        if fresh or not path.exists():
            header = [f"# {title}", f"# started: {datetime.now().isoformat(timespec='seconds')}"]
            header += [f"# {k}: {v}" for k, v in {**settings, **info}.items()]
            path.write_text("\n".join(header + ["\t".join(columns)]) + "\n")
            self.rows = []
            return
        lines = path.read_text().splitlines()
        found = dict(line[2:].split(": ", 1) for line in lines
                     if line.startswith("# ") and ": " in line)
        for k, v in settings.items():
            if found.get(k) != str(v):
                sys.exit(f"error: {path} has {k}: {found.get(k)}, not {v}. "
                         "Use the same options, or --fresh to start over.")
        self.rows = [line.split("\t") for line in lines if line and not line.startswith("#")][1:]

    def append(self, row):
        with self.path.open("a") as f:
            f.write("\t".join(str(x) for x in row) + "\n")


def benchmark(args, pairs, results, dim=None):
    """Run play (and Hera, if args.hera) on each pair (a, b, note).

    A pair counts as done when its row has A0 values and, with Hera
    columns, a Hera value; done pairs are skipped when resuming. ??? marks a
    pair skipped with Ctrl+N; it is not retried."""
    with_hera = dim is not None
    done = {(r[0], r[1]) for r in results.rows
            if len(r) > 2 and r[2] and (not with_hera or (len(r) > 7 and r[7]))}
    print(f"play: {args.play}  ({'AIDA' if args.aida else '--no-aida'})")
    if with_hera:
        print(f"hera: {args.hera or 'skipped'}")
    print(f"results: {results.path}")
    if done:
        print(f"resuming: {len(done)} pair(s) already done")
    for i, (a, b, note) in enumerate(pairs, 1):
        if (a.name, b.name) in done:
            continue
        print(f"[{i}/{len(pairs)}] {a.name}  vs  {b.name}{note}")
        hera = []
        if with_hera:
            hera = list(run_hera(args, a, b, dim)) if args.hera else ["skip", ""]
            if args.hera:
                print(f"  hera: {hera[0]}  ({hera[1]}s)")
        try:
            r = run_play(args, [a, b], dim)
        except KeyboardInterrupt:
            if args.hera and with_hera:   # keep Hera's result; play is redone on resume
                results.append([a.name, b.name, "", "", "", "", "", *hera])
            raise
        stats = a0_stats(r)
        print("  A0: " + outcome(r, stats))
        results.append([a.name, b.name, *a0_fields(r, stats), r.seconds, *hera])
    print(f"Done. Results: {results.path}")
    return 0


# --- the sets ---------------------------------------------------------------------

def cmd_small(args):
    cases = SMALL + ([] if args.fast else SMALL_SLOW)
    print(f"play: {args.play}  (--no-aida, timeout {args.timeout:g}s)")
    print("columns: A0 lines, vertices, edges, faces\n")
    failed = 0
    for name, expected, play_args in cases:
        print(f"{name:<42}", end="", flush=True)
        r = run_play(args, play_args)
        stats = a0_stats(r)
        got = " ".join(stats) if stats else ""
        if got == expected:
            print(f"ok  {got}  ({r.seconds}s)")
            continue
        failed += 1
        if stats:
            print(f"FAIL got '{got}' expected '{expected}'  ({r.seconds}s)")
        else:
            print(outcome(r, lines=20))
    print(f"\nAll {len(cases)} tests passed." if not failed
          else f"\n{len(cases) - failed} passed, {failed} failed.")
    return 1 if failed else 0


def cmd_medium(args):
    def path(entry, ext):
        return entry if "/" in entry else f"tests/{entry}.{ext}"

    print(f"play: {args.play}  ({'AIDA' if args.aida else '--no-aida'}, timeout {args.timeout:g}s)")
    if args.hera:
        print(f"hera: {args.hera}  (reference only, timeout {args.hera_timeout:g}s)\n")
    else:
        print("hera: skipped\n")
    table = ["name\tknown_dist\ta0_lines\tvertices\tedges\tfaces\tseconds\tstatus\thera\thera_seconds"]
    failed = 0
    try:
        for name, known, a, b in MEDIUM:
            print(f"{name:<40}", end="", flush=True)
            hera, hera_sec = "skip", ""
            hera_a, hera_b = path(a, "bif"), path(b, "bif")
            if args.hera and os.path.exists(hera_a) and os.path.exists(hera_b):
                hera, hera_sec = run_hera(args, hera_a, hera_b, "0")
            ext = "bif" if args.aida else "scc"
            r = run_play(args, [path(a, ext), path(b, ext)], "0")
            stats = a0_stats(r)
            ref = f"known {known}" + (f", hera {hera} ({hera_sec}s)" if hera != "skip" else "")
            if stats:
                print(f"ok  {stats[0]} lines, {stats[1]} vertices  ({r.seconds}s)  {ref}")
                status = "ok"
            else:
                failed += 1
                status = r.status
                print(f"{ref}\n{outcome(r)}")
            table.append("\t".join([name, known, *(stats or [""] * 4), r.seconds, status,
                                    hera, hera_sec]))
    finally:
        print("\n" + "\n".join(table))
    print(f"\nAll {len(MEDIUM)} medium cases built A0." if not failed
          else f"\n{len(MEDIUM) - failed} passed, {failed} failed.")
    return 1 if failed else 0


HERA_COLUMNS = ["file_a", "file_b", "a0_lines", "vertices", "edges", "faces", "seconds",
                "hera_distance", "hera_seconds"]


def cmd_finite(args):
    pairs = [(BIF_DIR / a, BIF_DIR / b, f"  (Hera: {h})") for a, b, h in FINITE]
    missing = [p for a, b, _ in pairs for p in (a, b) if not p.exists()]
    if missing:
        sys.exit("error: missing " + ", ".join(map(str, missing)))
    mode = "aida" if args.aida else "noaida"
    results = Results(args.results or ROOT / f"tests/finite_results_a0_{mode}.txt",
                      "A0 benchmark: the 18 finite pairs of tests/results.txt",
                      {"dim": args.dim, "preprocess": mode},
                      {"play": args.play, "hera": args.hera or "skipped"},
                      HERA_COLUMNS, args.fresh)
    return benchmark(args, pairs, results, dim=args.dim)


def cmd_random(args):
    files = sorted(BIF_DIR.glob("*.bif"))
    if len(files) < 2:
        sys.exit(f"error: need at least two .bif files in {BIF_DIR}")
    pairs = [(a, b, "") for i, a in enumerate(files) for b in files[i + 1:]]
    mode = "aida" if args.aida else "noaida"
    results = Results(args.results or ROOT / "tests/results_a0.txt",
                      "A0 benchmark: all pairs of tests/random_bifiltrations",
                      {"dim": args.dim, "preprocess": mode},
                      {"play": args.play, "hera": args.hera or "skipped",
                       "files": len(files), "pairs": len(pairs)},
                      HERA_COLUMNS, args.fresh)
    return benchmark(args, pairs, results, dim=args.dim)


def kerber_files(cap):
    """(distinct grades, suffix group, name) of each file with at most cap
    distinct bigrades, sorted. The group keeps homology degrees apart: the
    trailing index of .scc_k / .minpres_k is reflected homology."""
    def group_of(name):
        for tag in (".minpres_", ".scc_"):
            if tag in name:
                return tag[1:] + name.rsplit("_", 1)[-1]
        return "scc" if name.endswith(".scc") else "other"

    def distinct_grades(path):
        lines = path.read_text(errors="replace").splitlines()
        starts = [k for k, line in enumerate(lines) if line.strip() == "scc2020"]
        if not starts:
            return None
        grades = set()
        for line in lines[starts[0] + 3:]:
            line = line.split("#", 1)[0].strip()
            if line == "scc2020":
                break
            parts = line.split(";", 1)[0].split()
            if len(parts) >= 2:
                grades.add((parts[0], parts[1]))
        return len(grades)

    rows = []
    for path in sorted(KERBER_DIR.iterdir()):
        n = distinct_grades(path) if path.is_file() else None
        if n is not None and n <= cap:
            rows.append((n, group_of(path.name), path.name))
    return sorted(rows)


def cmd_kerber(args):
    if not KERBER_DIR.is_dir():
        sys.exit(f"error: missing directory {KERBER_DIR}")
    files = kerber_files(args.max_grades)
    print(f"Feasible files (<= {args.max_grades} distinct bigrades): {len(files)}")
    by_group = defaultdict(list)
    for n, group, name in files:
        print(f"  {n} grades  {group}  {name}")
        by_group[group].append((name, n))
    pairs = []
    for group, items in by_group.items():
        items.sort()
        pairs += [(na + nb, a, b, na, nb, group)
                  for i, (a, na) in enumerate(items) for b, nb in items[i + 1:]]
    pairs = [(KERBER_DIR / a, KERBER_DIR / b, f"  ({group}, grades {na} + {nb})")
             for _, a, b, na, nb, group in sorted(pairs)]
    if not pairs:
        sys.exit("error: the feasible files do not form a same-suffix pair")
    results = Results(args.results or ROOT / "tests/small_instance_results_a0.txt",
                      "A0 benchmark: data_benchmarks/small_instances, same suffix",
                      {"max_grades": args.max_grades,
                       "preprocess": "aida" if args.aida else "noaida"},
                      {"play": args.play, "files": len(files), "pairs": len(pairs)},
                      HERA_COLUMNS[:7], args.fresh)
    return benchmark(args, pairs, results)


# --- command line -----------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sets = parser.add_subparsers(dest="set", required=True, metavar="SET")

    def add(name, func, help, timeout=None, aida=None, hera_timeout=None, results=False,
            dim=False):
        s = sets.add_parser(name, help=help, description=help)
        # Defaults set here take precedence over those of add_argument.
        s.set_defaults(func=func, aida=False, use_hera=False, hera=None, hera_timeout=None)
        s.add_argument("--play", default=os.environ.get("PLAY_BIN"),
                       help="play binary (default: $PLAY_BIN or build/play)")
        s.add_argument("--timeout", type=float, default=timeout, metavar="S",
                       help=f"seconds per play run (default: {timeout or 'none'})")
        if aida is not None:
            s.set_defaults(aida=aida)
            g = s.add_mutually_exclusive_group()
            g.add_argument("--aida", dest="aida", action="store_true",
                           help="run mpfree and AIDA first" + (" (default)" if aida else ""))
            g.add_argument("--no-aida", dest="aida", action="store_false",
                           help="skip mpfree and AIDA" + ("" if aida else " (default)"))
        if dim:
            s.add_argument("--dim", default="0", help="homology dimension (default: 0)")
        if results:
            s.add_argument("--results", type=Path, metavar="PATH", help="results file")
            s.add_argument("--fresh", action="store_true",
                           help="overwrite the results file instead of resuming")
        if func in (cmd_medium, cmd_finite, cmd_random):
            s.set_defaults(use_hera=True, hera=os.environ.get("HERA_BIN"),
                           hera_timeout=hera_timeout)
            s.add_argument("--hera", metavar="PATH",
                           help="Hera's matching_dist (default: $HERA_BIN, else searched)")
            s.add_argument("--no-hera", dest="use_hera", action="store_false",
                           help="do not run Hera")
            s.add_argument("--hera-timeout", type=float, metavar="S",
                           help=f"seconds per Hera run (default: {hera_timeout or 'none'})")
        return s

    small = add("small", cmd_small, "exact A0 counts for small pairs (always --no-aida)",
                timeout=30)
    small.add_argument("--fast", action="store_true", help="skip the ~140-line cases")
    add("medium", cmd_medium, "larger pairs with known distances", timeout=60, aida=False,
        hera_timeout=30)
    add("finite", cmd_finite, "the 18 random pairs with a finite distance", aida=True,
        results=True, dim=True)
    add("random", cmd_random, "every pair of tests/random_bifiltrations", aida=True,
        results=True, dim=True)
    kerber = add("kerber", cmd_kerber, "feasible pairs of data_benchmarks/small_instances",
                 aida=True, results=True)
    kerber.add_argument("--max-grades", type=int, default=200, metavar="N",
                        help="largest number of distinct bigrades per file (default: 200)")

    args = parser.parse_args()
    # Paths given by the user are relative to where they called the script.
    args.play = str(Path(args.play).resolve()) if args.play else str(ROOT / "build/play")
    if args.hera:
        args.hera = str(Path(args.hera).resolve())
    if getattr(args, "results", None):
        args.results = args.results.resolve()
    args.hera = find_hera(args)
    return args


def main():
    sys.stdout.reconfigure(line_buffering=True)
    args = parse_args()
    if not os.access(args.play, os.X_OK):
        sys.exit(f"error: play not executable: {args.play}\n"
                 "Build it with: cmake -S . -B build && cmake --build build")
    os.chdir(ROOT)

    def stop(signum, frame):
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM, stop)
    listen_for_ctrl_n()
    try:
        return args.func(args)
    except KeyboardInterrupt:
        print("\nInterrupted. Completed rows are kept.", file=sys.stderr)
        return 130


if __name__ == "__main__":
    sys.exit(main())
