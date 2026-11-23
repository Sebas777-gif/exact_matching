#!/usr/bin/env python3
"""Regression tests and benchmarks for play, which builds the A0 arrangement
and evaluates F on its vertices (max F on A0 vertices, a lower bound on the
matching distance), and with --distance computes the matching distance.

  ./bench.py small [--fast]    exact A0 counts and max F for small pairs; exit 1 on a mismatch
  ./bench.py medium [--aida]   larger pairs with known distances; checks max F against them
  ./bench.py finite            the 18 random pairs with a finite matching distance
  ./bench.py random            every pair of tests/random_bifiltrations/*.bif
  ./bench.py kerber            feasible pairs of data_benchmarks/small_instances
  ./bench.py kerber-extended   repeats of the same construction in data_benchmarks/instances

finite, random, kerber, and kerber-extended append one row per pair to a results file and
resume where they stopped. --a0-only skips F (play --a0-only) in every set.
--incremental (with --samples K, --seed N, --bound B) uses the red/green
insertion instead of evaluating F at every A0 vertex; the result must be the
same. --distance adds the local phase and reports the matching distance.
./bench.py SET --help lists the options of a set.
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
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BIF_DIR = ROOT / "tests" / "random_bifiltrations"
KERBER_DIR = ROOT / "data_benchmarks" / "small_instances"
KERBER_ALL_DIR = ROOT / "data_benchmarks" / "instances"


# --- pairs ------------------------------------------------------------------

def v10(sample):
    return f"tests/random_bifiltrations/r_bif_v_10_dim_2_n_top_simplices_50_sample_{sample}.bif"


# small: name, expected "lines vertices edges faces", expected max F on the
# A0 vertices (exact, as play prints it), expected matching distance d
# (checked with --distance), play arguments. Always run with --no-aida. Lines are the flat strip's segments; the other three are the
# totals of the flat and steep strips.
#
# d is the matching distance each pair was built for. Max F is derived by
# hand unless the comment says otherwise. Since F <= d on every slice, a
# vertex with F = d shows max F = d. Pushes: max(y, s x + t) in the flat
# strip, max(x, s y + t) in the steep strip.
SMALL = [
    # Rank-1 intervals (d = 1). Grades (0,0), (1,1), (2,2): all three
    # duals meet at (s, t) = (1, 0), so the wall s = 1 is a single point
    # and is left out (3 segments + the wall s = 0). Each strip has
    # 4 vertices, 5 edges, 3 faces. Max F = 1: at the vertex (1, 0) the
    # bars are [0, 2) and [0, 1), and pairing them costs 1. The .bif
    # version adds one infinite bar at push(0, 0) to both sides (cost 0).
    ("interval scc", "4 8 10 6", "1", "1", ["tests/interval_long.scc", "tests/interval_short.scc"]),
    ("interval bif", "4 8 10 6", "1", "1", ["tests/interval_long.bif", "tests/interval_short.bif"]),
    # Two origin vertices merging at (51,0) vs a single infinite bar
    # (d = 25.5). The .scc version adds unused grades at x, y in
    # {10, 20, 30} that only cut A0; the .bif version does not (its second
    # vertex and edge at (30,0) give a bar of length 0). Max F = 51/2: the
    # infinite bars both start at push(0,0); the finite bar
    # [max(0, t), max(0, 51 s + t)) goes to the diagonal, which costs 51/2
    # at the vertex (1, 0) (end of the dual of (0,0)). Same in the steep strip.
    ("far_merge vs A0 anchors", "10 25 37 16", "51/2", "51/2",
     ["tests/far_merge.scc", "tests/a0_anchors.scc"]),
    ("far_merge bif", "4 10 12 6", "51/2", "51/2", ["tests/far_merge.bif", "tests/a0_anchors.bif"]),
    # Horizontal merge at (50,0) vs vertical merge at (0,50) (d = 25).
    # Max F = 25: at the vertex (0, 0), the horizontal slice y = 0, the
    # horizontal bar has length 0 and the vertical one is [0, 50).
    ("horizontal vs vertical merge", "5 10 12 6", "25", "25",
     ["tests/merge_horizontal.scc", "tests/merge_vertical.scc"]),
    # Different numbers of infinite bars (d = infinity): F is infinite
    # on every slice.
    ("unmatched infinite bars", "3 6 6 4", "inf", "inf", ["tests/free_one.scc", "tests/free_two.scc"]),
    # Raw H_2 presentations of two prism 2-skeleta (d = infinity): 8 vs 9
    # generators and no relations, so 8 vs 9 infinite bars.
    ("prism H_2", "7 14 20 10", "inf", "inf", ["tests/prism_1.bif", "tests/prism_2.bif", "--dim", "2"]),
    # L of two intervals vs one spanning interval (d = 1.5). The maximum of
    # F sits at (s, t) = (1, 1/2), on an A0 edge, not an A0 vertex. Max F
    # on the vertices is 1, checked at all 7 flat and 8 steep vertices; it
    # is reached at (1, 0) in both strips (span [0,3) pairs with [1,2) at
    # cost 1, the other L bar [2,3) goes to the diagonal at 1/2).
    ("ell vs span", "7 15 23 12", "1", "3/2", ["tests/ell.scc", "tests/span.scc"]),
]

# small, skipped with --fast.
SMALL_SLOW = [
    # far_merge with 140 unused relation grades (d = 25.5). The extra
    # relations are empty columns, so the barcodes and the argument are
    # those of "far_merge vs A0 anchors".
    ("far_merge vs 140 A0 anchors", "144 427 707 284", "51/2", "51/2",
     ["tests/far_merge.scc", "tests/a0_anchors_medium.scc"]),
    # Merge pair with ~140 distinct grades (d = 25), again only empty extra
    # relations: same argument as "horizontal vs vertical merge".
    ("horizontal vs vertical merge medium", "145 608 1072 468", "25", "25",
     ["tests/merge_horizontal_medium.scc", "tests/merge_vertical_medium.scc"]),
    # Real 10-vertex 1-skeleta, H_0 (d = 3 and 2.5). Max F not derived by
    # hand: recorded from play. It equals d, the largest value possible.
    # d is None: without AIDA each module is one large summand, and the
    # local phase takes minutes here, so --distance skips these two cases.
    ("v10 bif sample 1 vs 2", "74 1131 2338 1211", "3", None, [v10(1), v10(2), "--dim", "0"]),
    ("v10 bif sample 1 vs 5", "71 983 2032 1053", "5/2", None, [v10(1), v10(5), "--dim", "0"]),
]

# medium: name, known matching distance d, how max F on the A0 vertices is
# checked against d ("=" or "<="), the two inputs. A bare name is a
# synthetic pair: tests/<name>.scc (empty dummy relation columns that only
# add grades) without AIDA, and tests/<name>.bif with AIDA, because AIDA
# removes the dummy columns but keeps the extra short H_0 bars of the .bif
# 1-skeleta. Hera always reads the .bif version.
#
# The synthetic pairs reach d at an A0 vertex, by the arguments of the
# matching small cases ("far_merge vs A0 anchors", "horizontal vs vertical
# merge"); the extra short bars of the .bif versions cost at most 1/2. For
# the v10 pairs only max F <= d is known.
MEDIUM = [
    ("far_merge vs 140 A0 anchors", "25.5", "=", "far_merge", "a0_anchors_medium"),
    ("far_merge_41 vs 140 A0 anchors", "20.5", "=", "far_merge_41", "a0_anchors_medium"),
    ("far_merge vs 10x10 A0 grid", "25.5", "=", "far_merge", "a0_anchors_grid"),
    ("horizontal vs vertical merge L=50", "25", "=", "merge_horizontal_medium",
     "merge_vertical_medium"),
    ("horizontal vs vertical merge L=40", "20", "=", "merge_h40_medium", "merge_v40_medium"),
    ("horizontal vs vertical merge L=60", "30", "=", "merge_h60_medium", "merge_v60_medium"),
] + [
    (f"v10 n50 sample {i} vs {j}", d, "<=", v10(i), v10(j))
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
    if args.a0_only:
        cmd.append("--a0-only")
    elif args.incremental or args.distance:
        cmd += ["--incremental", "--samples", str(args.samples), "--seed", str(args.seed)]
        if args.no_skip:
            cmd.append("--no-skip")
        cmd += ["--bound", args.bound]
        if args.distance:
            cmd += ["--distance", "--max-candidates", str(args.max_candidates),
                    "--split-depth", str(args.split_depth),
                    "--level-lines", str(args.level_lines)]
    return run(cmd, args.timeout)


# The values read from play's output, in the order of the results columns.
#   max_F, max_F_decimal     max F on the A0 vertices, exact ("3/2", "inf") and decimal
#   evaluations              F evaluations for that maximum
#   pruned_pct               share of the A0 vertex evaluations avoided (incremental only)
#   distance, distance_decimal, distance_status
#                            the matching distance (--distance only); status "exact" or
#                            "lower bound" (faces over the candidate limit)
#   local_evaluations        F evaluations of the local phase (--distance only)
STATS = ["a0_lines", "vertices", "edges", "faces", "max_F", "max_F_decimal", "evaluations",
         "pruned_pct", "distance", "distance_decimal", "distance_status", "local_evaluations"]


def play_stats(r, args):
    """The STATS values from play's output, as a dict of strings ("" where
    play did not print a value), or None if a required line is missing."""
    def find(pattern):
        return re.search(pattern, r.stdout, re.M)

    lines = find(r"^A0 lines: (\d+)")
    arr = find(r"^Arrangement has (\d+) vertices, (\d+) edges, and (\d+) faces")
    max_f = find(r"^Max F on A0 vertices: (\S+)(?: \((\S+)\))?")
    dist = find(r"^Matching distance: (\S+?)(?: \(([^)]*)\))?, (exact|lower bound)")
    if r.status != "ok" or not (lines and arr) or not (max_f or args.a0_only):
        return None
    if args.distance and not args.a0_only and not dist:
        return None
    st = dict.fromkeys(STATS, "")
    st["a0_lines"] = lines[1]
    st["vertices"], st["edges"], st["faces"] = arr.groups()
    if max_f:
        st["max_F"], st["max_F_decimal"] = max_f[1], max_f[2] or max_f[1]
    if m := find(r"^F evaluations: (\d+)"):
        st["evaluations"] = m[1]
    if m := find(r"^Pruned: (-?[\d.]+)%"):
        st["pruned_pct"] = m[1]
    if dist:
        st["distance"], st["distance_decimal"], st["distance_status"] = (
            dist[1], dist[2] or dist[1], dist[3])
    if m := find(r"^Local phase: .*?, (\d+) evaluations,"):
        st["local_evaluations"] = m[1]
    return st


def as_number(value):
    """A distance string ("3/2", "25.5", "inf", "infinity") as a Fraction or
    float('inf'); None if it is not a number (error, timeout, ...)."""
    if value.lower().lstrip("+") in ("inf", "infinity"):
        return float("inf")
    try:
        return Fraction(value)
    except ValueError:
        return None


def stat_fields(r, stats):
    """The STATS columns of a results row."""
    if stats:
        return [stats[k] for k in STATS]
    if r.status == "skip":
        return ["???"] * len(STATS)
    return ["timeout" if r.status == "timeout" else "error"] + [""] * (len(STATS) - 1)


def summary(st):
    """The play values of one run, as one line of text."""
    text = f"{st['a0_lines']} lines, {st['vertices']} vertices, {st['edges']} edges, {st['faces']} faces"
    if st["max_F"]:
        text += f", max F {st['max_F']}"
        if st["max_F_decimal"] != st["max_F"]:
            text += f" ({st['max_F_decimal']})"
    if st["evaluations"]:
        text += f", {st['evaluations']} F evaluations"
    if st["pruned_pct"]:
        text += f" ({st['pruned_pct']}% pruned)"
    if st["distance"]:
        text += f", distance {st['distance']}"
        if st["distance_decimal"] != st["distance"]:
            text += f" ({st['distance_decimal']})"
        text += f" [{st['distance_status']}], {st['local_evaluations']} local evaluations"
    return text


def outcome(r, stats=None, lines=8):
    """One-line summary of a play run; on failure, the last lines of its output."""
    if stats:
        return summary(stats) + f"  ({r.seconds}s)"
    if r.status == "ok":
        return f"FAIL a line is missing in play's output  ({r.seconds}s)"
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

TAB = 4   # tab width that aligned result files assume (editors; less -x4)


class Results:
    """Tab-separated file: '# key: value' header lines, a column line, one row
    per pair. Rows are appended right away, so a crash loses nothing. When a
    benchmark finishes, align() pads the columns with tabs; a file read back
    is turned into plain tab-separated rows again before new rows are added."""

    def __init__(self, path, title, settings, info, columns, fresh):
        self.path = path
        if fresh or not path.exists():
            header = [f"# {title}", f"# started: {datetime.now().isoformat(timespec='seconds')}"]
            header += [f"# {k}: {v}" for k, v in {**settings, **info}.items()]
            path.write_text("\n".join(header + ["\t".join(columns)]) + "\n")
            self.header, self.columns, self.rows = "\t".join(columns), list(columns), []
            return
        lines = path.read_text().splitlines()
        found = dict(line[2:].split(": ", 1) for line in lines
                     if line.startswith("# ") and ": " in line)
        for k, v in settings.items():
            if found.get(k) != str(v):
                sys.exit(f"error: {path} has {k}: {found.get(k)}, not {v}. "
                         "Use the same options, or --fresh to start over.")
        comments = [line for line in lines if line.startswith("#")]
        table = [line for line in lines if line and not line.startswith("#")]
        self.header = table[0] if table else "\t".join(columns)
        self.columns = re.split(r"\t+", self.header)
        self.rows = [self.split(line) for line in table[1:]]
        if table and table[0] != "\t".join(self.columns):   # aligned: write it back plain
            path.write_text("\n".join(comments + ["\t".join(r) for r in [self.columns] + self.rows])
                            + "\n")

    def split(self, line):
        """Cells of a plain row, or of a row aligned by align(), which may
        have empty cells: those are found by the tab stops of the columns."""
        cells = line.split("\t")
        if len(cells) == len(self.columns):
            return cells
        starts = [m.start() for m in re.finditer(r"[^\t ]+", self.header.expandtabs(TAB))]
        wide = line.expandtabs(TAB)
        return [wide[a:b].strip() for a, b in zip(starts, starts[1:] + [None])]

    @staticmethod
    def widths(rows):
        """Width of each column in characters: whole tabs of TAB
        characters, at least one more than its longest cell."""
        return [(max(len(c) for c in col) // TAB + 1) * TAB for col in zip(*rows)]

    def append(self, row):
        with self.path.open("a") as f:
            f.write("\t".join(str(x) for x in row) + "\n")

    def align(self):
        """Pads every cell but the last of each row with tabs so that the
        columns line up. Header comment lines stay as they are."""
        lines = self.path.read_text().splitlines()
        comments = [line for line in lines if line.startswith("#")]
        table = [line for line in lines if line and not line.startswith("#")]
        if not table:
            return
        self.header = table[0]
        columns = self.columns = re.split(r"\t+", self.header)
        rows = [columns] + [r for r in map(self.split, table[1:])]
        rows = [r + [""] * (len(columns) - len(r)) for r in rows]
        widths = self.widths(rows)
        out = ["".join(c + "\t" * ((w - len(c) - 1) // TAB + 1) for c, w in zip(r[:-1], widths)) + r[-1]
               for r in rows]
        self.path.write_text("\n".join(comments + out) + "\n")


def benchmark(args, pairs, results, dim=None):
    """Run play (and Hera, if args.hera) on each pair (a, b, note).

    A pair counts as done when its row has A0 values and, with Hera
    columns, a Hera value; done pairs are skipped when resuming. ??? marks a
    pair skipped with Ctrl+N; it is not retried."""
    with_hera = dim is not None
    hera_col = 2 + len(STATS) + 1  # after file_a, file_b, STATS, seconds
    done = {(r[0], r[1]) for r in results.rows
            if len(r) > 2 and r[2] and (not with_hera or (len(r) > hera_col and r[hera_col]))}
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
                results.append([a.name, b.name, *[""] * (len(STATS) + 1), *hera])
            raise
        stats = play_stats(r, args)
        print("  play: " + outcome(r, stats))
        results.append([a.name, b.name, *stat_fields(r, stats), r.seconds, *hera])
    results.align()
    print(f"Done. Results: {results.path}")
    return 0


# --- the sets ---------------------------------------------------------------------

def cmd_small(args):
    cases = SMALL + ([] if args.fast else SMALL_SLOW)
    print(f"play: {args.play}  (--no-aida, timeout {args.timeout:g}s)")
    print(f"method: {method(args)}")
    print("columns: A0 lines, vertices, edges, faces" +
          ("" if args.a0_only else ", max F on A0 vertices") +
          (", d = matching distance" if args.distance and not args.a0_only else "") +
          ("" if args.a0_only else "; then the F evaluations") + "\n")
    failed = skipped = 0
    for name, expected_a0, expected_f, expected_d, play_args in cases:
        if args.distance and not args.a0_only and expected_d is None:
            print(f"{name:<42}skipped (too slow with --distance and --no-aida)")
            skipped += 1
            continue
        expected = expected_a0 if args.a0_only else f"{expected_a0} {expected_f}"
        if args.distance and not args.a0_only:
            expected += f" d={expected_d} exact"
        print(f"{name:<42}", end="", flush=True)
        r = run_play(args, play_args)
        stats = play_stats(r, args)
        got = ""
        if stats:
            got = " ".join(stats[k] for k in ["a0_lines", "vertices", "edges", "faces", "max_F"]
                           if stats[k])
            if stats["distance"]:
                got += f" d={stats['distance']} {stats['distance_status']}"
        if got == expected:
            evals = ""
            if stats["evaluations"]:
                evals = f"  [{stats['evaluations']} evaluations"
                if stats["pruned_pct"]:
                    evals += f", {stats['pruned_pct']}% pruned"
                if stats["local_evaluations"]:
                    evals += f", {stats['local_evaluations']} local"
                evals += "]"
            print(f"ok  {got}{evals}  ({r.seconds}s)")
            continue
        failed += 1
        if stats:
            print(f"FAIL got '{got}' expected '{expected}'  ({r.seconds}s)")
        else:
            print(outcome(r, lines=20))
    ran = len(cases) - skipped
    note = f" ({skipped} skipped)" if skipped else ""
    print(f"\nAll {ran} tests passed{note}." if not failed
          else f"\n{ran - failed} passed, {failed} failed{note}.")
    return 1 if failed else 0


def cmd_medium(args):
    def path(entry, ext):
        return entry if "/" in entry else f"tests/{entry}.{ext}"

    print(f"play: {args.play}  ({'AIDA' if args.aida else '--no-aida'}, timeout {args.timeout:g}s)")
    print(f"method: {method(args)}")
    if args.hera:
        print(f"hera: {args.hera}  (reference only, timeout {args.hera_timeout:g}s)\n")
    else:
        print("hera: skipped\n")
    table = ["\t".join(["name", "known_dist", *STATS, "seconds", "status", "hera", "hera_seconds"])]
    failed = skipped = 0
    try:
        for name, known, check, a, b in MEDIUM:
            print(f"{name:<40}", end="", flush=True)
            # Real pairs (paths) without AIDA: one large summand per module,
            # and the local phase takes minutes; use --aida with --distance.
            if args.distance and not args.aida and "/" in a:
                print("skipped (too slow with --distance and --no-aida; use --aida)")
                skipped += 1
                continue
            hera, hera_sec = "skip", ""
            hera_a, hera_b = path(a, "bif"), path(b, "bif")
            if args.hera and os.path.exists(hera_a) and os.path.exists(hera_b):
                hera, hera_sec = run_hera(args, hera_a, hera_b, "0")
            ext = "bif" if args.aida else "scc"
            r = run_play(args, [path(a, ext), path(b, ext)], "0")
            stats = play_stats(r, args)
            ref = f"known {known}" + (f", hera {hera} ({hera_sec}s)" if hera != "skip" else "")
            status = "ok" if stats else r.status
            if stats and not args.a0_only:
                f, d = as_number(stats["max_F"]), Fraction(known)
                if not (f == d if check == "=" else f <= d):
                    status = f"max F not {check} known"
                # The matching distance must equal the known distance; a lower
                # bound (faces over the candidate limit) must not exceed it.
                if args.distance:
                    got_d = as_number(stats["distance"])
                    if stats["distance_status"] == "exact" and got_d != d:
                        status = "distance not = known"
                    elif stats["distance_status"] != "exact" and not (f <= got_d <= d):
                        status = "lower bound not between max F and known"
            if stats:
                text = f"{stats['a0_lines']} lines, {stats['vertices']} vertices"
                if not args.a0_only:
                    text += f", max F {stats['max_F']} ({check} {known}), {stats['evaluations']} evaluations"
                    if stats["pruned_pct"]:
                        text += f" ({stats['pruned_pct']}% pruned)"
                if stats["distance"]:
                    text += (f", distance {stats['distance']} [{stats['distance_status']}],"
                             f" {stats['local_evaluations']} local evaluations")
                print(f"{'ok  ' if status == 'ok' else 'FAIL '}{text}  ({r.seconds}s)  {ref}")
                if status != "ok":
                    print(f"     {status}")
            else:
                print(f"{ref}\n{outcome(r)}")
            if status != "ok":
                failed += 1
            table.append("\t".join([name, known, *stat_fields(r, stats), r.seconds, status,
                                    hera, hera_sec]))
    finally:
        print("\n" + "\n".join(table))
    ran = len(MEDIUM) - skipped
    note = f" ({skipped} skipped)" if skipped else ""
    print(f"\nAll {ran} medium cases passed{note}." if not failed
          else f"\n{ran - failed} passed, {failed} failed{note}.")
    return 1 if failed else 0


HERA_COLUMNS = ["file_a", "file_b", *STATS, "seconds", "hera_distance", "hera_seconds"]


def evaluated(args):
    """Settings value: whether and how the max_F column is filled."""
    return "skipped (--a0-only)" if args.a0_only else f"on A0 vertices, {method(args)}"


def method(args):
    """How play computes max F: every A0 vertex, or the red/green insertion
    (with --distance followed by the local phase)."""
    if args.a0_only:
        return "A0 only (--a0-only)"
    if args.incremental or args.distance:
        skip = "no skipping" if args.no_skip else "skipping covered vertices"
        text = f"incremental, {args.samples} samples, seed {args.seed}, {skip}, {args.bound} bound"
        if args.distance:
            text += f", local phase (at most {args.max_candidates} candidates per face"
            if args.split_depth:
                text += f", split depth {args.split_depth}"
            if args.level_lines:
                text += f", level lines up to {args.level_lines}"
            text += ")"
        return text
    return "all A0 vertices"


def cmd_finite(args):
    pairs = [(BIF_DIR / a, BIF_DIR / b, f"  (Hera: {h})") for a, b, h in FINITE]
    missing = [p for a, b, _ in pairs for p in (a, b) if not p.exists()]
    if missing:
        sys.exit("error: missing " + ", ".join(map(str, missing)))
    mode = "aida" if args.aida else "noaida"
    results = Results(args.results or ROOT / f"tests/finite_results_a0_{mode}.txt",
                      "A0 benchmark: the 18 finite pairs of tests/results.txt",
                      {"dim": args.dim, "preprocess": mode, "max_F": evaluated(args)},
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
                      {"dim": args.dim, "preprocess": mode, "max_F": evaluated(args)},
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
                       "preprocess": "aida" if args.aida else "noaida",
                       "max_F": evaluated(args)},
                      {"play": args.play, "files": len(files), "pairs": len(pairs)},
                      HERA_COLUMNS[:-2], args.fresh)
    return benchmark(args, pairs, results)


# kerber-extended: data_benchmarks/README.txt names files
#   [instance_type]_[instance_size]_[instance_number].scc_[dimension],
# where instance numbers 1-5 are repeats of the same construction with
# different random seeds. Pairs of repeats (same type, size, and dimension)
# are natural test pairs: same distribution, same homology degree, and
# usually a finite nonzero matching distance. The full_rips, function_rips,
# and random files use other spellings of the same scheme.
CONSTRUCTION_PATTERNS = [
    re.compile(r"^(?P<group>.+_\d+)_(?P<instance>\d+)\.(?P<dim>scc_\d+)$"),
    re.compile(r"^(?P<group>full_rips_size_\d+)_instance_(?P<instance>\d+)\.(?P<dim>scc\.minpres_\d+)$"),
    re.compile(r"^(?P<group>.+)_(?P<instance>\d+)\.(?P<dim>scc\.minpres_\d+|scc)$"),
]


def construction(name):
    """(construction, instance number) of a benchmark file name, or None.
    The construction includes the dimension suffix, so only files of the same
    homology degree share one."""
    for pattern in CONSTRUCTION_PATTERNS:
        if m := pattern.match(name):
            return f"{m['group']}.{m['dim']}", m["instance"]
    return None


def count_grades(path, cap):
    """Distinct bigrades of the first scc2020 block of path, or None if there
    are more than cap. Stops reading as soon as the cap is exceeded, so large
    files cost little."""
    grades, started, header = set(), False, 0
    with path.open(errors="replace") as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            if not started:
                started = line == "scc2020"
                continue
            if header < 2:   # the number of parameters and the sizes
                header += 1
                continue
            if line == "scc2020":
                break
            parts = line.split(";", 1)[0].split()
            if len(parts) >= 2:
                grades.add((parts[0], parts[1]))
                if len(grades) > cap:
                    return None
    return len(grades) if started else None


def cmd_kerber_extended(args):
    if not KERBER_ALL_DIR.is_dir():
        sys.exit(f"error: missing directory {KERBER_ALL_DIR}")
    groups = defaultdict(list)
    names = sorted(p.name for p in KERBER_ALL_DIR.iterdir() if p.is_file())
    print(f"Scanning {len(names)} files for <= {args.max_grades} distinct bigrades ...", flush=True)
    for name in names:
        key = construction(name)
        if key is None:
            continue
        n = count_grades(KERBER_ALL_DIR / name, args.max_grades)
        if n is not None:
            groups[key[0]].append((name, n))
    pairs = []
    for group in sorted(groups):
        items = sorted(groups[group])
        if len(items) < 2:
            continue
        print(f"  {group}: {len(items)} files, grades {[n for _, n in items]}")
        pairs += [(na + nb, a, b, na, nb, group)
                  for i, (a, na) in enumerate(items) for b, nb in items[i + 1:]]
    pairs = [(KERBER_ALL_DIR / a, KERBER_ALL_DIR / b, f"  ({group}, grades {na} + {nb})")
             for _, a, b, na, nb, group in sorted(pairs)]
    if not pairs:
        sys.exit("error: no construction has two feasible repeats")
    files = len({p for a, b, _ in pairs for p in (a, b)})
    print(f"{files} files, {len(pairs)} pairs of repeats\n")
    results = Results(args.results or ROOT / "tests/kerber_extended_results.txt",
                      "Benchmark: repeats of the same construction in data_benchmarks/instances",
                      {"max_grades": args.max_grades,
                       "preprocess": "aida" if args.aida else "noaida",
                       "max_F": evaluated(args)},
                      {"play": args.play, "files": files, "pairs": len(pairs)},
                      HERA_COLUMNS[:-2], args.fresh)
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
        s.add_argument("--a0-only", action="store_true",
                       help="only build A0; do not evaluate F (play --a0-only)")
        s.add_argument("--incremental", action="store_true",
                       help="red/green insertion instead of evaluating F at every A0 vertex")
        s.add_argument("--samples", type=int, default=50, metavar="K",
                       help="with --incremental: random A0 vertices evaluated first (default: 50)")
        s.add_argument("--seed", type=int, default=1, metavar="N",
                       help="with --incremental: random seed (default: 1)")
        s.add_argument("--no-skip", action="store_true",
                       help="with --incremental: evaluate vertices covered by a sample too")
        s.add_argument("--bound", choices=["both", "barbox", "variation", "lipschitz"],
                       default="both", help="with --incremental: face bound (default: both)")
        s.add_argument("--distance", action="store_true",
                       help="also compute the matching distance (play --distance; implies "
                            "--incremental)")
        s.add_argument("--max-candidates", type=int, default=500, metavar="K",
                       help="with --distance: skip red faces with more candidate functions; "
                            "the distance is then only a lower bound (default: 500)")
        s.add_argument("--split-depth", type=int, default=0, metavar="D",
                       help="with --distance: split such faces in two up to D times first "
                            "(default: 0)")
        s.add_argument("--level-lines", type=int, default=0, metavar="L",
                       help="with --distance: first try to certify each red face via its "
                            "level lines e = Delta, if it has at most L of them (default: 0, off)")
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

    small = add("small", cmd_small,
                "exact A0 counts and max F for small pairs (always --no-aida)",
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
    extended = add("kerber-extended", cmd_kerber_extended,
                   "pairs of repeats (instance numbers) of the same construction in "
                   "data_benchmarks/instances", aida=True, results=True)
    extended.add_argument("--max-grades", type=int, default=400, metavar="N",
                          help="largest number of distinct bigrades per file (default: 400)")
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
