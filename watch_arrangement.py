#!/usr/bin/env python3
"""Draw the matching-distance arrangement while play refines it.

The red/green coloring exists only on the local overlay of the A0 face
currently being refined. play overwrites one JSON snapshot at a time:

    ./build/play a.scc b.scc --vis /tmp/arr.json
    python3 watch_arrangement.py /tmp/arr.json

The first frames are the uncolored A0 arrangement (gray), written when
maximize_F starts, so the picture is up during the A0 sample. Later
frames follow the face being refined: translucent red and green are the
local subfaces inside that A0 cell, clipped to a window that ignores far
intersections. face_index is the 0-based index among A0 faces that have
an interior sample with s > 0. Faces pruned before refinement leave no
frame, so the index can jump.

    python3 watch_arrangement.py /tmp/arr.json --once -o /tmp/arr.png

--once reads the current file and exits. -o writes a PNG and does not
open a window. Without either flag the script polls about five times a
second and redraws when the file changes. A half-written read is ignored.
"""

import argparse
import json
import os
import sys


RED = "#d95f5f"
GREEN = "#3f9d62"
A0_COLOR = "#8a8a8a"
GUIDE = "#444444"


def load_snapshot(path):
    try:
        with open(path, "r", encoding="utf-8") as handle:
            data = json.load(handle)
    except (OSError, json.JSONDecodeError, UnicodeError):
        return None
    if not isinstance(data, dict):
        return None
    window = data.get("window")
    if not isinstance(window, list) or len(window) != 4:
        return None
    try:
        coords = [float(v) for v in window]
    except (TypeError, ValueError):
        return None
    if any(c != c or c in (float("inf"), float("-inf")) for c in coords):
        return None
    if coords[0] >= coords[1] or coords[2] >= coords[3]:
        return None
    data["window"] = coords
    return data


def draw(ax, data):
    ax.clear()
    s0, s1, t0, t1 = data["window"]
    seen = set()
    for region in data.get("regions") or []:
        if not isinstance(region, dict):
            continue
        color_name = region.get("color")
        poly = region.get("poly") or []
        if color_name not in ("red", "green") or len(poly) < 3:
            continue
        try:
            xs = [float(p[0]) for p in poly]
            ys = [float(p[1]) for p in poly]
        except (TypeError, ValueError, IndexError):
            continue
        xs.append(xs[0])
        ys.append(ys[0])
        color = RED if color_name == "red" else GREEN
        ax.fill(xs, ys, color=color, alpha=0.45, linewidth=0, zorder=2)
        ax.plot(xs, ys, color=color, linewidth=0.6, solid_capstyle="round", zorder=3)
        seen.add(color_name)
    for seg in data.get("a0_segments") or []:
        if not isinstance(seg, list) or len(seg) != 4:
            continue
        try:
            x0, y0, x1, y1 = (float(v) for v in seg)
        except (TypeError, ValueError):
            continue
        ax.plot([x0, x1], [y0, y1], color=A0_COLOR, linewidth=0.7,
                solid_capstyle="round", zorder=1)
    for s in (0.0, 1.0):
        if s0 <= s <= s1:
            ax.axvline(s, color=GUIDE, linewidth=0.8, linestyle="--", zorder=4)
    ax.set_xlim(s0, s1)
    ax.set_ylim(t0, t1)
    ax.set_aspect("auto")
    ax.set_xlabel("s")
    ax.set_ylabel("t")
    phase = data.get("phase")
    a0_lines = data.get("a0_lines", 0)
    total = data.get("total_lines", 0)
    deleted = data.get("green_green_deletions", 0)
    if phase == "refine":
        face = int(data.get("face_index", -1)) + 1
        count = int(data.get("face_count", 0))
        title = "A0 face %d/%d    lines %s (A0 %s)    green-green %s" % (
            face, count, total, a0_lines, deleted)
    else:
        title = "A0 overview    A0 lines %s" % a0_lines
    ax.set_title(title)
    handles = []
    import matplotlib.patches as patches
    if "red" in seen:
        handles.append(patches.Patch(color=RED, alpha=0.45, label="red"))
    if "green" in seen:
        handles.append(patches.Patch(color=GREEN, alpha=0.45, label="green"))
    if data.get("a0_segments"):
        handles.append(plt_line_proxy(ax))
    if handles:
        ax.legend(handles=handles, loc="upper right", frameon=False)
    ax.figure.canvas.draw_idle()


def plt_line_proxy(ax):
    line, = ax.plot([], [], color=A0_COLOR, linewidth=0.7, label="A0")
    return line


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Draw play --vis snapshots of the (s, t) arrangement.")
    parser.add_argument("path", help="JSON file written by play --vis")
    parser.add_argument("--once", action="store_true",
                        help="Draw the current snapshot and exit")
    parser.add_argument("-o", "--output",
                        help="Write a PNG and exit (no window)")
    parser.add_argument("--interval", type=float, default=0.2,
                        help="Seconds between polls (default 0.2)")
    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    if args.output:
        import matplotlib
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if args.output or args.once:
        data = load_snapshot(args.path)
        if data is None:
            sys.stderr.write("No snapshot in %s\n" % args.path)
            return 1
        fig, ax = plt.subplots(figsize=(8, 6))
        draw(ax, data)
        if args.output:
            fig.savefig(args.output, dpi=140, bbox_inches="tight")
        else:
            plt.show()
        return 0

    plt.ion()
    fig, ax = plt.subplots(figsize=(8, 6))
    try:
        fig.canvas.manager.set_window_title("matching distance")
    except AttributeError:
        pass
    ax.set_xlabel("s")
    ax.set_ylabel("t")
    ax.set_title("waiting for %s" % args.path)
    seen_mtime = None
    while plt.fignum_exists(fig.number):
        mtime = None
        try:
            mtime = os.stat(args.path).st_mtime_ns
        except OSError:
            mtime = None
        if mtime is not None and mtime != seen_mtime:
            data = load_snapshot(args.path)
            if data is not None:
                draw(ax, data)
                seen_mtime = mtime
        try:
            plt.pause(max(args.interval, 0.05))
        except Exception:
            break
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
