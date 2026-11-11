# Sample presentations, tests, and benchmarks

`play` builds the coarse arrangement $A_0$ (flat and steep strip) and
evaluates $F$ at its vertices. The sets of `../bench.py` check or record its
output: the number of $A_0$ lines, the vertex, edge, and face counts, and
max $F$ on the $A_0$ vertices. The matching distances listed here are the
values each pair was built for. Since $F \le d_\mathrm{match}$ on every
slice, max $F$ on the vertices can never exceed them; it reaches them
unless the maximum of $F$ lies off the vertices.

The pair lists of all sets are tables at the top of `bench.py`.
`./bench.py SET --help` lists the options of a set. In every set, Ctrl+N
skips the running pair and Ctrl+C stops.

`interval_long.scc` and `interval_short.scc` are two interval modules with
matching distance **1**.

```bash
./play tests/interval_long.scc tests/interval_short.scc
./play tests/interval_long.bif tests/interval_short.bif
```

The `.bif` pair is the same comparison as a 1-critical bifiltration: two
vertices at (0, 0) and one edge at the death grade. The $H_0$ presentation
has two generators at the origin and one relation, so the diagrams also
contain an infinite bar (the connected component). Those infinite bars match
at cost 0; the finite bars still give matching distance 1.

## The modules

Each file is a rank-1 interval: one generator, killed by one relation.

| File | Generator | Relation | Support |
| --- | --- | --- | --- |
| `interval_long.scc` | (0, 0) | (2, 2) | $k$ at $x$ iff $(0,0)\le x$ and $(2,2)\not\le x$ |
| `interval_short.scc` | (0, 0) | (1, 1) | same, dying at (1, 1) |

On any affine line $\ell$ of positive slope the slice of an interval is a
single bar, from the push of the birth grade to the push of the death grade.

## Why the matching distance is 1

The matching distance is $\sup_\ell F(\ell)$, where $F(\ell)=w(\ell)\,d_B(M^\ell,N^\ell)$
and $w(\ell)$ is the weight from Kerber–Lesnick–Oudot (so that $F$ measures
distance in the original $L_\infty$ metric on the plane).

On every slice both modules contribute one finite bar, born at the same point
(the push of the origin) and dying at the pushes of (2, 2) and (1, 1). The
bottleneck cost of matching those two bars is the weighted distance between
the two death points. Unmatching both bars cannot cost more than that.

The weighted distance between two grades is at most their $L_\infty$ distance
in the plane, and $\|(2,2)-(1,1)\|_\infty=1$. So $F(\ell)\le 1$ for every
slice.

Equality holds on the diagonal $\ell\colon y=x$ (slope $s=1$, intercept $t=0$):
both death points lie on $\ell$, the weighted grades are 0 (birth), 1 and 2
(deaths), the bars are $(0,2)$ and $(0,1)$, and $d_B=1$. For $s=1$ the
weight times Euclidean factor cancel, so $F=1$.

The same value 1 is the interleaving distance of the two modules: you must
shift by 1 in both coordinates to push (1, 1) onto (2, 2).

## Small regression tests

`./bench.py small` runs play on each pair and compares the five numbers
(lines, vertices, edges, faces, max $F$ on the $A_0$ vertices) with the
expected ones. Max $F$ is compared exactly, as a fraction. Every case uses
`--no-aida`. It exits with 1 if a case fails, times out (30s), or is
skipped. With `--a0-only` only the four $A_0$ counts are checked. With
`--incremental` the same expected values apply: the red/green insertion
must find the same maximum, and the line also shows how many $F$
evaluations it needed.

```bash
./bench.py small                     # all cases
./bench.py small --fast              # tiny cases only
./bench.py small --play other/play   # or PLAY_BIN=other/play
```

| Pair | $A_0$ (L V E F) | max $F$ | $d_\mathrm{match}$ | Why the pair exists |
| --- | --- | --- | --- | --- |
| `interval_long` vs `interval_short` (`.scc` and `.bif`) | 4 8 10 6 | 1 | 1 | Rank-1 intervals. All three duals meet at $(1,0)$ on $s=1$ |
| `far_merge` vs `a0_anchors` | 10 25 37 16 | 51/2 | 25.5 | $A_0$ faces cut by extra dual segments at $x,y\in\{10,20,30\}$ |
| `far_merge.bif` vs `a0_anchors.bif` | 4 10 12 6 | 51/2 | 25.5 | Same distance as a 1-critical $H_0$ bifiltration (no extra grades) |
| `merge_horizontal` vs `merge_vertical` | 5 10 12 6 | 25 | 25 | Type-III slices; the maximum is at the horizontal slice $y=0$ ($s=0$) |
| `free_one` vs `free_two` | 3 6 6 4 | inf | infinity | Unmatched infinite bars, $s>0$ |
| `prism_1` vs `prism_2` `--dim 2` | 7 14 20 10 | inf | infinity | Raw $H_2$ presentations of two prism 2-skeleta (8 vs 9 triangles, no 3-simplices) |
| `ell` vs `span` | 7 15 23 12 | 1 | 1.5 | Maximum of $F$ on an $A_0$ edge, not at an $A_0$ vertex |
| `far_merge` vs `a0_anchors_medium` | 144 427 707 284 | 51/2 | 25.5 | Same as above with 140 unused grades |
| `merge_horizontal_medium` vs `merge_vertical_medium` | 145 608 1072 468 | 25 | 25 | Merge pair with ~140 distinct grades |
| `r_bif_v_10_…_sample_1` vs `sample_2` (`--dim 0`) | 74 1131 2338 1211 | 3 | 3 | Real $H_0$ 1-skeleton |
| `r_bif_v_10_…_sample_1` vs `sample_5` (`--dim 0`) | 71 983 2032 1053 | 5/2 | 2.5 | Same, different sample |

The last four rows are skipped with `--fast`. L is the number of segments
in the flat strip (distinct grades, plus the walls $s=0$ and $s=1$ unless
one is a single point). V, E, F are the totals of the flat and steep
strips. The counts were recorded from the current CGAL build. The small
vertex counts can be checked by hand, as in the interval case; the vertex
counts of the `.scc` rows were also recounted with exact fractions in
Python. Each strip satisfies Euler's formula $V - E + F = 2$.

Max $F$ is derived by hand for every row except the two `r_bif_v_10`
rows; the derivations are comments in the `SMALL` table of `bench.py`. In
short:

- **Intervals:** at the vertex $(1,0)$ (the slice $y=x$) the bars are
  $[0,2)$ and $[0,1)$; pairing them costs 1 $= d$.
- **`far_merge`:** both infinite bars start at the push of $(0,0)$; the
  finite bar goes to the diagonal at cost $51/2 = d$ at the vertex $(1,0)$.
  The extra grades of `a0_anchors` and `a0_anchors_medium` are empty
  relations, which change $A_0$ but not the barcodes.
- **Merges:** at the vertex $(0,0)$ of the flat strip (the horizontal slice
  $y=0$) the horizontal bar has length 0 and the vertical bar is $[0,50)$,
  so $F = 25 = d$.
- **`free`, `prism`:** different numbers of infinite bars, so $F=\infty$
  on every slice.
- **`ell` vs `span`:** evaluated by hand at all 7 flat and 8 steep
  vertices; the largest value is 1, at $(1,0)$ in both strips.

The two `r_bif_v_10` values were recorded from play. They equal the known
distances, which are the largest values possible.

`ell.scc` is two intervals forming an L: $(1,0)\to(1,2)$ and $(0,2)\to(3,2)$.
`span.scc` is one interval $(0,0)\to(3,2)$. The flat strip has 7 vertices; on
$s=1$ they sit at $t\in\{-1,0,1,2\}$. Between $t=0$ and $t=1$, $F$
rises to $3/2$ (span's bar has weighted length 3 and goes to the diagonal;
cross-matches with the L bars are no cheaper). That point $(s,t)=(1,1/2)$
is not a vertex of $A_0$, and $F$ at the $A_0$ vertices only reaches 1.
So evaluating $F$ only at $A_0$ vertices is not enough: the exact algorithm
needs the vertices of the refined arrangement $A$.

`far_merge` is two generators at the origin killed together at $(51,0)$: one
finite bar of length 51 and one infinite bar. `a0_anchors` is a single
infinite bar (plus unused relation grades that only cut $A_0$). The
infinite bars match at cost 0; the finite bar goes to the diagonal at 25.5.

Dummy relation columns in the `*_medium.scc` files are empty: they do not
change the slice barcodes, only the arrangement.

## Medium set

`./bench.py medium` runs play on a larger set of pairs and prints a TSV
summary (`name`, `known_dist`, `a0_lines`, `vertices`, `edges`, `faces`,
`max_F`, `max_F_decimal`, `seconds`, `status`, `hera`, `hera_seconds`). A case passes when
play exits cleanly and max $F$ on the $A_0$ vertices passes its check
against the known distance $d$:

- **Synthetic pairs** (`far_merge*`, `merge_*`): max $F = d$. The arguments
  are the ones of the corresponding small cases; the extra short bars of the
  `.bif` versions used with `--aida` cost at most 1/2.
- **`v10` pairs:** max $F \le d$, since $F \le d$ on every slice. Nothing
  more is known by hand.

The script exits with 1 if a case fails. Hera's `matching_dist` value (on
the matching `.bif` files; `--hera PATH` or `HERA_BIN` if set, otherwise
auto-detected) is a reference column only. Play takes well under a second
per case here without AIDA and up to about 3s with it; Hera is the slow
part (`--hera-timeout`, default 30s).

```bash
./bench.py medium                     # --no-aida
./bench.py medium --aida              # mpfree+AIDA
./bench.py medium --no-hera           # skip the Hera reference
```

| Pair | Known distance | Check |
| --- | --- | --- |
| `far_merge` vs `a0_anchors_medium` | 25.5 | = |
| `far_merge_41` vs `a0_anchors_medium` | 20.5 | = |
| `far_merge` vs `a0_anchors_grid` (10×10 dummy grades) | 25.5 | = |
| `merge_horizontal_medium` vs `merge_vertical_medium` (L=50) | 25 | = |
| `merge_h40_medium` vs `merge_v40_medium` | 20 | = |
| `merge_h60_medium` vs `merge_v60_medium` | 30 | = |
| v10 n50 samples 1 vs 2 | 3 | ≤ |
| v10 n50 samples 1 vs 5 | 2.5 | ≤ |
| v10 n50 samples 2 vs 8 | 3 | ≤ |
| v10 n50 samples 3 vs 7 | 4.5 | ≤ |
| v10 n50 samples 4 vs 9 | 6 | ≤ |
| v10 n50 samples 7 vs 10 | 3 | ≤ |
| v10 n50 samples 8 vs 10 | 3 | ≤ |

In the current version every `v10` pair reaches its known distance.

`--aida` drops `--no-aida` so play runs mpfree then AIDA. Dummy empty
relation columns in the `*_medium.scc` files disappear under AIDA, which
would shrink $A_0$ to a few lines. So `--aida` switches the synthetics to
the matching `.bif` files (`far_merge.bif`, `far_merge_41.bif`,
`a0_anchors_medium.bif`, `a0_anchors_grid.bif`, `merge_*medium.bif`). Those
1-skeleta put length-1 extra $H_0$ bars at grades around $x,y\ge 200$, far
from the length-40–60 features. The bars survive mpfree+AIDA, keep $A_0$
large, and cost 0.5 to unmatch, so the constructed distances stay the same.
The `--aida` grid stand-in is an 8×8 of real bars.

Extra $H_0$ bars are a **star into the origin**. Isolated vertex-edge
pairs would each add an extra infinite bar and send $d_\mathrm{match}$ to
infinity. With x-axis stars all extra duals pass through $(s,t)=(0,0)$,
so the flat strip has many segments but few vertices (144 segments, 144
vertices for `far_merge.bif` vs `a0_anchors_medium.bif`; the wall $s=0$ is a
single point and is left out).

### Unused AIDA files

`aida_xanc_{33,34,35}.bif`, `aida_yanc_{21,22}.bif`,
`aida_mh8_{40,50,60,70}.bif`, `aida_mv8_{40,50,60,70}.bif`, and
`far_merge_{61,71,81}.bif` were built to time the removed refinement step
(25 pairs at 30–60s each). No script uses them now. They are still valid
inputs with known distances: `far_merge_L` vs `aida_xanc_K` or
`aida_yanc_K`, and `aida_mh8_L` vs `aida_mv8_L`, all have
$d_\mathrm{match}=L/2$.

## Finite random pairs

`tests/results.txt` (dim 0, old `--no-refine`) has 18 finite-distance
pairs among the first 1577 unordered `random_bifiltrations` comparisons;
the rest returned infinity in 1–3s. `./bench.py finite` builds $A_0$ for
those 18 pairs and runs Hera on each as a reference. The pair list also
stores Hera's distance for each pair, which is printed next to the pair.

```bash
./bench.py finite                    # AIDA (play default)
./bench.py finite --no-aida
./bench.py finite --no-hera          # A0 only
```

The old `--no-refine` distances in `results.txt` came from an early, buggy
version of $F$; several are far above Hera's value (for example 85 vs
34.5).

## All random pairs

`./bench.py random` builds $A_0$ for every unordered pair of
`tests/random_bifiltrations/*.bif` (77,815 pairs) and runs Hera on each as
a reference. The old `tests/results.txt` is not touched. `--no-aida` and
`--no-hera` work as above.

## Results files

`finite`, `random`, and `kerber` write tab-separated files:

| Set | Default file |
| --- | --- |
| `finite` | `tests/finite_results_a0_<aida\|noaida>.txt` |
| `random` | `tests/results_a0.txt` |
| `kerber` | `tests/small_instance_results_a0.txt` |

Columns: `file_a`, `file_b`, `a0_lines`, `vertices`, `edges`, `faces`,
`max_F` (exact, as play prints it), `max_F_decimal` (play's decimal; the
exact fractions get long when the input grades are decimals, as in
Kerber's files), `evaluations` (number of $F$ evaluations), `seconds`, and
for `finite` and `random` also `hera_distance` and `hera_seconds`. These
three columns are empty with `--a0-only`. The method (all vertices, or
`--incremental` with its samples and seed) is part of the settings, so a
results file never mixes the two.

`--results PATH` writes elsewhere. The `#` lines at the top record the
settings, including whether `max_F` was evaluated. Each pair's row is appended as soon as the pair is done, so
Ctrl+C keeps every completed row. If Ctrl+C hits while play runs after
Hera finished, the row has Hera's value and empty $A_0$ columns.

A second run with the same settings skips the pairs that are already done
(`--fresh` starts over; different settings are an error). A pair skipped
with Ctrl+N gets `???` in its $A_0$ columns and is not retried. `error` or
`timeout` (with `--timeout S`) in `a0_lines` means play failed.

## Kerber `small_instances`

`data_benchmarks/small_instances` is small for mpfree/AIDA, not for $A_0$.
A file is feasible when it has at most 200 distinct bigrades (override with
`--max-grades`). That keeps a pair under roughly 330 $A_0$ lines. The
5k-point `cube2_density` files have about 14k grades each and are excluded.
Pairs are only within the same suffix (`.scc_k`, `.minpres_k`, or plain
`.scc`), because the trailing index is reflected homology.

```bash
./bench.py kerber                   # AIDA, max 200 grades
./bench.py kerber --no-aida
./bench.py kerber --max-grades 400
```

The results file has no Hera columns: Hera reads only `.bif` files. The old
distance logs `tests/small_instance_results.txt` and
`tests/small_instance_results_refine.txt` are not touched.
