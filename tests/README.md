# Sample presentations and test scripts

`play` currently builds the coarse arrangement $A_0$ and nothing more. The
scripts below check or record its output: the number of $A_0$ lines and the
vertex, edge, and face counts. The matching distances listed here are the
values each pair was built for. They are kept as reference values for when
$F$ is evaluated again. No script checks them right now.

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

`run_small.sh` builds $A_0$ for each pair and compares the four numbers
(lines, vertices, edges, faces) with the expected ones. Every case uses
`--no-aida`.

```bash
./tests/run_small.sh                 # all cases
./tests/run_small.sh --fast          # tiny cases only
PLAY_BIN=./build/play ./tests/run_small.sh
```

| Pair | $A_0$ (L V E F) | $d_\mathrm{match}$ | Why the pair exists |
| --- | --- | --- | --- |
| `interval_long` vs `interval_short` (`.scc` and `.bif`) | 5 4 15 12 | 1 | Rank-1 intervals. All three duals and $s=1$ meet at $(1,0)$ |
| `far_merge` vs `a0_anchors` | 10 24 62 39 | 25.5 | Unbounded $A_0$ faces cut by extra dual lines at $x,y\in\{10,20,30\}$ |
| `far_merge.bif` vs `a0_anchors.bif` | 5 4 15 12 | 25.5 | Same distance as a 1-critical $H_0$ bifiltration (no extra grades) |
| `merge_horizontal` vs `merge_vertical` | 5 6 18 13 | 25 | Type-III slices; the limit $s\to+\infty$ matters |
| `free_one` vs `free_two` | 4 3 11 9 | infinity | Unmatched infinite bars, $s>0$ |
| `prism_1` vs `prism_2` `--dim 2` | 7 12 34 23 | infinity | Raw $H_2$ presentations of two prism 2-skeleta (8 vs 9 triangles, no 3-simplices) |
| `ell` vs `span` | 7 9 29 21 | 1.5 | Maximum of $F$ on an $A_0$ edge, not at an $A_0$ vertex |
| `far_merge` vs `a0_anchors_medium` | 144 5183 10581 5399 | 25.5 | Same as above with 140 unused grades |
| `merge_horizontal_medium` vs `merge_vertical_medium` | 145 5493 11222 5730 | 25 | Merge pair with ~140 distinct grades |
| `r_bif_v_10_…_sample_1` vs `sample_2` (`--dim 0`) | 74 1966 4221 2256 | 3 | Real $H_0$ 1-skeleton |
| `r_bif_v_10_…_sample_1` vs `sample_5` (`--dim 0`) | 71 1806 3883 2078 | 2.5 | Same, different sample |

The last four rows are skipped with `--fast`. The expected counts were
recorded from the current CGAL build. The line counts follow from the
grades (distinct grades + 2). The small vertex counts can be checked by
hand, as in the interval case.

`ell.scc` is two intervals forming an L: $(1,0)\to(1,2)$ and $(0,2)\to(3,2)$.
`span.scc` is one interval $(0,0)\to(3,2)$. $A_0$ has nine vertices; on
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

## Medium suite

`run_medium.sh` builds $A_0$ for a larger set of pairs and prints a TSV
summary (`name`, `known_dist`, `a0_lines`, `vertices`, `edges`, `faces`,
`seconds`, `status`, `hera`, `hera_seconds`). A case passes when play exits
cleanly and prints the $A_0$ counts. The known distance and Hera's
`matching_dist` value (on the matching `.bif` files, `HERA_BIN` if set,
otherwise auto-detected) are reference columns only. Building $A_0$ takes
well under a second per case here; Hera is the slow part (30s timeout).

```bash
./tests/run_medium.sh                 # --no-aida
./tests/run_medium.sh --aida          # mpfree+AIDA
./tests/run_medium.sh --no-hera       # skip the Hera reference
PLAY_BIN=./build/play ./tests/run_medium.sh --aida
```

| Pair | Known distance |
| --- | --- |
| `far_merge` vs `a0_anchors_medium` | 25.5 |
| `far_merge_41` vs `a0_anchors_medium` | 20.5 |
| `far_merge` vs `a0_anchors_grid` (10×10 dummy grades) | 25.5 |
| `merge_horizontal_medium` vs `merge_vertical_medium` (L=50) | 25 |
| `merge_h40_medium` vs `merge_v40_medium` | 20 |
| `merge_h60_medium` vs `merge_v60_medium` | 30 |
| v10 n50 samples 1 vs 2 | 3 |
| v10 n50 samples 1 vs 5 | 2.5 |
| v10 n50 samples 2 vs 8 | 3 |
| v10 n50 samples 3 vs 7 | 4.5 |
| v10 n50 samples 4 vs 9 | 6 |
| v10 n50 samples 7 vs 10 | 3 |
| v10 n50 samples 8 vs 10 | 3 |

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
so $A_0$ has many lines but few vertices (145 lines, 144 vertices for
`far_merge.bif` vs `a0_anchors_medium.bif`).

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
the rest returned infinity in 1–3s. `finite_benchmarks.sh` builds $A_0$
for those 18 pairs and records Hera's distance next to it:

```bash
./finite_benchmarks.sh               # AIDA (play default)
./finite_benchmarks.sh --no-aida
./finite_benchmarks.sh --no-hera     # A0 only
```

Results go to `tests/finite_results_a0_<aida|noaida>.txt`. The script also
prints the old `--no-refine` distance from `results.txt` for each pair.
Several of those are far above Hera's value (for example 85 vs 34.5). They
came from an early, buggy version of $F$, so use Hera's column as the
reference. Ctrl+C keeps completed rows.

## All random pairs

`./benchmarks.sh` builds $A_0$ for every unordered pair of
`tests/random_bifiltrations/*.bif` (77,815 pairs) and runs Hera on each as
a reference. Results go to `tests/results_a0.txt`. The old
`tests/results.txt` is not touched. `--no-aida` and `--no-hera` work as above.

## Kerber `small_instances`

`data_benchmarks/small_instances` is small for mpfree/AIDA, not for $A_0$.
A file is feasible when it has at most 200 distinct bigrades (override with
`--max-grades`). That keeps a pair under roughly 330 $A_0$ lines. The
5k-point `cube2_density` files have about 14k grades each and are excluded.
Pairs are only within the same suffix (`.scc_k`, `.minpres_k`, or plain
`.scc`), because the trailing index is reflected homology.

```bash
./small_instance_benchmarks.sh              # AIDA, max 200 grades
./small_instance_benchmarks.sh --no-aida
./small_instance_benchmarks.sh --max-grades 400
```

Results go to `tests/small_instance_results_a0.txt` with columns `a0_lines`,
`vertices`, `edges`, `faces`, and `seconds`. The old distance logs
`tests/small_instance_results.txt` and
`tests/small_instance_results_refine.txt` are not touched. Ctrl+N skips the
pair that is running and writes `???` in its measurement columns (that pair
is not retried). Ctrl+C keeps completed rows.
