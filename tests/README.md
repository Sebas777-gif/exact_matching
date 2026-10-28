# Sample presentations

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

`run_small.sh` runs presentations that hit the same maximizer pitfalls
as the slow random `.bif` pairs:

```bash
./tests/run_small.sh                 # incremental refine (default)
./tests/run_small.sh --no-refine     # A0 sample only
./tests/run_small.sh --fast          # tiny cases only
PLAY_BIN=./build/play ./tests/run_small.sh
```

Every case uses `--no-aida`. Default is the incremental $A\setminus A_0$
path. `--no-refine` is the $A_0$ sample that used to overestimate and can
still miss a comparison vertex. The ell vs span miss is always checked.

| Pair | Expected | What it covers |
| --- | --- | --- |
| `interval_long` vs `interval_short` (`.scc` and `.bif`) | 1 | Rank-1 intervals, refine and not |
| `far_merge` vs `a0_anchors` | 25.5 | Unbounded $A_0$ faces with extra dual lines at $x,y\in\{10,20,30\}$. A sample that leaves the face mixed those partitions with slope $r=-30$ and reported 30 (or 51) instead of the diagonal cost of the length-51 bar |
| `far_merge.bif` vs `a0_anchors.bif` | 25.5 | Same distance as a 1-critical $H_0$ bifiltration |
| `merge_horizontal` vs `merge_vertical` | 25 | Type-III slices. The Lemma-10 closed form $c\|r+q_x\|$ on an outer $A_0$ ray is 50; the true $s\to+\infty$ bottleneck is 25 |
| `free_one` vs `free_two` | infinity | Unmatched infinite bars, $s>0$ |
| `prism_1` vs `prism_2` `--dim 2 --no-aida` | infinity | Raw $H_2$ presentations of the two prism 2-skeleta (8 vs 9 triangles, no 3-simplices) |
| `ell` vs `span` refine | 1.5 | Comparison vertex on the $A_0$ edge $s=1$, $t=1/2$. `--no-refine` reports ~1.057 (open-face sample); $A_0$ vertices only reach 1 |

`ell.scc` is two intervals forming an L: $(1,0)\to(1,2)$ and $(0,2)\to(3,2)$.
`span.scc` is one interval $(0,0)\to(3,2)$. $A_0$ has nine vertices; on
$s=1$ they sit at $t\in\{-1,0,1,2\}$. Between $t=0$ and $t=1$, $F$
rises to $3/2$ (span's bar has weighted length 3 and goes to the diagonal;
cross-matches with the L bars are no cheaper). That point is a vertex of
$A$ (comparison line $\cap$ the $A_0$ line $s=1$) and is not a
vertex of $A_0$, so `--no-refine` misses it.

```bash
./play --no-aida --no-refine tests/ell.scc tests/span.scc   # ~1.057
./play --no-aida tests/ell.scc tests/span.scc               # 1.5
```

`far_merge` is two generators at the origin killed together at $(51,0)$: one
finite bar of length 51 and one infinite bar. `a0_anchors` is a single
infinite bar (plus unused relation grades that only cut $A_0$). The
infinite bars match at cost 0; the finite bar goes to the diagonal at 25.5.

### Medium cases (about 10–30s each)

These keep the same expected distances, but enlarge $A_0$ (or use a real
10-vertex 1-skeleton). With `--no-refine` they run in about 10–30s; default
incremental refine can take longer. Skip them with `--fast`.

| Pair | Expected | Typical time | What it covers |
| --- | --- | --- | --- |
| `far_merge` vs `a0_anchors_medium` (140 unused grades) | 25.5 | ~19s | Same 25.5 vs 30 trap, with enough dual lines for unbounded faces to be thin |
| `merge_horizontal_medium` vs `merge_vertical_medium` | 25 | ~19s | Type-III pair with ~140 distinct grades |
| `r_bif_v_10_…_sample_1` vs `sample_2` (`--dim 0`) | 3 | ~16s | Real $H_0$ 1-skeleton, same pipeline as the slow random `.bif` pairs |
| `r_bif_v_10_…_sample_1` vs `sample_5` | 2.5 | ~13s | Same, different 10-vertex sample |

Dummy relation columns in the `*_medium.scc` files are empty: they do not
change the slice barcodes, only the arrangement.

## Medium benchmark suite

`run_medium.sh` is the 10–30s set for timing code tweaks. It includes the
four medium cases above plus more synthetics and 10-vertex bifiltration
pairs, and prints a TSV summary (`name`, `expected`, `got`, `seconds`,
`status`).

```bash
./tests/run_medium.sh                 # --no-aida, incremental refine
./tests/run_medium.sh --no-refine     # --no-aida, A0 sample only (10–30s)
./tests/run_medium.sh --aida          # mpfree+AIDA, incremental refine (~30–60s synthetics)
./tests/run_medium.sh --aida --no-refine
PLAY_BIN=./build/play ./tests/run_medium.sh --aida
TIMEOUT_SEC=90 ./tests/run_medium.sh --aida
```

`--aida` drops `--no-aida` so play runs mpfree then AIDA. Dummy empty
relation columns in the `*_medium.scc` files disappear under AIDA, which
would shrink $A_0$ to a few lines. The matching `*.bif` 1-skeleta put
length-1 extra $H_0$ bars at grades around $x,y\ge 200$ (far from the
length-40–60 features). Those bars survive mpfree+AIDA, keep $A_0$ large,
and cost 0.5 to unmatch, so the constructed distances stay the same.

| Pair | Expected | Typical time |
| --- | --- | --- |
| `far_merge` vs `a0_anchors_medium` | 25.5 | ~19s |
| `far_merge_41` vs `a0_anchors_medium` | 20.5 | ~19s |
| `far_merge` vs `a0_anchors_grid` (10×10 dummy grades) | 25.5 | ~30s |
| `merge_horizontal_medium` vs `merge_vertical_medium` (L=50) | 25 | ~19s |
| `merge_h40_medium` vs `merge_v40_medium` | 20 | ~29s |
| `merge_h60_medium` vs `merge_v60_medium` | 30 | ~28s |
| v10 n50 samples 1 vs 2 | 3 | ~16s |
| v10 n50 samples 1 vs 5 | 2.5 | ~13s |
| v10 n50 samples 2 vs 8 | 3 | ~24s |
| v10 n50 samples 3 vs 7 | 4.5 | ~26s |
| v10 n50 samples 4 vs 9 | 6 | ~25s |
| v10 n50 samples 7 vs 10 | 3 | ~26s |
| v10 n50 samples 8 vs 10 | 3 | ~27s |

Default is `--no-aida` with incremental refine on the `.scc` synthetics
(empty dummy columns) and the v10 `.bif` pairs. `--no-refine` is the
10–30s $A_0$-sample timing band. `--aida --no-refine` switches the
synthetics to the matching `.bif` files (`far_merge.bif`,
`far_merge_41.bif`, `a0_anchors_medium.bif`, `a0_anchors_grid.bif`,
`merge_*medium.bif`) and the same v10 pairs. Times in the table are from
one `--no-aida --no-refine` run; they move with CPU load but should stay
in the same band. AIDA adds preprocess time and may change maximizer time
because the presentation is minimal. The `--aida` grid stand-in is an 8×8
of real bars (the 10×10 dummy grid is too slow once the extra bars are in
the barcode).

`--aida` with incremental refine (the play default) does **not** use those
140-anchor / merge-70 files: they take minutes once comparison lines are
inserted. The suite switches to 25 constructed `.bif` pairs sized so
AIDA+refine lands around 30–60s each (default timeout 90s for slack).
v10 n50 pairs still run; with AIDA they finish in a few seconds.

Extra $H_0$ bars are a **star into the origin** (offset 200, step 3).
Isolated vertex-edge pairs each add an extra infinite bar, which sends
$d_\mathrm{match}$ to infinity and makes `maximize_F` return immediately.
Unmatch cost of a length-1 bar is 0.5, so the constructed distances stay
$L/2$. Changing $L$ keeps runtime and changes the distance.

| Pair | Expected | Typical AIDA+refine time |
| --- | --- | --- |
| `far_merge_L` vs `aida_xanc_K`, $L\in\{41,51,61,71,81\}$, $K\in\{33,34,35\}$ (15 pairs) | $L/2$ | ~30–45s |
| `far_merge_L` vs `aida_yanc_K`, $L\in\{41,51,61\}$, $K\in\{21,22\}$ (6 pairs) | $L/2$ | ~40s |
| `aida_mh8_L` vs `aida_mv8_L`, $L\in\{40,50,60,70\}$ (4 pairs) | $L/2$ | ~50–55s |

`far_merge.bif` is $L=51$; `far_merge_41.bif` already existed;
`far_merge_{61,71,81}.bif` are the same two-vertex merge at $(L,0)$.
x-axis stars produce concurrent duals (fewer vertices, refine on unbounded
wedges); y-axis stars produce more vertices, so a smaller $K$ hits the
same band. Merge $K=8$ extras are y-axis on the horizontal module and
x-axis on the vertical module.
