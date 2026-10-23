# Sample presentations

`interval_long.scc` and `interval_short.scc` are two interval modules with
matching distance **1**.

```bash
./play tests/interval_long.scc tests/interval_short.scc
./play tests/interval_long.bif tests/interval_short.bif
```

The `.bif` pair is the same comparison as a 1-critical bifiltration: two
vertices at (0, 0) and one edge at the death grade. The \(H_0\) presentation
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

`run_small.sh` runs tiny presentations that hit the same maximizer pitfalls
as the slow random `.bif` pairs, each in well under a minute:

```bash
./tests/run_small.sh
PLAY_BIN=./build/play ./tests/run_small.sh
```

Every case uses `--no-aida`. `--no-refine` is the path that used to
overestimate; refine is also checked on the smallest pairs.

| Pair | Expected | What it covers |
| --- | --- | --- |
| `interval_long` vs `interval_short` (`.scc` and `.bif`) | 1 | Rank-1 intervals, refine and not |
| `far_merge` vs `a0_anchors` | 25.5 | Unbounded \(A_0\) faces with extra dual lines at \(x,y\in\{10,20,30\}\). A sample that leaves the face mixed those partitions with slope \(r=-30\) and reported 30 (or 51) instead of the diagonal cost of the length-51 bar |
| `far_merge.bif` vs `a0_anchors.bif` | 25.5 | Same distance as a 1-critical \(H_0\) bifiltration |
| `merge_horizontal` vs `merge_vertical` | 25 | Type-III slices. The Lemma-10 closed form \(c\|r+q_x\|\) on an outer \(A_0\) ray is 50; the true \(s\to+\infty\) bottleneck is 25 |
| `free_one` vs `free_two` | infinity | Unmatched infinite bars, \(s>0\) |
| `prism_1` vs `prism_2` `--dim 2 --no-aida` | infinity | Raw \(H_2\) presentations of the two prism 2-skeleta (8 vs 9 triangles, no 3-simplices) |

`far_merge` is two generators at the origin killed together at \((51,0)\): one
finite bar of length 51 and one infinite bar. `a0_anchors` is a single
infinite bar (plus unused relation grades that only cut \(A_0\)). The
infinite bars match at cost 0; the finite bar goes to the diagonal at 25.5.
