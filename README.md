# Exact matching distance

Goal: compute the matching distance between two 2-parameter persistence
modules from finite presentations, following Kerber, Lesnick, and Oudot,
[*Exact Computation of the Matching Distance on 2-Parameter Persistence Modules*](https://arxiv.org/abs/1812.09085).

**Current state:** the program only builds the coarse arrangement $A_0$,
in two strips of the $(s,t)$ parameter plane:

- the **flat strip**: slices $y = sx + t$ with $0 \le s \le 1$;
- the **steep strip**: slices $x = sy + t$ with $0 \le s \le 1$, built the
  same way from the modules with $x$ and $y$ swapped.

Together they cover every slice of positive slope, including the horizontal
($s=0$ in the flat strip) and vertical ($s=0$ in the steep strip) limits,
and every crossing appears in exactly one strip (crossings on $s=1$ in both).
In each strip, $A_0$ has one dual segment $t = v_y - s\,v_x$,
$0 \le s \le 1$, per distinct input grade $v$ of either module, plus the
walls $s=0$ and $s=1$. Each wall runs from the lowest to the highest segment
endpoint on it, and is left out if those coincide. The program does **not**
evaluate $F$, refine $A_0$ to $A$, or print a matching distance yet.

An earlier version also evaluated $F$ on $A_0$ and refined it towards $A$
per face, using only the comparison lines of the current realizing pair.
That refinement was not correct, and the code was removed (see the git
history up to commit `322a0c8`). There is no `--no-refine` flag any more:
building $A_0$ without refinement is the only mode.

## Dependencies

- A C++14 compiler
- [CGAL](https://www.cgal.org/) (arrangement of lines; Core is used if present)
- [PHAT](https://bitbucket.org/phat-code/phat) headers
- [mpp_utils](https://bitbucket.org/mkerber/mpfree) headers (`mpp_utils_mod`)
- [scc2020 parser](https://bitbucket.org/mkerber/chain_complex_format) headers

Default include paths in `CMakeLists.txt` are:

| CMake variable | Default |
| --- | --- |
| `PHAT_INCLUDE_DIR` | `~/phat/include` |
| `MPP_UTILS_INCLUDE_DIR` | `~/mpfree/mpp_utils_mod/include` |
| `SCC_INCLUDE_DIR` | `~/chain_complex_format/include` |
| `CGAL_DIR` | `~/CGAL-4.14` |

Override them on the CMake command line if the libraries live elsewhere, for
example `-DPHAT_INCLUDE_DIR=/path/to/phat/include`. System CGAL is used when
CMake finds it.

The stripped-down code only includes CGAL headers. PHAT, mpp_utils, and the
scc2020 parser are still on the include path in `CMakeLists.txt` but are
not used until slice barcodes come back. The `mpfree` and `aida` binaries
are needed at run time unless `--no-aida` is given.

## Build

```bash
mkdir -p build
cd build
cmake ..
cmake --build .
```

The executable is `build/play`.

## Input

Each module is either a **2-parameter presentation** in [scc2020](https://bitbucket.org/mkerber/chain_complex_format)
format, or a **1-critical bifiltration** in Hera's `.bif` format
(`bifiltration_phat_like`). A presentation block looks like:

```plain
scc2020
2
<n_relations> <n_generators> 0
x y ; i j k
...
x y ;
```

The first numeric line after `2` is the sizes of the generating sets. Each of
the next `n_relations` lines is a relation: bigrade `x y`, then generator
indices (0-based) after `;`. Then `n_generators` lines of generator bigrades.

Example (`ex1.scc`):

```plain
scc2020
2
3 2 0
2.0 2.0; 0 1
3.0 0.0; 0
0.0 3.0; 1
1.0 0.0;
0.0 1.0;
```

`.scc` inputs must be presentations (one matrix: relations vs generators), not
raw chain complexes with more than two non-zero degrees. Minimal presentations
from mpfree are the expected `.scc` input.

Kerber's min-pres datasets use `foo.scc_k` or `foo.scc.minpres_k`. Those are
already presentations; `k` is the chain-complex position, reflected homology
(`.scc_2` is $H_0$, `.scc_1` is $H_1$). Play treats them as `.scc`. Do not
pass that `k` as `--dim` (`--dim` only selects homology when reading a `.bif`).

By default each file is one module and is preprocessed with **mpfree** (for
`.bif`) and **AIDA** (for `.bif` and `.scc`). AIDA's indecomposable summands
become the presentations of that module. $A_0$ is built from the union of
the distinct grades of all summands of both modules. The summands are kept
separately for when slice barcodes come back (the matching distance does not
factor over summands, but each slice barcode does).

A `.sccsum` file is already decomposed and is read as-is. `--no-aida` skips
mpfree and AIDA: a `.bif` is converted in-process to a (generally non-minimal)
presentation of $H_d$ by taking $d$-simplices as generators and
$(d+1)$-simplices as relations.

$A_0$ is built for **two modules**, given as exactly two files (one module
each). The two files may mix `.scc`, `.scc_k`, `.minpres_k`, `.sccsum`, and
`.bif`. A file with several `scc2020` blocks is one module whose summands are
those blocks; it is not split into two modules.

Set `$MPFREE` and `$AIDA` (or `$AIDA_WITHOUT_STATS`) if the binaries are not
on `PATH`. Their input and output files go to a fresh folder
`/tmp/exact_matching_XXXXXX` per input file, which is deleted once the file
has been loaded (also on errors). A process that is killed while mpfree or
AIDA runs can leave that folder behind.

## Usage

Two files, one module each:

```bash
./play module_a.scc module_b.scc
./play a.bif b.bif
./play a.bif b.scc --dim 0
```

### Flags

| Flag | Effect |
| --- | --- |
| `--dim d` | Homology dimension used when reading `.bif` files (default 0). Ignored for presentations (`.scc`, `.scc_k`, `.minpres_k`). |
| `--no-aida` | Do not run mpfree or AIDA. Use the file as a raw presentation (or the in-process $H_d$ presentation of a `.bif`). |

Flags may appear before or after the file names. Any other `--flag`
(including the removed `--no-refine` and `--no-cgal`) is an error.

```bash
./play a.scc b.scc
./play a.bif b.bif --dim 0
./play --no-aida tests/interval_long.scc tests/interval_short.scc
```

### Tests and benchmarks

`bench.py` (Python 3, standard library only) runs play on fixed sets of
pairs. Every set checks or records what play prints now: the number of
$A_0$ lines and the vertex, edge, and face counts. See `tests/README.md`.

| Command | What it does |
| --- | --- |
| `./bench.py small` | Regression: exact $A_0$ counts for small pairs, exit 1 on a mismatch (`--fast` skips the ~140-line cases) |
| `./bench.py medium` | Timing set, `--no-aida` by default, `--aida` for mpfree+AIDA. Prints the known distance and Hera's value as a reference (`--no-hera` skips Hera) |
| `./bench.py finite` | $A_0$ for the 18 finite pairs of the old `tests/results.txt`, AIDA by default |
| `./bench.py kerber` | $A_0$ for feasible pairs of Kerber's `data_benchmarks/small_instances` |
| `./bench.py random` | $A_0$ for every pair of `tests/random_bifiltrations` |

`./bench.py SET --help` lists the options of a set. `finite`, `kerber`, and
`random` append one row per pair to a `*_a0.txt` file under `tests/` and
resume where they stopped. They leave the old distance logs (`results.txt`,
`finite_results_new.txt`, `small_instance_results*.txt`,
`refined_results.txt`) alone. In every set, Ctrl+N skips the running pair
and Ctrl+C stops.

## Output

Typical run:

```plain
$ ./play --no-aida tests/interval_long.scc tests/interval_short.scc
Loaded 2 module(s) with 1 summand(s) and 1 summand(s)
A0 lines: 4
Flat strip:  4 vertices, 5 edges, 3 faces
Steep strip: 4 vertices, 5 edges, 3 faces
Arrangement has 8 vertices, 10 edges, and 6 faces
```

- `Loaded` reports the number of summands per module (AIDA's indecomposables,
  or 1 with `--no-aida`).
- `A0 lines` is the number of segments inserted into the flat strip: one
  dual segment per distinct grade of both modules, plus the walls that are
  not a single point.
- `Flat strip` and `Steep strip` are CGAL's counts for each strip.
  `Arrangement has` is their sum (the line `bench.py` reads). Each strip is
  connected and the counts include its unbounded outer face, so Euler's
  formula $V - E + F = 2$ holds per strip and is a quick sanity check.

In the example the grades are $(0,0)$, $(1,1)$, and $(2,2)$. Their three
dual segments run from $(0,0)$, $(0,1)$, $(0,2)$ and all end at
$(s,t)=(1,0)$. So the wall $s=1$ is a single point and is left out, and
each strip has 4 vertices, the wall $s=0$ split in two plus the three
segments as 5 edges, and two triangles plus the outer face as 3 faces.

Two dual lines meet at the slice through both grades. Within $0 < s < 1$
of the flat strip that slice has slope $s$; a crossing with slope $s > 1$ is
in the steep strip at slope $1/s$, and a crossing with $s < 0$ (incomparable
grades) is in neither. Clipping to the two strips gives about half the
vertices of the old full arrangement of lines (up to commit
`626481d`).

If not exactly two input files are given, the program exits with an error.

## Complexity

The coarse arrangement $A_0$ has $O(n)$ segments and $O(n^2)$ vertices per
strip, where $n$ is the number of distinct grades. Each strip is built with
exact arithmetic in one `CGAL::insert` (a sweep). For two 100-vertex random
bifiltrations (`r_bif_v_100_…_sample_1` vs `sample_2`, `--dim 0`) that is
324 segments and about 23k vertices over both strips in about 1.1s with
AIDA (preprocessing included), or 419 segments and about 43k vertices in
about 3s with `--no-aida`. The full arrangement of lines had 42k and 72k
vertices.
