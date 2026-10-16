# Exact matching distance

Computes the matching distance between two 2-parameter persistence modules from
finite presentations, following Kerber, Lesnick, and Oudot,
[*Exact Computation of the Matching Distance on 2-Parameter Persistence Modules*](https://arxiv.org/abs/1812.09085).

The program builds the coarse slice-space arrangement $A_0$ (one line per
input grade, plus $s=0$ and $s=1$), evaluates the weighted bottleneck
distance $F$ at every vertex (Section 5), and takes the type-III limits along
unbounded outer regions. Comparison lines that refine $A_0$ to $A$ are
handled per face: only the $O(n^2)$ lines of the current realizing pair are
intersected, so the huge global arrangement of $O(n^4)$ lines is never stored.
The largest of those values is the matching distance.

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

## Build

```bash
mkdir -p build
cd build
cmake ..
cmake --build .
```

The executable is `build/play`.

## Input

Each module is a **2-parameter presentation** in [scc2020](https://bitbucket.org/mkerber/chain_complex_format)
format. A presentation block looks like:

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

Inputs must be presentations (one matrix: relations vs generators), not raw
chain complexes with more than two non-zero degrees. Minimal presentations
from mpfree are the expected input.

The matching distance is taken between the **first two** presentations that
are loaded.

## Usage

Two presentations in one file (two `scc2020` blocks):

```bash
./play modules.scc
```

Two files, one presentation each:

```bash
./play module_a.scc module_b.scc
```

### Flags

| Flag | Effect |
| --- | --- |
| `--no-cgal` | Only count generator/relation joins. Do not build the arrangement or compute the distance. |
| `--no-refine` | Maximize $F$ only on vertices and outer cells of $A_0$. Skip per-face comparison-line intersections. Faster, and can miss the exact maximum. |

Flags may appear before or after the file names.

```bash
./play a.scc b.scc --no-refine
./play --no-cgal modules.scc
```

## Output

Typical run:

```plain
In total, found 2 generator joins and 3 relation joins
Sum: 5
Loaded 2 presentation(s)
Arrangement has 18 vertices, 27 edges, and 11 faces
Matching distance: 1.5
```

- Join counts come from incomparable generator grades and incomparable
  relation grades of every loaded presentation.
- The arrangement statistics are for the slice-space subdivision used in the
  maximization.
- `Matching distance` is the computed value of $d_{\mathrm{match}}$, printed
  as a double.

If fewer than two presentations are loaded, the program exits after the join
report (or with an error if an arrangement was requested).

## Complexity

The coarse arrangement $A_0$ has $O(n)$ lines and $O(n^2)$ vertices, where
$n$ is the number of distinct grades. Per-face refinement intersects
$O(n^2)$ lines of one realizing pair at a time (never the global $O(n^4)$
arrangement $A$). Each candidate point still requires two slice persistence
computations, so large presentations remain slow. Use `--no-refine` to skip
those extra intersection points.
