# TVM

[![Tests](https://github.com/ZhangTao-SJTU/tvm/actions/workflows/rust.yml/badge.svg)](https://github.com/ZhangTao-SJTU/tvm/actions/workflows/rust.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE.md)

TVM is a three-dimensional vertex model for simulating the mechanics and rearrangements of confluent cellular tissues. Cells are represented as polyhedra whose vertices evolve under volume and surface-area forces. I–H and H–I reconnections allow the network to change its neighbor relationships.

The simulation engine is written in safe Rust and uses only the Rust standard library. It supports a rectangular simulation box with periodic boundaries in all three directions. Python tools are included for generating initial Voronoi networks and processing simulation output.

## Quick start

Download the Linux or Windows executable package from [tvm-1.0](https://github.com/ZhangTao-SJTU/tvm/releases/tag/tvm-1.0), or build from source using the steps below. Executable packages include a ready-to-run example and platform-specific instructions.

Install a stable Rust toolchain with Cargo and a linker for your platform. On Windows, the MSVC toolchain uses the Visual Studio C++ build tools; the GNU toolchain uses a compatible MinGW linker.

```sh
git clone https://github.com/ZhangTao-SJTU/tvm.git
cd tvm
cargo build --release --locked
```

The executable is `target/release/tvm` on Linux and `target/release/tvm.exe` on Windows. A ready-to-run, 64-cell example is included; running it requires no Python packages.

**Linux or macOS shell:**

```sh
mkdir -p runs/example
cp tests/fixtures/periodic64/conf tests/fixtures/periodic64/sample.topo runs/example/
cargo run --release --locked -- --dir runs/example
```

**Windows PowerShell:**

```powershell
New-Item -ItemType Directory -Path runs/example
Copy-Item tests/fixtures/periodic64/conf, tests/fixtures/periodic64/sample.topo runs/example
cargo run --release --locked -- --dir runs/example
```

TVM reads `conf` and `sample.topo` from the directory selected by `--dir` and writes results into that directory. With no arguments, it uses the current working directory. Run `tvm --help` for the command-line syntax.

Use a fresh directory for each simulation. TVM refuses to overwrite existing result files, so choose another directory name when repeating the example.

## Physical model

The dimensionless energy is

$$
E = k_v\sum_c (V_c - 1)^2 + \sum_c (A_c - s_0)^2,
$$

where $V_c$ and $A_c$ are the volume and surface area of cell $c$. The target cell volume is one, $s_0$ is the target surface area, and $k_v$ sets the relative volume stiffness. The reported cell shape index is $A_c/V_c^{2/3}$.

Forces are the negative energy gradients and include the derivatives of the perimeter-weighted face centers. Overdamped dynamics uses Euler integration at zero temperature and Euler–Maruyama integration at finite temperature, with mobility and Boltzmann's constant set to one. The mean deterministic vertex velocity is removed before each position update.

Short edges and triangular faces are candidates for I–H and H–I reconnections. Each proposed change is checked for valid geometry and topology before being applied. The periodic bulk network has four edges and four cells incident to each vertex, three faces incident to each edge, and two cells incident to each face. Faces and cells must form connected, closed structures with consistent incidence.

## Configuration

The `conf` file contains one parameter per line. Blank lines and comments beginning with `#` are allowed.

```text
time 0 0.1 0.001
dtr 0.01
dump vtk 0.1
log 0.01
s0 5.6
Lth 0.02
T 0
kv 10
box 4 4 4 p p p
seed 20261028
validate_topology on
```

| Parameter | Meaning |
| --- | --- |
| `time start end dt` | Initial time, final time, and integration step. The last step is shortened when needed to reach `end` exactly. |
| `dtr interval` | Interval between reconnection sweeps. Optional; defaults to `10 * dt`. |
| `dump vtk interval` | Save network and cell measurements, plus VTK geometry, at this interval. Use `dump topo interval` for text outputs without VTK. |
| `log interval` | Interval for console output of time, volume, energies, and reconnection counts. |
| `s0 value` | Target cell surface area. |
| `Lth length` | Reconnection length threshold. Zero disables reconnection. Add `verbose` to save detailed event diagnostics. |
| `T value` | Temperature controlling the amplitude of thermal noise. Zero gives deterministic dynamics. |
| `kv value` | Volume-energy coefficient $k_v$. |
| `box Lx Ly Lz p p p` | Positive box lengths and periodic boundaries in all three directions. |
| `seed integer` | Unsigned 64-bit random seed. Optional; defaults to `5489`. |
| `validate_topology on` | Enable topology checks after reconnection sweeps. Optional; defaults to `on`. Input validation always runs. |

All entries except `dtr`, `seed`, and `validate_topology` are required. The time step must be positive; the final time must be at least the initial time. Output and reconnection intervals must be at least `dt`. The parameters `s0`, `Lth`, `T`, and `kv` must be nonnegative. Unknown keys, duplicate entries, and nonfinite numbers are rejected.

## Initial topology

`sample.topo` describes the network in four sections, in this order:

| Section header | Each following row |
| --- | --- |
| `vertices N` | Vertex ID followed by three coordinates: `id x y z`. |
| `edges N` | Edge ID followed by two endpoint IDs: `id v0 v1`. |
| `polygons N` | Face ID followed by its edge IDs. |
| `cells N` | Cell ID followed by its face IDs. |

Each header specifies the number of rows in that section. IDs are nonnegative and need not be contiguous. TVM checks references, closed face cycles, cell shells, periodic bulk incidence, and cell volumes when loading the network.

To generate a new periodic Voronoi network, use `scripts/tvm/main.py` in a Python environment with `numpy` and `pyvoro`. For example, after creating a fresh directory and placing a matching `conf` file there:

```sh
python scripts/tvm/main.py --seed 5489 --box 4 4 4 --cells 64 --output runs/custom/sample.topo
```

The output directory must already exist. The generator writes `sample.topo` and a companion JSON file with its settings. Its box lengths must match the `box` entry in `conf`. The generation seed controls the initial network; the `seed` in `conf` controls thermal noise during simulation. With unit target cell volume, choosing the number of cells equal to the box volume gives a mean cell volume of one.

## Simulation outputs

| File | Contents |
| --- | --- |
| `topo.txt` | Network snapshots with simulation times, vertex coordinates, and connectivity. |
| `cellVolume.txt` | Per-cell volumes and total volume. |
| `cellCenter.txt` | Per-cell center coordinates. |
| `cellShapeIndex.txt` | Per-cell shape indices and their mean. |
| `*.sample.vtk` | Cell geometry for visualization in ParaView when `dump vtk` is selected. |
| `reconnections.txt` | Reconnection events and their simulation times. |
| `IH_*` / `HI_*` diagnostic files | Before-and-after geometry for reconnection events when `Lth ... verbose` is enabled. |
| `run_metadata.txt` | Model settings, random seed, and initial network size. |
| `run_summary.txt` | Successful completion, final time, energy, total volume, and rejected-event counts. |
| `rng_state.txt` | Final random-generator state. |

Negative start times can be used for relaxation: intermediate snapshots before zero are omitted, while the final state is always saved.

## Reproducibility and testing

Thermal noise uses MT19937-64 with Marsaglia polar Gaussian sampling. Preserve the configuration, input topology, executable, and execution environment when reproducing a trajectory. Identical seeds reproduce the random sequence within the same implementation; floating-point results can differ across platforms and math libraries. The random-state file alone is not a complete simulation restart checkpoint.

Run the automated tests with:

```sh
cargo test --release --locked
```

The tests cover analytical force gradients, periodic geometry, topology constraints, reference numerical values, I–H/H–I round trips, invalid-input handling, time scheduling, stochastic reproducibility, VTK output, and protection of existing result files. See [Testing](docs/TESTING.md) for details. GitHub Actions builds and tests TVM on Linux and Windows.

For a scientific application, check time-step convergence, reconnection-parameter sensitivity, and simulation duration for the quantities being measured.

## Code organization

| Module | Responsibility |
| --- | --- |
| `model` | Vertices, edges, faces, cells, periodic geometry, and topology checks. |
| `mechanics` | Volume and surface-area energies and analytical forces. |
| `reconnection` | I–H/H–I transformations and event diagnostics. |
| `io` | Configuration and topology parsing, topology and VTK serialization. |
| `simulation` | Time integration, scheduling, measurements, and output. |
| `rng` | Random-number and Gaussian sampling. |

## Authors and reference

- **Tao Zhang**, Shanghai Jiao Tong University — zhangtao.scholar@sjtu.edu.cn
- **Jennifer Schwarz**, Syracuse University — jschwarz@physics.syr.edu

The three-dimensional vertex-model framework follows Okuda et al. (2013), [DOI: 10.1007/s10237-012-0430-7](https://doi.org/10.1007/s10237-012-0430-7).

## License

TVM is distributed under the **GNU General Public License, version 3 or later**. See [LICENSE.md](LICENSE.md) for the full license.
