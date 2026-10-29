# Testing TVM

Run the complete test suite from the repository root:

```sh
cargo test --release --locked
```

The tests use the bundled periodic 64-cell configuration and small constructed meshes. No Python packages are needed for the Rust test suite. GitHub Actions runs the same tests and builds the executable on Linux and Windows.

## Coverage

| Area | Checks |
| --- | --- |
| Mechanics | Analytical forces against finite differences of the energy; planar, nonplanar, wrapped and large-ID meshes; translation invariance |
| Geometry | Periodic distance conventions, positive cell volumes, degenerate faces and inverted cells |
| Topology | Connected polygon cycles and cell shells, incidence counts, duplicate edges and Euler relations |
| Reconnection | I-H-I transitions, periodic translations, threshold lengths, rejected proposals, ID limits and diagnostic output |
| Runtime | Input validation, shortened final steps, negative and large starting times, unstable geometry, CLI and output preservation |
| Random numbers | MT19937-64 reference sequence, cached Gaussian samples, seed reproducibility and zero-temperature behavior |
| Numerical references | All 445 initial vertex forces, 64 cell volumes, both energy terms, total volume and Gaussian samples |

Individual groups can be run with `cargo test --release --test mechanics`, `--test reconnection`, `--test runtime`, or `--test reference`.

## Test data

`tests/fixtures/periodic64/sample.topo` describes 64 cells in a periodic box of size `4 × 4 × 4`. `conf` provides a short simulation configuration. The accompanying `provenance.json` records the geometry generator, seed, dimensions and file checksum. `reference.tsv` stores numerical values for the reference test at `kv = 10`, `s0 = 5.6` and random seed `20261028`.

Tests create isolated temporary output directories. The fixture files are never modified by a simulation test. Runtime tests check that existing outputs are preserved, that successful runs reach their requested final time, and that invalid states return an error without writing a completion summary.

## Assessing a simulation

Enable `validate_topology on` when preparing a new configuration. Inspect the final time, energy, total volume and rejected-event counts in `run_summary.txt`; use VTK output to examine cell geometry. Choose a timestep and reconnection interval appropriate to the mesh, and check convergence for the quantities being measured.

The test suite covers defined meshes and numerical cases. Long-time statistics, arbitrary-mesh nonintersection and application-specific convergence need to be assessed for each study. Reproducible stochastic trajectories require the same executable, input configuration, topology and seed; bitwise equality across platforms is not guaranteed.
