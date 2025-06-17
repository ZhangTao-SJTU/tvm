# TVM [![Build Status][1]][2] [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

[1]: https://travis-ci.com/ZhangTao-SJTU/tvm.svg?token=YPqm5yHsQT7PD3VM6WG5&branch=main
[2]: https://travis-ci.com/ZhangTao-SJTU/tvm

## Authors
- Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn
- Shabeeb Ameen @ Syracuse University, mameen@syr.edu
- Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu

## Description
This code simulates tissue dynamics based on a 3D vertex model, inspired by the work of Okuda et al. (2013) [https://link.springer.com/article/10.1007/s10237-012-0430-7](https://link.springer.com/article/10.1007/s10237-012-0430-7). The simulation focuses on cellular structures and their interactions, including energy calculations, topological changes (reconnections), and interactions with an Extracellular Matrix (ECM).

The network topology adheres to the following conditions:
1. Two edges never share two vertices simultaneously.
2. Two polygonal faces never share two or more edges simultaneously.
3. **EXTRA RULE** Two polyhedral cells never share two or more polygonal faces simultaneously.

The codebase is organized into several modules:
- `Cell`: Manages individual cell properties (e.g., volume) and their constituent polygons.
- `Edge`: Defines the edges connecting vertices, including their length and type.
- `Polygon`: Represents the polygonal faces of cells, managing their area, tension, and constituent edges/vertices.
- `Vertex`: Represents the points where edges meet, handling forces and movement.
- `Run`: Orchestrates the simulation, loading configurations, managing time steps, and coordinating updates across components like vertices, edges, polygons, and cells. It also handles dumping simulation data.
- `Energy`: Calculates various energy contributions, including volume constraints (`Volume.cpp`, `Volume.h`), interface tensions (`Interface.cpp`, `Interface.h`), and ECM interactions (e.g., fiber elasticity in `FiberElasticity.cpp`, `FiberElasticity.h`, and fiber-link interactions in `FiberLink.cpp`, `FiberLink.h`).
- `Reconnection`: Handles topological changes in the cellular network, such as T1 transitions (I_H and H_I processes detailed in `Reconnection.cpp` and `Reconnection.h`).
- `ECM`: Models the Extracellular Matrix, including nodes (`Node.cpp`, `Node.h`), fibers (`Fiber.cpp`, `Fiber.h`), and links (`Link.cpp`, `Link.h`), and their mechanical properties.

Initial configurations can be generated using scripts in the `scripts/tvm/` folder (e.g., `main.py`), which may utilize libraries like `voro++` for Voronoi tessellation to create `sample.topo`.

The current version of the code is primarily designed for bulk systems with periodic boundary conditions, managed by the `Box` component within the `Run` module.

Should you encounter any warning prompts, please reach out to the authors as we work towards making the code more robust for more general geometries, deformations, and energy functionals.

## Quick Start
#### Compiling
The project uses CMake for building.
```bash
mkdir build
cd build
cmake ../
make
```

#### Usage
1.  **Generate Initial Configuration (Optional):**
    Use scripts in the `scripts/tvm/` directory (e.g., `main.py`) to generate an initial `sample.topo` file if needed. These scripts often use `voro++` for generating Voronoi-based initial cell structures.

2.  **Create Configuration File:**
    Assuming the working directory is `build`, create a configuration file named `conf`.
    ```bash
    touch conf
    # Edit conf with a text editor, e.g., vim conf or notepad conf
    ```
    The `conf` file defines simulation parameters. An example (`sample_conf` in the root directory) is shown below:
    ```
    time 0 25000 0.005       # start_time end_time timestep
    dump vtk 500.0          # dump_format dump_period (vtk for Paraview)
    log 100.0               # logging_period_to_screen
    s0 5.2 0.25             # target_cell_area parameters (specific to energy model)
    Lth 0.02                # threshold_length_for_reconnection
    T 1e-4                  # temperature_or_noise_level
    kv 10.                  # volume_constraint_stiffness
    box 64. 64. 64. p p p   # box_dimensions (Lx Ly Lz) and boundary_conditions (p for periodic)
    fiber 10.0 2.8284271 0.001 # ECM fiber parameters (if ECM is used)
    link 100 10. 1.5 0.2 0.004 1000. # ECM link parameters (if ECM is used)
    ```
    Parameter details:
    - `time`: Specifies start time, end time, and timestep for the simulation.
    - `dump`: Defines the output file format (e.g., `.vtk` for visualization with ParaView) and the time interval for writing output files.
    - `log`: Sets the time interval for logging simulation information to the console.
    - `s0`: Specifies parameters related to the target cell area or interface energy.
    - `Lth`: Defines the threshold length for edge reconnection events.
    - `T`: Represents a temperature or noise level in the system.
    - `kv`: Stiffness constant for the cell volume constraint.
    - `box`: Defines the simulation box dimensions (Length_x, Length_y, Length_z) and boundary conditions for each axis ('p' for periodic).
    - `fiber`, `link`: Parameters for the Extracellular Matrix (ECM) components, if modeled.

3.  **Run Simulation:**
    Execute the compiled program from the `build` directory.
    ```bash
    ./tvm
    ```
    The simulation will read `sample.topo` (for initial structure) and `conf` (for parameters).

## Output
The simulation typically outputs `.vtk` files that can be visualized using software like [ParaView](https://www.paraview.org/).

## License

[GNU GPL v3 License](./LICENSE.md)

Copyright 2021-2023 Tao Zhang, Shabeeb Ameen
