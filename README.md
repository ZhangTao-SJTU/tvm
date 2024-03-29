# TVM [![Build Status][1]][2] [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

[1]: https://travis-ci.com/ZhangTao-SJTU/tvm.svg?token=YPqm5yHsQT7PD3VM6WG5&branch=main
[2]: https://travis-ci.com/ZhangTao-SJTU/tvm

A 3D vertex model code

Author:

Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn

Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu

## Description
This code is based on Okuda 2013 paper: https://link.springer.com/article/10.1007/s10237-012-0430-7.

- The network topology satisfies the following conditions:
  1. Two edges never share two vertices simultaneously.
  2. Two polygonal faces never share two or more edges simultaneously.
  3. **EXTRA RULE** Two polyhedral cells never share two or more polygonal faces simultaneously

The scripts are in the folder "scripts" on the main branch.

Current version of the code works ONLY for bulk system with periodic boundary condition. 

Should you receive any warning prompts, please reach out to us as we work towards making the code more robust for more general geometries, deformations, and energy functionals. 

## Quick Start
#### Compiling
```bash
mkdir build
cd build
cmake ../
make
```

#### Usage
Assuming the working directory is `build`, first create a configuration file `conf`
```bash
touch conf
vim conf
```
The `conf` file defines values of simulation parameters. See an example as below:
```
time -10000.0 1000000.0 0.005
dump vtk 500.0
log 50.0
s0 5.60
Lth 0.02
T 1.0e-04
kv 10.
box 8. 8. 8. p p p
```
Here `time` specifies start time, end time, timestep;
`dump` specifies the dumped file format (the `.vtk` file format is for visualization using `paraview`), the dumping time period;
`log` specifies the time period logging simulation information to screen;
`s0` specifies the target cell area;
`Lth` specifies the threshold edge length to trigger reconnection;
`T` specifies the temperature;
`kv` specifies the strength of the volume elasticity;
`box` specifies the simulation box sizes and boundary conditions.

Use python script `scripts/tvm/main.py` to generate the initial configuration file `sample.topo`,
which contains the information of vertices' coordinates and the system topology. Please install pyvoro in advance. 
```bash
pip install pyvoro-mmalahe
python ../scripts/tvm/main.py
```

Run `tvm`
```bash
./tvm
```

## License

[GNU GPL v3 License](./LICENSE.md)

Copyright 2021-2023 Tao Zhang
