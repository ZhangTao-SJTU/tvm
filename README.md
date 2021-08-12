# TVM [![Build Status][1]][2] [![MIT licensed][3]][4]

[1]: https://travis-ci.com/ZhangTao-SJTU/tvm.svg?token=YPqm5yHsQT7PD3VM6WG5&branch=main
[2]: https://travis-ci.com/ZhangTao-SJTU/tvm
[3]: https://img.shields.io/badge/license-MIT-blue.svg
[4]: LICENSE

A 3D vertex model code 

Author: Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn

Corresponding author: Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu

## Description
This code is based on Okuda 2013 paper: https://link.springer.com/article/10.1007/s10237-012-0430-7. 

- The network topology satisfies the following conditions:
  1. Two edges never share two vertices simultaneously. 
  2. Two polygonal faces never share two or more edges simultaneously. 
  3. **EXTRA RULE** Two polyhedral cells never share two or more polygonal faces simultaneously

The scripts are in the folder "scripts" on the main branch.

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
log 5.0
s0 5.30
Lth 0.020
T 1.0e-05
```
Here `time` specifies start time, end time, timestep;
`dump` specifies the dumped file format (the `.vtk` file format is for visualization using `paraview`), the dumping time period;
`log` specifies the time period logging simulation information to screen;
`s0` specifies the target cell area;
`Lth` specifies the threshold edge length to trigger reconnection;
`T` specifies the temperature. 

Use python script `scripts/tvm/main.py` to generate the initial configuration file `sample.topo`, 
which contains the information of vertices' coordinates and the system topology. 
```bash
python ../scripts/tvm/main.py
```

Run `tvm`
```bash
./tvm
```

## License
[MIT License](./LICENSE)