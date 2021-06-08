// Copyright June 2021 Tao Zhang @ Shanghai Jiao Tong University.  All Rights Reserved.
// Author: Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn
// Corresponding author: Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu

#include <cstdio>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <unordered_map>

#include "Vertex.h"

using namespace std;

Vertex::Vertex(Run * run, long int id) {
    run_ = run;
    id_ = id;
    for (int i = 0; i < 3; i++) {
        position_[i] = 0.;
        velocity_[i] = 0.;
    }
}
