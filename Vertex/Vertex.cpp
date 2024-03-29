/* ---------------------------------------------------------------------------------
 * Copyright 2021-2023 Tao Zhang
 *
 * This file is part of TVM.
 *
 * TVM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation,
 * either version 3 of the License,
 * or (at your option) any later version.
 *
 * TVM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with TVM. If not, see <https://www.gnu.org/licenses/>.
 *
 * Author: Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn
 * Coauthor: Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu
 * ---------------------------------------------------------------------------------
 */

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
        volumeForce_[i] = 0.;
        interfaceForce_[i] = 0.;
        pullingForce_[i] = 0.;
        velocity_[i] = 0.;
    }
}

int Vertex::logCells(std::string name) {
    printf("%s %ld\n",name.c_str(), id_);
    printf("%ld",cells_.size());
    for (auto cell : cells_) {
        printf(" %ld",cell->id_);
    }
    printf("\n");

    return 0;
}

int Vertex::logEdges(std::string name) {
    printf("%s %ld\n",name.c_str(), id_);
    printf("%ld",edges_.size());
    for (auto edge : edges_) {
        printf(" %ld",edge->id_);
    }
    printf("\n");

    return 0;
}
