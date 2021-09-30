// Copyright June 2021 Tao Zhang @ Shanghai Jiao Tong University.  All Rights Reserved.
// Author: Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn
// Corresponding author: Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu

#include <cstdio>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <unordered_map>
#include <random>
#include "Box.h"

using namespace std;

Box::Box(Run * run) {
    run_ = run;
}

int     Box::resetPosition(double * r) {
    if (fabs(r[0]) > 1e6) {
        printf("%e\n",r[0]);
    }
    if (fabs(r[1]) > 1e6) {
        printf("%e\n",r[1]);
    }
    if (fabs(r[2]) > 1e6) {
        printf("%e\n",r[2]);
    }
    for (int i = 0; i < 3; i++) {
        if (boundaryCondition_[i]) {
            r[i] = r[i] - size_[i] * floor(r[i] / size_[i]);
        }
    }

    return 0;
}

int     Box::resetDistance(double * dx) {
    for (int i = 0; i < 3; i++) {
        if (boundaryCondition_[i]) {
            while (dx[i] > size_[i] / 2.0) {
                dx[i] = dx[i] - size_[i];
            }
            while (dx[i] < (-1.0) * size_[i] / 2.0) {
                dx[i] = dx[i] + size_[i];
            }
        }
    }

    return 0;
}