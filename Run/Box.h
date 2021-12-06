// Copyright June 2021 Tao Zhang @ Shanghai Jiao Tong University.  All Rights Reserved.
// Author: Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn
// Corresponding author: Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu

#ifndef BOX_H_INCLUDED
#define BOX_H_INCLUDED

class Box;
#include <vector>
#include "Run.h"

class Box {
  public:
    Run * run_;
    double size_[3];
    bool boundaryCondition_[3];

    Box(Run * run);
    int     resetPosition(double *);
    int     resetDistance(double *);
};

#endif
