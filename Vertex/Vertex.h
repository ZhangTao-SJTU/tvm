// Copyright June 2021 Tao Zhang @ Shanghai Jiao Tong University.  All Rights Reserved.
// Author: Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn
// Corresponding author: Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu

#ifndef VERTEX_H_INCLUDED
#define VERTEX_H_INCLUDED

class Vertex;
#include "../Run/Run.h"
#include "../Edge/Edge.h"

class Vertex {
public:
    long int id_;
    double position_[3];
    double volumeForce_[3];
    double velocity_[3];
    std::vector<Edge *> edges_;
    explicit Vertex(Run *, long int);
private:
    Run * run_;
};

#endif
