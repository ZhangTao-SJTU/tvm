// Copyright June 2021 Tao Zhang @ Shanghai Jiao Tong University.  All Rights Reserved.
// Author: Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn
// Corresponding author: Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu

#ifndef EDGE_H_INCLUDED
#define EDGE_H_INCLUDED

class Edge;
#include "../Run/Run.h"

class Edge {
public:
    long int id_;
    double vv_[3];  // the vector pointing from vertex 0 to 1
    double center_[3];
    double length_;
    std::vector<Vertex *> vertices_;
    explicit Edge(Run *, long int);

    bool crossBoundary();
    int update();
private:
    Run * run_;
};

#endif
