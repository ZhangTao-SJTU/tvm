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
    int reconnection_candidate_;    // reconnection candidate edge type
    // -1: length larger than Lth, not considered
    // 0: length less than Lth, and connected to no triangle
    // 1: length less than Lth, and connected to 1 triangle
    // 2: length less than Lth, and connected to 2 triangles
    // 3: length less than Lth, and connected to 3 triangles

    std::vector<Vertex *> vertices_;
    explicit Edge(Run *, long int);

    bool crossBoundary();
    int update();
private:
    Run * run_;
};

#endif
