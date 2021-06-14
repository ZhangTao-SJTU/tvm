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
    bool candidate_;    // reconnection candidate edge, length shorter than Lth
    int triangle_count_;
    bool markToDelete_;
    // 0: connected to no triangle
    // 1: connected to 1 triangle
    // 2: connected to 2 triangles
    // 3: connected to 3 triangles

    std::vector<Vertex *> vertices_;
    explicit Edge(Run *, long int);

    bool crossBoundary();
    int update();
    bool checkI();
    bool checkH();
    Vertex * otherVertex(Vertex *);
private:
    Run * run_;
};

#endif
