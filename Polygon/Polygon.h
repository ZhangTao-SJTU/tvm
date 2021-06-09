// Copyright June 2021 Tao Zhang @ Shanghai Jiao Tong University.  All Rights Reserved.
// Author: Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn
// Corresponding author: Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu

#ifndef POLYGON_H_INCLUDED
#define POLYGON_H_INCLUDED

class Polygon;
#include "../Run/Run.h"
#include "../Edge/Edge.h"

class Polygon {
public:
    long int id_;
    std::vector<Edge *> edges_;
    std::vector<Vertex *> vertices_;
    double center_[3];
    double volumeForce_[3];
    explicit Polygon(Run *, long int);

    int updateVertices();
    int updateCenter();
    bool crossBoundary();
private:
    Run * run_;
};

#endif
