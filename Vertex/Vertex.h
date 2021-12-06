// Copyright June 2021 Tao Zhang @ Shanghai Jiao Tong University.  All Rights Reserved.
// Author: Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn
// Corresponding author: Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu

#ifndef VERTEX_H_INCLUDED
#define VERTEX_H_INCLUDED

class Vertex;
#include "../Run/Run.h"
#include "../Edge/Edge.h"
#include "../Cell/Cell.h"

class Vertex {
public:
    long int id_;
    long int dumpID_;
    double position_[3];
    double volumeForce_[3];
    double interfaceForce_[3];
    double pullingForce_[3];
    double velocity_[3];
    bool pull_;
    std::vector<Edge *> edges_;
    std::vector<Cell *> cells_;
    explicit Vertex(Run *, long int);

    int logCells(std::string);
    int logEdges(std::string);
private:
    Run * run_;
};

#endif
