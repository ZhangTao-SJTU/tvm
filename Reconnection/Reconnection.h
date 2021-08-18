// Copyright June 2021 Tao Zhang @ Shanghai Jiao Tong University.  All Rights Reserved.
// Author: Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn
// Corresponding author: Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu

#ifndef RECONNECTION_H_INCLUDED
#define RECONNECTION_H_INCLUDED

class Reconnection;
#include "../Run/Run.h"

class Reconnection {
public:
    double  Lth_;   // threshold length of network reconnection
    long int count_IH_;
    long int count_HI_;
    bool verbose_;

    explicit Reconnection(Run *);

    int start();
    int I_H(Edge *);
    int H_I(Polygon *);
    Polygon * commonPolygon(Cell *, Cell *);
    Edge * commonEdge(Polygon *, Polygon *);
    int computeDirection(double *, double *, double *);
    int computeDistance(double *, double *, double *);
    int dumpVertices(std::vector<Vertex *>);
    int dumpCells(bool, bool, std::vector<Cell *>);
private:
    Run * run_;
};

#endif
