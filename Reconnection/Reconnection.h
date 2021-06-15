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

    explicit Reconnection(Run *);

    int start();
    int I_H(Edge *, bool verbose);
    int H_I(Polygon *, bool verbose);
    Polygon * commonPolygon(Cell *, Cell *);
    Edge * commonEdge(Polygon *, Polygon *);
    int computeDirection(double *, double *, double *);
    int dumpVtk(std::vector<Polygon *>, bool, bool);
private:
    Run * run_;
};

#endif
