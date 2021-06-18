// Copyright June 2021 Tao Zhang @ Shanghai Jiao Tong University.  All Rights Reserved.
// Author: Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn
// Corresponding author: Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu

#ifndef CELL_H_INCLUDED
#define CELL_H_INCLUDED

class Cell;
#include "../Run/Run.h"
#include "../Polygon/Polygon.h"

class Cell {
public:
    long int id_;
    bool growing_;
    double volume_;
    double pressure_;
    std::vector<Polygon *> polygons_;
    std::unordered_map<long int, bool> polygonDirections_;
    explicit Cell(Run *, long int);

    int updatePolygonDirections();
    int updateVolume();
    int logPolygons(std::string);
private:
    Run * run_;
};

#endif
