// Copyright June 2021 Tao Zhang @ Shanghai Jiao Tong University.  All Rights Reserved.
// Author: Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn
// Corresponding author: Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu

#ifndef VOLUME_H_INCLUDED
#define VOLUME_H_INCLUDED

class Volume;
#include "../Run/Run.h"

class Volume {
public:
    double kv_;
    double totalVolume_;
    double energy_;

    explicit Volume(Run *);

    int updateForces();
    int updatePolygonDirections();
    int updateVolume();
    int updatePressure();
    int updatePolygonForces(Cell *, Polygon *);
    int updateEnergy();
private:
    Run * run_;
};

#endif
