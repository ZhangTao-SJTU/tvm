// Copyright June 2021 Tao Zhang @ Shanghai Jiao Tong University.  All Rights Reserved.
// Author: Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn
// Corresponding author: Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu

#ifndef VOLUME_H_INCLUDED
#define VOLUME_H_INCLUDED

class Volume;
#include "../Run/Run.h"

class Volume {
public:
    explicit Volume(Run *);
private:
    Run * run_;
};

#endif
