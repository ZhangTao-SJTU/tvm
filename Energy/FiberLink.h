/* ---------------------------------------------------------------------------------
 * Copyright 2021-2023 Tao Zhang
 *
 * This file is part of TVM.
 *
 * TVM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation,
 * either version 3 of the License,
 * or (at your option) any later version.
 *
 * TVM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with TVM. If not, see <https://www.gnu.org/licenses/>.
 *
 * Author: Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn
 * Coauthor: Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu
 * ---------------------------------------------------------------------------------
 */

#ifndef FIBERLINK_H_INCLUDED
#define FIBERLINK_H_INCLUDED

class FiberLink;
#include "../Run/Run.h"

class FiberLink {
public:
    double ks_;
    double l0_;
    double l1_;
    int N_;
    double stretchEnergy_;
    double minArea_;
    double minEdgeLength_;
    double maxEdgeLength_;
    long int nFormedLinks_;
    double shrinkSpeed_;
    double shrinkStartTime_;

    explicit FiberLink(Run *);

    int updateForces();
    int updateEnergy();
private:
    Run * run_;
};

#endif
