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

#ifndef FIBERELASTICITY_H_INCLUDED
#define FIBERELASTICITY_H_INCLUDED

class FiberElasticity;
#include "../Run/Run.h"

class FiberElasticity {
public:
    double ks_;
    double l0_;
    double kb_;
    double stretchEnergy_;
    double bendEnergy_;

    explicit FiberElasticity(Run *);

    int updateForces();
    int updateStretchForces(Node *, Node *);
    int updateBendForces(Node *, Node *, Node *);
    int updateEnergy();
private:
    Run * run_;
};

#endif
