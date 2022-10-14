/* ---------------------------------------------------------------------------------
 * Copyright 2021-2022 Tao Zhang
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

#ifndef BOX_H_INCLUDED
#define BOX_H_INCLUDED

class Box;
#include <vector>
#include "Run.h"

class Box {
  public:
    Run * run_;
    double size_[3];
    bool boundaryCondition_[3];

    Box(Run * run);
    int     resetPosition(double *);
    int     resetDistance(double *);
};

#endif
