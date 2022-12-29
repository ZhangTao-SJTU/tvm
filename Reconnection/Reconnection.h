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
    int dumpVtk(bool, bool, std::vector<Cell *>);
    int dumpEdgesVtk(bool, bool, std::vector<Edge *>);
private:
    Run * run_;
};

#endif
