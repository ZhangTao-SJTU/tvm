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

#ifndef POLYGON_H_INCLUDED
#define POLYGON_H_INCLUDED

class Polygon;
#include "../Run/Run.h"
#include "../Edge/Edge.h"
#include "../Cell/Cell.h"

class Polygon {
public:
    long int id_;
    double center_[3];
    double area_;
    double tension_;
    double volumeForce_[3];
    double interfaceForce_[3];
    std::vector<Edge *> edges_;
    std::vector<Vertex *> vertices_;
    std::vector<Cell *> cells_;
    explicit Polygon(Run *, long int);

    int updateVertices();
    int updateCenter();
    int updateArea();
    bool crossBoundary();
    bool checkH();
    int shrink(Edge *);
    int expand(Edge *);
    int logEdges(std::string);
private:
    Run * run_;
};

#endif
