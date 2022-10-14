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

#ifndef CELL_H_INCLUDED
#define CELL_H_INCLUDED

class Cell;
#include "../Run/Run.h"
#include "../Polygon/Polygon.h"

class Cell {
public:
    long int id_;
    double volume_;
    double pressure_;
    double shapeIndex_;
    std::vector<Polygon *> polygons_;
    std::vector<Vertex *> vertices_;
    std::unordered_map<long int, bool> polygonDirections_;
    explicit Cell(Run *, long int);

    int updatePolygonDirections();
    int updateSurfacePolygonDirections();
    int updateVolume();
    int logPolygons(std::string);
private:
    Run * run_;
};

#endif
