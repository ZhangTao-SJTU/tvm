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

#ifndef RUN_H_INCLUDED
#define RUN_H_INCLUDED

class Run;
#include <vector>
#include "../Vertex/Vertex.h"
#include "../Edge/Edge.h"
#include "../Polygon/Polygon.h"
#include "../Cell/Cell.h"
#include "../Energy/Volume.h"
#include "../Energy/Interface.h"
#include "../Reconnection/Reconnection.h"
#include "Box.h"

class Run {
  public:
    double  dt_;    // integration time step
    double  dtr_;   // time interval of network reconnection
    double  mu_;   // inverse damping coefficient of vertex
    double  kB_;
    double  temperature_;
    int     NCell_;
    double  simulation_time_;
    double   t_start_;
    double   t_end_;
    long int count_dump_;
    double   dump_period_;
    long int count_log_;
    double   log_period_;
    long int count_reconnect_;
    long int count_vertices_;
    long int count_edges_;
    long int count_polygons_;
    long int count_cells_;
    Volume * volume_;
    Interface * interface_;
    Reconnection * reconnection_;
    Box * box_;
    std::stringstream verboseReconnection_;

    std::vector<Vertex *> vertices_;
    std::vector<Edge *> edges_;
    std::vector<Polygon *> polygons_;
    std::vector<Cell *> cells_;
    std::vector<Cell *> emptyCells_;

    Run();
    int     start();
    int     updatePolygonVertices();
    int     updatePolygonCells();
    int     updateCellVertices();
    int     updateCellShapeIndex();
    int     updateVertexEdges();
    int     updateVertexCells();
    int     updateGeoinfo();
    int     updateVerticesVelocity();
    int     updateVerticesPosition();
    int     deleteVertex(Vertex *);
    int     deleteEdge(Edge *);
    int     deletePolygon(Polygon *);
    Edge *  addEdge(Vertex *, Vertex *);
    int     dumpConfigurationVtk();
    int     dumpCellCenter();
    int     dumpCellShapeIndex();
    int     dumpCellVolume();
    int     dumpTopo();
    int     dumpReconnection();
};

#endif
