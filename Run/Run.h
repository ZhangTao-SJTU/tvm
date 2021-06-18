// Copyright June 2021 Tao Zhang @ Shanghai Jiao Tong University.  All Rights Reserved.
// Author: Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn
// Corresponding author: Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu

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

class Run {
  public:
    double  dt_;    // integration time step
    double  dtr_;   // time interval of network reconnection
    double  eta_;   // friction coefficient of vertex
    double  Lx_;
    double  Ly_;
    int     NCell_;
    double  Aic_;
    double  rho_growth_;
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
    Cell * cellTop_;    // virtual cell on the top to set the boundary
    Cell * cellBottom_; // virtual cell on the bottom to set the boundary
    Volume * volume_;
    Interface * interface_;
    Reconnection * reconnection_;

    std::vector<Vertex *> vertices_;
    std::vector<Edge *> edges_;
    std::vector<Polygon *> polygons_;
    std::vector<Cell *> cells_;

    Run();
    int     start();
    double  computeD(double *, double *);
    double  computeD(double , double, double *);
    int     updatePolygonVertices();
    int     updatePolygonCells();
    int     updatePolygonType();
    int     updatePolygonDumpType();
    int     updateVertexEdges();
    int     updateVertexCells();
    int     updateGeoinfo();
    int     updateVerticesVelocity();
    int     updateVerticesPosition();
    int     deleteVertex(Vertex *);
    int     deleteEdge(Edge *);
    int     deletePolygon(Polygon *);
    int     resetPosition(double *);
    Edge *  addEdge(Vertex *, Vertex *);
    int     dumpConfigurationVtk();
};

#endif
