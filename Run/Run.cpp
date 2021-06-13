// Copyright June 2021 Tao Zhang @ Shanghai Jiao Tong University.  All Rights Reserved.
// Author: Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn
// Corresponding author: Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu

#include <cstdio>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <unordered_map>
#include "Run.h"

using namespace std;

int DumpConfigurationVtk(double simulation_time, Run *);
int DumpITypeEdgesVtk(double simulation_time, Run *);

Run::Run() {
    dt_ = 1.0e-4;
    dtr_ = 1.0e-3;
    eta_ = 1.0;
    Lx_ = 13.2;
    Ly_ = 200.0/13.2;
    NCell_ = 400;
    Aic_ = 0.5;
    t_start_ = 0.;
    t_end_ = 5.;
    dump_period_ = 1.;
    log_period_ = 0.01;
}

int Run::start() {
    double simulation_time;
    double t_roundError = 0.01*dt_;

    count_reconnect_ = 0;
    count_dump_ = 0;
    count_log_ = 0;
    simulation_time = t_start_;
    auto start = chrono::steady_clock::now();

    printf("\nSimulation Start ...\n");
    printf("Real time elapsed: Rte\n");
    printf("Time        ");
    printf("Rte         ");
    printf("Volume      ");
    printf("E_volume    ");
    printf("E_interface ");
    printf("Energy      \n");

    while (simulation_time < t_end_ + t_roundError) {
        // update geometry information
        updateGeoinfo();
        // update volumeForces
        volume_->updateForces();
        // update interfaceForces
        interface_->updateForces();
        // update velocities
        updateVerticesVelocity();

        // log to screen
        if (simulation_time - t_start_ + t_roundError  > count_log_ * log_period_) {
            volume_->updateEnergy();
            double sum_volume = 0.;
            for (long int i = 0; i < cells_.size(); i++) {
                sum_volume += cells_[i]->volume_;
            }
            interface_->updateEnergy();
            printf("%-12.2f%-12.3f%-12.3f%-12.6f%-12.6f%-12.6f\n", simulation_time,
                   (chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - start).count())/1.0e6,
                   sum_volume,
                   volume_->energy_,
                   interface_->energy_,
                   volume_->energy_+interface_->energy_);
            start = chrono::steady_clock::now();
            count_log_++;
        }
        // dump
        if (simulation_time - t_start_ + t_roundError > count_dump_ * dump_period_) {
            DumpConfigurationVtk(t_start_ + count_dump_ * dump_period_, this);
            DumpITypeEdgesVtk(t_start_ + count_dump_ * dump_period_, this);
            count_dump_++;
        }

        // Euler dynamics
        updateVerticesPosition();

        // reconnect
        if (simulation_time - t_start_ + t_roundError > count_reconnect_ * dtr_) {
            reconnection_->start();
            count_reconnect_++;
        }

        simulation_time += dt_;
    }

//    for (long int i = 0; i < cells_.size(); i++) {
//        printf("%f\n", cells_[i]->volume_);
//    }

    return 0;
}

int     Run::updateVerticesVelocity() {
    for (long int i = 0; i < vertices_.size(); i++) {
        for (int m = 0; m < 3; m++) {
            vertices_[i]->velocity_[m] = 1.0/eta_ * (vertices_[i]->volumeForce_[m] + vertices_[i]->interfaceForce_[m]);
        }
    }

    return 0;
}

int     Run::updateVerticesPosition() {
    for (long int i = 0; i < vertices_.size(); i++) {
        for (int m = 0; m < 3; m++) {
            vertices_[i]->position_[m] = vertices_[i]->position_[m] + vertices_[i]->velocity_[m] * dt_;
        }
    }

    return 0;
}

double Run::computeD(double* v1, double* v2) {
    double dx = fabs(v2[0] - v1[0]);
    double dy = fabs(v2[1] - v1[1]);
    while (dx > Lx_/2.0) {
        dx -= Lx_;
    }
    while (dy > Ly_/2.0) {
        dy -= Ly_;
    }

    return sqrt(dx*dx + dy*dy);
}

double Run::computeD(double cx, double cy, double* v) {
    double dx = fabs(v[0] - cx);
    double dy = fabs(v[1] - cy);
    while (dx > Lx_/2.0) {
        dx -= Lx_;
    }
    while (dy > Ly_/2.0) {
        dy -= Ly_;
    }

    return sqrt(dx*dx + dy*dy);
}

int     Run::updatePolygonVertices() {
    // update vertices in polygon
    for (long int i = 0; i < polygons_.size(); i++) {
        polygons_[i]->updateVertices();
    }

    return 0;
}

int     Run::updateVertexEdges() {
    for (long int i = 0; i < vertices_.size(); i++) {
        vertices_[i]->edges_.clear();
    }
    for (long int i = 0; i < edges_.size(); i++) {
        edges_[i]->vertices_[0]->edges_.push_back(edges_[i]);
        edges_[i]->vertices_[1]->edges_.push_back(edges_[i]);
    }

    return 0;
}

int     Run::updateVertexCells() {
    for (long int i = 0; i < vertices_.size(); i++) {
        vertices_[i]->cells_.clear();
    }
    std::vector<Cell *> tpm_cells = cells_;
    tpm_cells.push_back(cellBottom_);
    tpm_cells.push_back(cellTop_);
    for (long int i = 0; i < tpm_cells.size(); i++) {
        Cell * cell = tpm_cells[i];
        for (long int j = 0; j < cell->polygons_.size(); j++) {
            for (long int k = 0; k < cell->polygons_[j]->edges_.size(); k++) {
                for (long int l = 0; l < 2; l++) {
                    Vertex * vertex = cell->polygons_[j]->edges_[k]->vertices_[l];
                    if (std::find(vertex->cells_.begin(), vertex->cells_.end(), cell) == vertex->cells_.end()) {
                        // new cell to be added
                        vertex->cells_.push_back(cell);
                    }
                }
            }
        }
    }

//    for (long int i = 0; i < vertices_.size(); i++) {
//        printf("%d\n", vertices_[i]->cells_.size());
//    }

    return 0;
}

int     Run::updatePolygonCells() {
    for (long int i = 0; i < polygons_.size(); i++) {
        polygons_[i]->cells_.clear();
    }
    for (long int i = 0; i < cells_.size(); i++) {
        for (int j = 0; j < cells_[i]->polygons_.size(); j++) {
            cells_[i]->polygons_[j]->cells_.push_back(cells_[i]);
        }
    }

    return 0;
}

int     Run::updatePolygonType() {
    updatePolygonCells();
    for (long int i = 0; i < polygons_.size(); i++) {
        if (polygons_[i]->cells_.size() > 1) {
            polygons_[i]->cell_cell = true;
        } else {
            polygons_[i]->cell_cell = false;
        }
    }

    return 0;
}

int     Run::updateGeoinfo() {
    // update edge midpoint and length
    for (long int i = 0; i < edges_.size(); i++) {
        edges_[i]->update();
//        printf("%6f\n", run->edges_[i]->length_);
    }
    // update polygon center position
    for (long int i = 0; i < polygons_.size(); i++) {
        polygons_[i]->updateCenter();
    }
    // update cell center position
    for (long int i = 0; i < cells_.size(); i++) {
        cells_[i]->updateCenter();
    }
    // update polygon type
    updatePolygonType();

    return 0;
}

int     Run::deleteVertex(Vertex * vertex) {
    auto it = find(vertices_.begin(), vertices_.end(), vertex);
    if (it != vertices_.end()) {
//        int index = it - vertices_.begin();
        vertices_.erase(it);
    } else {
        printf("vertex %d not found in vertices_\n", vertex->id_);
        exit(1);
    }
    delete vertex;

    return 0;
}

int     Run::resetPosition(double * r) {
    while (r[0] > Lx_) {
        r[0] = r[0] - Lx_;
    }
    while (r[0] < 0.) {
        r[0] = r[0] + Lx_;
    }
    while (r[1] > Ly_) {
        r[1] = r[1] - Ly_;
    }
    while (r[1] < 0.) {
        r[1] = r[1] + Ly_;
    }

    return 0;
}

Edge *  Run::addEdge(Vertex * v0, Vertex * v1) {
    Edge * edge = new Edge(this, count_edges_);
    count_edges_ += 1;
    edges_.push_back(edge);
    if (v0->id_ < v1->id_) {
        edge->vertices_.push_back(v0);
        edge->vertices_.push_back(v1);
    } else {
        edge->vertices_.push_back(v1);
        edge->vertices_.push_back(v0);
    }

    return edge;
}