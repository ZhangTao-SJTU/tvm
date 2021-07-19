// Copyright June 2021 Tao Zhang @ Shanghai Jiao Tong University.  All Rights Reserved.
// Author: Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn
// Corresponding author: Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu

#include <cstdio>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <unordered_map>

#include "Volume.h"

using namespace std;

Volume::Volume(Run * run) {
    run_ = run;
    kv_ = 10.;
    totalVolume_ = 0.;
    energy_ = 0.;
}

int     Volume::updateForces() {
    // reset all volumeForce values in vertices
    for (long int i = 0; i < run_->vertices_.size(); i++) {
        for (int j = 0; j < 3; j++) {
            run_->vertices_[i]->volumeForce_[j] = 0.;
        }
    }

    // update volume of each cell, and direction of polygons in each cell
    updateVolume();

    // update pressure in each cell
    updatePressure();

    // update volumeForce values
    for (long int i = 0; i < run_->cells_.size(); i++) {
        for (int j = 0; j < run_->cells_[i]->polygons_.size(); j++) {
            updatePolygonForces(run_->cells_[i], run_->cells_[i]->polygons_[j]);
        }
    }

    return 0;
}

int Volume::updatePolygonDirections() {
    run_->updatePolygonVertices();
    for (auto cell : run_->cells_) {
        cell->updatePolygonDirections();
    }

    return 0;
}

int Volume::updateVolume() {
    // update cell volume
    totalVolume_ = 0.;
    for (auto cell : run_->cells_) {
        cell->updateVolume();
//        if (run_->simulation_time_ < run_->t_start_+0.01*run_->dt_) {
//            printf("%6f\n", cell->volume_);
//        }
        totalVolume_ += cell->volume_;
    }

    return 0;
}

int Volume::updatePressure() {
    for (auto cell : run_->cells_) {
        cell->pressure_ = (-1.0)*2.0*kv_*(cell->volume_-1.0);
    }

    return 0;
}

int Volume::updatePolygonForces(Cell *cell, Polygon *polygon) {
    double pressure = cell->pressure_;
    bool correctDirection = cell->polygonDirections_[polygon->id_];
    int Nv = polygon->vertices_.size();
    // reset volumeForce values of the polygon center
    for (int m = 0; m < 3; m++) {
        polygon->volumeForce_[m] = 0.;
    }

    double cv[Nv][3];   // the vectors pointing from polygon center to edge vertices
    for (int i = 0; i < Nv; i++) {
        for (int m = 0; m < 3; m++) {
            cv[i][m] = polygon->vertices_[i]->position_[m] - polygon->center_[m];
        }
        while (cv[i][0] > run_->Lx_/2.0) {
            cv[i][0] = cv[i][0] - run_->Lx_;
        }
        while (cv[i][0] < (-1.0)*run_->Lx_/2.0) {
            cv[i][0] = cv[i][0] + run_->Lx_;
        }
        while (cv[i][1] > run_->Ly_/2.0) {
            cv[i][1] = cv[i][1] - run_->Ly_;
        }
        while (cv[i][1] < (-1.0)*run_->Ly_/2.0) {
            cv[i][1] = cv[i][1] + run_->Ly_;
        }
        while (cv[i][2] > run_->Lz_/2.0) {
            cv[i][2] = cv[i][2] - run_->Lz_;
        }
        while (cv[i][2] < (-1.0)*run_->Lz_/2.0) {
            cv[i][2] = cv[i][2] + run_->Lz_;
        }
    }
    for (int i = 0; i < Nv; i++) {
        // compute the vector of the triangle interface formed by polygon center, and edge vertices
        int j = (i + 1)%Nv;
        double interface[3];
        interface[0] = 0.5*(cv[i][1]*cv[j][2] - cv[j][1]*cv[i][2]);
        interface[1] = 0.5*(cv[j][0]*cv[i][2] - cv[i][0]*cv[j][2]);
        interface[2] = 0.5*(cv[i][0]*cv[j][1] - cv[j][0]*cv[i][1]);
        // make the interface vector pointing outwards
        if (!correctDirection) {
            for (int m = 0; m < 3; m++) {
                interface[m] = (-1.0)*interface[m];
            }
        }
        // update volumeForces
        Vertex * v0 = polygon->vertices_[i];
        Vertex * v1 = polygon->vertices_[j];
        for (int m = 0; m < 3; m++) {
            v0->volumeForce_[m] = v0->volumeForce_[m] + 1.0/3.0*pressure*interface[m];
            v1->volumeForce_[m] = v1->volumeForce_[m] + 1.0/3.0*pressure*interface[m];
            polygon->volumeForce_[m] = polygon->volumeForce_[m] + 1.0/3.0*pressure*interface[m];
        }
    }

    // redistribute polygon center volumeForces back to vertices
    double sum_l = 0.;
    for (int i = 0; i < polygon->edges_.size(); i++) {
        sum_l += polygon->edges_[i]->length_;
    }
    for (int i = 0; i < polygon->edges_.size(); i++) {
        double weight = polygon->edges_[i]->length_/sum_l;
        for (int k = 0; k < 2; k++) {
            Vertex *vertex = polygon->edges_[i]->vertices_[k];
            for (int m = 0; m < 3; m++) {
                vertex->volumeForce_[m] = vertex->volumeForce_[m] + 0.5*weight*polygon->volumeForce_[m];
            }
        }
    }

    return 0;
}

int Volume::updateEnergy() {
    energy_ = 0.;
    for (auto cell : run_->cells_) {
        energy_ += kv_*pow(cell->volume_-1.0, 2.0);
    }

    return 0;
}