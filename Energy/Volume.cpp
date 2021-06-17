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
    kcv_ = 14.;
    vu0_ = 1.;
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

int Volume::updateVolume() {
    run_->updatePolygonVertices();
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
    for (long int i = 0; i < run_->cells_.size(); i++) {
        run_->cells_[i]->pressure_ = (-1.0)*kcv_/vu0_*(run_->cells_[i]->volume_/vu0_-1.0);
    }

    return 0;
}

int Volume::updatePolygonForces(Cell *cell, Polygon *polygon) {
    double pressure = cell->pressure_;
    // reset volumeForce values of the polygon center
    for (int m = 0; m < 3; m++) {
        polygon->volumeForce_[m] = 0.;
    }

//    // the polygon center is the reference point
//    double cc[3];   // the vector pointing from polygon center to cell center
//                    // be used to determine the outer direction of triangle face
//    for (int m = 0; m < 3; m++) {
//        cc[m] = cell->center_[m] - polygon->center_[m];
//    }
//    while (cc[0] > run_->Lx_/2.0) {
//        cc[0] = cc[0] - run_->Lx_;
//    }
//    while (cc[0] < (-1.0)*run_->Lx_/2.0) {
//        cc[0] = cc[0] + run_->Lx_;
//    }
//    while (cc[1] > run_->Ly_/2.0) {
//        cc[1] = cc[1] - run_->Ly_;
//    }
//    while (cc[1] < (-1.0)*run_->Ly_/2.0) {
//        cc[1] = cc[1] + run_->Ly_;
//    }
//
//    for (int i = 0; i < polygon->edges_.size(); i++) {
//        Edge * edge = polygon->edges_[i];
//        double cv[2][3];   // the vectors pointing from polygon center to edge vertices
//        for (int k = 0; k < 2; k++) {
//            Vertex * vertex = edge->vertices_[k];
//            for (int m = 0; m < 3; m++) {
//                cv[k][m] = vertex->position_[m] - polygon->center_[m];
//            }
//            while (cv[k][0] > run_->Lx_/2.0) {
//                cv[k][0] = cv[k][0] - run_->Lx_;
//            }
//            while (cv[k][0] < (-1.0)*run_->Lx_/2.0) {
//                cv[k][0] = cv[k][0] + run_->Lx_;
//            }
//            while (cv[k][1] > run_->Ly_/2.0) {
//                cv[k][1] = cv[k][1] - run_->Ly_;
//            }
//            while (cv[k][1] < (-1.0)*run_->Ly_/2.0) {
//                cv[k][1] = cv[k][1] + run_->Ly_;
//            }
//        }
//        // compute the vector of the triangle interface formed by polygon center, and edge vertices
//        double interface[3];
//        interface[0] = 0.5*(cv[0][1]*cv[1][2] - cv[1][1]*cv[0][2]);
//        interface[1] = 0.5*(cv[1][0]*cv[0][2] - cv[0][0]*cv[1][2]);
//        interface[2] = 0.5*(cv[0][0]*cv[1][1] - cv[1][0]*cv[0][1]);
//        double dP = 0.;
//        for (int m = 0; m < 3; m++) {
//            dP += cc[m]*interface[m];
//        }
//        // make the interface vector pointing outwards
//        if (dP > 0.) {
//            for (int m = 0; m < 3; m++) {
//                interface[m] = (-1.0)*interface[m];
//            }
//        }
//        // update volumeForces
//        for (int m = 0; m < 3; m++) {
//            edge->vertices_[0]->volumeForce_[m] = edge->vertices_[0]->volumeForce_[m] + 1.0/3.0*pressure*interface[m];
//            edge->vertices_[1]->volumeForce_[m] = edge->vertices_[1]->volumeForce_[m] + 1.0/3.0*pressure*interface[m];
//            polygon->volumeForce_[m] = polygon->volumeForce_[m] + 1.0/3.0*pressure*interface[m];
//        }
//    }

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
    for (long int i = 0; i < run_->cells_.size(); i++) {
        energy_ += 0.5*kcv_*pow(run_->cells_[i]->volume_/vu0_-1.0, 2.0);
    }

    return 0;
}