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

#include "Interface.h"

using namespace std;

Interface::Interface(Run * run) {
    run_ = run;
    epsilon_cc_ = 1.;
    epsilon_co_ = 0.7;
    energy_ = 0.;
}

int     Interface::updateForces() {
    // reset all interfaceForce values in vertices
    for (long int i = 0; i < run_->vertices_.size(); i++) {
        for (int j = 0; j < 3; j++) {
            run_->vertices_[i]->interfaceForce_[j] = 0.;
        }
    }

    // update interfaceForce values
    for (long int i = 0; i < run_->polygons_.size(); i++) {
        updatePolygonForces(run_->polygons_[i]);
    }

    return 0;
}

int Interface::updatePolygonForces(Polygon *polygon) {
    // reset interfaceForce values of the polygon center
    for (int m = 0; m < 3; m++) {
        polygon->interfaceForce_[m] = 0.;
    }

    double epsilon;
    if (polygon->cell_cell) {
        epsilon = epsilon_cc_;
    } else {
        epsilon = epsilon_co_;
    }

    // reset polygon area
    polygon->area_ = 0.;

    // the polygon center is the reference point
    for (int i = 0; i < polygon->edges_.size(); i++) {
        Edge * edge = polygon->edges_[i];
        // the vectors pointing from polygon center to edge vertices
        double cv[2][3];
        for (int k = 0; k < 2; k++) {
            Vertex * vertex = edge->vertices_[k];
            for (int m = 0; m < 3; m++) {
                cv[k][m] = vertex->position_[m] - polygon->center_[m];
            }
            while (cv[k][0] > run_->Lx_/2.0) {
                cv[k][0] = cv[k][0] - run_->Lx_;
            }
            while (cv[k][0] < (-1.0)*run_->Lx_/2.0) {
                cv[k][0] = cv[k][0] + run_->Lx_;
            }
            while (cv[k][1] > run_->Ly_/2.0) {
                cv[k][1] = cv[k][1] - run_->Ly_;
            }
            while (cv[k][1] < (-1.0)*run_->Ly_/2.0) {
                cv[k][1] = cv[k][1] + run_->Ly_;
            }
        }
        // the edge vector
        double vv[3];
        for (int m = 0; m < 3; m++) {
            vv[m] = cv[1][m] - cv[0][m];
            // vv[m] = edge->vv_[m];
        }
        // compute the normal vector of the triangle interface formed by polygon center, and edge vertices
        double nv[3];
        nv[0] = cv[0][1]*cv[1][2] - cv[1][1]*cv[0][2];
        nv[1] = cv[1][0]*cv[0][2] - cv[0][0]*cv[1][2];
        nv[2] = cv[0][0]*cv[1][1] - cv[1][0]*cv[0][1];
        double norm_nv = sqrt(nv[0]*nv[0] + nv[1]*nv[1] + nv[2]*nv[2]);
        polygon->area_ = polygon->area_ + 0.5*norm_nv;
        nv[0] = nv[0] / norm_nv;
        nv[1] = nv[1] / norm_nv;
        nv[2] = nv[2] / norm_nv;
        // compute forces on triangle edges
        double Fcv0[3];
        double Fcv1[3];
        double Fvv[3];
        Fvv[0] = epsilon*(nv[1]*vv[2] - vv[1]*nv[2]);
        Fvv[1] = epsilon*(vv[0]*nv[2] - nv[0]*vv[2]);
        Fvv[2] = epsilon*(nv[0]*vv[1] - vv[0]*nv[1]);
        Fcv0[0] = epsilon*(nv[1]*cv[0][2] - cv[0][1]*nv[2]);
        Fcv0[1] = epsilon*(cv[0][0]*nv[2] - nv[0]*cv[0][2]);
        Fcv0[2] = epsilon*(nv[0]*cv[0][1] - cv[0][0]*nv[1]);
        Fcv1[0] = epsilon*(cv[1][1]*nv[2] - nv[1]*cv[1][2]);
        Fcv1[1] = epsilon*(nv[0]*cv[1][2] - cv[1][0]*nv[2]);
        Fcv1[2] = epsilon*(cv[1][0]*nv[1] - nv[0]*cv[1][1]);
        // update interfaceForces
        for (int m = 0; m < 3; m++) {
            edge->vertices_[0]->interfaceForce_[m] = edge->vertices_[0]->interfaceForce_[m] + 0.5*(Fcv0[m]+Fvv[m]);
            edge->vertices_[1]->interfaceForce_[m] = edge->vertices_[1]->interfaceForce_[m] + 0.5*(Fcv1[m]+Fvv[m]);
            polygon->interfaceForce_[m] = polygon->interfaceForce_[m] + 0.5*(Fcv0[m]+Fcv1[m]);
        }
    }

    // redistribute polygon center interfaceForces back to vertices
    double sum_l = 0.;
    for (int i = 0; i < polygon->edges_.size(); i++) {
        sum_l += polygon->edges_[i]->length_;
    }
    for (int i = 0; i < polygon->edges_.size(); i++) {
        double weight = polygon->edges_[i]->length_/sum_l;
        for (int k = 0; k < 2; k++) {
            Vertex *vertex = polygon->edges_[i]->vertices_[k];
            for (int m = 0; m < 3; m++) {
                vertex->interfaceForce_[m] = vertex->interfaceForce_[m] + 0.5*weight*polygon->interfaceForce_[m];
            }
        }
    }

    return 0;
}

int Interface::updateEnergy() {
    energy_ = 0.;
    for (long int i = 0; i < run_->polygons_.size(); i++) {
        double epsilon;
        if (run_->polygons_[i]->cell_cell) {
            epsilon = epsilon_cc_;
        } else {
            epsilon = epsilon_co_;
        }
        energy_ += epsilon*run_->polygons_[i]->area_;
    }

    return 0;
}
