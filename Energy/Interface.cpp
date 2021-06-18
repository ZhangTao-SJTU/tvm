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

    double cv[polygon->vertices_.size()][3];
    double nv[polygon->vertices_.size()][3];
    double vv[polygon->vertices_.size()][3];
    bool edgeDirections[polygon->vertices_.size()];
    // the polygon center is the reference point
    for (int i = 0; i < polygon->vertices_.size(); i++) {
        for (int m = 0; m < 3; m++) {
            cv[i][m] = polygon->vertices_[i]->position_[m] - polygon->center_[m];
        }
        while (cv[i][0] > run_->Lx_ / 2.0) {
            cv[i][0] = cv[i][0] - run_->Lx_;
        }
        while (cv[i][0] < (-1.0) * run_->Lx_ / 2.0) {
            cv[i][0] = cv[i][0] + run_->Lx_;
        }
        while (cv[i][1] > run_->Ly_ / 2.0) {
            cv[i][1] = cv[i][1] - run_->Ly_;
        }
        while (cv[i][1] < (-1.0) * run_->Ly_ / 2.0) {
            cv[i][1] = cv[i][1] + run_->Ly_;
        }
    }
    for (int i = 0; i < polygon->vertices_.size(); i++) {
        int j = (i + 1) % polygon->vertices_.size();
        // the edge vector
        for (int m = 0; m < 3; m++) {
            vv[i][m] = cv[j][m] - cv[i][m];
            // vv[m] = edge->vv_[m];
        }
        // compute the normal vector of the triangle interface formed by polygon center, and edge vertices
        nv[i][0] = cv[i][1] * cv[j][2] - cv[j][1] * cv[i][2];
        nv[i][1] = cv[j][0] * cv[i][2] - cv[i][0] * cv[j][2];
        nv[i][2] = cv[i][0] * cv[j][1] - cv[j][0] * cv[i][1];
        double norm_nv = sqrt(nv[i][0] * nv[i][0] + nv[i][1] * nv[i][1] + nv[i][2] * nv[i][2]);
        nv[i][0] = nv[i][0] / norm_nv;
        nv[i][1] = nv[i][1] / norm_nv;
        nv[i][2] = nv[i][2] / norm_nv;
        if (i == 0) {
            edgeDirections[i] = true;
            polygon->area_ = polygon->area_ + 0.5 * norm_nv;
            continue;
        }
        // check direction of nv with nv0
        double dP = 0.;
        for (int m = 0; m < 3; m++) {
            dP += nv[i][m]*nv[0][m];
        }
        if (dP > 0.) {
            edgeDirections[i] = true;
            polygon->area_ = polygon->area_ + 0.5 * norm_nv;
        } else {
            edgeDirections[i] = false;
            polygon->area_ = polygon->area_ - 0.5 * norm_nv;
        }
    }

    if (polygon->area_ < 0.) {
        polygon->area_ = fabs(polygon->area_);
        for (int i = 0; i < polygon->vertices_.size(); i++) {
            edgeDirections[i] = (!edgeDirections[i]);
        }
    }

    for (int i = 0; i < polygon->vertices_.size(); i++) {
        int j = (i + 1) % polygon->vertices_.size();
        // compute forces on triangle edges
        double Fcv0[3];
        double Fcv1[3];
        double Fvv[3];
        Fvv[0] = epsilon*(nv[i][1]*vv[i][2] - vv[i][1]*nv[i][2]);
        Fvv[1] = epsilon*(vv[i][0]*nv[i][2] - nv[i][0]*vv[i][2]);
        Fvv[2] = epsilon*(nv[i][0]*vv[i][1] - vv[i][0]*nv[i][1]);
        Fcv0[0] = epsilon*(nv[i][1]*cv[i][2] - cv[i][1]*nv[i][2]);
        Fcv0[1] = epsilon*(cv[i][0]*nv[i][2] - nv[i][0]*cv[i][2]);
        Fcv0[2] = epsilon*(nv[i][0]*cv[i][1] - cv[i][0]*nv[i][1]);
        Fcv1[0] = epsilon*(cv[j][1]*nv[i][2] - nv[i][1]*cv[j][2]);
        Fcv1[1] = epsilon*(nv[i][0]*cv[j][2] - cv[j][0]*nv[i][2]);
        Fcv1[2] = epsilon*(cv[j][0]*nv[i][1] - nv[i][0]*cv[j][1]);
        // update interfaceForces
        double sign = 1.0;
        if (!edgeDirections[i]) {
            sign = (-1.0);
        }
        for (int m = 0; m < 3; m++) {
            polygon->vertices_[i]->interfaceForce_[m] = polygon->vertices_[i]->interfaceForce_[m] + sign*0.5*(Fcv0[m]+Fvv[m]);
            polygon->vertices_[j]->interfaceForce_[m] = polygon->vertices_[j]->interfaceForce_[m] + sign*0.5*(Fcv1[m]+Fvv[m]);
            polygon->interfaceForce_[m] = polygon->interfaceForce_[m] + sign*0.5*(Fcv0[m]+Fcv1[m]);
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
