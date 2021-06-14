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

#include "Cell.h"

using namespace std;

Cell::Cell(Run * run, long int id) {
    run_ = run;
    id_ = id;
    growing_ = false;
    volume_ = 0.;
    pressure_ = 0.;
    for (int i = 0; i < 3; i++) {
        center_[i] = 0.;
    }
}

int Cell::updateCenter() {
    // set reference point
    double tmp_origin[3];
    for (int i = 0; i < 3; i++) {
        tmp_origin[i] = polygons_[0]->edges_[0]->vertices_[0]->position_[i];
    }

    double sum_x = 0.;
    double sum_y = 0.;
    double sum_z = 0.;
    for (int i = 0; i < polygons_.size(); i++) {
        double dx = polygons_[i]->center_[0] - tmp_origin[0];
        double dy = polygons_[i]->center_[1] - tmp_origin[1];
        double dz = polygons_[i]->center_[2] - tmp_origin[2];
        while (dx > run_->Lx_/2.0) {
            dx -= run_->Lx_;
        }
        while (dx < (-1.0)*run_->Lx_/2.0) {
            dx += run_->Lx_;
        }
        while (dy > run_->Ly_/2.0) {
            dy -= run_->Ly_;
        }
        while (dy < (-1.0)*run_->Ly_/2.0) {
            dy += run_->Ly_;
        }
        sum_x += dx;
        sum_y += dy;
        sum_z += dz;
    }
    center_[0] = sum_x/polygons_.size() + tmp_origin[0];
    center_[1] = sum_y/polygons_.size() + tmp_origin[1];
    center_[2] = sum_z/polygons_.size() + tmp_origin[2];

    return 0;
}

int Cell::updateVolume() {
    volume_ = 0.;
    for (int i = 0; i < polygons_.size(); i++) {
        // the polygon center is the reference point
        double cc[3];   // the vector pointing from polygon center to cell center
        for (int m = 0; m < 3; m++) {
            cc[m] = center_[m] - polygons_[i]->center_[m];
        }
        while (cc[0] > run_->Lx_/2.0) {
            cc[0] = cc[0] - run_->Lx_;
        }
        while (cc[0] < (-1.0)*run_->Lx_/2.0) {
            cc[0] = cc[0] + run_->Lx_;
        }
        while (cc[1] > run_->Ly_/2.0) {
            cc[1] = cc[1] - run_->Ly_;
        }
        while (cc[1] < (-1.0)*run_->Ly_/2.0) {
            cc[1] = cc[1] + run_->Ly_;
        }

        for (int j = 0; j < polygons_[i]->edges_.size(); j++) {
            double cv[2][3];   // the vectors pointing from polygon center to edge vertices
            for (int k = 0; k < 2; k++) {
                Vertex * vertex = polygons_[i]->edges_[j]->vertices_[k];
                for (int m = 0; m < 3; m++) {
                    cv[k][m] = vertex->position_[m] - polygons_[i]->center_[m];
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
            // compute the volume of the tetrahedron formed by cell center, polygon center, and edge vertices
            double cP[3];
            double dP = 0.;
            cP[0] = cv[0][1]*cv[1][2] - cv[1][1]*cv[0][2];
            cP[1] = cv[1][0]*cv[0][2] - cv[0][0]*cv[1][2];
            cP[2] = cv[0][0]*cv[1][1] - cv[1][0]*cv[0][1];
            for (int m = 0; m < 3; m++) {
                dP += cc[m]*cP[m];
            }
            volume_ += 1.0/6.0*fabs(dP);
        }
    }

    return 0;
}

int Cell::logPolygons(std::string name) {
    printf("%s %ld\n",name.c_str(), id_);
    printf("%ld",polygons_.size());
    for (auto polygon : polygons_) {
        printf(" %ld",polygon->id_);
    }
    printf("\n");
}