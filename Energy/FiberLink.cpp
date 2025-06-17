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

#include "FiberLink.h"

using namespace std;

FiberLink::FiberLink(Run * run) {
    run_ = run;
//    ks_ = 10.;
//    l0_ = 1.5;
//    N_ = 10;
    minArea_ = 0.1;
    minEdgeLength_ = 0.02;
    maxEdgeLength_ = 1.5;
    nFormedLinks_ = 0;
}

int     FiberLink::updateForces() {
    // reset all linkForce values in vertices and nodes
    for (auto vertex : run_->vertices_) {
        for (int m = 0; m < 3; m++) {
            vertex->linkForce_[m] = 0.;
        }
    }
    for (auto node : run_->nodes_) {
        for (int m = 0; m < 3; m++) {
            node->linkForce_[m] = 0.;
        }
    }

    // update stretchForce values
    for (auto link : run_->links_) {
        Polygon * polygon = link->polygon_;
        Node * node = link->node_;
//        cout << run_->simulation_time_ << " " << polygon->id_ << endl;

        double dx[3];
        dx[0] = node->position_[0] - polygon->center_[0];
        dx[1] = node->position_[1] - polygon->center_[1];
        dx[2] = node->position_[2] - polygon->center_[2];
        run_->box_->resetDistance(dx);
        double l = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
        double Fabs = ks_ * (-1. + l0_ / l);
        double FPolygon[3];
        FPolygon[0] = (-1.0) * Fabs * dx[0];
        FPolygon[1] = (-1.0) * Fabs * dx[1];
        FPolygon[2] = (-1.0) * Fabs * dx[2];
        node->linkForce_[0] = node->linkForce_[0] + Fabs * dx[0];
        node->linkForce_[1] = node->linkForce_[1] + Fabs * dx[1];
        node->linkForce_[2] = node->linkForce_[2] + Fabs * dx[2];

        // redistribute polygon center interfaceForces back to vertices
        double sum_l = 0.;
        for (auto edge : polygon->edges_) {
            sum_l += edge->length_;
        }
        for (auto edge : polygon->edges_) {
            double weight = edge->length_/sum_l;
            for (auto vertex : edge->vertices_) {
                for (int m = 0; m < 3; m++) {
                    vertex->linkForce_[m] = vertex->linkForce_[m] + 0.5*weight*FPolygon[m];
                }
            }
        }
    }

    return 0;
}

int FiberLink::updateEnergy() {
    stretchEnergy_ = 0.;

    for (auto link : run_->links_) {
        Polygon * polygon = link->polygon_;
        Node * node = link->node_;

        double dx[3];
        dx[0] = node->position_[0] - polygon->center_[0];
        dx[1] = node->position_[1] - polygon->center_[1];
        dx[2] = node->position_[2] - polygon->center_[2];
        run_->box_->resetDistance(dx);
        double l = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
        stretchEnergy_ += 0.5 * ks_ * pow(l - l0_, 2.0);
    }

    return 0;
}