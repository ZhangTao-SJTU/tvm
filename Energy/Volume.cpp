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

#include "Volume.h"

using namespace std;

Volume::Volume(Run * run) {
    run_ = run;
    kv_ = 10.;  // 0.1, 1, 10, 100, 1000
    totalVolume_ = 0.;
    energy_ = 0.;
}
/* The original function to update vertex  -> volumeForces_,*/

// This works by distributing the pressure

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

/* Updated Volume::updateForces() function

// This calcuates the exact volume force on each vertex
int     Volume::updateForces() {
    // initialize volumeForce values in all vertices
    for (auto vertex : run_->vertices_) {
        for (int m = 0; m < 3; m++) {
            vertex->volumeForce_[m] = 0.;
        }
    }
    updateVolume();
    // Main algorithm
    for (auto cell : run_->cells_) {
        for (auto polygon : cell->polygons_) {
            //Calculate sum_term = 1/n sum_mu(r_mu cross r_{mu+1})
            // Step 1: Set origin to the first vertex of the polygon
            double origin[3];
            for (int m = 0; m < 3; m++) {
                origin[m] = polygon->vertices_[0]->position_[m];
            }
            // Step 2: Calculate the sum_term = 1/n sum_mu(r_mu cross r_{mu+1})
            double sum_term[3];
            for (int m = 0; m < 3; m++) {
                sum_term[m] = 0.;
            }
            for (int i = 0; i < polygon->vertices_.size(); i++) {
                int j = (i + 1)%polygon->vertices_.size();
                double r_mu[3];
                double r_mu_plus_1[3];
                for (int m = 0; m < 3; m++) {
                    r_mu[m] = polygon->vertices_[i]->position_[m] - origin[m];
                    r_mu_plus_1[m] = polygon->vertices_[j]->position_[m] - origin[m];
                }
                // adjust for periodic boundary conditions
                run_->box_->resetDistance(r_mu);
                run_->box_->resetDistance(r_mu_plus_1);
                // Calculate the cross product of r_mu and r_{mu+1}
                double cross_product[3];
                cross_product[0] = r_mu[1]*r_mu_plus_1[2] - r_mu[2]*r_mu_plus_1[1];
                cross_product[1] = r_mu[2]*r_mu_plus_1[0] - r_mu[0]*r_mu_plus_1[2];
                cross_product[2] = r_mu[0]*r_mu_plus_1[1] - r_mu[1]*r_mu_plus_1[0];
                // Add the cross product to sum_term
                for (int m = 0; m < 3; m++) {
                    sum_term[m] += cross_product[m];
                }
            }
            for (int m = 0; m < 3; m++) {
                sum_term[m] = sum_term[m]/polygon->vertices_.size();
            }
            // Step 3: For each vertex on the polygon:
            //(a) Add sum_term to vertex->volumeForce_
            //(b) Add (r_mu+1 -r_mu-1) cross poly_center) to vertex->volumeForce_
            // Add these with a factor of -1 if the polygon direction is incorrect
            for (int i = 0; i < polygon->vertices_.size(); i++) {
                int j = (i + 1)%polygon->vertices_.size();
                int k = (i + polygon->vertices_.size() - 1)%polygon->vertices_.size();
                double r_mu_plus_1[3];
                double r_mu_minus_1[3];
                double poly_center[3];
                for (int m = 0; m < 3; m++) {
                    r_mu_plus_1[m] = polygon->vertices_[j]->position_[m] - origin[m];
                    r_mu_minus_1[m] = polygon->vertices_[k]->position_[m] - origin[m];
                    poly_center[m] = polygon->center_[m] - origin[m];
                }
                // adjust for periodic boundary conditions
                run_->box_->resetDistance(r_mu_plus_1);
                run_->box_->resetDistance(r_mu_minus_1);
                run_->box_->resetDistance(poly_center);
                // Calculate the cross product of (r_mu+1 - r_mu-1) and poly_center
                double cross_product[3];
                cross_product[0] = (r_mu_plus_1[1] - r_mu_minus_1[1])*poly_center[2] - (r_mu_plus_1[2] - r_mu_minus_1[2])*poly_center[1];
                cross_product[1] = (r_mu_plus_1[2] - r_mu_minus_1[2])*poly_center[0] - (r_mu_plus_1[0] - r_mu_minus_1[0])*poly_center[2];
                cross_product[2] = (r_mu_plus_1[0] - r_mu_minus_1[0])*poly_center[1] - (r_mu_plus_1[1] - r_mu_minus_1[1])*poly_center[0];
                // Add (or subtract) to vertex->volumeForce_ based on polygon directions
                for (int m = 0; m < 3; m++) {
                    if (cell->polygonDirections_[polygon->id_]) {
                        polygon->vertices_[i]->volumeForce_[m] -= 2*kv_*(cell->volume_-1)*(sum_term[m]+cross_product[m]);
                    } else {
                        polygon->vertices_[i]->volumeForce_[m] += 2*kv_*(cell->volume_-1)*(sum_term[m]+cross_product[m]);
                    }
                }
           }
        }
    }

    return 0;
}
*/
int Volume::updatePolygonDirections() {
    run_->updatePolygonVertices();
    run_->updatePolygonCells();
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
        run_->box_->resetDistance(cv[i]);
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