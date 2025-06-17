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
#include <deque>

#include "Volume.h"

using namespace std;

Volume::Volume(Run * run) {
    run_ = run;
    kv_ = 10.;  // 0.1, 1, 10, 100, 1000
    totalVolume_ = 0.;
    energy_ = 0.;
}

int     Volume::updateForces() {
    // reset all volumeForce values in vertices
    for (auto vertex : run_->vertices_) {
        for (int m = 0; m < 3; m++) {
            vertex->volumeForce_[m] = 0.;
        }
    }

    // update volume of each cell, and direction of polygons in each cell
    updateVolume();

    // update pressure in each cell
    updatePressure();

    // update volumeForce values
    for (auto cell : run_->cells_) {
        if (cell->type_ < 0) {
            continue;
        }
        for (auto polygon : cell->polygons_) {
            updatePolygonForces(cell, polygon);
        }
    }

    return 0;
}

int Volume::updatePolygonDirections() {
    run_->updatePolygonVertices();
    run_->updatePolygonCells();
    for (auto cell : run_->cells_) {
        if (cell->type_ < 0) {
            continue;
        }
        cell->updatePolygonDirections();
    }
    // adjust polygon directions of each cell, so that one polygon's directions in two cells are opposite
    std::deque<Cell *> queneCells;
    std::unordered_map<long int, bool> visited;
    for (auto cell : run_->cells_) {
        if (cell->type_ < 0) {
            visited[cell->id_] = true;;
        } else {
            visited[cell->id_] = false;
        }
    }
    visited[run_->cells_[0]->id_] = true;
    queneCells.push_back(run_->cells_[0]);
    while (queneCells.size() > 0) {
        Cell * cell = queneCells[0];
        queneCells.pop_front();
        for (auto polygon : cell->polygons_) {
            Cell * nextCell = polygon->cells_[0];
            if (nextCell->id_ == cell->id_) {
                nextCell = polygon->cells_[1];
                if (nextCell->id_ == cell->id_) {
                    cout << "bullshit" <<endl;
                    exit(1);
                }
            }
            if (!visited[nextCell->id_]) {
                visited[nextCell->id_] = true;
                queneCells.push_back(nextCell);
                if (cell->polygonDirections_[polygon->id_] == nextCell->polygonDirections_[polygon->id_]) {
                    for (auto p: nextCell->polygons_) {
                        nextCell->polygonDirections_[p->id_] = (!nextCell->polygonDirections_[p->id_]);
                    }
                }
            }
        }
    }
    // depend on total volume of the system, flip all cell polygon directions
    updateVolume();
    if (totalVolume_ < 0) {
        for (auto cell : run_->cells_) {
            if (cell->type_ < 0) {
                continue;
            }
            for (auto polygon: cell->polygons_) {
                cell->polygonDirections_[polygon->id_] = (!cell->polygonDirections_[polygon->id_]);
            }
        }
//        cout << "flipped" << endl;
    }

    // update polygonDirections of emptySpace_
    for (auto polygon: run_->emptySpace_->polygons_) {
        Cell * cell = polygon->cells_[0];
        if (cell->type_ < 0) {
            cell = polygon->cells_[1];
        }
        run_->emptySpace_->polygonDirections_[polygon->id_] = (!cell->polygonDirections_[polygon->id_]);
    }

    return 0;
}

int Volume::updateVolume() {
    // update cell volume
    totalVolume_ = 0.;
    for (auto cell : run_->cells_) {
        if (cell->type_ < 0) {
            continue;
        }
        cell->updateVolume();
        totalVolume_ += cell->volume_;
    }

    return 0;
}

int Volume::updatePressure() {
    for (auto cell : run_->cells_) {
        if (cell->type_ < 0) {
            continue;
        }
        //**EMPTY CELL CHECK**//
        if (cell->type_ == 0) {
            cell->pressure_ = (-1.0) * 2.0 * kv_ * (cell->volume_ - 1.0);
        } else {
            cell->pressure_ = (-1.0) * 2.0 * kv_ * (cell->volume_ - 1.0);
        }
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
        //**EMPTY CELL CHECK**//
        if (cell->type_ == 0) {
            if (v0->type_ < 3) {
                for (int m = 0; m < 3; m++) {
                    v0->volumeForce_[m] = v0->volumeForce_[m] + 1.0/3.0*pressure*interface[m];
                }
            }
            if (v1->type_ < 3) {
                for (int m = 0; m < 3; m++) {
                    v1->volumeForce_[m] = v1->volumeForce_[m] + 1.0/3.0*pressure*interface[m];
                }
            }
            for (int m = 0; m < 3; m++) {
                polygon->volumeForce_[m] = polygon->volumeForce_[m] + 1.0/3.0*pressure*interface[m];
            }
        } else {
            for (int m = 0; m < 3; m++) {
                v0->volumeForce_[m] = v0->volumeForce_[m] + 1.0 / 3.0 * pressure * interface[m];
                v1->volumeForce_[m] = v1->volumeForce_[m] + 1.0 / 3.0 * pressure * interface[m];
                polygon->volumeForce_[m] = polygon->volumeForce_[m] + 1.0 / 3.0 * pressure * interface[m];
            }
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
            //**EMPTY CELL CHECK**//
            if (cell->type_ == 0 && vertex->type_ >= 3) {
                continue;
            }
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
        //**EMPTY CELL CHECK**//
        if (cell->type_ <= 0) {
            continue;
        }
        energy_ += kv_*pow(cell->volume_-1.0, 2.0);
    }

    return 0;
}