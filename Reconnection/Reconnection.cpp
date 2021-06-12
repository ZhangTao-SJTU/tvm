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

#include "Reconnection.h"

using namespace std;

Reconnection::Reconnection(Run * run) {
    run_ = run;
    Lth_ = 1.0e-3;
}

int     Reconnection::start() {
    run_->updateGeoinfo();
    // set edge candidate
    for (long int i = 0; i < run_->edges_.size(); i++) {
        if (run_->edges_[i]->length_ < Lth_) {
            run_->edges_[i]->candidate_ = true;
        } else {
            run_->edges_[i]->candidate_ = false;
        }
        run_->edges_[i]->triangle_count_ = 0;
    }
    // set edge triangle_count_
    for (long int i = 0; i < run_->polygons_.size(); i++) {
        if (run_->polygons_[i]->edges_.size() == 3) {
            for (int j = 0; j < 3; j++) {
                run_->polygons_[i]->edges_[j]->triangle_count_ = run_->polygons_[i]->edges_[j]->triangle_count_ + 1;
            }
        }
    }
    // look for H type triangle candidates
    std::vector<Polygon *> triangleCandidates;
    for (long int i = 0; i < run_->polygons_.size(); i++) {
        if (run_->polygons_[i]->checkH()) {
            triangleCandidates.push_back(run_->polygons_[i]);
        }
    }

    // I -> H reconnection
    for (long int i = 0; i < run_->edges_.size(); i++) {
        if (!run_->edges_[i]->checkI()) {
            continue;
        }
        I_H(run_->edges_[i]);
    }

    // H -> I reconnection
    for (long int i = 0; i < triangleCandidates.size(); i++) {
        if (!triangleCandidates[i]->checkH()) {
            continue;
        }
        H_I(triangleCandidates[i]);
    }

    // update topology and geometry information
    run_->updateVertexCells();
    run_->updateGeoinfo();

    return 0;
}

int Reconnection::I_H(Edge * edge) {
    // locate two vertices: 10 and 11
    Vertex * v10 = edge->vertices_[0];
    Vertex * v11 = edge->vertices_[1];

    // locate five cells: 10-123, 11-456, 1011-1245, 1011-2356, 1011-1346
    if (v10->cells_.size() != 4) {
        printf("Topology Error: vertex %d has %d neighboring cells\n", v10->id_, v10->cells_.size());
        exit(1);
    }
    if (v11->cells_.size() != 4) {
        printf("Topology Error: vertex %d has %d neighboring cells\n", v11->id_, v11->cells_.size());
        exit(1);
    }
    Cell * c123 = NULL;
    Cell * c456 = NULL;
    Cell * c1245 = NULL;
    Cell * c2356 = NULL;
    Cell * c1346 = NULL;
    std::vector<Cell *> sideCells;
    for (auto cell10 : v10->cells_) {
        if (std::find(v11->cells_.begin(), v11->cells_.end(), cell10) != v11->cells_.end()) {
            sideCells.push_back(cell10);
        } else {
            c123 = cell10;
        }
    }
    if (sideCells.size() != 3) {
        printf("Topology Error: edge %d has %d neighboring side cells\n", edge->id_, sideCells.size());
        exit(1);
    }
    c1245 = sideCells[0];
    c2356 = sideCells[1];
    c1346 = sideCells[2];
    for (auto cell11 : v11->cells_) {
        if (std::find(v10->cells_.begin(), v10->cells_.end(), cell11) == v10->cells_.end()) {
            c456 = cell11;
            break;
        }
    }
    if (c456 == NULL) {
        printf("Topology Error: edge %d has 0 neighboring bottom cells\n", edge->id_);
        exit(1);
    }

}

int Reconnection::H_I(Polygon * polygon) {

}