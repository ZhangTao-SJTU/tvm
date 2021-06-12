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
    // set triangle candidate
    for (long int i = 0; i < run_->polygons_.size(); i++) {
        Polygon * polygon = run_->polygons_[i];
        if (polygon->edges_.size() == 3) {
            polygon->candidate_ = true;
            for (int j = 0; j < 3; j++) {
                if (!polygon->edges_[j]->candidate_) {
                    polygon->candidate_ = false;
                    break;
                }
                if (polygon->edges_[j]->triangle_count_ != 1) {
                    polygon->candidate_ = false;
                    break;
                }
            }
        } else {
            polygon->candidate_ = false;
        }
    }

    // I -> H reconnection
    for (long int i = 0; i < run_->edges_.size(); i++) {
        if (!run_->edges_[i]->candidate_) {
            continue;
        }
        if (run_->edges_[i]->triangle_count_ > 0) {
            continue;
        }
        // TODO
    }

    // H -> I reconnection
    for (long int i = 0; i < run_->polygons_.size(); i++) {
        Polygon * polygon = run_->polygons_[i];
        bool good_candidate = true;
        if (!polygon->candidate_) {
            good_candidate = false;
        } else if (polygon->edges_.size() != 3) {
            good_candidate = false;
        } else {
            for (int j = 0; j < 3; j++) {
                if (!polygon->edges_[j]->candidate_) {
                    good_candidate = false;
                    break;
                }
                if (polygon->edges_[j]->triangle_count_ != 1) {
                    good_candidate = false;
                    break;
                }
            }
        }
        if (!good_candidate) {
            continue;
        }
        // TODO
    }

    // update topology and geometry information
    run_->updateVertexCells();
    run_->updateGeoinfo();

    return 0;
}
