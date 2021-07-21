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

#include "Edge.h"

using namespace std;

Edge::Edge(Run * run, long int id) {
    run_ = run;
    id_ = id;
    for (int i = 0; i < 3; i++) {
        vv_[i] = 0.;
        center_[i] = 0.;
    }
    length_ = 0.;
    candidate_ = false;
    triangle_count_ = 0;
    markToDelete_ = false;
}

bool Edge::crossBoundary() {
    if (fabs(vertices_[1]->position_[0] - vertices_[0]->position_[0]) > run_->Lx_/2.0) {
        return true;
    }
    if (fabs(vertices_[1]->position_[1] - vertices_[0]->position_[1]) > run_->Ly_/2.0) {
        return true;
    }
    if (fabs(vertices_[1]->position_[2] - vertices_[0]->position_[2]) > run_->Lz_/2.0) {
        return true;
    }

    return false;
}

int Edge::update() {
    double dx = vertices_[1]->position_[0] - vertices_[0]->position_[0];
    double dy = vertices_[1]->position_[1] - vertices_[0]->position_[1];
    double dz = vertices_[1]->position_[2] - vertices_[0]->position_[2];
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
    while (dz > run_->Lz_/2.0) {
        dz -= run_->Lz_;
    }
    while (dz < (-1.0)*run_->Lz_/2.0) {
        dz += run_->Lz_;
    }
//    dx = dx - run_->Lx_ * floor((dx + run_->Lx_/2.0) / run_->Lx_);
//    dy = dy - run_->Ly_ * floor((dy + run_->Ly_/2.0) / run_->Ly_);
    vv_[0] = dx;
    vv_[1] = dy;
    vv_[2] = dz;
    length_ = sqrt(dx*dx + dy*dy + dz*dz);
    center_[0] = vertices_[0]->position_[0] + dx/2.0;
    center_[1] = vertices_[0]->position_[1] + dy/2.0;
    center_[2] = vertices_[0]->position_[2] + dz/2.0;

    return 0;
}

bool Edge::checkI() {
    if (!candidate_) {
        return false;
    }
    if (triangle_count_ > 0) {
        return false;
    }

    return true;
}

bool Edge::checkH() {
    if (!candidate_) {
        return false;
    }
    if (triangle_count_ != 1) {
        return false;
    }

    return true;
}

Vertex * Edge::otherVertex(Vertex * v0) {
    if (v0 == vertices_[0]) {
        return vertices_[1];
    } else if (v0 == vertices_[1]) {
        return vertices_[0];
    } else {
        return NULL;
    }
}