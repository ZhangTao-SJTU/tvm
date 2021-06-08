// Copyright June 2021 Tao Zhang @ Shanghai Jiao Tong University.  All Rights Reserved.
// Author: Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn
// Corresponding author: Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu

#include <cstdio>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <unordered_map>
#include "Run.h"

using namespace std;

Run::Run() {
    Lx_ = 13.2;
    Ly_ = 200.0/13.2;
    NCell_ = 400;
}

int Run::start() {

    return 0;
}

double Run::computeD(double* v1, double* v2) {
    double dx = fabs(v2[0] - v1[0]);
    double dy = fabs(v2[1] - v1[1]);
    while (dx > Lx_/2.0) {
        dx -= Lx_;
    }
    while (dy > Ly_/2.0) {
        dy -= Ly_;
    }

    return sqrt(dx*dx + dy*dy);
}

double Run::computeD(double cx, double cy, double* v) {
    double dx = fabs(v[0] - cx);
    double dy = fabs(v[1] - cy);
    while (dx > Lx_/2.0) {
        dx -= Lx_;
    }
    while (dy > Ly_/2.0) {
        dy -= Ly_;
    }

    return sqrt(dx*dx + dy*dy);
}

int     Run::updatePolygonVertices() {
    // update vertices in polygon
    for (long int i = 0; i < polygons_.size(); i++) {
        polygons_[i]->updateVertices();
    }

    return 0;
}

int     Run::updateVertexEdges() {
    for (long int i = 0; i < vertices_.size(); i++) {
        vertices_[i]->edges_.clear();
    }
    for (long int i = 0; i < edges_.size(); i++) {
        edges_[i]->vertices_[0]->edges_.push_back(edges_[i]);
        edges_[i]->vertices_[1]->edges_.push_back(edges_[i]);
    }

    return 0;
}