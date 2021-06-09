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

int DumpConfigurationVtk(double simulation_time, Run *);

Run::Run() {
    dt_ = 1.0e-4;
    dtr_ = 1.0e-3;
    Lth_ = 1.0e-3;
    Lx_ = 13.2;
    Ly_ = 200.0/13.2;
    NCell_ = 400;
    t_start_ = 0.;
    t_end_ = 40.;
    dump_period_ = 1.;
    log_period_ = 0.01;
}

int Run::start() {
    double simulation_time;

    count_reconnect_ = 0;
    count_dump_ = 0;
    count_log_ = 0;
    simulation_time = t_start_;
    auto start = chrono::steady_clock::now();

    printf("Real time elapsed: Rte\n");
    printf("Start Simulation\n");
    printf("      Time    Rte\n");

    while (simulation_time <= t_end_) {
        // update energy, force, velocity
        updateGeoinfo();
        // TODO updateVolumeForces();
        // TODO updateInterfaceForces();

        // reconnect
        if (simulation_time - t_start_ + 1e-6 > count_reconnect_ * dtr_) {
            // TODO reconnect
            count_reconnect_++;
        }

        // log to screen
        if (simulation_time - t_start_ + 1e-6  > count_log_ * log_period_) {
            printf("%12.2f %.3f\n", simulation_time,
                   (chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - start).count())/1.0e6);
            start = chrono::steady_clock::now();
            count_log_++;
        }
        // dump
        if (simulation_time - t_start_ + 1e-6 > count_dump_ * dump_period_) {
            DumpConfigurationVtk(t_start_ + count_dump_ * dump_period_, this);
            count_dump_++;
        }

        // Euler dynamics
        // TODO UpdateXYZ();  // Updates nodal positions
        simulation_time += dt_;
    }

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

int     Run::updateGeoinfo() {
    for (long int i = 0; i < edges_.size(); i++) {
        edges_[i]->update();
//        printf("%6f\n", run->edges_[i]->length_);
    }
    for (long int i = 0; i < polygons_.size(); i++) {
        polygons_[i]->updateCenter();
    }
    for (long int i = 0; i < cells_.size(); i++) {
        cells_[i]->updateCenter();
    }
    // update cell volumes
    for (long int i = 0; i < cells_.size(); i++) {
        cells_[i]->updateVolume();
//        printf("%6f\n", run->cells_[i]->volume_);
    }
}