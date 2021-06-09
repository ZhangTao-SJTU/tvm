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
    eta_ = 1.0;
    Lth_ = 1.0e-3;
    Lx_ = 13.2;
    Ly_ = 200.0/13.2;
    NCell_ = 400;
    t_start_ = 0.;
    t_end_ = 5.;
    dump_period_ = 1.;
    log_period_ = 0.01;
}

int Run::start() {
    double simulation_time;
    double t_roundError = 0.01*dt_;

    count_reconnect_ = 0;
    count_dump_ = 0;
    count_log_ = 0;
    simulation_time = t_start_;
    auto start = chrono::steady_clock::now();

    printf("Real time elapsed: Rte\n");
    printf("Start Simulation\n");
    printf("Time        ");
    printf("Rte         ");
    printf("Volume      ");
    printf("Energy      \n");

    while (simulation_time <= t_end_ + t_roundError) {
        // update energy, force, velocity
        updateGeoinfo();
        // update volumeForces
        volume_->updateForces();
        // update interfaceForces
        // TODO interface_->updateForces();
        // update velocities
        updateVerticesVelocity();

        // log to screen
        if (simulation_time - t_start_ + t_roundError  > count_log_ * log_period_) {
            volume_->updateEnergy();
            double sum_volume = 0.;
            for (long int i = 0; i < cells_.size(); i++) {
                sum_volume += cells_[i]->volume_;
            }
            printf("%-12.2f%-12.3f%-12.3f%-12.6f\n", simulation_time,
                   (chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - start).count())/1.0e6,
                   sum_volume,
                   volume_->energy_);
            start = chrono::steady_clock::now();
            count_log_++;
        }
        // dump
        if (simulation_time - t_start_ + t_roundError > count_dump_ * dump_period_) {
            DumpConfigurationVtk(t_start_ + count_dump_ * dump_period_, this);
            count_dump_++;
        }

        // Euler dynamics
        updateVerticesPosition();

        // reconnect
        if (simulation_time - t_start_ + t_roundError > count_reconnect_ * dtr_) {
            // TODO reconnect
            count_reconnect_++;
        }

        simulation_time += dt_;
    }

    for (long int i = 0; i < cells_.size(); i++) {
        printf("%f\n", cells_[i]->volume_);
    }

    return 0;
}

int     Run::updateVerticesVelocity() {
    for (long int i = 0; i < vertices_.size(); i++) {
        for (int m = 0; m < 3; m++) {
            vertices_[i]->velocity_[m] = 1.0/eta_ * (vertices_[i]->volumeForce_[m]);
        }
    }

    return 0;
}

int     Run::updateVerticesPosition() {
    for (long int i = 0; i < vertices_.size(); i++) {
        for (int m = 0; m < 3; m++) {
            vertices_[i]->position_[m] = vertices_[i]->position_[m] + vertices_[i]->velocity_[m] * dt_;
        }
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
    // update edge midpoint and length
    for (long int i = 0; i < edges_.size(); i++) {
        edges_[i]->update();
//        printf("%6f\n", run->edges_[i]->length_);
    }
    // update polygon center position
    for (long int i = 0; i < polygons_.size(); i++) {
        polygons_[i]->updateCenter();
    }
    // update cell center position
    for (long int i = 0; i < cells_.size(); i++) {
        cells_[i]->updateCenter();
    }
    // update cell volume
    for (long int i = 0; i < cells_.size(); i++) {
        cells_[i]->updateVolume();
//        printf("%6f\n", run->cells_[i]->volume_);
    }
}