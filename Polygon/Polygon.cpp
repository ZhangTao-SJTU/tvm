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

#include "Polygon.h"

using namespace std;

Polygon::Polygon(Run * run, long int id) {
    run_ = run;
    id_ = id;
    for (int i = 0; i < 3; i++) {
        center_[i] = 0.;
        volumeForce_[i]  = 0.;
        interfaceForce_[i]  = 0.;
    }
    area_ = 0.;
    cell_cell = false;
    dumpType = 0;
    dumpVolumeRatio = 1.;
}

int Polygon::updateVertices() {
    vertices_.clear();
    std::vector<Edge *> tmp_edges = edges_;
    Vertex * currentVertex = NULL;
    for (int i = edges_.size()-1; i > 0; i--) {
        if (i == edges_.size()-1) {
            vertices_.push_back(tmp_edges[tmp_edges.size()-1]->vertices_[0]);
            currentVertex = tmp_edges[tmp_edges.size()-1]->vertices_[1];
            vertices_.push_back(currentVertex);
            tmp_edges.pop_back();
            continue;
        }
        for (int j = 0; j < tmp_edges.size(); j++) {
            if (currentVertex == tmp_edges[j]->vertices_[0]) {
                currentVertex = tmp_edges[j]->vertices_[1];
                vertices_.push_back(currentVertex);
                if (j < tmp_edges.size() - 1) {
                    tmp_edges[j] = tmp_edges[tmp_edges.size()-1];
                }
                tmp_edges.pop_back();
                break;
            }
            if (currentVertex == tmp_edges[j]->vertices_[1]) {
                currentVertex = tmp_edges[j]->vertices_[0];
                vertices_.push_back(currentVertex);
                if (j < tmp_edges.size() - 1) {
                    tmp_edges[j] = tmp_edges[tmp_edges.size()-1];
                }
                tmp_edges.pop_back();
                break;
            }
        }
    }

    return 0;
}

int Polygon::updateCenter() {
    // set reference point
    double tmp_origin[3];
    for (int i = 0; i < 3; i++) {
        tmp_origin[i] = edges_[0]->vertices_[0]->position_[i];
    }

    double sum_lx = 0.;
    double sum_ly = 0.;
    double sum_lz = 0.;
    double sum_l = 0.;
    for (int i = 0; i < edges_.size(); i++) {
        double length = edges_[i]->length_;
        double dx = edges_[i]->center_[0] - tmp_origin[0];
        double dy = edges_[i]->center_[1] - tmp_origin[1];
        double dz = edges_[i]->center_[2] - tmp_origin[2];
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
        sum_lx += length*dx;
        sum_ly += length*dy;
        sum_lz += length*dz;
        sum_l += length;
    }
    center_[0] = sum_lx/sum_l + tmp_origin[0];
    center_[1] = sum_ly/sum_l + tmp_origin[1];
    center_[2] = sum_lz/sum_l + tmp_origin[2];

    return 0;
}

bool Polygon::crossBoundary() {
    for (int i = 0; i < edges_.size(); i++) {
        if (edges_[i]->crossBoundary()) {
            return true;
        }
    }

    return false;
}

bool Polygon::checkH() {
    if (edges_.size() != 3) {
        return false;
    }
    for (int j = 0; j < 3; j++) {
        if (!edges_[j]->checkH()) {
            return false;
        }
    }

    return true;
}

int Polygon::shrink(Edge* edge) {
    auto it = find(edges_.begin(), edges_.end(), edge);
    if (it != edges_.end()) {
        edges_.erase(it);
    } else {
        printf("edge %ld not found in polygon %ld\n", edge->id_, id_);
        exit(1);
    }
    if (edges_.size() == 3) {
        for (auto e : edges_) {
            e->triangle_count_ = e->triangle_count_ + 1;
        }
    }
    if (edges_.size() == 2) {
        for (auto e : edges_) {
            e->triangle_count_ = e->triangle_count_ - 1;
        }
    }

    return 0;
}

int Polygon::expand(Edge* edge) {
    edges_.push_back(edge);
    if (edges_.size() == 3) {
        for (auto e : edges_) {
            e->triangle_count_ = e->triangle_count_ + 1;
        }
    }
    if (edges_.size() == 4) {
        for (auto e : edges_) {
            e->triangle_count_ = e->triangle_count_ - 1;
        }
    }

    return 0;
}

int Polygon::logEdges(std::string name) {
    printf("%s %ld\n",name.c_str(), id_);
    printf("%ld",edges_.size());
    for (auto edge : edges_) {
        printf(" %ld",edge->id_);
    }
    printf("\n");

    return 0;
}