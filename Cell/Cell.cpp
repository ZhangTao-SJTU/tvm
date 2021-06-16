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

#include "Cell.h"

using namespace std;

Cell::Cell(Run * run, long int id) {
    run_ = run;
    id_ = id;
    growing_ = false;
    volume_ = 0.;
    pressure_ = 0.;
}

int Cell::updateVolume() {
    // compute direction of polygons
    polygonDirections_.clear();

    // update direction of each polygon
    std::unordered_map<long int, bool> edgeDirections;
    std::vector<Polygon *> tmp_polygons;
    while (tmp_polygons.size() < polygons_.size()) {
        for (auto polygon : polygons_) {
            if (std::find(tmp_polygons.begin(), tmp_polygons.end(), polygon) != tmp_polygons.end()) {
                continue;
            }
            if (tmp_polygons.size() == 0) {
                tmp_polygons.push_back(polygon);
                // no need to adjust the direction of the first polygon
                polygonDirections_[polygon->id_] = true;
                // record the direction of edges of the first polygon
                for (int i = 0; i < polygon->vertices_.size(); i++) {
                    Vertex * v0 = polygon->vertices_[i];
                    Vertex * v1 = polygon->vertices_[(i + 1)%polygon->vertices_.size()];
                    if (v0->id_ < v1->id_) {
                        long int edgeID = v0->id_*run_->count_vertices_ + v1->id_;
                        edgeDirections[edgeID] = true;
                        // true: v0->v1 is from lower vertex id to larger vertex id
                    } else {
                        long int edgeID = v1->id_*run_->count_vertices_ + v0->id_;
                        edgeDirections[edgeID] = false;
                        // false: v0->v1 is from larger vertex id to lower vertex id
                    }
                }
                continue;
            }
            // check if this polygon has any registered edge
            bool registered = false;
            bool oppositeDirection = true;
            for (int i = 0; i < polygon->vertices_.size(); i++) {
                Vertex * v0 = polygon->vertices_[i];
                Vertex * v1 = polygon->vertices_[(i + 1)%polygon->vertices_.size()];
                long int edgeID;
                bool edgeDirection;
                if (v0->id_ < v1->id_) {
                    edgeID = v0->id_*run_->count_vertices_ + v1->id_;
                    edgeDirection = true;
                } else {
                    edgeID = v1->id_*run_->count_vertices_ + v0->id_;
                    edgeDirection = false;
                }
                if (edgeDirections.find(edgeID) != edgeDirections.end()) {
                    registered = true;
                    if (edgeDirection != edgeDirections[edgeID]) {
                        oppositeDirection = true;
                    } else {
                        oppositeDirection = false;
                    }
                    break;
                }
            }
            if (!registered) {
                continue;
            }
            tmp_polygons.push_back(polygon);
            // adjust the direction of the polygon
            polygonDirections_[polygon->id_] = oppositeDirection;
            // record the direction of extra edges of the polygon
            for (int i = 0; i < polygon->vertices_.size(); i++) {
                Vertex * v0 = polygon->vertices_[i];
                Vertex * v1 = polygon->vertices_[(i + 1)%polygon->vertices_.size()];
                long int edgeID;
                bool edgeDirection;
                if (v0->id_ < v1->id_) {
                    edgeID = v0->id_*run_->count_vertices_ + v1->id_;
                    edgeDirection = true;
                } else {
                    edgeID = v1->id_*run_->count_vertices_ + v0->id_;
                    edgeDirection = false;
                }
                if (edgeDirections.find(edgeID) == edgeDirections.end()) {
                    if (!oppositeDirection) {
                        edgeDirection = (!edgeDirection);
                    }
                    edgeDirections[edgeID] = edgeDirection;
                }
            }
        }
        printf("Cell::updateVolume Error: cannot find next polygon in polygons_\n");
        exit(1);
    }
    edgeDirections.clear();

    // compute cell volume
    volume_ = 0.;
    for (auto polygon : polygons_) {
        // the origin is the reference point
        double cc[3];   // the vector pointing from origin to polygon center
        for (int m = 0; m < 3; m++) {
            cc[m] = polygon->center_[m];
        }
        while (cc[0] > run_->Lx_/2.0) {
            cc[0] = cc[0] - run_->Lx_;
        }
        while (cc[0] < (-1.0)*run_->Lx_/2.0) {
            cc[0] = cc[0] + run_->Lx_;
        }
        while (cc[1] > run_->Ly_/2.0) {
            cc[1] = cc[1] - run_->Ly_;
        }
        while (cc[1] < (-1.0)*run_->Ly_/2.0) {
            cc[1] = cc[1] + run_->Ly_;
        }
        for (int i = 0; i < polygon->vertices_.size(); i++) {
            double cv[2][3];   // the vectors pointing from polygon center to edge vertices
            for (int k = 0; k < 2; k++) {
                Vertex * vertex = polygon->vertices_[(i + k)%polygon->vertices_.size()];
                for (int m = 0; m < 3; m++) {
                    cv[k][m] = vertex->position_[m] - polygons_[i]->center_[m];
                }
                while (cv[k][0] > run_->Lx_/2.0) {
                    cv[k][0] = cv[k][0] - run_->Lx_;
                }
                while (cv[k][0] < (-1.0)*run_->Lx_/2.0) {
                    cv[k][0] = cv[k][0] + run_->Lx_;
                }
                while (cv[k][1] > run_->Ly_/2.0) {
                    cv[k][1] = cv[k][1] - run_->Ly_;
                }
                while (cv[k][1] < (-1.0)*run_->Ly_/2.0) {
                    cv[k][1] = cv[k][1] + run_->Ly_;
                }
            }
            // compute the volume of the tetrahedron formed by cell center, polygon center, and edge vertices
            double cP[3];
            double dP = 0.;
            cP[0] = cv[0][1]*cv[1][2] - cv[1][1]*cv[0][2];
            cP[1] = cv[1][0]*cv[0][2] - cv[0][0]*cv[1][2];
            cP[2] = cv[0][0]*cv[1][1] - cv[1][0]*cv[0][1];
            for (int m = 0; m < 3; m++) {
                dP += cc[m]*cP[m];
            }
            if (polygonDirections_[polygon->id_]) {
                volume_ += 1.0/6.0*dP;
            } else {
                volume_ -= 1.0/6.0*dP;
            }
        }
    }

    if (volume_ < 0.) {
        volume_ = fabs(volume_);
        for (auto polygon : polygons_) {
            polygonDirections_[polygon->id_] = (!polygonDirections_[polygon->id_]);
        }
    }

    return 0;
}

int Cell::logPolygons(std::string name) {
    printf("%s %ld\n",name.c_str(), id_);
    printf("%ld",polygons_.size());
    for (auto polygon : polygons_) {
        printf(" %ld",polygon->id_);
    }
    printf("\n");
}