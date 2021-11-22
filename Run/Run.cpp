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
#include <random>
#include "Run.h"

using namespace std;

Run::Run() {
//    dt_ = 0.001;
//    dtr_ = 10*dt_;
//    dump_period_ = 10000*dt_;
//    log_period_ = 100*dt_;
//    t_start_ = 0.;
//    t_end_ = 10000.;
    mu_ = 1.0;
    kB_ = 1.0;
//    temperature_ = 1.0e-5;
    NCell_ = 512;
}

int Run::start() {
    count_reconnect_ = 0;
    count_dump_ = 0;
    count_log_ = 0;
    simulation_time_ = t_start_;
    double t_roundError = 0.01*dt_;
    auto start = chrono::steady_clock::now();

    printf("\nSimulation Start ...\n");
    printf("Real time elapsed: Rte\n");
    printf("Time        ");
    printf("Rte         ");
    printf("Volume      ");
    printf("I->H        ");
    printf("H->I        ");
    printf("E_volume    ");
    printf("E_interface ");
    printf("Energy      \n");

    while (simulation_time_ < t_end_ + t_roundError) {
        // update geometry information
        updateGeoinfo();
        // update volumeForces
        volume_->updateForces();
        // update interfaceForces
        interface_->updateForces();
        // update radialForces
        if (simulation_time_ - t_start_ + t_roundError > 500.) {
            if (simulation_time_ - t_start_ < 500. + t_roundError) {
                assignPullingPolygons();
            }
            updatePullingForces();
        }
        // update velocities
        updateVerticesVelocity();

        // log to screen
        if (simulation_time_ - t_start_ + t_roundError  > count_log_ * log_period_) {
            volume_->updateEnergy();
            interface_->updateEnergy();
            printf("%-12.2f%-12.3f%-12.3f%-12ld%-12ld%-12.6f%-12.6f%-12.6f\n", simulation_time_,
                   (chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - start).count())/1.0e6,
                   volume_->totalVolume_,
                   reconnection_->count_IH_,
                   reconnection_->count_HI_,
                   volume_->energy_,
                   interface_->energy_,
                   volume_->energy_+interface_->energy_);
            start = chrono::steady_clock::now();
            reconnection_->count_IH_ = 0;
            reconnection_->count_HI_ = 0;
            count_log_++;
        }
        // dump
        if (simulation_time_ - t_start_ + t_roundError > count_dump_ * dump_period_) {
            if (simulation_time_ > (-0.01)*dt_) {
//                dumpTopo();
//                dumpCellCenter();
//                dumpCellShapeIndex();
//                dumpReconnection();
                dumpConfigurationVtk();
            }
//            dumpCellCenter();
//            dumpCellShapeIndex();
            count_dump_++;
        }

        // Euler dynamics
        updateVerticesPosition();

        // reconnect
        if (simulation_time_ - t_start_ + t_roundError > count_reconnect_ * dtr_) {
            reconnection_->start();
            count_reconnect_++;
        }

        simulation_time_ += dt_;
    }

//    for (long int i = 0; i < cells_.size(); i++) {
//        printf("%f\n", cells_[i]->volume_);
//    }

    return 0;
}

int     Run::assignPullingPolygons() {
    for (auto cell : emptyCells_) {
        for (auto polygon: cell->polygons_) {
            double xa = 0.;
            double ya = 0.;
            double za = 0.;
            for (auto vertex: polygon->vertices_) {
                xa += vertex->position_[0];
                ya += vertex->position_[1];
                za += vertex->position_[2];
            }
            xa /= polygon->vertices_.size();
            ya /= polygon->vertices_.size();
            za /= polygon->vertices_.size();
            if (fabs(xa) > 4.0 && sqrt(ya*ya + za*za) < 2.0) {
                polygon->pull_ = true;
            }
        }
    }

    for (auto polygon : polygons_) {
        if (polygon->pull_) {
            for (auto vertex: polygon->vertices_) {
                vertex->pull_ = true;
            }
        }
    }

    return 0;
}

int     Run::updatePullingForces() {
    // reset all volumeForce values in vertices
    for (auto vertex : vertices_) {
        for (int m = 0; m < 3; m++) {
            vertex->pullingForce_[m] = 0.;
        }
    }

    // update vertices on surface
    std::vector<Vertex *> verticesSurface;
    for (auto polygon : polygons_) {
        if (polygon->pull_) {
            for (auto vertex: polygon->vertices_) {
                if (std::find(verticesSurface.begin(), verticesSurface.end(), vertex) == verticesSurface.end()) {
                    // new vertex to be added
                    verticesSurface.push_back(vertex);
                }
            }
        }
    }

    // update radialForces for each vertex on surface
    for (auto vertex : verticesSurface) {
        double x = vertex->position_[0];
        if (fabs(x) >= pullxMax_) {
            vertex->pullingForce_[0] = 0.;
            vertex->volumeForce_[0] = 0.;
            vertex->volumeForce_[1] = 0.;
            vertex->volumeForce_[2] = 0.;
            vertex->interfaceForce_[0] = 0.;
            vertex->interfaceForce_[1] = 0.;
            vertex->interfaceForce_[2] = 0.;
            continue;
        }
        if (x < 0) {
            vertex->pullingForce_[0] = (-1.0)*pullForce_;
        } else {
            vertex->pullingForce_[0] = pullForce_;
        }
        vertex->volumeForce_[0] = 0.;
        vertex->volumeForce_[1] = 0.;
        vertex->volumeForce_[2] = 0.;
        vertex->interfaceForce_[0] = 0.;
        vertex->interfaceForce_[1] = 0.;
        vertex->interfaceForce_[2] = 0.;
    }

    return 0;
}

int     Run::updateVerticesVelocity() {
    for (auto vertex : vertices_) {
        for (int m = 0; m < 3; m++) {
            vertex->velocity_[m] = mu_ * (vertex->volumeForce_[m] + vertex->interfaceForce_[m] + vertex->pullingForce_[m]);
        }
    }
    // remove drift velocity
    double averageVelocity[3] = {0., 0., 0.};
    for (auto vertex : vertices_) {
        for (int m = 0; m < 3; m++) {
            averageVelocity[m] = averageVelocity[m] + vertex->velocity_[m];
        }
    }
    for (int m = 0; m < 3; m++) {
        averageVelocity[m] = averageVelocity[m] / vertices_.size();
    }
    for (auto vertex : vertices_) {
        for (int m = 0; m < 3; m++) {
            vertex->velocity_[m] = vertex->velocity_[m] - averageVelocity[m];
        }
    }

    return 0;
}

int     Run::updateVerticesPosition() {
    std::default_random_engine generator(std::random_device{}());
    std::normal_distribution<double> ndist(0., 1.);
    double cR = sqrt(2.0*mu_*kB_*temperature_*dt_);
    for (long int i = 0; i < vertices_.size(); i++) {
        if (vertices_[i]->pull_) {
            for (int m = 0; m < 3; m++) {
                vertices_[i]->position_[m] = vertices_[i]->position_[m] + vertices_[i]->velocity_[m] * dt_;
            }
        } else {
            for (int m = 0; m < 3; m++) {
                vertices_[i]->position_[m] = vertices_[i]->position_[m] + vertices_[i]->velocity_[m] * dt_ + cR*ndist(generator);
            }
        }
        box_->resetPosition(vertices_[i]->position_);
    }

    return 0;
}

int     Run::updatePolygonVertices() {
    // update vertices in polygon
    for (long int i = 0; i < polygons_.size(); i++) {
        polygons_[i]->updateVertices();
    }

    return 0;
}

int     Run::updateCellVertices() {
    // update vertices in cell
    updatePolygonVertices();
    for (auto cell : cells_) {
        cell->vertices_.clear();
        for (auto polygon : cell->polygons_) {
            for (auto vertex : polygon->vertices_) {
                if (std::find(cell->vertices_.begin(), cell->vertices_.end(), vertex) == cell->vertices_.end()) {
                    // new vertex to be added
                    cell->vertices_.push_back(vertex);
                }
            }
        }
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

int     Run::updateVertexCells() {
    for (auto vertex : vertices_) {
        vertex->cells_.clear();
    }
    std::vector<Cell *> tmp_cells = cells_;
    for (auto cell : emptyCells_) {
        tmp_cells.push_back(cell);
    }
    for (auto cell : tmp_cells) {
        for (auto polygon : cell->polygons_) {
            for (auto edge : polygon->edges_) {
                for (auto vertex : edge->vertices_) {
                    if (std::find(vertex->cells_.begin(), vertex->cells_.end(), cell) == vertex->cells_.end()) {
                        // new cell to be added
                        vertex->cells_.push_back(cell);
                    }
                }
            }
        }
    }

//    for (long int i = 0; i < vertices_.size(); i++) {
//        printf("%d\n", vertices_[i]->cells_.size());
//    }

    return 0;
}

int     Run::updatePolygonCells() {
    for (long int i = 0; i < polygons_.size(); i++) {
        polygons_[i]->cells_.clear();
    }
    for (long int i = 0; i < cells_.size(); i++) {
        for (int j = 0; j < cells_[i]->polygons_.size(); j++) {
            cells_[i]->polygons_[j]->cells_.push_back(cells_[i]);
        }
    }
    for (long int i = 0; i < emptyCells_.size(); i++) {
        for (int j = 0; j < emptyCells_[i]->polygons_.size(); j++) {
            emptyCells_[i]->polygons_[j]->cells_.push_back(emptyCells_[i]);
        }
    }

    return 0;
}

int     Run::updateCellShapeIndex() {
    for (auto cell : cells_) {
        double area = 0.;
        for (auto polygon : cell->polygons_) {
            area += polygon->area_;
        }
        cell->shapeIndex_ = area * pow(cell->volume_, (-1.0)*2.0/3.0);
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

    return 0;
}

int     Run::deleteVertex(Vertex * vertex) {
    auto it = find(vertices_.begin(), vertices_.end(), vertex);
    if (it != vertices_.end()) {
//        int index = it - vertices_.begin();
        vertices_.erase(it);
    } else {
        printf("vertex %ld not found in vertices_\n", vertex->id_);
        exit(1);
    }
    delete vertex;

    return 0;
}

int     Run::deleteEdge(Edge * edge) {
    auto it = find(edges_.begin(), edges_.end(), edge);
    if (it != edges_.end()) {
        edges_[it-edges_.begin()]->markToDelete_ = true;
//        edges_.erase(it);
    } else {
        printf("edge %ld not found in edges_\n", edge->id_);
        exit(1);
    }
//    delete edge;

    return 0;
}

int     Run::deletePolygon(Polygon * polygon) {
    auto it = find(polygons_.begin(), polygons_.end(), polygon);
    if (it != polygons_.end()) {
//        int index = it - vertices_.begin();
        polygons_.erase(it);
    } else {
        printf("polygon %ld not found in polygons_\n", polygon->id_);
        exit(1);
    }
    delete polygon;

    return 0;
}

Edge *  Run::addEdge(Vertex * v0, Vertex * v1) {
    Edge * edge = new Edge(this, count_edges_);
    count_edges_ += 1;
    edges_.push_back(edge);
    if (v0->id_ < v1->id_) {
        edge->vertices_.push_back(v0);
        edge->vertices_.push_back(v1);
    } else {
        edge->vertices_.push_back(v1);
        edge->vertices_.push_back(v0);
    }
    edge->update();
    edge->candidate_ = false;

    return edge;
}

int Run::dumpConfigurationVtk() {
    //////////////////////////////////////////////////////////////////////////////////////
    stringstream filename;
    filename << setw(7) << setfill('0') << (long int)(floor(simulation_time_+0.01*dt_)) << ".sample.vtk";
    ofstream out(filename.str().c_str());
    if (!out.is_open()) {
        cout << "Error opening output file " << filename.str().c_str() << endl;
        exit(1);
    }
    out << "# vtk DataFile Version 2.0" << endl;
    out << "polydata" << endl;
    out << "ASCII" << endl;
    out << "DATASET POLYDATA" << endl;
    out << "POINTS " << vertices_.size() << " double" << endl;
    for (long int i = 0; i < vertices_.size(); i++) {
        // reset vertex id for dumping polygons
        vertices_[i]->dumpID_ = i;
        out << right << setw(12) << scientific << setprecision(5) << vertices_[i]->position_[0];
        out << " " << right << setw(12) << scientific << setprecision(5) << vertices_[i]->position_[1];
        out << " " << right << setw(12) << scientific << setprecision(5) << vertices_[i]->position_[2];
        out << endl;
    }
    out << endl;

//    long int Nedges = 0;
//    for (long int i = 0; i < run->edges_.size(); i++) {
//        if (!run->edges_[i]->crossBoundary()) {
//            Nedges++;
//        }
//    }
//    out << "LINES " << Nedges << " " << 3*Nedges << endl;
//    for (long int i = 0; i < run->edges_.size(); i++) {
//        if (!run->edges_[i]->crossBoundary()) {
//            out << left << setw(6) << 2;
//            for (int j = 0; j < run->edges_[i]->vertices_.size(); j++) {
//                out << " " << left << setw(6) << run->edges_[i]->vertices_[j]->id_;
//            }
//            out << endl;
//        }
//    }
//    out << endl;

    updatePolygonVertices();
    long int Npolygons = 0;
    long int NpolygonVertices = 0;
    for (long int i = 0; i < polygons_.size(); i++) {
        if (!polygons_[i]->crossBoundary()) {
            Npolygons++;
            NpolygonVertices += polygons_[i]->vertices_.size();
        }
    }
    out << "POLYGONS " << Npolygons << " " << Npolygons + NpolygonVertices << endl;
    for (long int i = 0; i < polygons_.size(); i++) {
        if (!polygons_[i]->crossBoundary()) {
            out << left << setw(6) << polygons_[i]->vertices_.size();
            for (int j = 0; j < polygons_[i]->vertices_.size(); j++) {
                out << " " << left << setw(6) << polygons_[i]->vertices_[j]->dumpID_;
            }
            out << endl;
        }
    }
    out << endl;

//    updatePolygonVolumeRatio();
//    out << "CELL_DATA " << Npolygons << endl;
//    out << "SCALARS volumeRatio double 1" << endl;
//    out << "LOOKUP_TABLE default" << endl;
//    for (long int i = 0; i < polygons_.size(); i++) {
//        if (!polygons_[i]->crossBoundary()) {
//            out << left << setw(6) << polygons_[i]->dumpVolumeRatio_ << endl;
//        }
//    }
//    out << endl;

    out.close();

    return 0;
}

int     Run::dumpCellCenter() {
    updateCellVertices();
    stringstream filename;
    filename << "cellCenter.txt";
    ofstream out(filename.str().c_str(), std::ios_base::app);
    if (!out.is_open()) {
        cout << "Error opening output file " << filename.str().c_str() << endl;
        exit(1);
    }
    out << "time ";
    out << left << setw(12) << simulation_time_;
    out << endl;

    for (auto cell : cells_) {
        double center[3] = {0., 0., 0.};
        double reference[3];
        for (int m = 0; m < 3; m++) {
            reference[m] = cell->vertices_[0]->position_[m];
        }
        for (auto vertex : cell->vertices_) {
            double dx[3];
            for (int m = 0; m < 3; m++) {
                dx[m] = (vertex->position_[m] - reference[m]);
            }
            box_->resetDistance(dx);
            for (int m = 0; m < 3; m++) {
                center[m] = center[m] + dx[m];
            }
        }
        for (int m = 0; m < 3; m++) {
            center[m] = center[m]/cell->vertices_.size() + reference[m];
        }
        box_->resetPosition(center);
        out << left << setw(6) << cell->id_;
        out << " " << right << setw(12) << scientific << setprecision(5) << center[0];
        out << " " << right << setw(12) << scientific << setprecision(5) << center[1];
        out << " " << right << setw(12) << scientific << setprecision(5) << center[2];
        out << endl;
    }
    out << endl;

    out.close();

    return 0;
}

int     Run::dumpCellShapeIndex() {
    updateCellShapeIndex();
    stringstream filename;
    filename << "cellShapeIndex.txt";
    ofstream out(filename.str().c_str(), std::ios_base::app);
    if (!out.is_open()) {
        cout << "Error opening output file " << filename.str().c_str() << endl;
        exit(1);
    }
    out << "time ";
    out << left << setw(12) << simulation_time_;
    double averageShapeIndex = 0.;
    for (auto cell : cells_) {
        averageShapeIndex += cell->shapeIndex_;
    }
    averageShapeIndex /= cells_.size();
    out << " " << left << setw(12) << averageShapeIndex;
    out << endl;

    for (auto cell : cells_) {
        out << left << setw(6) << cell->id_;
        out << " " << cell->shapeIndex_;
        out << endl;
    }
    out << endl;

    out.close();

    return 0;
}

int     Run::dumpTopo() {
    stringstream filename;
    filename << "topo.txt";
//    ofstream out(filename.str().c_str(), std::ios::binary | std::ios_base::app);
    ofstream out(filename.str().c_str(), std::ios_base::app);
    if (!out.is_open()) {
        cout << "Error opening output file " << filename.str().c_str() << endl;
        exit(1);
    }
    out << "time ";
    out << left << setw(12) << simulation_time_;
    out << endl;

    out << "vertices ";
    out << left << setw(12) << vertices_.size();
    out << endl;
    for (auto vertex : vertices_) {
        out << left << setw(6) << vertex->id_;
        out << " " << right << setw(12) << scientific << setprecision(5) << vertex->position_[0];
        out << " " << right << setw(12) << scientific << setprecision(5) << vertex->position_[1];
        out << " " << right << setw(12) << scientific << setprecision(5) << vertex->position_[2];
        out << endl;
    }

    out << "edges ";
    out << left << setw(12) << edges_.size();
    out << endl;
    for (auto edge : edges_) {
        out << left << setw(6) << edge->id_;
        for (auto vertex : edge->vertices_) {
            out << " " << right << setw(12) << scientific << setprecision(5) << vertex->id_;
        }
        out << endl;
    }

    out << "polygons ";
    out << left << setw(12) << polygons_.size();
    out << endl;
    for (auto polygon : polygons_) {
        out << left << setw(6) << polygon->id_;
        for (auto edge : polygon->edges_) {
            out << " " << right << setw(12) << scientific << setprecision(5) << edge->id_;
        }
        out << endl;
    }

    out << "cells ";
    out << left << setw(12) << cells_.size();
    out << endl;
    for (auto cell : cells_) {
        out << left << setw(6) << cell->id_;
        for (auto polygon : cell->polygons_) {
            out << " " << right << setw(12) << scientific << setprecision(5) << polygon->id_;
        }
        out << endl;
    }

    out << endl;
    out.close();

    return 0;
}

int     Run::dumpReconnection() {
    stringstream filename;
    filename << "reconnections.txt";
    ofstream out(filename.str().c_str(), std::ios_base::app);
    if (!out.is_open()) {
        cout << "Error opening output file " << filename.str().c_str() << endl;
        exit(1);
    }

    out << verboseReconnection_.str();
    verboseReconnection_.str("");

    out.close();

    return 0;
}