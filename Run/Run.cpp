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
    long int simulation_time_counter = 0;
    double t_roundError = 0.01*dt_;
    auto start = chrono::steady_clock::now();

    printf("\nSimulation Start ...\n");
    printf("Real time elapsed: Rte\n");
    printf("Time        ");
    printf("Rte   ");
    printf("Volume   ");
    printf("I->H     ");
    printf("H->I     ");
    printf("Linker     ");
//    printf("NTS       ");
//    printf("ETE       ");
    printf("E_V      ");
    printf("E_S      ");
    printf("ECM_str  ");
    printf("ECM_bend ");
    printf("Link_str ");
    printf("Energy      \n");

    while (simulation_time_ < t_end_ + t_roundError) {
        // update geometry information
        updateGeoinfo();
        // update volumeForces
        volume_->updateForces();
        // update interfaceForces
        interface_->updateForces();
        // update fiberForces
        fiberElasticity_->updateForces();
        // update linkForces
        fiberLink_->updateForces();
//        // update contactForces
//        contact_->updateForces();

        // update velocities
        updateVerticesVelocity();
        updateNodesVelocity();

        if (simulation_time_ > 500. - t_roundError && simulation_time_ < 500. + t_roundError) {
            initializeLinks();
        }


        // log to screen
        if (simulation_time_ - t_start_ + t_roundError  > count_log_ * log_period_) {
            volume_->updateEnergy();
            interface_->updateEnergy();
            fiberElasticity_->updateEnergy();
            fiberLink_->updateEnergy();
            printf("%-12.2f%-6.1f%-9.1f%-9ld%-9ld%-5d+%-5ld%-9.1f%-9.1f%-9.1f%-9.1f%-9.1f%-9.1f\n", simulation_time_,
                   (chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - start).count())/1.0e6,
                   volume_->totalVolume_,
                   reconnection_->count_IH_,
                   reconnection_->count_HI_,
                   fiberLink_->N_, fiberLink_->nFormedLinks_,
                   volume_->energy_,
                   interface_->energy_,
                   fiberElasticity_->stretchEnergy_,
                   fiberElasticity_->bendEnergy_,
                   fiberLink_->stretchEnergy_,
                   volume_->energy_+interface_->energy_+fiberElasticity_->stretchEnergy_+fiberElasticity_->bendEnergy_+fiberLink_->stretchEnergy_);
            start = chrono::steady_clock::now();
            reconnection_->count_IH_ = 0;
            reconnection_->count_HI_ = 0;
            count_log_++;
        }
        // dump
        if (simulation_time_ - t_start_ + t_roundError > count_dump_ * dump_period_) {
            if (simulation_time_ > (-1.)*t_roundError) {
                dumpTopo();
                dumpCellCenter();
                dumpCellShapeIndex();
                dumpCellVolume();
//                dumpReconnection();
                dumpConfigurationVtk();
                dumpLinkInfo();
                dumpSpheroidShape();
                // dumpVertexForce();
            }
//            dumpCellCenter();
//            dumpCellShapeIndex();
            count_dump_++;
        }

        // Euler dynamics
        updateVerticesPosition();
        updateNodesPosition();

        // reconnect
        if (simulation_time_ - t_start_ + t_roundError > count_reconnect_ * dtr_) {
            reconnection_->start();
            count_reconnect_++;

            if (simulation_time_ > 500. + t_roundError) {

//                cout << "check 1" << endl;
                checkLinks();
//                cout << "check 2" << endl;
//                updateLinks();
            }
        }

        // decrease the rest length of link spring from 2.0 to 0.2 in time interval 5000, every 200 time intervals.
        if (simulation_time_ > fiberLink_->shrinkStartTime_ - t_roundError  && fiberLink_->l0_ > fiberLink_->l1_) {
            if ((simulation_time_ - fiberLink_->shrinkStartTime_ + t_roundError ) > simulation_time_counter * 10.){
                fiberLink_->l0_ = fiberLink_->l0_ - fiberLink_->shrinkSpeed_;
                simulation_time_counter = simulation_time_counter + 1;
            }
        }

        simulation_time_ += dt_;
    }

//    for (long int i = 0; i < cells_.size(); i++) {
//        printf("%f\n", cells_[i]->volume_);
//    }

    return 0;
}

int     Run::updateVerticesVelocity() {
    for (auto vertex : vertices_) {
        for (int m = 0; m < 3; m++) {
            vertex->velocity_[m] = mu_ * (vertex->volumeForce_[m] + vertex->interfaceForce_[m] + vertex->linkForce_[m]);
        }
    }
    // remove drift velocity
//    if (true) {
//        double averageVelocity[3] = {0., 0., 0.};
//        long int nRealVertices = 0;
//        for (auto vertex : vertices_) {
//            if (vertex->type_ > 0) {
//                for (int m = 0; m < 3; m++) {
//                    averageVelocity[m] = averageVelocity[m] + vertex->velocity_[m];
//                }
//                nRealVertices++;
//            }
//        }
//        for (int m = 0; m < 3; m++) {
//            averageVelocity[m] = averageVelocity[m] / nRealVertices;
//        }
//        for (auto vertex : vertices_) {
//            if (vertex->type_ > 0) {
//                for (int m = 0; m < 3; m++) {
//                    vertex->velocity_[m] = vertex->velocity_[m] - averageVelocity[m];
//                }
//            }
//        }
//    }

    return 0;
}

int     Run::updateVerticesPosition() {
    std::default_random_engine generator(std::random_device{}());
    std::normal_distribution<double> ndist(0., 1.);
    double cR = sqrt(2.0*mu_*kB_*temperature_*dt_);
    double cRe = sqrt(2.0*mu_*kB_*temperature_*dte_);
    for (auto vertex : vertices_) {
        if (vertex->type_ > 0) {
            if (vertex->pull_) {
                for (int m = 0; m < 3; m++) {
                    vertex->position_[m] = vertex->position_[m] + vertex->velocity_[m] * dt_;
                }
            } else {
                for (int m = 0; m < 3; m++) {
                    vertex->position_[m] =
                            vertex->position_[m] + vertex->velocity_[m] * dt_ + cR * ndist(generator);
                }
            }
        } else {
            for (int m = 0; m < 3; m++) {
                if (vertex->velocity_[m] * dte_ > 0.1) {
//                    cout << vertex->velocity_[m] * dte_ << endl;
                    vertex->velocity_[m] = 0.1/dte_;
                }
                if (vertex->velocity_[m] * dte_ < -0.1) {
//                    cout << vertex->velocity_[m] * dte_ << endl;
                    vertex->velocity_[m] = -0.1/dte_;
                }
                vertex->position_[m] =
                        vertex->position_[m] + vertex->velocity_[m] * dte_ + cRe * ndist(generator);
            }
        }
        box_->resetPosition(vertex->position_);
    }

    return 0;
}

int     Run::updateNodesVelocity() {
    for (auto node : nodes_) {
        for (int m = 0; m < 3; m++) {
            node->velocity_[m] = mu_ * (node->stretchForce_[m] + node->bendForce_[m] + node->linkForce_[m] + node->contactForce_[m]);
        }
    }

    return 0;
}

int     Run::updateNodesPosition() {
    for (auto node : nodes_) {
        for (int m = 0; m < 3; m++) {
            node->position_[m] = node->position_[m] + node->velocity_[m] * dt_;
        }
        box_->resetPosition(node->position_);
    }

    return 0;
}

int     Run::initializeLinks() {
    updateGeoinfo();
    updatePolygonCells();
    links_.clear();

    boundaryPolygons_.clear();
    for (auto polygon : polygons_) {
        if (polygon->type_ == 2) {
            boundaryPolygons_.push_back(polygon);
        }
    }
    std::vector<Node *> boundaryNodes_;
    std::vector<Node *> tmpNodes_;
    for (auto node : nodes_) {
        if (node->link_) {
            boundaryNodes_.push_back(node);
        }
    }
    if (boundaryNodes_.size() < fiberLink_->N_) {
        fiberLink_->N_ = boundaryNodes_.size();
    }
    std::default_random_engine generator(std::random_device{}());
    std::shuffle(boundaryNodes_.begin(), boundaryNodes_.end(), generator);

    // set pairs of polygon and node
    int countNode = 0;
    for (auto node : boundaryNodes_) {
        Polygon * pairPolygon = NULL;
        double minDist = 1000.;
        for (auto polygon : boundaryPolygons_) {
            if (polygon->link_) {
                continue;
            }
            bool skipFlag = false;
            for (auto cell : polygon->cells_) {
                if (cell->type_ > 0 && cell->link_) {
                    skipFlag = true;
                    break;
                }
            }
            if (skipFlag) {
                continue;
            }
            polygon->updateCenter();
            polygon->updateArea();
            if (polygon->area_ < fiberLink_->minArea_) {
                continue;
            }
            double dx[3];
            dx[0] = node->position_[0] - polygon->center_[0];
            dx[1] = node->position_[1] - polygon->center_[1];
            dx[2] = node->position_[2] - polygon->center_[2];
//            box_->resetDistance(dx);
            double dist = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
            if (dist > 3.) {
                continue;
            }
            if (dist < minDist) {
                minDist = dist;
                pairPolygon = polygon;
            }
        }
        if (pairPolygon != NULL) {
            Link * link = new Link(this, pairPolygon, node);
            links_.push_back(link);
            pairPolygon->link_ = true;
            for (auto cell : pairPolygon->cells_) {
                if (cell->type_ > 0) {
                    cell->link_ = true;
                }
            }
            countNode++;
        }
        if(countNode == fiberLink_->N_) {
            break;
        }
    }

//    int count = 0;
//    for (auto link : links_) {
//        cout << count << " " << link->polygon_->id_ << " " << link->node_->id_ << endl;
//        count++;
//    }

    return 0;
}

int     Run::checkLinks() {
    updatePolygonCells();
    std::vector<Link *> tmpLinks = links_;
    links_.clear();
    std::vector<Node *> freeNodes;

    for (auto link : tmpLinks) {
        Polygon * polygon = link->polygon_;
        Node * node = link->node_;
//        cout << run_->simulation_time_ << " " << polygon->id_ << endl;

        if (polygon->id_ > 0) {
            if (polygon->checkShape()) {
                links_.push_back(link);
            } else {
                polygon->link_ = false;
                for (auto cell : polygon->cells_) {
                    if (cell->type_ > 0) {
                        cell->link_ = false;
                    }
                }
                delete link;
                freeNodes.push_back(node);
            }
        } else {
            delete link;
            freeNodes.push_back(node);
        }
    }

    boundaryPolygons_.clear();
    for (auto polygon : polygons_) {
        if (polygon->type_ == 2) {
            boundaryPolygons_.push_back(polygon);
        }
    }

    // add pairs of polygon and free node
    for (auto node : freeNodes) {
        Polygon * pairPolygon = NULL;
        double minDiff = 1000.;
        for (auto polygon : boundaryPolygons_) {
            if (polygon->link_) {
                continue;
            }
            bool skipFlag = false;
            for (auto cell : polygon->cells_) {
                if (cell->type_ > 0 && cell->link_) {
                    skipFlag = true;
                    break;
                }
            }
            if (skipFlag) {
                continue;
            }
            if (!polygon->checkShape()) {
                continue;
            }
            double dx[3];
            dx[0] = node->position_[0] - polygon->center_[0];
            dx[1] = node->position_[1] - polygon->center_[1];
            dx[2] = node->position_[2] - polygon->center_[2];
            box_->resetDistance(dx);
            double diff = fabs(sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]) - fiberLink_->l0_);
//            if (diff > 1.0) {
//                continue;
//            }
            if (diff < minDiff) {
                minDiff = diff;
                pairPolygon = polygon;
            }
        }
        if (pairPolygon != NULL) {
            Link * link = new Link(this, pairPolygon, node);
            links_.push_back(link);
            pairPolygon->link_ = true;
            for (auto cell : pairPolygon->cells_) {
                if (cell->type_ > 0) {
                    cell->link_ = true;
                }
            }
            fiberLink_->nFormedLinks_++;
        }
    }

    return 0;
}

int     Run::updateLinks() {
    std::vector<Link *> tmpLinks = links_;
    links_.clear();
    std::vector<Node *> freeNodes;

    for (auto link : tmpLinks) {
        Polygon * polygon = link->polygon_;
        Node * node = link->node_;
//        cout << run_->simulation_time_ << " " << polygon->id_ << endl;

        double dx[3];
        dx[0] = node->position_[0] - polygon->center_[0];
        dx[1] = node->position_[1] - polygon->center_[1];
        dx[2] = node->position_[2] - polygon->center_[2];
        box_->resetDistance(dx);
//        double l = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
        if (polygon->checkShape()) {
            links_.push_back(link);
        } else {
            delete link;
            polygon->link_ = false;
            freeNodes.push_back(node);
        }
    }

    boundaryPolygons_.clear();
    for (auto polygon : polygons_) {
        if (polygon->type_ == 2) {
            boundaryPolygons_.push_back(polygon);
        }
    }

    // add pairs of polygon and free node
    for (auto node : freeNodes) {
        Polygon * pairPolygon = NULL;
        double minDiff = 1000.;
        for (auto polygon : boundaryPolygons_) {
            if (polygon->link_) {
                continue;
            }
            if (!polygon->checkShape()) {
                continue;
            }
            double dx[3];
            dx[0] = node->position_[0] - polygon->center_[0];
            dx[1] = node->position_[1] - polygon->center_[1];
            dx[2] = node->position_[2] - polygon->center_[2];
            box_->resetDistance(dx);
            double diff = fabs(sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]) - fiberLink_->l0_);
//            if (diff > 1.0) {
//                continue;
//            }
            if (diff < minDiff) {
                minDiff = diff;
                pairPolygon = polygon;
            }
        }
        if (pairPolygon != NULL) {
            Link * link = new Link(this, pairPolygon, node);
            links_.push_back(link);
            pairPolygon->link_ = true;
        }
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
    for (auto cell : cells_) {
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

int     Run::updateEdgeCells() {
    for (auto edge : edges_) {
        edge->cells_.clear();
    }
    for (auto cell : cells_) {
        for (auto polygon : cell->polygons_) {
            for (auto edge : polygon->edges_) {
                if (std::find(edge->cells_.begin(), edge->cells_.end(), cell) == edge->cells_.end()) {
                    // new cell to be added
                    edge->cells_.push_back(cell);
                }
            }
        }
    }

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

int     Run::updateEmptyCells() {
    for (auto polygon : polygons_) {
        polygon->type_ = 0;
    }
    for (auto vertex : vertices_) {
        vertex->type_ = 0;
    }
    for (auto edge : edges_) {
        edge->type_ = 0;
    }
    for (auto cell : cells_) {
        if (cell->type_ == 1) {
            for (auto polygon : cell->polygons_) {
                polygon->type_ = 1;
            }
        }
    }
    for (auto polygon : polygons_) {
        if (polygon->type_ == 1) {
            if (polygon->cells_[0]->type_ == 0 || polygon->cells_[1]->type_ == 0) {
                polygon->type_ = 2;
            }
//            if (polygon->cells_[0]->type_ < 0 || polygon->cells_[1]->type_ < 0) {
//                polygon->type_ = 2;
//            }
        }
    }
    for (auto edge : edges_) {
        int n = 0;
        for (auto cell : edge->cells_) {
            if (cell->type_ > 0) {
                n++;
            }
        }
        edge->type_ = n;
    }
    for (auto vertex : vertices_) {
        int n = 0;
        for (auto cell : vertex->cells_) {
            if (cell->type_ > 0) {
                n++;
            }
        }
        vertex->type_ = n;
    }

    // check isolated cell
    for (auto cell : cells_) {
        if (cell->type_ == 1) {
            bool isolated = true;
            for (auto polygon : cell->polygons_) {
                if (polygon->type_ != 2) {
                    isolated = false;
                    break;
                }
            }
            if (isolated) {
                cell->type_ = 0;
                for (auto polygon : cell->polygons_) {
                    polygon->type_ = 0;
                    for (auto edge: polygon->edges_) {
                        edge->type_ = 0;
                        for (auto vertex : edge->vertices_) {
                            vertex->type_ = 0;
                        }
                    }
                }
            }
        }
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
//    delete polygon;
    polygon->id_ = -1;

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
    updatePolygonVertices();
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

    long int Npolygons = 0;
    long int NpolygonVertices = 0;
    for (long int i = 0; i < polygons_.size(); i++) {
        if (polygons_[i]->type_ > 0 && !polygons_[i]->crossBoundary()) {
            Npolygons++;
            NpolygonVertices += polygons_[i]->vertices_.size();
        }
    }
    out << "POLYGONS " << Npolygons << " " << Npolygons + NpolygonVertices << endl;
    for (long int i = 0; i < polygons_.size(); i++) {
        if (polygons_[i]->type_ > 0 && !polygons_[i]->crossBoundary()) {
            out << left << setw(6) << polygons_[i]->vertices_.size();
            for (int j = 0; j < polygons_[i]->vertices_.size(); j++) {
                out << " " << left << setw(6) << polygons_[i]->vertices_[j]->dumpID_;
            }
            out << endl;
        }
    }
    out << endl;

    out << "CELL_DATA " << Npolygons << endl;
    out << "SCALARS type int 1" << endl;
    out << "LOOKUP_TABLE default" << endl;
    for (long int i = 0; i < polygons_.size(); i++) {
        if (polygons_[i]->type_ > 0 && !polygons_[i]->crossBoundary()) {
            out << left << setw(6) << polygons_[i]->type_ << endl;
        }
    }
    out << endl;

    out.close();

    //////////////////////////////////////////////////////////////////////////////////////
    stringstream filename2;
    filename2 << setw(7) << setfill('0') << (long int)(floor(simulation_time_+0.01*dt_)) << ".extra.vtk";
    ofstream out2(filename2.str().c_str());
    if (!out2.is_open()) {
        cout << "Error opening output file " << filename2.str().c_str() << endl;
        exit(1);
    }
    out2 << "# vtk DataFile Version 2.0" << endl;
    out2 << "polydata" << endl;
    out2 << "ASCII" << endl;
    out2 << "DATASET POLYDATA" << endl;

    long int Nvertices = 0;
    for (auto vertex : vertices_) {
        if (vertex->type_ > 0) {
            // reset vertex id for dumping polygons
            vertex->dumpID_ = Nvertices;
            Nvertices++;
        }
    }
    out2 << "POINTS " << Nvertices << " double" << endl;
    for (auto vertex : vertices_) {
        if (vertex->type_ > 0) {
            out2 << right << setw(12) << scientific << setprecision(5) << vertex->position_[0];
            out2 << " " << right << setw(12) << scientific << setprecision(5) << vertex->position_[1];
            out2 << " " << right << setw(12) << scientific << setprecision(5) << vertex->position_[2];
            out2 << endl;
        }
    }
    out2 << endl;

    long int Nedges = 0;
    long int NedgeVertices = 0;
    for (auto edge : edges_) {
        if (edge->type_ >= 2 && !edge->crossBoundary()) {
            Nedges++;
            NedgeVertices += 2;
        }
    }
    out2 << "LINES " << Nedges << " " << Nedges + NedgeVertices << endl;
    for (auto edge : edges_) {
        if (edge->type_ >= 2 && !edge->crossBoundary()) {
            out2 << left << setw(6) << 2;
            out2 << " " << left << setw(6) << edge->vertices_[0]->dumpID_;
            out2 << " " << left << setw(6) << edge->vertices_[1]->dumpID_;
            out2 << endl;
        }
    }
    out2 << endl;

    out2 << "POINT_DATA " << Nvertices << endl;
    out2 << "SCALARS type int 1" << endl;
    out2 << "LOOKUP_TABLE default" << endl;
    for (auto vertex : vertices_) {
        if (vertex->type_ > 0) {
            out2 << left << setw(6) << vertex->type_ + 10 << endl;
        }
    }
    out2 << endl;

    out2.close();

    //////////////////////////////////////////////////////////////////////////////////////
    stringstream filename1;
    filename1 << setw(7) << setfill('0') << (long int)(floor(simulation_time_+0.01*dt_)) << ".empty.vtk";
    ofstream out1(filename1.str().c_str());
    if (!out1.is_open()) {
        cout << "Error opening output file " << filename1.str().c_str() << endl;
        exit(1);
    }
    out1 << "# vtk DataFile Version 2.0" << endl;
    out1 << "polydata" << endl;
    out1 << "ASCII" << endl;
    out1 << "DATASET POLYDATA" << endl;
    out1 << "POINTS " << vertices_.size() << " double" << endl;
    for (long int i = 0; i < vertices_.size(); i++) {
        // reset vertex id for dumping polygons
        vertices_[i]->dumpID_ = i;
        out1 << right << setw(12) << scientific << setprecision(5) << vertices_[i]->position_[0];
        out1 << " " << right << setw(12) << scientific << setprecision(5) << vertices_[i]->position_[1];
        out1 << " " << right << setw(12) << scientific << setprecision(5) << vertices_[i]->position_[2];
        out1 << endl;
    }
    out1 << endl;

    Npolygons = 0;
    NpolygonVertices = 0;
    for (long int i = 0; i < polygons_.size(); i++) {
        if (polygons_[i]->type_ != 1 && !polygons_[i]->crossBoundary()) {
            Npolygons++;
            NpolygonVertices += polygons_[i]->vertices_.size();
        }
    }
    out1 << "POLYGONS " << Npolygons << " " << Npolygons + NpolygonVertices << endl;
    for (long int i = 0; i < polygons_.size(); i++) {
        if (polygons_[i]->type_ != 1 && !polygons_[i]->crossBoundary()) {
            out1 << left << setw(6) << polygons_[i]->vertices_.size();
            for (int j = 0; j < polygons_[i]->vertices_.size(); j++) {
                out1 << " " << left << setw(6) << polygons_[i]->vertices_[j]->dumpID_;
            }
            out1 << endl;
        }
    }
    out1 << endl;

    out1 << "CELL_DATA " << Npolygons << endl;
    out1 << "SCALARS type int 1" << endl;
    out1 << "LOOKUP_TABLE default" << endl;
    for (long int i = 0; i < polygons_.size(); i++) {
        if (polygons_[i]->type_ != 1 && !polygons_[i]->crossBoundary()) {
            out1 << left << setw(6) << polygons_[i]->type_ << endl;
        }
    }
    out1 << endl;

    out1.close();

    //////////////////////////////////////////////////////////////////////////////////////
    stringstream ECMfilename;
    ECMfilename << setw(7) << setfill('0') << (long int)(floor(simulation_time_+0.01*dt_)) << ".ECM.vtk";
    ofstream ECMout(ECMfilename.str().c_str());
    if (!ECMout.is_open()) {
        cout << "Error opening output file " << ECMfilename.str().c_str() << endl;
        exit(1);
    }
    ECMout << "# vtk DataFile Version 2.0" << endl;
    ECMout << "polydata" << endl;
    ECMout << "ASCII" << endl;
    ECMout << "DATASET POLYDATA" << endl;
    ECMout << "POINTS " << nodes_.size() << " double" << endl;
    for (long int i = 0; i < nodes_.size(); i++) {
        // reset vertex id for dumping polygons
        nodes_[i]->dumpID_ = i;
        ECMout << right << setw(12) << scientific << setprecision(5) << nodes_[i]->position_[0];
        ECMout << " " << right << setw(12) << scientific << setprecision(5) << nodes_[i]->position_[1];
        ECMout << " " << right << setw(12) << scientific << setprecision(5) << nodes_[i]->position_[2];
        ECMout << endl;
    }
    ECMout << endl;

    long int Nfibers = 0;
    for (auto fiber : fibers_) {
        for (int i = 0; i < fiber->nodes_.size() - 1; i++) {
            Node * ni = fiber->nodes_[i];
            Node * nj = fiber->nodes_[i + 1];
            if (!box_->crossBoundary(ni->position_, nj->position_)) {
                Nfibers++;
            }
        }
    }
    ECMout << "LINES " << Nfibers << " " << 3*Nfibers << endl;
    for (auto fiber : fibers_) {
        for (int i = 0; i < fiber->nodes_.size() - 1; i++) {
            Node * ni = fiber->nodes_[i];
            Node * nj = fiber->nodes_[i + 1];
            if (!box_->crossBoundary(ni->position_, nj->position_)) {
                ECMout << left << setw(6) << 2;
                ECMout << " " << left << setw(6) << ni->dumpID_;
                ECMout << " " << left << setw(6) << nj->dumpID_;
                ECMout << endl;
            }
        }
    }
    ECMout << endl;

    ECMout << "CELL_DATA " << Nfibers << endl;
    ECMout << "SCALARS strain double" << endl;
    ECMout << "LOOKUP_TABLE default" << endl;
    for (auto fiber : fibers_) {
        for (int i = 0; i < fiber->nodes_.size() - 1; i++) {
            Node * ni = fiber->nodes_[i];
            Node * nj = fiber->nodes_[i + 1];
            if (!box_->crossBoundary(ni->position_, nj->position_)) {
                double lt = 0.;
                for (int m = 0; m < 3; m++) {
                    lt = lt+pow ((ni->position_[m]-nj->position_[m]), 2.0) ;
                }
                double strain = (sqrt(lt) - fiberElasticity_->l0_) / fiberElasticity_->l0_;
                ECMout << left << scientific << setprecision(7) << strain << endl;
            }
        }
    }

    ECMout << endl;

    ECMout.close();

    //////////////////////////////////////////////////////////////////////////////////////
    stringstream Linkfilename;
    Linkfilename << setw(7) << setfill('0') << (long int)(floor(simulation_time_+0.01*dt_)) << ".link.vtk";
    ofstream Linkout(Linkfilename.str().c_str());
    if (!Linkout.is_open()) {
        cout << "Error opening output file " << Linkfilename.str().c_str() << endl;
        exit(1);
    }
    Linkout << "# vtk DataFile Version 2.0" << endl;
    Linkout << "polydata" << endl;
    Linkout << "ASCII" << endl;
    Linkout << "DATASET POLYDATA" << endl;
    Linkout << "POINTS " << links_.size()*2 << " double" << endl;
    for (long int i = 0; i < links_.size(); i++) {
        Link * link = links_[i];
        // reset vertex id for dumping polygons
        Linkout << right << setw(12) << scientific << setprecision(5) << link->polygon_->center_[0];
        Linkout << " " << right << setw(12) << scientific << setprecision(5) << link->polygon_->center_[1];
        Linkout << " " << right << setw(12) << scientific << setprecision(5) << link->polygon_->center_[2];
        Linkout << endl;
        Linkout << right << setw(12) << scientific << setprecision(5) << link->node_->position_[0];
        Linkout << " " << right << setw(12) << scientific << setprecision(5) << link->node_->position_[1];
        Linkout << " " << right << setw(12) << scientific << setprecision(5) << link->node_->position_[2];
        Linkout << endl;
    }
    Linkout << endl;

    Linkout << "LINES " << links_.size() << " " << 3*links_.size() << endl;
    for (long int i = 0; i < links_.size(); i++) {
        Link * link = links_[i];
        Linkout << left << setw(6) << 2;
        Linkout << " " << left << setw(6) << 2*i;
        Linkout << " " << left << setw(6) << 2*i + 1;
        Linkout << endl;
    }
    Linkout << endl;

    Linkout << "CELL_DATA " << links_.size() << endl;
    Linkout << "SCALARS strain double" << endl;
    Linkout << "LOOKUP_TABLE default" << endl;
    for (long int i = 0; i < links_.size(); i++) {
        Link * link = links_[i];
        double linklength = 0.;
        for (int m = 0; m < 3; m++) {
            linklength = linklength + pow ((link->polygon_->center_[m]-link->node_->position_[m]), 2.0) ;
        }
        double currentL0 = max(fiberLink_->l0_, fiberLink_->l1_);
        double strain = (sqrt(linklength) - currentL0) / currentL0;
        Linkout << left << scientific << setprecision(7) << strain << endl;}
    Linkout << endl;

    Linkout.close();

    return 0;
}

int     Run::dumpLinkInfo(){
    stringstream Linkfilename;
    Linkfilename << setw(7) << setfill('0') << (long int)(floor(simulation_time_+0.01*dt_)) << ".linkInfo.txt";
    ofstream Linkout(Linkfilename.str().c_str());
    if (!Linkout.is_open()) {
        cout << "Error opening output file " << Linkfilename.str().c_str() << endl;
        exit(1);
    }
    Linkout << left << setw(12) <<"PolygonID";
    Linkout << left << setw(12) <<"NodeID"<<endl;

    for (long int i = 0; i < links_.size(); i++) {
        Link * link = links_[i];
        Linkout << left << setw(12) << link->polygon_->id_;
        Linkout << left << setw(12) << link->node_->id_;
        Linkout << endl;
    }
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
        if (cell->type_ <= 0) {
            continue;
        }
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
    int countCells = 0;
    for (auto cell : cells_) {
        if (cell->type_ <= 0) {
            continue;
        }
        averageShapeIndex += cell->shapeIndex_;
        countCells++;
    }
    averageShapeIndex /= countCells;
    out << " " << left << setw(12) << averageShapeIndex;
    out << endl;

    for (auto cell : cells_) {
        if (cell->type_ <= 0) {
            continue;
        }
        out << left << setw(6) << cell->id_;
        out << " " << cell->shapeIndex_;
        out << endl;
    }
    out << endl;

    out.close();

    return 0;
}

int     Run::dumpCellVolume() {
    stringstream filename;
    filename << "cellVolume.txt";
    ofstream out(filename.str().c_str(), std::ios_base::app);
    if (!out.is_open()) {
        cout << "Error opening output file " << filename.str().c_str() << endl;
        exit(1);
    }
    out << "time ";
    out << left << setw(12) << simulation_time_;
    double totalVolume = 0.;
    int countCells = 0;
    for (auto cell : cells_) {
        if (cell->type_ <= 0) {
            continue;
        }
        totalVolume += cell->volume_;
        countCells++;
    }
    out << " " << left << setw(12) << totalVolume;
    out << endl;

    for (auto cell : cells_) {
        if (cell->type_ <= 0) {
            continue;
        }
        out << left << setw(6) << cell->id_;
        out << " " << cell->volume_;
        out << endl;
    }
    out << endl;

    out.close();

    return 0;
}

int     Run::dumpVertexForce() {
    stringstream filename;
    filename << setw(7) << setfill('0') << (long int)(floor(simulation_time_+0.01*dt_)) << ".vertexForce.txt";
    //filename << "vertexForce.txt";
    ofstream out(filename.str().c_str(), std::ios_base::app);
    if (!out.is_open()) {
        cout << "Error opening output file " << filename.str().c_str() << endl;
        exit(1);
    }
    //out << "time ";
    //out << left << setw(12) << simulation_time_;
    //out << endl;


    for (auto vertex : vertices_) {
        out << left << setw(6) << vertex->id_;
        for (int m = 0; m < 3; m++) {
            out << " " << right << setw(12) << scientific << setprecision(5) << (vertex->volumeForce_[m] + vertex->interfaceForce_[m] + vertex->linkForce_[m]);
        }
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
    out << " box";
    out << " " << left << setw(12) << box_->size_[0];
    out << " " << left << setw(12) << box_->size_[1];
    out << " " << left << setw(12) << box_->size_[2];
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
    out << left << setw(12) << cells_.size() - 1;
    out << endl;
    for (auto cell : cells_) {
        if (cell->type_ < 0) {
            continue;
        }
        out << left << setw(6) << cell->id_;
        for (auto polygon : cell->polygons_) {
            out << " " << right << setw(12) << scientific << setprecision(5) << polygon->id_;
        }
        out << " " << right << setw(6) << cell->type_;
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

int     Run::dumpSpheroidShape() {
    stringstream filename;
    filename << "spheroidShape.txt";

    ofstream out(filename.str().c_str(), std::ios_base::app);
    if (!out.is_open()) {
        cout << "Error opening output file " << filename.str().c_str() << endl;
        exit(1);
    }

    out << left << setw(12) << simulation_time_;
    double area = 0.;
    for (auto polygon : polygons_) {
        if (polygon->type_ == 2) {
            area += polygon->area_;
        }
    }
    out << " " << left << setw(12) << area;        
    
    double volume = 0.;
    for (auto cell : cells_) {
        if (cell->type_ <= 0) {
            continue;
        }
        volume += cell->volume_;
    }
    out << " " << left << setw(12) << volume;
    out << " " << left << setw(12) << area * pow(volume, (-1.0)*2.0/3.0);
    out << endl;
    
    out.close();

    return 0;
}