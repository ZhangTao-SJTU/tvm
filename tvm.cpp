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
#include "Run/Run.h"

using namespace std;

int     InitializeAll(Run *);
//int     LoadConf(string filename, Run *);
int     DumpConfigurationVtk(double simulation_time, Run *);
int     DumpCentersVtk(double simulation_time, Run *);

int main(int argc, char *argv[]) {
    Run * run = new Run();
    InitializeAll(run);
    run->start();
//    DumpConfigurationVtk(0., run);
//    DumpCentersVtk(0., run);
    return 0;
}

int InitializeAll(Run * run) {
//    LoadConf("conf", run);
    printf("Initialization start ...\n");

    // generate initial configuration
    double areaCell = run->Lx_*run->Ly_/run->NCell_;
    double unitL = sqrt(areaCell/(3.0/2.0*sqrt(3.0)));
    printf("cell surface area: %6f, unit length: %6f\n", areaCell, unitL);
    // generate center of hexagons
    std::vector< double > centerVerticesX;
    std::vector< double > centerVerticesY;
    int Nx = ceil(run->Lx_/(3.0/2.0*unitL))+1;
    int Ny = ceil(run->Ly_/(sqrt(3.0)*unitL))+1;
    for (int i = 0; i < Nx; i++) {
        for (int j = (-Ny); j < Ny; j++) {
            double x = (0.1 + i*(3.0/2.0))*unitL;
            double y = (0.1 + i*(sqrt(3.0)/2.0) + j*(sqrt(3.0)))*unitL;
            if (x > 0 && x < run->Lx_ && y > 0 && y < run->Ly_) {
                centerVerticesX.push_back(x);
                centerVerticesY.push_back(y);
            }
        }
    }
    printf("%ld hexagon centers generated on bottom surface\n", centerVerticesX.size());
    // generate vertices, each center generates two vertices in directions 0 and 1
    for (long int i = 0; i < centerVerticesX.size(); i++) {
        double cx = centerVerticesX[i];
        double cy = centerVerticesY[i];
        Vertex* v1 = new Vertex(run, run->vertices_.size());
        run->vertices_.push_back(v1);
        Vertex* v2 = new Vertex(run, run->vertices_.size());
        run->vertices_.push_back(v2);
        v1->position_[0] = cx + unitL;
        v1->position_[1] = cy;
        v2->position_[0] = cx + 0.5*unitL;
        v2->position_[1] = cy + sqrt(3.0)/2.0*unitL;
        if (v1->position_[0] > run->Lx_) {
            v1->position_[0] = v1->position_[0] - run->Lx_;
        }
        if (v1->position_[1] > run->Ly_) {
            v1->position_[1] = v1->position_[1] - run->Ly_;
        }
        if (v2->position_[0] > run->Lx_) {
            v2->position_[0] = v2->position_[0] - run->Lx_;
        }
        if (v2->position_[1] > run->Ly_) {
            v2->position_[1] = v2->position_[1] - run->Ly_;
        }
    }
    printf("%ld vertices generated on bottom surface\n", run->vertices_.size());
    // generate edges on bottom surface
    for (long int i = 0; i < run->vertices_.size(); i++) {
        for (long int j = i + 1; j < run->vertices_.size(); j++) {
            double dij = run->computeD(run->vertices_[i]->position_,run->vertices_[j]->position_);
//            printf("%6f ",dij);
            // check the distance between any two vertices on bottom surface
            if (fabs(dij - unitL) < 0.1*unitL) {
                Edge * edge = new Edge(run, run->edges_.size());
                run->edges_.push_back(edge);
                edge->vertices_.push_back(run->vertices_[i]);
                edge->vertices_.push_back(run->vertices_[j]);
            }
        }
    }
    printf("%ld edges generated on bottom surface\n", run->edges_.size());
    // generate polygons on bottom surface
    for (long int i = 0; i < centerVerticesX.size(); i++) {
        double cx = centerVerticesX[i];
        double cy = centerVerticesY[i];
        Polygon * polygon = new Polygon(run, run->polygons_.size());
        run->polygons_.push_back(polygon);
        for (long int j = 0; j < run->edges_.size(); j++) {
            Vertex * v1 = run->edges_[j]->vertices_[0];
            Vertex * v2 = run->edges_[j]->vertices_[1];
            // check the distance between the center and two edge vertices
            if (fabs(run->computeD(cx, cy, v1->position_)-unitL) < 0.1*unitL) {
                if (fabs(run->computeD(cx, cy, v2->position_)-unitL) < 0.1*unitL) {
                    polygon->edges_.push_back(run->edges_[j]);
                }
            }
        }
        if (polygon->edges_.size() != 6) {
            printf("Number of edges in polygon %ld: %ld not equal to 6\n",polygon->id_, polygon->edges_.size());
            exit(1);
        }
    }
    printf("%ld polygons generated on bottom surface\n", run->polygons_.size());
    // generate vertices on top surface and store edges in z direction
    std::vector <Edge *> tmp_edges;
    long int NverticesBottom = run->vertices_.size();
    for (long int i = 0; i < NverticesBottom; i++) {
        Vertex * vertex = new Vertex(run, run->vertices_.size());
        run->vertices_.push_back(vertex);
        vertex->position_[0] = run->vertices_[i]->position_[0];
        vertex->position_[1] = run->vertices_[i]->position_[1];
        vertex->position_[2] = run->vertices_[i]->position_[2] + 2.0;
        Edge * edge = new Edge(run, -1);
        edge->vertices_.push_back(run->vertices_[i]);
        edge->vertices_.push_back(vertex);
        tmp_edges.push_back(edge);
    }
    printf("%ld vertices generated on top surface\n", run->vertices_.size()-NverticesBottom);
    // generate edges on top surface
    long int NedgesBottom = run->edges_.size();
    for (long int i = 0; i < NedgesBottom; i++) {
        Edge * edge = new Edge(run, run->edges_.size());
        Vertex * v1 = run->vertices_[run->edges_[i]->vertices_[0]->id_+NverticesBottom];
        Vertex * v2 = run->vertices_[run->edges_[i]->vertices_[1]->id_+NverticesBottom];
        edge->vertices_.push_back(v1);
        edge->vertices_.push_back(v2);
        run->edges_.push_back(edge);
    }
    printf("%ld edges generated on top surface\n", run->edges_.size()-NedgesBottom);
    // add edges in z direction
    for (long int i = 0; i < tmp_edges.size(); i++) {
        tmp_edges[i]->id_ = run->edges_.size();
        run->edges_.push_back(tmp_edges[i]);
    }
    printf("%ld edges generated in z direction\n", tmp_edges.size());
    // generate polygons on top surface
    long int NpolygonsBottom = run->polygons_.size();
    for (long int i = 0; i < NpolygonsBottom; i++) {
        Polygon * polygon = new Polygon(run, run->polygons_.size());
        run->polygons_.push_back(polygon);
        for (int j = 0; j < run->polygons_[i]->edges_.size(); j++) {
            long int edgeID = run->polygons_[i]->edges_[j]->id_ + NedgesBottom;
            polygon->edges_.push_back(run->edges_[edgeID]);
        }
    }
    printf("%ld polygons generated on top surface\n", run->polygons_.size() - NpolygonsBottom);
    // generate polygons in z direction
    run->updateVertexEdges();
    for (long int i = 0; i < NedgesBottom; i++) {
        Polygon * polygon = new Polygon(run, run->polygons_.size());
        run->polygons_.push_back(polygon);
        polygon->edges_.push_back(run->edges_[i]);
        polygon->edges_.push_back(run->edges_[NedgesBottom+i]);
        for (int j = 0; j < 2; j++) {
            Vertex * vertex = polygon->edges_[0]->vertices_[j];
            for (int k = 0; k < vertex->edges_.size(); k++) {
                Edge * edge = vertex->edges_[k];
                if (fabs(edge->vertices_[1]->position_[2] - edge->vertices_[0]->position_[2]) > 1.0) {
                    polygon->edges_.push_back(edge);
                    break;
                }
            }
        }
    }
    printf("%ld polygons generated in z direction\n", run->polygons_.size() - 2*NpolygonsBottom);
    // generate cells
    for (long int i = 0; i < NpolygonsBottom; i++) {
        Cell * cell = new Cell(run, run->cells_.size());
        run->cells_.push_back(cell);
        cell->polygons_.push_back(run->polygons_[i]);
        cell->polygons_.push_back(run->polygons_[NpolygonsBottom+i]);
        for (int j = 0; j < run->polygons_[i]->edges_.size(); j++) {
            long int edgeID = run->polygons_[i]->edges_[j]->id_;
            cell->polygons_.push_back(run->polygons_[2*NpolygonsBottom+edgeID]);
        }
    }
    printf("%ld cells generated\n", run->cells_.size());

//    // add initial curved surface
//    for (long int i = 0; i < run->vertices_.size(); i++) {
//        Vertex * vertex = run->vertices_[i];
//        double x = vertex->position_[0];
//        double y = vertex->position_[1];
//        double z = vertex->position_[2];
//        vertex->position_[2] = z + run->Aic_*cos(2.0*M_PI/run->Lx_*x)*cos(2.0*M_PI/run->Ly_*y);
//    }

    run->updatePolygonVertices();
    run->updateVertexEdges();

    run->updateGeoinfo();
    double sum_volume = 0.;
    for (long int i = 0; i < run->cells_.size(); i++) {
        sum_volume += run->cells_[i]->volume_;
    }
    printf("Total volume: %6f\n", sum_volume);

    // initialize volume object
    run->volume_ = new Volume(run);
    // initialize interface object
    run->interface_ = new Interface(run);

    return 0;
}

int DumpConfigurationVtk(double simulation_time, Run * run) {
    //////////////////////////////////////////////////////////////////////////////////////
    stringstream filename;
    filename << setw(10) << setfill('0') << (long int)(floor(simulation_time)) << ".sample.vtk";
    ofstream out(filename.str().c_str());
    if (!out.is_open()) {
        cout << "Error opening output file " << filename.str().c_str() << endl;
        exit(1);
    }
    out << "# vtk DataFile Version 2.0" << endl;
    out << "polydata" << endl;
    out << "ASCII" << endl;
    out << "DATASET POLYDATA" << endl;
    out << "POINTS " << run->vertices_.size() << " double" << endl;
    for (long int i = 0; i < run->vertices_.size(); i++) {
            out << right << setw(12) << scientific << setprecision(5) << run->vertices_[i]->position_[0];
            out << " " << right << setw(12) << scientific << setprecision(5) << run->vertices_[i]->position_[1];
            out << " " << right << setw(12) << scientific << setprecision(5) << run->vertices_[i]->position_[2];
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

    long int Npolygons = 0;
    long int NpolygonVertices = 0;
    for (long int i = 0; i < run->polygons_.size(); i++) {
        if (!run->polygons_[i]->crossBoundary()) {
            Npolygons++;
            NpolygonVertices += run->polygons_[i]->vertices_.size();
        }
    }
    out << "POLYGONS " << Npolygons << " " << Npolygons + NpolygonVertices << endl;
    for (long int i = 0; i < run->polygons_.size(); i++) {
        if (!run->polygons_[i]->crossBoundary()) {
            out << left << setw(6) << run->polygons_[i]->vertices_.size();
            for (int j = 0; j < run->polygons_[i]->vertices_.size(); j++) {
                out << " " << left << setw(6) << run->polygons_[i]->vertices_[j]->id_;
            }
            out << endl;
        }
    }
    out << endl;

    out.close();

    return 0;
}

int DumpCentersVtk(double simulation_time, Run * run) {
    //////////////////////////////////////////////////////////////////////////////////////
    stringstream filename;
    filename << setw(10) << setfill('0') << (long int)(floor(simulation_time)) << ".centers.vtk";
    ofstream out(filename.str().c_str());
    if (!out.is_open()) {
        cout << "Error opening output file " << filename.str().c_str() << endl;
        exit(1);
    }
    out << "# vtk DataFile Version 2.0" << endl;
    out << "polydata" << endl;
    out << "ASCII" << endl;
    out << "DATASET POLYDATA" << endl;
    out << "POINTS " << run->cells_.size() << " double" << endl;
    for (long int i = 0; i < run->cells_.size(); i++) {
        out << right << setw(12) << scientific << setprecision(5) << run->cells_[i]->center_[0];
        out << " " << right << setw(12) << scientific << setprecision(5) << run->cells_[i]->center_[1];
        out << " " << right << setw(12) << scientific << setprecision(5) << run->cells_[i]->center_[2];
        out << endl;
    }
    out << endl;

    out.close();

    return 0;
}

//int LoadConf(string filename, Run * run) {
//    ifstream conf(filename.c_str());
//
//    if (!conf.is_open()) {
//        cout << "Error opening conf file" << endl;
//        exit(1);
//    }
//
//    string buffer;
//    string delimiter = " ";
//    size_t pos = 0;
//    vector<string> tokens;
//    vector<vector<string>> lines;
//
//    int time_written = 0;
//    int dx_written = 0;
//    int rc_written = 0;
//    int dump_written = 0;
//    int log_screen_written = 0;
//    int temperature_written = 0;
//    int fiber_written = 0;
//
//    while (getline(conf, buffer))
//    {
//        pos = buffer.find((char)13);
//        if (pos != string::npos) {
//            buffer = buffer.substr(0, pos);
//        }
//        if ((buffer.length() == 0) || (buffer[0] == '#')) continue;
//
//        tokens.clear();
//        while ((pos = buffer.find(delimiter)) != string::npos) {
//            string token = buffer.substr(0, pos);
//            if (token.length() > 0) {
//                tokens.push_back(token);
//            }
//            buffer.erase(0, pos + delimiter.length());
//        }
//        if (buffer.length() > 0) {
//            tokens.push_back(buffer);
//        }
//        lines.push_back(tokens);
//    }
//
//    for (int i = 0; i < lines.size(); i++) {
//        tokens = lines[i];
//        if (tokens[0] == "time") {
//            if (tokens.size() != 3) {
//                cerr << "conf file error: ";
//                for (int j = 0; j < tokens.size(); j++) {
//                    cerr << tokens[j] << " ";
//                }
//                cerr << endl;
//                exit(1);
//            }
//            run->t_start_ = atof(tokens[1].c_str());
//            run->t_end_ = atof(tokens[2].c_str());
//            time_written = 1;
//            cout << "time: " << run->t_start_ << " ~ " << run->t_end_ << endl;
//        }
//        else if (tokens[0] == "dx") {
//            if (tokens.size() != 2) {
//                cerr << "conf file error: ";
//                for (int j = 0; j < tokens.size(); j++) {
//                    cerr << tokens[j] << " ";
//                }
//                cerr << endl;
//                exit(1);
//            }
//            run->dx_ = atof(tokens[1].c_str());
//            run->updatePairsOnSurfaceClose_period_ = 1.0;
//            dx_written = 1;
//            cout << "dx: " << run->dx_ << endl;
//        }
//        else if (tokens[0] == "rc") {
//            if (tokens.size() != 5) {
//                cerr << "conf file error: ";
//                for (int j = 0; j < tokens.size(); j++) {
//                    cerr << tokens[j] << " ";
//                }
//                cerr << endl;
//                exit(1);
//            }
//            run->contact_->Rc_ = atof(tokens[1].c_str());
//            run->contact_->rc_ = atof(tokens[2].c_str());
//            run->contact_->D_ = atof(tokens[3].c_str());
//            run->contact_->a_ = atof(tokens[4].c_str());
//            rc_written = 1;
//            cout << "rc: " << run->contact_->Rc_ << " " << run->contact_->rc_ << endl;
//        }
//        else if (tokens[0] == "dump") {
//            if (tokens.size() != 3) {
//                cerr << "conf file error: ";
//                for (int j = 0; j < tokens.size(); j++) {
//                    cerr << tokens[j] << " ";
//                }
//                cerr << endl;
//                exit(1);
//            }
//            run->dump_period_ = atof(tokens[2].c_str());
//            dump_written = 1;
//            cout << "dump: " << tokens[1] << " " << run->dump_period_ << endl;
//        }
//        else if (tokens[0] == "log") {
//            if (tokens.size() != 3) {
//                cerr << "conf file error: ";
//                for (int j = 0; j < tokens.size(); j++) {
//                    cerr << tokens[j] << " ";
//                }
//                cerr << endl;
//                exit(1);
//            }
//            run->log_period_ = atof(tokens[2].c_str());
//            log_screen_written = 1;
//            cout << "log: " << tokens[1] << " " << run->log_period_ << endl;
//        }
//        else if (tokens[0] == "temperature") {
//            if (tokens.size() != 2) {
//                cerr << "conf file error: ";
//                for (int j = 0; j < tokens.size(); j++) {
//                    cerr << tokens[j] << " ";
//                }
//                cerr << endl;
//                exit(1);
//            }
//            run->temperature_ = atof(tokens[1].c_str());
//            temperature_written = 1;
//            cout << "temperature: " << run->temperature_ << endl;
//        }
//        else if (tokens[0] == "fiber") {
//            if (tokens.size() != 4) {
//                cerr << "conf file error: ";
//                for (int j = 0; j < tokens.size(); j++) {
//                    cerr << tokens[j] << " ";
//                }
//                cerr << endl;
//                exit(1);
//            }
//            run->Bcilia_ = atof(tokens[1].c_str());
//            run->Kbend_ = atof(tokens[2].c_str());
//            run->fiber_equi_temperature_ = atof(tokens[3].c_str());
//            fiber_written = 1;
//            cout << "fiber: " << run->Bcilia_ << " " << run->Kbend_ << " " << run->fiber_equi_temperature_ << endl;
//        }
//        else
//        {
//            cerr << "conf file error: ";
//            for (int j = 0; j < tokens.size(); j++) {
//                cerr << tokens[j] << " ";
//            }
//            cerr << endl;
//            exit(1);
//        }
//    }
//
//    if (conf.bad() || !conf.eof()) {
//        cout << "Error reading file [" << filename << "]" << endl;
//        exit(1);
//    }
//
//    if (time_written == 0) {
//        cout << "conf file error: missing time" << endl;
//        exit(1);
//    }
//
//    if (dx_written == 0) {
//        cout << "conf file error: dx" << endl;
//        exit(1);
//    }
//
//    if (rc_written == 0) {
//        cout << "conf file error: rc" << endl;
//        exit(1);
//    }
//
//    if (dump_written == 0) {
//        cout << "conf file error: missing dump" << endl;
//        exit(1);
//    }
//
//    if (log_screen_written == 0) {
//        cout << "conf file error: missing log screen" << endl;
//        exit(1);
//    }
//
//    if (temperature_written == 0) {
//        cout << "conf file error: missing temperature" << endl;
//        exit(1);
//    }
//
//    if (fiber_written == 0) {
//        cout << "conf file error: missing fiber" << endl;
//        exit(1);
//    }
//
//    conf.close();
//
//    return 0;
//}
