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
#include <random>
#include "Run/Run.h"

using namespace std;

int     InitializeAll(Run *);
int     LoadConf(string filename, Run *);

int main(int argc, char *argv[]) {
    Run * run = new Run();
    InitializeAll(run);
    run->updatePolygonVertices();
//    run->dumpConfigurationVtk();

//    for (auto cell : run->cells_) {
//        cell->updateVolume();
//        printf("%f\n", cell->volume_);
//    }
//    run->cells_[200]->updateVolume();
//    printf("%f\n", run->cells_[200]->volume_);
//    run->reconnection_->Lth_ = 0.5;
//    run->reconnection_->I_H(run->edges_[2700], true);
//    run->reconnection_->Lth_ = 2.0;
//    run->reconnection_->H_I(run->polygons_[run->polygons_.size()-1], true);

    run->start();

    return 0;
}

int InitializeAll(Run * run) {
    printf("Initialization start ...\n");

    // load initial configuration
    ifstream topofile("sample.topo");
    if (!topofile.is_open()) {
        cout << "Error opening sample topo file" << endl;
        exit(1);
    }

    char buffer[256];
    long int tmp_long;
    long int tmp_id;
    int tmp_int;

//    topofile.getline(buffer, 100);

    // initialize vertices_ objects
    topofile >> buffer;
    cout << "Reading " << buffer << endl;
    topofile >> tmp_long;
    cout << "Number of vertices: " << tmp_long << endl;
    for (long int i = 0; i < tmp_long; i++) {
        Vertex * vertex = new Vertex(run, i);
        topofile >> tmp_id;
        for (int j = 0; j < 3; j++) {
            topofile >> vertex->position_[j];
        }
        run->vertices_.push_back(vertex);
    }

    // initialize edges_ objects
    topofile.getline(buffer, 100);
    topofile >> buffer;
    cout << "Reading " << buffer << endl;
    topofile >> tmp_long;
    cout << "Number of edges: " << tmp_long << endl;
    for (long int i = 0; i < tmp_long; i++) {
        Edge * edge = new Edge(run, i);
        topofile >> tmp_int;
        for (int j = 0; j < tmp_int; j++) {
            topofile >> tmp_id;
            edge->vertices_.push_back(run->vertices_[tmp_id]);
        }
        run->edges_.push_back(edge);
    }

    // initialize polygons_ objects
    topofile.getline(buffer, 100);
    topofile >> buffer;
    cout << "Reading " << buffer << endl;
    topofile >> tmp_long;
    cout << "Number of polygons: " << tmp_long << endl;
    for (long int i = 0; i < tmp_long; i++) {
        Polygon * polygon = new Polygon(run, i);
        topofile >> tmp_int;
        for (int j = 0; j < tmp_int; j++) {
            topofile >> tmp_id;
            polygon->edges_.push_back(run->edges_[tmp_id]);
        }
        run->polygons_.push_back(polygon);
    }

    // initialize cells_ objects
    topofile.getline(buffer, 100);
    topofile >> buffer;
    cout << "Reading " << buffer << endl;
    topofile >> tmp_long;
    cout << "Number of cells: " << tmp_long << endl;
    for (long int i = 0; i < tmp_long; i++) {
        Cell * cell = new Cell(run, i);
        topofile >> tmp_int;
        for (int j = 0; j < tmp_int; j++) {
            topofile >> tmp_id;
            cell->polygons_.push_back(run->polygons_[tmp_id]);
        }
        run->cells_.push_back(cell);
    }

    run->count_vertices_ = run->vertices_.size();
    run->count_edges_ = run->edges_.size();
    run->count_polygons_ = run->polygons_.size();
    run->count_cells_ = run->cells_.size();

    // initialize volume object
    run->volume_ = new Volume(run);
    // initialize interface object
    run->interface_ = new Interface(run);
    // initialize reconnection object
    run->reconnection_ = new Reconnection(run);

    LoadConf("conf", run);

    // update geometry and topology information
    run->updateGeoinfo();
    run->updateVertexCells();
    run->volume_->updatePolygonDirections();

    return 0;
}

//int DumpCentersVtk(double simulation_time, Run * run) {
//    //////////////////////////////////////////////////////////////////////////////////////
//    stringstream filename;
//    filename << setw(10) << setfill('0') << (long int)(floor(simulation_time)) << ".centers.vtk";
//    ofstream out(filename.str().c_str());
//    if (!out.is_open()) {
//        cout << "Error opening output file " << filename.str().c_str() << endl;
//        exit(1);
//    }
//    out << "# vtk DataFile Version 2.0" << endl;
//    out << "polydata" << endl;
//    out << "ASCII" << endl;
//    out << "DATASET POLYDATA" << endl;
//    out << "POINTS " << run->cells_.size() << " double" << endl;
//    for (long int i = 0; i < run->cells_.size(); i++) {
//        out << right << setw(12) << scientific << setprecision(5) << run->cells_[i]->center_[0];
//        out << " " << right << setw(12) << scientific << setprecision(5) << run->cells_[i]->center_[1];
//        out << " " << right << setw(12) << scientific << setprecision(5) << run->cells_[i]->center_[2];
//        out << endl;
//    }
//    out << endl;
//
//    out.close();
//
//    return 0;
//}

int LoadConf(string filename, Run * run) {
    ifstream conf(filename.c_str());

    if (!conf.is_open()) {
        cout << "Error opening conf file" << endl;
        exit(1);
    }

    string buffer;
    string delimiter = " ";
    size_t pos = 0;
    vector<string> tokens;
    vector<vector<string>> lines;

    int time_written = 0;
    int dump_written = 0;
    int log_screen_written = 0;
    int s0_written = 0;
    int Lth_written = 0;

    while (getline(conf, buffer))
    {
        pos = buffer.find((char)13);
        if (pos != string::npos) {
            buffer = buffer.substr(0, pos);
        }
        if ((buffer.length() == 0) || (buffer[0] == '#')) continue;

        tokens.clear();
        while ((pos = buffer.find(delimiter)) != string::npos) {
            string token = buffer.substr(0, pos);
            if (token.length() > 0) {
                tokens.push_back(token);
            }
            buffer.erase(0, pos + delimiter.length());
        }
        if (buffer.length() > 0) {
            tokens.push_back(buffer);
        }
        lines.push_back(tokens);
    }

    for (int i = 0; i < lines.size(); i++) {
        tokens = lines[i];
        if (tokens[0] == "time") {
            if (tokens.size() != 4) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->t_start_ = atof(tokens[1].c_str());
            run->t_end_ = atof(tokens[2].c_str());
            run->dt_ = atof(tokens[3].c_str());
            run->dtr_ = 10*run->dt_;
            time_written = 1;
            cout << "time: " << run->t_start_ << " ~ " << run->t_end_ << " ~ " << run->dt_ << " ~ " << run->dtr_ << endl;
        }
        else if (tokens[0] == "dump") {
            if (tokens.size() != 3) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            if (tokens[1] == "vtk") {
                run->dump_period_ = atof(tokens[2].c_str());
                dump_written = 1;
                cout << "dump: " << tokens[1] << " " << run->dump_period_ << endl;
            }
        }
        else if (tokens[0] == "log") {
            if (tokens.size() != 2) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->log_period_ = atof(tokens[1].c_str());
            log_screen_written = 1;
            cout << "log: " << run->log_period_ << endl;
        }
        else if (tokens[0] == "s0") {
            if (tokens.size() != 2) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->interface_->s0_ = atof(tokens[1].c_str());
            s0_written = 1;
            cout << "s0: " << run->interface_->s0_ << endl;
        }
        else if (tokens[0] == "Lth") {
            if (tokens.size() != 2) {
                cerr << "conf file error: ";
                for (int j = 0; j < tokens.size(); j++) {
                    cerr << tokens[j] << " ";
                }
                cerr << endl;
                exit(1);
            }
            run->reconnection_->Lth_ = atof(tokens[1].c_str());
            Lth_written = 1;
            cout << "Lth: " << run->reconnection_->Lth_ << endl;
        }
        else {
            cerr << "conf file error: ";
            for (int j = 0; j < tokens.size(); j++) {
                cerr << tokens[j] << " ";
            }
            cerr << endl;
            exit(1);
        }
    }

    if (conf.bad() || !conf.eof()) {
        cout << "Error reading file [" << filename << "]" << endl;
        exit(1);
    }

    if (time_written == 0) {
        cout << "conf file error: missing time" << endl;
        exit(1);
    }

    if (dump_written == 0) {
        cout << "conf file error: missing dump" << endl;
        exit(1);
    }

    if (log_screen_written == 0) {
        cout << "conf file error: missing log screen" << endl;
        exit(1);
    }

    if (s0_written == 0) {
        cout << "conf file error: s0" << endl;
        exit(1);
    }

    if (Lth_written == 0) {
        cout << "conf file error: Lth" << endl;
        exit(1);
    }

    conf.close();

    return 0;
}
