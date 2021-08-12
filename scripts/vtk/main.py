#!/usr/bin/env python3
'''
generate vtk files from simulation output topo.txt file
author: Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn, May 2021
corresponding author: Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu
'''

import os
import numpy as np

class Vertex:
    def __init__(self, id, xyz):
        self.id_ = id
        self.dumpID_ = id
        self.position_ = xyz

class Edge:
    def __init__(self, vertices):
        self.vertices_ = vertices
        self.dump_ = True

class Polygon:
    def __init__(self, edges):
        self.edges_ = edges
        self.vertices_ = []
        self.dump_ = True

def main(begin_num, end_num):
    for runid in range(begin_num, end_num+1):
        run(runid)

def run(runid):
    L = 8.0
    runDir = "{:06d}".format(runid)
    readFilePath = os.path.join(runDir, "topoFixed.txt")
    if not os.path.exists(readFilePath):
        readFilePath = os.path.join(runDir, "topo.txt")
    with open(readFilePath, "r") as file:
        verticesFlag = False
        edgesFlag = False
        polygonsFlag = False
        for line in file:
            if len(line) <= 1:
                continue
            lineSplit = line.split()
            if lineSplit[0] == "time":
                timestamp = float(lineSplit[1])
                vertices = {}
                edges = {}
                polygons = {}
                continue
            if lineSplit[0] == "vertices":
                verticesFlag = True
                continue
            if lineSplit[0] == "edges":
                verticesFlag = False
                edgesFlag = True
                continue
            if lineSplit[0] == "polygons":
                edgesFlag = False
                polygonsFlag = True
                continue
            if lineSplit[0] == "cells":
                polygonsFlag = False
                # set polygon.vertices_
                for tmpPolygonID in polygons:
                    for edge in polygons[tmpPolygonID].edges_:
                        for vertex in edge.vertices_:
                            if vertex not in polygons[tmpPolygonID].vertices_:
                                polygons[tmpPolygonID].vertices_.append(vertex)
                # set vertex.dumpID_
                vertexIDs_List = sorted(list(vertices.keys()))
                for i in range(len(vertexIDs_List)):
                    vertices[vertexIDs_List[i]].dumpID_ = i

                # set edge.dump_ and polygon.dump_
                for tmpEdgeID in edges:
                    edge = edges[tmpEdgeID]
                    dx = edge.vertices_[1].position_[0] - edge.vertices_[0].position_[0]
                    dy = edge.vertices_[1].position_[1] - edge.vertices_[0].position_[1]
                    dz = edge.vertices_[1].position_[2] - edge.vertices_[0].position_[2]
                    if abs(dx) > L/2.0 or abs(dy) > L/2.0 or abs(dz) > L/2.0:
                        edge.dump_ = False
                for tmpPolygonID in polygons:
                    polygon = polygons[tmpPolygonID]
                    for edge in polygon.edges_:
                        if not edge.dump_:
                            polygon.dump_ = False
                            break
                # dump vtk file
                if not os.path.exists("../vtk/{:06d}".format(runid)):
                    os.mkdir("../vtk/{:06d}".format(runid))
                dumpVTK(runid, timestamp, vertices, edges, polygons)
                print("dumped {:06d} {:07d}".format(runid, int(timestamp)))
                continue
            if verticesFlag:
                vertexID = int(lineSplit[0])
                x = float(lineSplit[1])
                y = float(lineSplit[2])
                z = float(lineSplit[3])
                vertex = Vertex(vertexID, (x, y, z))
                vertices[vertexID] = vertex
            if edgesFlag:
                edgeID = int(lineSplit[0])
                verticesList = []
                for i in range(1, len(lineSplit)):
                    verticesList.append(vertices[int(lineSplit[i])])
                edge = Edge(verticesList)
                edges[edgeID] = edge
            if polygonsFlag:
                polygonID = int(lineSplit[0])
                edgesList = []
                for i in range(1, len(lineSplit)):
                    edgesList.append(edges[int(lineSplit[i])])
                polygon = Polygon(edgesList)
                polygons[polygonID] = polygon

def dumpVTK(runid, timestamp, vertices, edges, polygons):
    with open(os.path.join("../vtk/{:06d}".format(runid), "{:07d}.sample.vtk".format(int(timestamp))), "w") as file:
        file.write("# vtk DataFile Version 2.0\npolydata\nASCII\nDATASET POLYDATA\n")
        file.write("POINTS {:d} double\n".format(len(vertices)))
        # vertices = dict(sorted(vertices.items(), key=lambda item: item[1].id_))
        for key in vertices:
            vertex = vertices[key]
            x = vertex.position_[0]
            y = vertex.position_[1]
            z = vertex.position_[2]
            file.write("{:12.5e} {:12.5e} {:12.5e}\n".format(x, y, z))

        Npolygons = 0
        NpolygonVertices = 0
        for key in polygons:
            polygon = polygons[key]
            if polygon.dump_:
                Npolygons += 1
                NpolygonVertices += len(polygon.vertices_)

        file.write("\nPOLYGONS {:d} {:d}\n".format(Npolygons, Npolygons + NpolygonVertices))
        for key in polygons:
            polygon = polygons[key]
            if polygons[key].dump_:
                file.write("{:d}".format(len(polygon.vertices_)))
                for vertex in polygon.vertices_:
                    file.write(" {:6d}".format(vertex.dumpID_))
                file.write("\n")
        file.write("\n")

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--begin_num", action="store", type=int, help="the id of the first job to be run")
    parser.add_argument("-e", "--end_num", action="store", type=int, default=-1,
                        help="the id of the last job to be run")
    parser.add_argument("-n", "--num_runs", action="store", type=int, default=1,
                        help="number of runs")
    args = parser.parse_args()

    if args.end_num > 0:
        end_num = args.end_num
    else:
        end_num = args.begin_num + args.num_runs - 1

    main(args.begin_num, end_num)