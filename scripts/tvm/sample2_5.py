#!/usr/bin/env python3
'''
generate the initial configuration for the code TVM where Nc cells are generated using the voro++ library
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
'''

import pyvoro
import numpy as np

class Vertex:
    def __init__(self, id):
        self.id_ = id
        self.position_ = [0., 0., 0.]
        self.edges_ = []
        self.polygons_ = []
        self.cells_ = []
        self.on_ = False
        self.dumpID = id
    def setPosition(self, position):
        self.position_ = list(position)
    def addCell(self, cell):
        self.cells_.append(cell)

class Edge:
    def __init__(self, id):
        self.id_ = id
        self.dumpOn_ = True
        self.vertices_ = []
        self.polygons_ = []
        self.cells_ = []
        self.on_ = False
    def addVertex(self, vertex):
        self.vertices_.append(vertex)

class Polygon:
    def __init__(self, id):
        self.id_ = id
        self.dumpOn_ = True
        self.edges_ = []
        self.vertices_ = []
        self.cells_ = []
        self.on_ = False
        self.type_ = 0
    def addEdge(self, edge):
        self.edges_.append(edge)
    def addVertex(self, vertex):
        self.vertices_.append(vertex)

class Cell:
    def __init__(self, id):
        self.id_ = id
        self.polygons_ = []
        self.edges_ = []
        self.vertices_ = []
        self.type_ = -1
    def addPolygon(self, polygon):
        self.polygons_.append(polygon)

def main():
    Lx, Ly, Lz = (20, 20, 20)
    points = generatePoints(Lx, Ly, Lz)
    # points = [[1.0, 2.0, 3.0], [4.0, 5.5, 6.0]]
    voroDict = pyvoro.compute_voronoi(
        points,  # point positions
        [[-Lx/2.0, Lx/2.0], [-Ly/2.0, Ly/2.0], [-Lz/2.0, Lz/2.0]],  # limits
        1.0,  # block size
        periodic = [False, False, False]    # periodic boundary condition
    )
    vertices, edges, polygons, cells = topologyVoro(voroDict, (Lx, Ly, Lz))
    print("box before cut:")
    print("number of vertices:", len(vertices))
    print("number of edges:", len(edges))
    print("number of polygons:", len(polygons))
    print("number of cells:", len(cells))

    # cut sphere
    R = 4.0
    vertices, edges, polygons, cells = cutSphere(vertices, edges, polygons, cells, R)
    print("cylinder:")
    print("number of vertices:", len(vertices))
    print("number of edges:", len(edges))
    print("number of polygons:", len(polygons))
    print("number of cells:", len(cells))
    # topologyCheck(vertices, edges, polygons, cells)

    for key in cells:
        cell = cells[key]
        if cell.type_ > 0:
            for polygon in cell.polygons_:
                polygon.type_ = 1

    dumpVTK(vertices, polygons)
    dumpEmptyVTK(vertices, polygons)

    dumpTopo(vertices, edges, polygons, cells)

def topologyVoro(voroDict, Lxyz):
    Lx, Ly, Lz = Lxyz
    vertices = {}
    edges = {}
    polygons = {}
    cells = {}
    for rcell in voroDict:
        roriginal = rcell['original']
        rvolume = rcell['volume']
        rvertices = rcell['vertices']
        radjacency = rcell['adjacency']
        rfaces = rcell['faces']
        # print("original:", roriginal)
        # print("volume:", rvolume)
        # print("vertices:", rvertices)
        # print("adjacency:", radjacency)
        # print("faces:", rfaces)
        # print("")

        # generate cell
        cell = Cell(len(cells))
        cells[cell.id_] = cell

        # generate vertices
        vertexKeys = []
        for vPosition in rvertices:
            # vPosition[0] = vPosition[0] % Lx
            # vPosition[1] = vPosition[1] % Ly
            # vPosition[2] = (vPosition[2] + Lz/2.0) % Lz - Lz/2.0
            key = "%.6f,%.6f,%.6f"%(vPosition[0], vPosition[1], vPosition[2])
            vertexKeys.append(key)
            if key not in vertices:
                vertex = Vertex(len(vertices))
                vertex.setPosition(vPosition)
                vertices[key] = vertex
            vertices[key].addCell(cell)

        # generate polygons and edges
        for face in rfaces:
            tmpVertices = face['vertices']
            edgeKeys = []
            for i in range(len(tmpVertices)):
                v0 = vertices[vertexKeys[tmpVertices[i]]]
                v1 = vertices[vertexKeys[tmpVertices[(i + 1)%len(tmpVertices)]]]
                edgekey = (min(v0.id_, v1.id_), max(v0.id_, v1.id_))
                edgeKeys.append(edgekey)
                if edgekey not in edges:
                    edge = Edge(len(edges))
                    edge.addVertex(v0)
                    edge.addVertex(v1)
                    edges[edgekey] = edge
            edgeIDs = []
            for edgeKey in edgeKeys:
                edgeIDs.append(edges[edgeKey].id_)
            edgeIDs.sort()
            polygonKey = ""
            for edgeID in edgeIDs:
                polygonKey = polygonKey + "%d "%(edgeID)
            if polygonKey not in polygons:
                polygon = Polygon(len(polygons))
                for edgeKey in edgeKeys:
                    polygon.addEdge(edges[edgeKey])
                for i in range(len(tmpVertices)):
                    v0 = vertices[vertexKeys[tmpVertices[i]]]
                    polygon.addVertex(v0)
                polygons[polygonKey] = polygon
            cell.addPolygon(polygons[polygonKey])

    return vertices, edges, polygons, cells

def generatePoints(Lx, Ly, Lz):
    # from datetime import datetime
    # np.random.seed(int(datetime.utcnow().timestamp()))
    # np.random.seed(2161133)
    Nvertices = int(Lx*Ly*Lz)
    points = []
    for i in range(Nvertices):
        point = []
        point.append(np.random.uniform(-Lx/2.0, Lx/2.0))
        point.append(np.random.uniform(-Ly/2.0, Ly/2.0))
        point.append(np.random.uniform(-Lz/2.0, Lz/2.0))
        points.append(point)

    return points

def cutSphere(vertices, edges, polygons, cells, R):
    # look for vertices in each cell
    for key in cells:
        cell = cells[key]
        for polygon in cell.polygons_:
            for vertex in polygon.vertices_:
                if vertex not in cell.vertices_:
                    cell.vertices_.append(vertex)
    # compute center x, y coordinates of each cell
    for key in cells:
        cell = cells[key]
        cx = [0., 0., 0.]
        for vertex in cell.vertices_:
            cx[0] += vertex.position_[0]
            cx[1] += vertex.position_[1]
            cx[2] += vertex.position_[2]
        cx[0] /= len(cell.vertices_)
        cx[1] /= len(cell.vertices_)
        cx[2] /= len(cell.vertices_)
        # cut cylinder with center (0, 0), and radius R
        dx = cx[0]
        dy = cx[1]
        dz = cx[2]
        dd = np.sqrt(dx*dx + dy*dy + dz*dz)
        if dd < R:
            cell.type_ = 1
            # print(cx[0], cx[1], dx, dy, dd)

    # locate two layers of empty cells surrounding the spheroid
    # locate the first layer of empty cells
    for key in cells:
        cell = cells[key]
        if cell.type_ > 0:
            for polygon in cell.polygons_:
                polygon.on_ = True
    for key in cells:
        cell = cells[key]
        if cell.type_ < 0:
            addCell = False
            for polygon in cell.polygons_:
                if polygon.on_:
                    addCell = True
                    break
            if addCell:
                cell.type_ = 0
    # locate the second layer of empty cells
    for key in cells:
        cell = cells[key]
        if cell.type_ == 0:
            for polygon in cell.polygons_:
                polygon.on_ = True
    for key in cells:
        cell = cells[key]
        if cell.type_ < 0:
            addCell = False
            for polygon in cell.polygons_:
                if polygon.on_:
                    addCell = True
                    break
            if addCell:
                cell.type_ = 0

    for key in cells:
        cell = cells[key]
        if cell.type_ >= 0:
            for polygon in cell.polygons_:
                polygon.on_ = True
    for key in polygons:
        polygon = polygons[key]
        if polygon.on_:
            for edge in polygon.edges_:
                edge.on_ = True
    for key in edges:
        edge = edges[key]
        if edge.on_:
            for vertex in edge.vertices_:
                vertex.on_ = True

    verticesCut = {}
    edgesCut = {}
    polygonsCut = {}
    cellsCut = {}
    for key in vertices:
        if vertices[key].on_:
            verticesCut[key] = vertices[key]
    for key in edges:
        if edges[key].on_:
            edgesCut[key] = edges[key]
    for key in polygons:
        if polygons[key].on_:
            polygonsCut[key] = polygons[key]
    for key in cells:
        if cells[key].type_ >= 0:
            cellsCut[key] = cells[key]

    return verticesCut, edgesCut, polygonsCut, cellsCut

# def topologyCheck(vertices, edges, polygons, cells):
#     # initialization
#     for key in vertices:
#         vertex = vertices[key]
#         vertex.cells_ = []
#         vertex.edges_ = []
#     for key in edges:
#         edge = edges[key]
#         edge.polygons_ = []
#         edge.cells_ = []
#     for key in polygons:
#         polygon = polygons[key]
#         polygon.cells_ = []
#     for key in cells:
#         cell = cells[key]
#         for polygon in cell.polygons_:
#             if cell not in polygon.cells_:
#                 polygon.cells_.append(cell)
#             for edge in polygon.edges_:
#                 if polygon not in edge.polygons_:
#                     edge.polygons_.append(polygon)
#                 if cell not in edge.cells_:
#                     edge.cells_.append(cell)
#                 for vertex in edge.vertices_:
#                     if edge not in vertex.edges_:
#                         vertex.edges_.append(edge)
#                     if cell not in vertex.cells_:
#                         vertex.cells_.append(cell)
#
#     # check number of neighboring cells of each vertex
#     n1 = 0
#     n2 = 0
#     n3 = 0
#     n4 = 0
#     for key in vertices:
#         vertex = vertices[key]
#         if len(vertex.cells_) == 1:
#             n1 += 1
#         elif len(vertex.cells_) == 2:
#             n2 += 1
#         elif len(vertex.cells_) == 3:
#             n3 += 1
#         elif len(vertex.cells_) == 4:
#             n4 += 1
#         else:
#             print("found vertex with", len(vertex.cells_), "neighboring cells")
#     print("the number of vertices with 1/2/3/4 neighboring cells is ", n1, n2, n3, n4)
#
#     # check number of neighboring edges of each vertex
#     n1 = 0
#     n2 = 0
#     n3 = 0
#     n4 = 0
#     for key in vertices:
#         vertex = vertices[key]
#         if len(vertex.edges_) == 1:
#             n1 += 1
#         elif len(vertex.edges_) == 2:
#             n2 += 1
#         elif len(vertex.edges_) == 3:
#             n3 += 1
#         elif len(vertex.edges_) == 4:
#             n4 += 1
#         else:
#             print("found vertex with", len(vertex.edges_), "neighboring edges")
#     print("the number of vertices with 1/2/3/4 neighboring edges is ", n1, n2, n3, n4)
#
#     # check number of neighboring cells of each edge
#     n1 = 0
#     n2 = 0
#     n3 = 0
#     n4 = 0
#     for key in edges:
#         edge = edges[key]
#         if len(edge.cells_) == 1:
#             n1 += 1
#         elif len(edge.cells_) == 2:
#             n2 += 1
#         elif len(edge.cells_) == 3:
#             n3 += 1
#         elif len(edge.cells_) == 4:
#             n4 += 1
#         else:
#             print("found edge with", len(edge.cells_), "neighboring cells")
#     print("the number of edges with 1/2/3/4 neighboring cells is ", n1, n2, n3, n4)
#
#     # check number of neighboring polygons of each edge
#     n1 = 0
#     n2 = 0
#     n3 = 0
#     n4 = 0
#     for key in edges:
#         edge = edges[key]
#         if len(edge.polygons_) == 1:
#             n1 += 1
#         elif len(edge.polygons_) == 2:
#             n2 += 1
#         elif len(edge.polygons_) == 3:
#             n3 += 1
#         elif len(edge.polygons_) == 4:
#             n4 += 1
#         else:
#             print("found edge with", len(edge.polygons_), "neighboring polygons")
#     print("the number of edges with 1/2/3/4 neighboring polygons is ", n1, n2, n3, n4)
#
#     # check if two polygons have more than one common edge
#     commonEdgeLog = {}
#     for key in edges:
#         edge = edges[key]
#         for i in range(len(edge.polygons_)-1):
#             for j in range(i+1, len(edge.polygons_)):
#                 p1 = edge.polygons_[i]
#                 p2 = edge.polygons_[j]
#                 p12 = (min(p1.id_, p2.id_), max(p1.id_, p2.id_))
#                 if p12 in commonEdgeLog:
#                     commonEdgeLog[p12] += 1
#                 else:
#                     commonEdgeLog[p12] = 1
#     n1 = 0
#     n2 = 0
#     for p12 in commonEdgeLog:
#         if commonEdgeLog[p12] == 1:
#             n1 += 1
#         elif commonEdgeLog[p12] == 2:
#             n2 += 1
#         else:
#             print("found polygons", p12, "with more than 2 common edges")
#     print("the number of pairs of polygons with 1/2 common edges is ", n1, n2)
#
#     # check if two cells have more than one common polygon
#     commonPolygonLog = {}
#     for key in polygons:
#         polygon = polygons[key]
#         for i in range(len(polygon.cells_) - 1):
#             for j in range(i + 1, len(polygon.cells_)):
#                 c1 = polygon.cells_[i]
#                 c2 = polygon.cells_[j]
#                 c12 = (min(c1.id_, c2.id_), max(c1.id_, c2.id_))
#                 if c12 in commonPolygonLog:
#                     commonPolygonLog[c12] += 1
#                 else:
#                     commonPolygonLog[c12] = 1
#     n1 = 0
#     n2 = 0
#     for c12 in commonPolygonLog:
#         if commonPolygonLog[c12] == 1:
#             n1 += 1
#         elif commonPolygonLog[c12] == 2:
#             n2 += 1
#         else:
#             print("found cells", c12, "with more than 2 common polygons")
#     print("the number of pairs of cells with 1/2 common polygons is ", n1, n2)

def dumpVTK(vertices, polygons):
    with open("sample.vtk", "w") as file:
        file.write("# vtk DataFile Version 2.0\npolydata\nASCII\nDATASET POLYDATA\n")
        file.write("POINTS {:d} double\n".format(len(vertices)))
        vertices = dict(sorted(vertices.items(), key=lambda item: item[1].id_))
        i = 0
        for key in vertices:
            vertex = vertices[key]
            x = vertex.position_[0]
            y = vertex.position_[1]
            z = vertex.position_[2]
            file.write("{:12.5e} {:12.5e} {:12.5e}\n".format(x, y, z))
            vertex.dumpID = i
            i += 1

        Npolygons = 0
        NpolygonVertices = 0
        for key in polygons:
            polygon = polygons[key]
            if polygon.type_ > 0:
                Npolygons += 1
                NpolygonVertices += len(polygon.vertices_)

        file.write("\nPOLYGONS {:d} {:d}\n".format(Npolygons, Npolygons + NpolygonVertices))
        for key in polygons:
            polygon = polygons[key]
            if polygons[key].type_ > 0:
                file.write("{:d}".format(len(polygon.vertices_)))
                for vertex in polygon.vertices_:
                    file.write(" {:6d}".format(vertex.dumpID))
                file.write("\n")
        file.write("\n")

def dumpEmptyVTK(vertices, polygons):
    with open("empty.vtk", "w") as file:
        file.write("# vtk DataFile Version 2.0\npolydata\nASCII\nDATASET POLYDATA\n")
        file.write("POINTS {:d} double\n".format(len(vertices)))
        vertices = dict(sorted(vertices.items(), key=lambda item: item[1].id_))
        i = 0
        for key in vertices:
            vertex = vertices[key]
            x = vertex.position_[0]
            y = vertex.position_[1]
            z = vertex.position_[2]
            file.write("{:12.5e} {:12.5e} {:12.5e}\n".format(x, y, z))
            vertex.dumpID = i
            i += 1

        Npolygons = 0
        NpolygonVertices = 0
        for key in polygons:
            polygon = polygons[key]
            if polygon.type_ == 0:
                Npolygons += 1
                NpolygonVertices += len(polygon.vertices_)

        file.write("\nPOLYGONS {:d} {:d}\n".format(Npolygons, Npolygons + NpolygonVertices))
        for key in polygons:
            polygon = polygons[key]
            if polygons[key].type_ == 0:
                file.write("{:d}".format(len(polygon.vertices_)))
                for vertex in polygon.vertices_:
                    file.write(" {:6d}".format(vertex.dumpID))
                file.write("\n")
        file.write("\n")

def dumpTopo(vertices, edges, polygons, cells):
    with open("sample.topo", "w") as file:
        file.write("vertices {:d}\n".format(len(vertices)))
        for key in vertices:
            vertex = vertices[key]
            id = vertex.id_
            x = vertex.position_[0]
            y = vertex.position_[1]
            z = vertex.position_[2]
            file.write("{:6d} {:12.5e} {:12.5e} {:12.5e}\n".format(id, x, y, z))

        file.write("edges {:d}\n".format(len(edges)))
        for key in edges:
            edge = edges[key]
            file.write("{:d}".format(edge.id_))
            for vertex in edge.vertices_:
                file.write(" {:6d}".format(vertex.id_))
            file.write("\n")

        file.write("polygons {:d}\n".format(len(polygons)))
        for key in polygons:
            polygon = polygons[key]
            file.write("{:d}".format(polygon.id_))
            for edge in polygon.edges_:
                file.write(" {:6d}".format(edge.id_))
            file.write("\n")

        file.write("cells {:d}\n".format(len(cells)))
        for key in cells:
            cell = cells[key]
            file.write("{:d}".format(cell.id_))
            for polygon in cell.polygons_:
                file.write(" {:6d}".format(polygon.id_))
            file.write(" {:d}".format(cell.type_))
            file.write("\n")

if __name__ == '__main__':
    main()
