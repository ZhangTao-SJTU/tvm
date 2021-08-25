#!/usr/bin/env python3
'''
generate the initial configuration for the code TVM where Nc cells are generated using the voro++ library
author: Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn, May 2021
corresponding author: Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu
'''

import pyvoro
import numpy as np

class Vertex:
    def __init__(self, id):
        self.id_ = id
        self.position_ = [0., 0., 0.]
        self.cells_ = []
    def setPosition(self, position):
        self.position_ = list(position)
    def addCell(self, cell):
        self.cells_.append(cell)

class Edge:
    def __init__(self, id):
        self.id_ = id
        self.dumpOn_ = True
        self.vertices_ = []
    def addVertex(self, vertex):
        self.vertices_.append(vertex)

class Polygon:
    def __init__(self, id):
        self.id_ = id
        self.dumpOn_ = True
        self.edges_ = []
        self.vertices_ = []
    def addEdge(self, edge):
        self.edges_.append(edge)
    def addVertex(self, vertex):
        self.vertices_.append(vertex)

class Cell:
    def __init__(self, id):
        self.id_ = id
        self.polygons_ = []
    def addPolygon(self, polygon):
        self.polygons_.append(polygon)

def main():
    Lx, Ly, Lz = (8, 8, 8)
    points = generatePoints(Lx, Ly, Lz)
    # points = [[1.0, 2.0, 3.0], [4.0, 5.5, 6.0]]
    voroDict = pyvoro.compute_voronoi(
        points,  # point positions
        [[0.0, Lx], [0.0, Ly], [0.0, Lz]],  # limits
        2.0,  # block size
        periodic = [True, True, True]    # periodic boundary condition
    )
    vertices, edges, polygons, cells = topologyVoro(voroDict, (Lx, Ly, Lz))
    print("number of vertices:", len(vertices))
    print("number of edges:", len(edges))
    print("number of polygons:", len(polygons))
    print("number of cells:", len(cells))
    # for key0 in vertices:
    #     v0 = vertices[key0]
    #     for key1 in vertices:
    #         v1 = vertices[key1]
    #         if v0.id_ == v1.id_:
    #             continue
    #         dx = v1.position_[0] - v0.position_[0]
    #         dy = v1.position_[1] - v0.position_[1]
    #         dz = v1.position_[2] - v0.position_[2]
    #         tolerance = 1e-5
    #         if abs(np.sqrt(dx**2+dy**2+dz**2)) < tolerance:
    #             print("error: redundant vertices "+"{:d} ".format(v0.id_)+key0+" {:d} ".format(v1.id_)+key1)
    #             exit(1)

    checkVertexCells(vertices, edges, polygons, cells)

    # dumpVTK(vertices, edges, polygons, cells, points, (Lx, Ly, Lz))
    dumpTopo(vertices, edges, polygons, cells, points, (Lx, Ly, Lz))

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
            vPosition[0] = vPosition[0] % Lx
            vPosition[1] = vPosition[1] % Ly
            vPosition[2] = vPosition[2] % Lz
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
    from datetime import datetime
    np.random.seed(int(datetime.utcnow().timestamp()))
    # np.random.seed(2161133)
    Nvertices = int(Lx*Ly*Lz)
    points = []
    for i in range(Nvertices):
        point = []
        point.append(np.random.uniform(0, Lx))
        point.append(np.random.uniform(0, Ly))
        point.append(np.random.uniform(0, Lz))
        points.append(point)

    return points

def checkVertexCells(vertices, edges, polygons, cells):
    for key in cells:
        cell = cells[key]
        for polygon in cell.polygons_:
            for edge in polygon.edges_:
                for vertex in edge.vertices_:
                    if cell not in vertex.cells_:
                        vertex.addCell(cell)
    for key in vertices:
        vertex = vertices[key]
        if len(vertex.cells_) != 4:
            print("the number of neighboring cells of vertex "+key+" is not 4")
            exit(1)

# def dumpVTK(vertices, edges, polygons, cells, points, Lxyz):
#     Lx, Ly, Lz = Lxyz
#     with open("sample.vtk", "w") as file:
#         file.write("# vtk DataFile Version 2.0\npolydata\nASCII\nDATASET POLYDATA\n")
#         file.write("POINTS {:d} double\n".format(len(vertices)))
#         vertices = dict(sorted(vertices.items(), key=lambda item: item[1].id_))
#         for key in vertices:
#             vertex = vertices[key]
#             x = vertex.position_[0]
#             y = vertex.position_[1]
#             z = vertex.position_[2]
#             file.write("{:12.5e} {:12.5e} {:12.5e}\n".format(x, y, z))
#
#         Npolygons = 0
#         NpolygonVertices = 0
#         for key in edges:
#             v0 = edges[key].vertices_[0]
#             v1 = edges[key].vertices_[1]
#             dx = v1.position_[0] - v0.position_[0]
#             dy = v1.position_[1] - v0.position_[1]
#             dz = v1.position_[2] - v0.position_[2]
#             if abs(dx) > 0.5 * Lx or abs(dy) > 0.5 * Ly or abs(dz) > 0.5 * Lz:
#                 edges[key].dumpOn_ = False
#             else:
#                 edges[key].dumpOn_ = True
#         for key in polygons:
#             polygon = polygons[key]
#             polygon.dumpOn_ = True
#             for edge in polygon.edges_:
#                 if not edge.dumpOn_:
#                     polygon.dumpOn_ = False
#                     break
#             if polygon.dumpOn_:
#                 Npolygons += 1
#                 NpolygonVertices += len(polygon.vertices_)
#
#         file.write("\nPOLYGONS {:d} {:d}\n".format(Npolygons, Npolygons + NpolygonVertices))
#         for key in polygons:
#             polygon = polygons[key]
#             if polygons[key].dumpOn_:
#                 file.write("{:d}".format(len(polygon.vertices_)))
#                 for vertex in polygon.vertices_:
#                     file.write(" {:6d}".format(vertex.id_))
#                 file.write("\n")
#         file.write("\n")
#
#     with open("points.sample.vtk", "w") as file:
#         file.write("# vtk DataFile Version 2.0\npolydata\nASCII\nDATASET POLYDATA\n")
#         file.write("POINTS {:d} double\n".format(len(points)))
#         for point in points:
#             x = point[0]
#             y = point[1]
#             z = point[2]
#             file.write("{:12.5e} {:12.5e} {:12.5e}\n".format(x, y, z))

def dumpTopo(vertices, edges, polygons, cells, points, Lxyz):
    Lx, Ly, Lz = Lxyz
    with open("sample.topo", "w") as file:
        file.write("vertices {:d}\n".format(len(vertices)))
        # vertices = dict(sorted(vertices.items(), key=lambda item: item[1].id_))
        # count = 0
        for key in vertices:
            vertex = vertices[key]
            id = vertex.id_
            x = vertex.position_[0]
            y = vertex.position_[1]
            z = vertex.position_[2]
            file.write("{:6d} {:12.5e} {:12.5e} {:12.5e}\n".format(id, x, y, z))
            # if count != id:
            #     print("vertices dict disordered {:d} {:d}\n".format(count, id))
            #     exit(1)
            # count += 1
        # file.write("\n")

        # edges = dict(sorted(edges.items(), key=lambda item: item[1].id_))
        file.write("edges {:d}\n".format(len(edges)))
        # count = 0
        for key in edges:
            edge = edges[key]
            file.write("{:d}".format(edge.id_))
            for vertex in edge.vertices_:
                file.write(" {:6d}".format(vertex.id_))
            file.write("\n")
            # if count != edge.id_:
            #     print("edges dict disordered {:d} {:d}\n".format(count, edge.id_))
            #     exit(1)
            # count += 1
        # file.write("\n")

        # polygons = dict(sorted(polygons.items(), key=lambda item: item[1].id_))
        file.write("polygons {:d}\n".format(len(polygons)))
        # count = 0
        for key in polygons:
            polygon = polygons[key]
            file.write("{:d}".format(polygon.id_))
            for edge in polygon.edges_:
                file.write(" {:6d}".format(edge.id_))
            file.write("\n")
            # if count != polygon.id_:
            #     print("polygons dict disordered {:d} {:d}\n".format(count, polygon.id_))
            #     exit(1)
            # count += 1
        # file.write("\n")

        # cells = dict(sorted(cells.items(), key=lambda item: item[1].id_))
        file.write("cells {:d}\n".format(len(cells)))
        # count = 0
        for key in cells:
            cell = cells[key]
            file.write("{:d}".format(cell.id_))
            for polygon in cell.polygons_:
                file.write(" {:6d}".format(polygon.id_))
            file.write("\n")
            # if count != cell.id_:
            #     print("cells dict disordered {:d} {:d}\n".format(count, cell.id_))
            #     exit(1)
            # count += 1
        # file.write("\n")

if __name__ == '__main__':
    main()
