#!/usr/bin/env python3
'''
generate fcc lattice for ECM, and cut it with cubic/spherical/parallelepiped shell
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

import math
import numpy as np

class Vertex:
    def __init__(self, position, ijk, id):
        if len(position) == 3:
            self.position = list(position)
        else:
            print("wrong position when initiating Vertex")
            exit(1)
        if len(ijk) == 3:
            self.ijk = tuple(ijk)
        else:
            print("wrong ijk when initiating Vertex")
            exit(1)
        self.id = id
        self.edges = [[], [], [], [], [], []]
        self.cubes = []
        self.phantomVertices = []
        self.on = True
    def addEdge(self, edge, i):
        self.edges[i].append(edge)
    def addCube(self, cube):
        self.cubes.append(cube)

class PhantomVertex:
    def __init__(self, position, id):
        if len(position) == 3:
            self.position = list(position)
        else:
            print("wrong position when initiating Vertex")
            exit(1)
        self.id = id
        self.edges = [[], []]
        self.on = True
        self.type = 0
    def addEdge(self, edge):
        self.edges.append(edge)

class Edge:
    def __init__(self, vertices):
        if len(vertices) == 2:
            self.vertices = vertices
        else:
            print("wrong vertices when initiating Edge")
            exit(1)
        # edge id is a tuple
        v0 = vertices[0].id
        v1 = vertices[1].id
        self.id = (min(v0, v1), max(v0, v1))
        self.on = True
        self.orientation = 0
        self.end = False
        self.phantomVertices = []

class Cube:
    def __init__(self, vertices):
        if len(vertices) == 14:
            self.vertices = list(vertices)
        else:
            print("wrong vertices when initiating Cube")
            exit(1)
        self.edges = []
        self.on = True
        # save cube to vertex
        for vertex in self.vertices:
            vertex.addCube(self)
    def addEdge(self, edge):
        self.edges.append(edge)

def main():
    Nx, Ny, Nz = (19, 19, 19) #(16+3)
    l0 = 4.0 #halve this
    lk = 16. #->double this
    if Nx%2 == 0 or Ny%2 == 0 or Nz%2 == 0:
        print("Nx, Ny, Nz must be odd numbers")
        exit(1)
    vertices, edges, cubes = generateFCC(Nx, Ny, Nz, l0, lk)

    phantomVertices = generatePhantomNetwork(vertices, edges, cubes)

    print("Before cutting and diluting, number of phantom vertices:", len(phantomVertices))
    print("Before cutting and diluting, number of edges:", len(edges))

    # cut hole
    # for vertex in vertices:
    #     x, y, z = vertex.position
    #     if np.sqrt(x*x+y*y+z*z) < 5.1:
    #         vertex.on = False
    #         for pVertex in vertex.phantomVertices:
    #             pVertex.on = False
    #     # if x < 0:
    #     #     vertex.on = False
    # for vertex in vertices:
    #     if not vertex.on:
    #         for cube in vertex.cubes:
    #             cube.on = False
    # for edge in edges:
    #     edge.on = False
    # for cube in cubes:
    #     if cube.on:
    #         for edge in cube.edges:
    #             edge.on = True

    for pVertex in phantomVertices:
        x, y, z = pVertex.position
        if np.sqrt(x*x+y*y+z*z) < 5.5:
            pVertex.on = False
            for i in range(2):
                for edge in pVertex.edges[i]:
                    edge.on = False

    # find pVertices on the inner hole boundary
    for pVertex in phantomVertices:
        if pVertex.on:
            offFlag = False
            for i in range(2):
                for edge in pVertex.edges[i]:
                    if not edge.on:
                        offFlag = True
                        break
            if offFlag:
                pVertex.type = 1

    nEdgesCut = 0
    for edge in edges:
        if edge.on:
            nEdgesCut += 1

    # dilute edges
    p = 0.85
    for edge in edges:
        if edge.on:
            if np.random.uniform() < p:
                edge.on = True
            else:
                edge.on = False
    for pVertex in phantomVertices:
        offFlag = True
        for i in range(2):
            for edge in pVertex.edges[i]:
                if edge.on:
                    offFlag = False
                    break
        if offFlag:
            pVertex.on = False

    nPhantomVertices = 0
    nEdges = 0
    for pVertex in phantomVertices:
        if pVertex.on:
            nPhantomVertices += 1
    for edge in edges:
        if edge.on:
            nEdges += 1
    print("Number of phantom vertices:", nPhantomVertices)
    print("Number of edges:", nEdges)

    fibers = fibersClassification(phantomVertices, edges)

    zs = 0
    zn = 0
    for pVertex in phantomVertices:
        if pVertex.on:
            nEdgesOn = 0
            for i in range(2):
                for edge in pVertex.edges[i]:
                    if edge.on:
                        nEdgesOn += 1
            # if nEdgesOn < 4:
            #     print(pVertex.id, pVertex.on, nEdgesOn)
            zs += nEdgesOn
            zn += 1
    z = zs / zn
    print("connectivity: {:.4f}".format(z))

    nEdges = 0
    for edge in edges:
        if edge.on:
            nEdges += 1
    p = nEdges/nEdgesCut
    print("p: {:.4f}".format(p))

    nEndNodes = 0
    for pVertex in phantomVertices:
        if pVertex.type > 0:
            nEndNodes += 1
    print("number of fiber end nodes: {:d}".format(nEndNodes))

    dumpTopo(phantomVertices, fibers)
    dumpVTK(phantomVertices, fibers, l0, lk)

def generateFCC(Nx, Ny, Nz, l0, lk):
    vertices = []
    edges = []
    cubes = []
    verticesDict = {}
    edgesDict = {}

    lk = int(np.floor(lk))  # lk must be an integer

    # initiate vertices
    ddijk = ((0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,1), (1,0,1), (1,1,1), (0,1,1),
            (0.5,0.5,0), (0.5,0.5,1), (0.5,0,0.5), (0.5,1,0.5), (0,0.5,0.5), (1,0.5,0.5))
    orientations = ((0.0,  0.5,  0.5), (0.0,  0.5, -0.5), ( 0.5, 0.0,  0.5), ( 0.5, 0.0, -0.5), ( 0.5,  0.5, 0.0), ( 0.5, -0.5, 0.0),
                    (0.0, -0.5, -0.5), (0.0, -0.5,  0.5), (-0.5, 0.0, -0.5), (-0.5, 0.0,  0.5), (-0.5, -0.5, 0.0), (-0.5,  0.5, 0.0))
    countVertices = 0
    for ck in range(-(Nz // 2), Nz // 2):
        for cj in range(-(Ny // 2), Ny // 2):
            for ci in range(-(Nx // 2), Nx // 2):
                for dijk in ddijk:
                    i = ci + dijk[0]
                    j = cj + dijk[1]
                    k = ck + dijk[2]
                    # reset k into the range [-lk/2.0, lk/2.0]
                    i = (i + lk / 2.0) % lk - lk / 2.0
                    j = (j + lk / 2.0) % lk - lk / 2.0
                    k = (k + lk / 2.0) % lk - lk / 2.0
                    if (i, j, k) not in verticesDict:
                        x = l0 * i + 0.25 * l0
                        y = l0 * j + 0.25 * l0
                        z = l0 * k + 0.25 * l0
                        vertex = Vertex((x, y, z), (i, j, k), countVertices)
                        vertices.append(vertex)
                        verticesDict[(i, j, k)] = vertex
                        countVertices += 1

    # initiate cubes and edges
    for ck in range(-(Nz // 2), Nz // 2):
        for cj in range(-(Ny // 2), Ny // 2):
            for ci in range(-(Nx // 2), Nx // 2):
                # initiate cube
                verticesCube = []
                for dijk in ddijk:
                    i = ci + dijk[0]
                    j = cj + dijk[1]
                    k = ck + dijk[2]
                    # reset k into the range [-lk/2.0, lk/2.0]
                    i = (i + lk / 2.0) % lk - lk / 2.0
                    j = (j + lk / 2.0) % lk - lk / 2.0
                    k = (k + lk / 2.0) % lk - lk / 2.0
                    verticesCube.append(verticesDict[(i, j, k)])
                cube = Cube(verticesCube)
                cubes.append(cube)
                # initiate edges
                for vi in range(len(cube.vertices)):
                    for vj in range(vi + 1, len(cube.vertices)):
                        # check the distance between two vertices
                        v1 = cube.vertices[vi]
                        v2 = cube.vertices[vj]
                        di = v2.ijk[0] - v1.ijk[0]
                        dj = v2.ijk[1] - v1.ijk[1]
                        dk = v2.ijk[2] - v1.ijk[2]
                        # reset dk into the range [-lk/2.0, lk/2.0]
                        di = (di + lk / 2.0) % lk - lk / 2.0
                        dj = (dj + lk / 2.0) % lk - lk / 2.0
                        dk = (dk + lk / 2.0) % lk - lk / 2.0
                        for orientationID, orientation in enumerate(orientations):
                            if di == orientation[0] and dj == orientation[1] and dk == orientation[2]:
                                edgeID = (min(v1.id, v2.id), max(v1.id, v2.id))
                                if edgeID not in edgesDict:
                                    edge = Edge((v1, v2))
                                    edge.orientation = orientationID%6
                                    edges.append(edge)
                                    # print((di, dj, dk), orientationID%6)
                                    edgesDict[edgeID] = edge
                                    v1.edges[orientationID % 6].append(edge)
                                    v2.edges[orientationID % 6].append(edge)
                                else:
                                    edge = edgesDict[edgeID]
                                # save edge to cube
                                cube.addEdge(edge)
                                continue

    return vertices, edges, cubes

def generatePhantomNetwork(vertices, edges, cubes):
    countVertices = 0
    for vertex in vertices:
        # print(vertex.id)
        # copy three phantom vertices
        pVertex0 = PhantomVertex(vertex.position, countVertices)
        countVertices += 1
        pVertex1 = PhantomVertex(vertex.position, countVertices)
        countVertices += 1
        pVertex2 = PhantomVertex(vertex.position, countVertices)
        countVertices += 1
        vertex.phantomVertices.append(pVertex0)
        vertex.phantomVertices.append(pVertex1)
        vertex.phantomVertices.append(pVertex2)
        edgesOrder = np.random.permutation(6)
        # print(edgesOrder)
        for i in [0, 1]:
            orientation = edgesOrder[i]
            # print(orientation, vertex.edges[orientation][0].id, vertex.edges[orientation][1].id)
            for edge in vertex.edges[orientation]:
                pVertex0.edges[i].append(edge)
                edge.phantomVertices.append(pVertex0)
        for i in [0, 1]:
            orientation = edgesOrder[i + 2]
            # print(orientation, vertex.edges[orientation][0].id, vertex.edges[orientation][1].id)
            for edge in vertex.edges[orientation]:
                pVertex1.edges[i].append(edge)
                edge.phantomVertices.append(pVertex1)
        for i in [0, 1]:
            orientation = edgesOrder[i + 4]
            # print(orientation, vertex.edges[orientation][0].id, vertex.edges[orientation][1].id)
            for edge in vertex.edges[orientation]:
                pVertex2.edges[i].append(edge)
                edge.phantomVertices.append(pVertex2)

    phantomVertices = []
    for vertex in vertices:
        for pVertex in vertex.phantomVertices:
            # print(pVertex.id, pVertex.edges[0][0].id, pVertex.edges[0][1].id, pVertex.edges[1][0].id, pVertex.edges[1][1].id)
            phantomVertices.append(pVertex)
    # for edge in edges:
    #     print(edge.vertices[0].id, edge.vertices[1].id, edge.phantomVertices[0].id, edge.phantomVertices[1].id)

    return phantomVertices

def fibersClassification(phantomVertices, edges):
    import networkx as nx

    G = nx.Graph()
    for edge in edges:
        if edge.on:
            G.add_node(edge)
    for pVertex in phantomVertices:
        for pairEdges in pVertex.edges:
            if pairEdges[0].on and pairEdges[1].on:
                G.add_edge(pairEdges[0], pairEdges[1])
            if pairEdges[0].on and (not pairEdges[1].on):
                pairEdges[0].end = True
                # print(pairEdges[0].on, pairEdges[1].on)
            if (not pairEdges[0].on) and pairEdges[1].on:
                pairEdges[1].end = True
                # print(pairEdges[0].on, pairEdges[1].on)
    print("graph edges:", G.number_of_nodes())
    print("graph linkers:", G.number_of_edges())

    # remove loops
    for component in nx.connected_components(G):
        if len(component) >= 2:
            ends = []
            for edge in component:
                if edge.end:
                    ends.append(edge)
            if len(ends) == 0:
                list(component)[np.random.randint(len(component))].on = False
    G = nx.Graph()
    for edge in edges:
        if edge.on:
            G.add_node(edge)
    for pVertex in phantomVertices:
        for pairEdges in pVertex.edges:
            if pairEdges[0].on and pairEdges[1].on:
                G.add_edge(pairEdges[0], pairEdges[1])
            if pairEdges[0].on and (not pairEdges[1].on):
                pairEdges[0].end = True
                # print(pairEdges[0].on, pairEdges[1].on)
            if (not pairEdges[0].on) and pairEdges[1].on:
                pairEdges[1].end = True
                # print(pairEdges[0].on, pairEdges[1].on)
    print("graph edges without loops:", G.number_of_nodes())
    print("graph linkers without loops:", G.number_of_edges())

    # # remove single edge fibers
    # for component in nx.connected_components(G):
    #     if len(component) == 1:
    #         for edge in component:
    #             edge.on = False
    # for pVertex in phantomVertices:
    #     offFlag = True
    #     for i in range(2):
    #         for edge in pVertex.edges[i]:
    #             if edge.on:
    #                 offFlag = False
    #                 break
    #     if offFlag:
    #         pVertex.on = False
    # G = nx.Graph()
    # for edge in edges:
    #     if edge.on:
    #         G.add_node(edge)
    # for pVertex in phantomVertices:
    #     for pairEdges in pVertex.edges:
    #         if pairEdges[0].on and pairEdges[1].on:
    #             G.add_edge(pairEdges[0], pairEdges[1])
    #         if pairEdges[0].on and (not pairEdges[1].on):
    #             pairEdges[0].end = True
    #             # print(pairEdges[0].on, pairEdges[1].on)
    #         if (not pairEdges[0].on) and pairEdges[1].on:
    #             pairEdges[1].end = True
    #             # print(pairEdges[0].on, pairEdges[1].on)
    # print("graph edges without loops or single edge fibers:", G.number_of_nodes())
    # print("graph linkers without loops or single edge fibers:", G.number_of_edges())

    nEdges = 0
    fiberLengths = []
    fibers = []
    for component in nx.connected_components(G):
        nEdges += len(component)
        fiberLengths.append(len(component))
        if len(component) == 1:
            path = list(component)
            fibers.append(path)
        else:
            ends = []
            for edge in component:
                if edge.end:
                    ends.append(edge)
            path = list(nx.shortest_simple_paths(G, ends[0], ends[1]))[0]
            # for edge in path:
            #     print("({:d},{:d})".format(edge.phantomVertices[0].id, edge.phantomVertices[1].id), end = "")
            # print("")
            fibers.append(path)
        # for edge in component:
        #     print(edge.id, end="")
        # print("*")
        # for edge in path:
        #     print(edge.id, end="")
        # print("#")

    print("fiber edges:", nEdges)
    # print(sorted(fiberLengths))

    return fibers

def dumpVTK(phantomVertices, fibers, l0, lk):
    edges = []
    for fiber in fibers:
        for edge in fiber:
            edges.append(edge)
    for edge in edges:
        v0 = edge.vertices[0]
        v1 = edge.vertices[1]
        dx = v1.position[0] - v0.position[0]
        dy = v1.position[1] - v0.position[1]
        dz = v1.position[2] - v0.position[2]
        if np.fabs(dx) > l0*lk/2.0:
            edge.on = False
        if np.fabs(dy) > l0*lk/2.0:
            edge.on = False
        if np.fabs(dz) > l0*lk/2.0:
            edge.on = False
    with open("ECM.vtk", "w") as file:
        file.write("# vtk DataFile Version 2.0\npolydata\nASCII\nDATASET POLYDATA\n")
        file.write("POINTS {:d} double\n".format(len(phantomVertices)))
        for pVertex in phantomVertices:
            x = pVertex.position[0]
            y = pVertex.position[1]
            z = pVertex.position[2]
            file.write("{:12.5e} {:12.5e} {:12.5e}\n".format(x, y, z))

        num_edges = 0
        for edge in edges:
            if edge.on:
                num_edges += 1
        file.write("\nLINES {:d} {:d}\n".format(num_edges, 3*num_edges))
        for edge in edges:
            if edge.on:
                file.write("2 {:6d} {:6d}\n".format(edge.phantomVertices[0].id, edge.phantomVertices[1].id))
        file.write("\n")

        file.write("\nCELL_DATA {:d}\n".format(num_edges))
        file.write("SCALARS fiber int\n")
        file.write("LOOKUP_TABLE default\n")
        nFiber = 0
        for fiber in fibers:
            for edge in fiber:
                if edge.on:
                    file.write("{:d}\n".format(nFiber%10))
            nFiber += 1

        file.write("\nPOINT_DATA {:d}\n".format(len(phantomVertices)))
        file.write("SCALARS vertex int\n")
        file.write("LOOKUP_TABLE default\n")
        for pVertex in phantomVertices:
            file.write("{:d}\n".format(pVertex.type))
        file.write("\n")

def dumpTopo(phantomVertices, fibers):
    nPhantomVertices = 0
    for pVertex in phantomVertices:
        if pVertex.on:
            nPhantomVertices += 1
    with open("ECM.topo", "w") as file:
        file.write("nodes {:d}\n".format(nPhantomVertices))
        # vertices = dict(sorted(vertices.items(), key=lambda item: item[1].id_))
        # count = 0
        for pVertex in phantomVertices:
            if pVertex.on:
                id = pVertex.id
                x = pVertex.position[0]
                y = pVertex.position[1]
                z = pVertex.position[2]
                type = pVertex.type
                file.write("{:6d} {:12.5e} {:12.5e} {:12.5e} {:d}\n".format(id, x, y, z, type))

        file.write("fibers {:d}\n".format(len(fibers)))
        count = 0
        for fiber in fibers:
            if len(fiber) == 1:
                file.write("{:d} {:6d} {:6d}\n".format(count, fiber[0].phantomVertices[0].id, fiber[0].phantomVertices[1].id))
                count += 1
                continue
            fiberEdges = []
            for edge in fiber:
                fiberEdges.append(edge)
            file.write("{:d}".format(count))
            if fiberEdges[0].phantomVertices[0] in fiberEdges[1].phantomVertices:
                file.write(" {:6d}".format(fiberEdges[0].phantomVertices[1].id))
            else:
                file.write(" {:6d}".format(fiberEdges[0].phantomVertices[0].id))
            for i in range(len(fiberEdges) - 1):
                edge0 = fiberEdges[i]
                edge1 = fiberEdges[i + 1]
                if edge0.phantomVertices[0] in edge1.phantomVertices:
                    file.write(" {:6d}".format(edge0.phantomVertices[0].id))
                else:
                    file.write(" {:6d}".format(edge0.phantomVertices[1].id))
            if fiberEdges[-1].phantomVertices[0] in fiberEdges[-2].phantomVertices:
                file.write(" {:6d}".format(fiberEdges[-1].phantomVertices[1].id))
            else:
                file.write(" {:6d}".format(fiberEdges[-1].phantomVertices[0].id))
            file.write("\n")
            # if count != edge.id_:
            #     print("edges dict disordered {:d} {:d}\n".format(count, edge.id_))
            #     exit(1)
            count += 1
        # file.write("\n")

if __name__ == '__main__':
    main()
