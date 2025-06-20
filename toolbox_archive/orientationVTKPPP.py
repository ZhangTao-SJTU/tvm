#!/usr/bin/env python3
'''
generate vtk files from simulation output topo.txt file
author: Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn, May 2021
corresponding author: Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu
'''

from scipy.spatial import ConvexHull, convex_hull_plot_2d
import numpy as np
import os
from numpy.linalg import eig, inv

class Vertex:
    def __init__(self, id, xyz):
        self.id_ = id
        self.position_ = xyz

class Edge:
    def __init__(self, vertices):
        self.vertices_ = vertices

class Polygon:
    def __init__(self, edges):
        self.edges_ = edges
        self.cells_ = []

class Cell:
    def __init__(self, polygons):
        self.polygons_ = polygons
        self.vertices_ = []
        self.surface_ = False

def main(begin_num, end_num, remote_disk, timesteps):
    for runid in range(begin_num, end_num+1):
        run(runid, remote_disk, timesteps)

def ls_ellipsoid(P, tolerance):
    # https://www.mathworks.com/matlabcentral/fileexchange/9542-minimum-volume-enclosing-ellipsoid
    # https://stackoverflow.com/questions/1768197/bounding-ellipse/1768440#1768440

    # Dimension of the points
    d = 3
    # Number of points
    N = P.shape[1]

    # Add a row of 1s to the 3xN matrix P - so Q is 4xN now.
    Q = np.vstack((P,np.ones(N)))

    # Initialize
    count = 1
    err = 1
    #u is an Nx1 vector where each element is 1/N
    u = (1/N) * np.ones((N,1))

    # Khachiyan Algorithm
    while err > tolerance:
        # Matrix multiplication:
        # diag(u) : if u is a vector, places the elements of u in the diagonal of an NxN matrix of zeros
        X = np.dot(np.dot(Q, np.diag(u.flatten())), Q.transpose())

        # inv(X) returns the matrix inverse of X
        # diag(M) when M is a matrix returns the diagonal vector of M
        M = np.diag(np.dot(np.dot(Q.transpose(), np.linalg.inv(X)), Q))

        # Find the value and location of the maximum element in the vector M
        maximum = np.amax(M)
        j = np.argmax(M)

        # Calculate the step size for the ascent
        step_size = (maximum - d -1)/((d+1)*(maximum-1))

        # Calculate the new_u:
        # Take the vector u, and multiply all the elements in it by (1-step_size)
        new_u = (1 - step_size)*u

        # Increment the jth element of new_u by step_size
        new_u[j] = new_u[j] + step_size

        # Store the error by taking finding the square root of the SSD
        # between new_u and u
        # The SSD or sum-of-square-differences, takes two vectors
        # of the same size, creates a new vector by finding the
        # difference between corresponding elements, squaring
        # each difference and adding them all together.

        # So if the vectors were: a = [1 2 3] and b = [5 4 6], then:
        # SSD = (1-5)^2 + (2-4)^2 + (3-6)^2;
        # And the norm(a-b) = sqrt(SSD);
        err = np.sum((new_u - u) ** 2)

        # Increment count and replace u
        count = count + 1
        u = new_u

    # Put the elements of the vector u into the diagonal of a matrix
    # U with the rest of the elements as 0
    U = np.diag(u.flatten())

    # Compute the A-matrix
    Pu = np.dot(P, u)
    A = (1/d) * np.linalg.inv(np.dot(np.dot(P, U), P.transpose()) - np.dot(Pu, Pu.transpose()))

    # And the center
    c = np.dot(P, u)

    # compute the orientation
    UQV = np.linalg.svd(A)
    # print(c)
    # print(UQV[0])
    # print(1.0/np.sqrt(UQV[1]))
    # print(UQV[2])

    return (c.flatten(), 1.0/np.sqrt(UQV[1][2]), UQV[2][2,:])

def run(runid, remote_disk, timesteps):
    if remote_disk:
        runDir = "/mnt/hgfs/E/SJTU/VertexModel/run/{:06d}".format(runid)
    else:
        runDir = "{:06d}".format(runid)
    readFilePath = os.path.join(runDir, "topo.txt")
    cells_frames = {}
    with open(readFilePath, "r") as file:
        verticesFlag = False
        edgesFlag = False
        polygonsFlag = False
        cellsFlag = False
        for line in file:
            if len(line) <= 1:
                continue
            lineSplit = line.split()
            if lineSplit[0] == "time":
                timestamp = float(lineSplit[1])
                cells_frames[timestamp] = {}
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
                cellsFlag = True
                cellsCount = int(lineSplit[1])
                continue

            if verticesFlag:
                if lineSplit[-1] == "nan" or lineSplit[-1] == "-nan":
                    if remote_disk:
                        dumpDir = "/mnt/hgfs/E/SJTU/VertexModel/run/{:06d}".format(runid)
                    else:
                        dumpDir = "{:06d}".format(runid)
                    os.system("rm " + os.path.join(dumpDir, "cellOrientation.txt"))
                    return 1
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
            if cellsFlag:
                type = int(lineSplit[-1])
                if type > 0:
                    cellID = int(lineSplit[0])
                    polygonsList = []
                    for i in range(1, len(lineSplit) - 1):
                        polygonsList.append(polygons[int(lineSplit[i])])
                    cell = Cell(polygonsList)
                    cells_frames[timestamp][cellID] = cell
                cellsCount -= 1
                # end of current timestamp
                if cellsCount == 0:
                    for cellID in cells_frames[timestamp]:
                        cell = cells_frames[timestamp][cellID]
                        for polygon in cell.polygons_:
                            for edge in polygon.edges_:
                                for vertex in edge.vertices_:
                                    if vertex not in cell.vertices_:
                                        cell.vertices_.append(vertex)
                    print("{:06d} loaded timestamp: {:f}".format(runid, timestamp))
                    cellsFlag = False
                    # update polgyon.cells_ and cell.surface_
                    for cellID in cells_frames[timestamp]:
                        cell = cells_frames[timestamp][cellID]
                        for polygon in cell.polygons_:
                            if cell not in polygon.cells_:
                                polygon.cells_.append(cell)
                    for polygonID in polygons:
                        polygon = polygons[polygonID]
                        if len(polygon.cells_) < 2:
                            for cell in polygon.cells_:
                                cell.surface_ = True
                    # clean cell.polygons_
                    for cellID in cells_frames[timestamp]:
                        cell = cells_frames[timestamp][cellID]
                        cell.polygons_ = []
                    # if timestamp > 0:
                    #     break

    for timestamp in cells_frames:
        cellIDs = []
        centers = []
        lengths = []
        orientations = []
        cellSurfaces = []
        for cellID in cells_frames[timestamp]:
            cell = cells_frames[timestamp][cellID]
            # let us assume some definition of x, y and z
            xList = []
            yList = []
            zList = []
            for vertex in cell.vertices_:
                xList.append(vertex.position_[0])
                yList.append(vertex.position_[1])
                zList.append(vertex.position_[2])
            x = np.array(xList)
            y = np.array(yList)
            z = np.array(zList)
            P = np.vstack((x,y,z))

            center, l, orientation = ls_ellipsoid(P, 0.001)  # get ellipsoid polynomial coefficients
            cellIDs.append(cellID)
            centers.append(center)
            lengths.append(l)
            orientations.append(orientation)
            cellSurfaces.append(cell.surface_)

        dumpAxisVTK(runid, timestamp, centers, lengths, orientations, remote_disk)
        print("{:06d} processed timestamp: {:f}".format(runid, timestamp))

def dumpVTK(runid, timestamp, centers, orientations, remote_disk):
    if remote_disk:
        dumpDir = "/mnt/hgfs/E/SJTU/VertexModel/vtk/{:06d}".format(runid)
    else:
        dumpDir = "../vtk/{:06d}".format(runid)
    if not os.path.exists(dumpDir):
        os.mkdir(dumpDir)
    with open(os.path.join(dumpDir, "{:07d}.orientation.vtk".format(int(timestamp))), "w") as file:
        file.write("# vtk DataFile Version 2.0\npolydata\nASCII\nDATASET POLYDATA\n")
        file.write("POINTS {:d} double\n".format(len(centers)))
        # vertices = dict(sorted(vertices.items(), key=lambda item: item[1].id_))
        for center in centers:
            file.write("{:12.5e} {:12.5e} {:12.5e}\n".format(center[0], center[1], center[2]))

        file.write("\nPOINT_DATA {:d}\n".format(len(centers)))
        file.write("VECTORS orientation float\n")
        for orientation in orientations:
                file.write("{:12.5e} {:12.5e} {:12.5e}\n".format(orientation[0], orientation[1], orientation[2]))
        file.write("\n")

def dumpAxisVTK(runid, timestamp, centers, lengths, orientations, remote_disk):
    if remote_disk:
        dumpDir = "/mnt/hgfs/E/SJTU/VertexModel/vtk/{:06d}".format(runid)
    else:
        dumpDir = "../vtk/{:06d}".format(runid)
    if not os.path.exists(dumpDir):
        os.mkdir(dumpDir)
    with open(os.path.join(dumpDir, "{:07d}.orientation.vtk".format(int(timestamp))), "w") as file:
        vertices = []
        edges = []
        countVertex = 0
        for i in range(len(centers)):
            center = centers[i]
            l = lengths[i]/2.0
            orientation = orientations[i]
            x1 = center[0] + l * orientation[0]
            y1 = center[1] + l * orientation[1]
            z1 = center[2] + l * orientation[2]
            x2 = center[0] - l * orientation[0]
            y2 = center[1] - l * orientation[1]
            z2 = center[2] - l * orientation[2]
            v1 = Vertex(countVertex, (x1, y1, z1))
            countVertex += 1
            v2 = Vertex(countVertex, (x2, y2, z2))
            countVertex += 1
            vertices.append(v1)
            vertices.append(v2)
            edge = Edge((v1, v2))
            edges.append(edge)

        file.write("# vtk DataFile Version 2.0\npolydata\nASCII\nDATASET POLYDATA\n")
        file.write("POINTS {:d} double\n".format(len(vertices)))
        for vertex in vertices:
            x = vertex.position_[0]
            y = vertex.position_[1]
            z = vertex.position_[2]
            file.write("{:12.5e} {:12.5e} {:12.5e}\n".format(x, y, z))

        file.write("\nLINES {:d} {:d}\n".format(len(edges), 3*len(edges)))
        for edge in edges:
            file.write("2 {:6d} {:6d}\n".format(edge.vertices_[0].id_, edge.vertices_[1].id_))
        file.write("\n")

#text file with timed orientation of cells
def dumpCellOrientation(runid, timestamp, cellIDs, centers, orientations, cellSurfaces, remote_disk):
    if remote_disk:
        dumpDir = "/mnt/hgfs/E/SJTU/VertexModel/run/{:06d}".format(runid)
    else:
        dumpDir = "{:06d}".format(runid)
    if not os.path.exists(dumpDir):
        os.mkdir(dumpDir)
    if timestamp < 500.1:
        os.system("rm "+os.path.join(dumpDir, "cellOrientation.txt"))
    with open(os.path.join(dumpDir, "cellOrientation.txt"), "a") as file:
        file.write("time {:.0f}\n".format(timestamp))
        for i in range(len(centers)):
            file.write("{:<12d} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:d}\n".format(cellIDs[i],
                                                             centers[i][0], centers[i][1], centers[i][2],
                                                             orientations[i][0], orientations[i][1], orientations[i][2],
                                                             cellSurfaces[i]))

        file.write("\n")

def test():
    xc = np.array([ [np.random.uniform(0., 20.)],
                    [np.random.uniform(0., 20.)],
                    [np.random.uniform(0., 20.)]    ])
    l1 = np.random.uniform(0., 10.)
    l2 = np.random.uniform(0., 10.)
    l3 = np.random.uniform(0., 10.)
    alpha = np.random.uniform(0., 2.0 * np.pi)
    beta = np.random.uniform(0., 2.0 * np.pi)
    gamma = np.random.uniform(0., 2.0 * np.pi)
    Rx = np.array([[1.0, 0., 0.], [0., np.cos(alpha), -np.sin(alpha)], [0., np.sin(alpha), np.cos(alpha)]])
    Ry = np.array([[np.cos(beta), 0., np.sin(beta)], [0., 1., 0.], [-np.sin(beta), 0., np.cos(beta)]])
    Rz = np.array([[np.cos(gamma), -np.sin(gamma), 0.], [np.sin(gamma), np.cos(gamma), 0.], [0., 0., 1.]])
    R = np.dot(np.dot(Rz, Ry), Rx)
    N = 20
    P = np.zeros((3, N))
    print("center", xc)
    print("rotation'", R.transpose())
    print("axis", l1, l2, l3)
    for i in range(N):
        theta = np.random.uniform(0., 2.0 * np.pi)
        phi = np.random.uniform(0., np.pi)
        x = np.array([  [l1 * np.cos(theta) * np.sin(phi)],
                        [l2 * np.sin(theta) * np.sin(phi)],
                        [l3 * np.cos(phi)]  ])
        x = np.dot(R, x) + xc
        # x = x + xc
        P[:,i] = x.flatten()

    print("fitting:")
    center, orientation = ls_ellipsoid(P, 0.000001)  # get ellipsoid polynomial coefficients

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--begin_num", action="store", type=int, help="the id of the first job to be run")
    parser.add_argument("-e", "--end_num", action="store", type=int, default=-1,
                        help="the id of the last job to be run")
    parser.add_argument("-n", "--num_runs", action="store", type=int, default=1,
                        help="number of runs")
    parser.add_argument("-y", "--remote_disk", action="store", type=int, default=1,
                        help="read data on remote disk")
    parser.add_argument("-t", "--timesteps", action="store", type=int, nargs='+', default=[500, 65000], help="the timesteps to compute")
    args = parser.parse_args()

    if args.end_num > 0:
        end_num = args.end_num
    else:
        end_num = args.begin_num + args.num_runs - 1

    main(args.begin_num, end_num, args.remote_disk, args.timesteps)
    # test()