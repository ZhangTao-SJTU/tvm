import numpy as np
import copy
import os
from toolbox import topology

def levi_civita(i,j,k):
    if i==j or j==k or k==i:
        return 0
    elif (i,j,k) in [(0,1,2),(1,2,0),(2,0,1)]:
        return 1
    elif (i,j,k) in [(2,1,0),(0,2,1),(1,0,2)]: 
        return -1
    
def mapmaker(dictionary):
    result = {}
    for i, key in enumerate(dictionary):
        result[key] = i
    return result

# Creates a dictionary of type {time: [start_line, end_line]} from a .txt file
# Purpose: This is used to read the lines from topo.txt, cellCenter.txt, etc and
# obtain the tissue topology at a given time.

# We can deprecate the use of this function if the c++ code itself outputs 
# files like 0000000.time.txt, 000000.cellInfo.txt, etc.
# Alternatively, we can extract such files ourselves.

def timevals_to_lines_dict(textfile:str) -> dict[int,list[int]]: 
    timevals_dict={}
    last_line_number=0
    with open(textfile, 'r') as file:
        for i, line in enumerate(file):
            last_line_number+=1
            if len(line)==1: continue
            if line.split()[0]=="time": 
                timevals_dict[int(float(line.split()[1]))]=[]
                timevals_dict[int(float(line.split()[1]))].append(i)
    startvals=list(timevals_dict.values())
    for i,key in enumerate(timevals_dict.keys()):
        if i==len(startvals)-1:
            timevals_dict[key].append(last_line_number-1)
            break
        else:
            timevals_dict[key].append(startvals[i+1][0]-1)    
    return timevals_dict

def calculateOrigin(dir,time):
    if not dir.endswith("/"):
        raise ValueError("dir must end with a '/'")
    if not os.path.isfile("{}{:07d}.cellInfo.txt".format(dir,time)):
        writeTimeCellInfo(dir,time)

    center = np.array([0,0,0])
    counter = 0

    with open("{}{:07d}.cellInfo.txt".format(dir,time),"r") as file:
        lines = file.readlines()
        for i in range(1,len(lines)):
        
            if len(lines[i].split()):
                center = np.add(center,
                                np.array([float(lines[i].split()[1]),
                                          float(lines[i].split()[2]),
                                          float(lines[i].split()[3])]))
                counter += 1

    center = np.multiply(1 / counter,center)
    return center


# Given cellShapeIndex.txt, cellCenter.txt, and cellVolume.txt, this function
# returns a text file "{:07d}".format(time)"+".cellInfo.txt", 
# combining the information in the following format:
# id centerX centerY centerZ volume shapeIndex
def writeTimeCellInfo(output_dir,time):
    # equip id_to_cellInfo dictionary with the information from the three files
    id_to_cellInfo = {}
    time_Found = False
    with open(output_dir + "cellCenter.txt", "r") as file:
        time_Flag = False
        for line in file.readlines():
            # pass any blank lines
            if not len(line.split()): continue
            
            if (line.split()[0] == "time" 
                and int(float(line.split()[1])) == int(time)): 
                    time_Flag = True
                    time_Found = True
                    continue
            if time_Flag:
                # stop at the next time
                if line.split()[0] == "time": break
        
                else: 
                    id = int(line.split()[0])
                    # if the id is not already in the dictionary, add it
                    if id not in id_to_cellInfo.keys():
                        id_to_cellInfo[id] = {}
                    # add the center coordinates    
                    id_to_cellInfo[id]["centerX"] = float(line.split()[1])
                    id_to_cellInfo[id]["centerY"] = float(line.split()[2])
                    id_to_cellInfo[id]["centerZ"] = float(line.split()[3])
        file.close()
    
    if not time_Found:
        raise ValueError("Time not found in cellCenter.txt")
    with open(output_dir + "cellVolume.txt", "r") as file:
        time_Flag = False
        for line in file.readlines():
            # pass any blank lines
            if not len(line.split()): continue
            
            if (line.split()[0] == "time" 
                and int(float(line.split()[1])) == int(time)): 
                time_Flag = True
                continue
            if time_Flag:
                # stop at the next time
                if line.split()[0] == "time": break
                else:
                    id = int(line.split()[0])
                    # add the volume
                    id_to_cellInfo[id]["volume"] = float(line.split()[1])
        file.close()
    with open(output_dir+"cellShapeIndex.txt", "r") as file:
        time_Flag = False
        for line in file.readlines():
            # pass any blank lines
            if not len(line.split()): continue
            
            if (line.split()[0] == "time" 
                and int(float(line.split()[1])) == int(time)): 
                time_Flag = True
                continue
            if time_Flag:
                # stop at the next time
                if line.split()[0] == "time": break
                else:
                    id = int(line.split()[0])
                    # add the shape index
                    id_to_cellInfo[id]["shapeIndex"] = float(line.split()[1])
        file.close()
    # write the dictionary to the cellInfo.txt file
    with open(output_dir+"{:07d}".format(time)+".cellInfo.txt", "w") as file:
        file.write("{:6} {:12} {:12} {:12} {:12} {:12}\n".format("id", 
                                              "centerX", 
                                              "centerY", 
                                              "centerZ", 
                                              "volume", 
                                              "shapeIndex"))
        for key in id_to_cellInfo.keys():
            file.write("{:<6d} {:>12.5e} {:>12.5e} {:>12.5e} {:<12.6f} {:<12.6f}\n".format(
                key,
                id_to_cellInfo[key]["centerX"],
                id_to_cellInfo[key]["centerY"],
                id_to_cellInfo[key]["centerZ"],
                id_to_cellInfo[key]["volume"],
                id_to_cellInfo[key]["shapeIndex"]))
        file.close()
    return

# Given topo.txt, this function extracts the topology at a given time
def make_time_topo(output_dir,time):
    time_Flag = False
    time_Found = False
    extracted_lines = []
    with open(output_dir + "topo.txt", "r") as file:
        for line in file.readlines():
            # pass any blank lines
            if not len(line.split()): continue

            if (line.split()[0] == "time" 
                and int(float(line.split()[1])) == int(time)): 
                time_Flag = True
                time_Found = True
                continue
            if time_Flag:
                # stop at the next time
                if line.split()[0] == "time": break
                else: extracted_lines.append(line)
        file.close()
    if not time_Found:
        raise ValueError("Time not found in topo.txt")
    with open(output_dir + "{:07d}".format(time)+".topo.txt", "w") as file:
        for line in extracted_lines:
            file.write(line)
        file.close()
    return

def make_sample_dot_topo(output_dir,vertices,edges,polygons,cells):
    ### make_sample_dot_topo(): make a new sample.topo file from the currently stored class variables.
    with open(output_dir+"sample.topo", "w") as file:
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
            for vertexID in edge.vertices_:
                file.write(" {:6d}".format(vertexID))
            file.write("\n")
        file.write("polygons {:d}\n".format(len(polygons)))
        for key in polygons:
            polygon = polygons[key]
            file.write("{:d}".format(polygon.id_))
            for edgeID in polygon.edges_:
                file.write(" {:6d}".format(edgeID))
            file.write("\n")
        file.write("cells {:d}\n".format(len(cells)))
        for key in cells:
            cell = cells[key]
            file.write("{:d}".format(cell.id_))
            for polygonID in cell.polygons_:
                file.write(" {:6d}".format(polygonID))
            file.write("\n")
    return

def linker_springs_dict(dir,time):
    with open(dir+"{:07d}.link.vtk".format(time),"r") as f:
        lines=f.readlines()
        POINTSFlag=False#flag for coordinates
        LINESFlag=False#flag for edges
        coordinates=[]
        edges=[]
        for line in lines:
            if len(line.split()):
                if line.split()[0]=="POINTS":
                    POINTSFlag=True
                    continue
                if line.split()[0]=="LINES":
                    POINTSFlag=False
                    LINESFlag=True
                    continue
                if line.split()[0]=="CELL_DATA":
                    LINESFlag=False
                    continue
                if POINTSFlag:
                    coordinates.append([float(line.split()[0]),
                                        float(line.split()[1]),
                                        float(line.split()[2])])
                if LINESFlag:
                    edges.append([int(line.split()[1]),
                                  int(line.split()[2])])
    coordinates_dict = {i:np.array(coordinates[i]) 
                        for i in range(len(coordinates))}
    #print("Time to make coordinates dict: ",time.time()-start)
    edges_dict = {i: edge
                  for i,edge in enumerate(edges)}
    #print("Time to make coordinates and edges dict: ",time.time()-start)
    return coordinates_dict,edges_dict

def write_configuration(sample,filename = "sample.topo"):
    with open(filename, "w") as file:
        file.write("vertices {:d}\n".format(len(sample.vertices_)))
        for key,vertex in sample.vertices_.items():
            id = vertex.id_
            x = vertex.position_[0]
            y = vertex.position_[1]
            z = vertex.position_[2]
            file.write("{:6d} {:.14f} {:.14f} {:.14f}\n".format(id, x, y, z))
        file.write("edges {:d}\n".format(len(sample.edges_)))
        for key,edge in sample.edges_.items():
            file.write("{:d}".format(edge.id_))
            for vertexID in edge.vertices_:
                file.write(" {:6d}".format(vertexID))
            file.write("\n")
        file.write("polygons {:d}\n".format(len(sample.polygons_)))
        for key, polygon in sample.polygons_.items():
            file.write("{:d}".format(polygon.id_))
            for edgeID in polygon.edges_:
                file.write(" {:6d}".format(edgeID))
            file.write("\n")
        file.write("cells {:d}\n".format(len(sample.cells_)))
        for key, cell in sample.cells_.items():
            file.write("{:d}".format(cell.id_))
            for polygonID in cell.polygons_:
                file.write(" {:6d}".format(polygonID))
            file.write("\n")

# P is a 3XN matrix, structured (with numpy arrays) as:
# [[x1,...xn], [y1,...yN],[z1,...zN]]
def bounding_ellipsoid(sample, cellID, tolerance = 1e-6):
    # Algorithm from:
    # https://www.mathworks.com/matlabcentral/fileexchange/9542-minimum-volume-enclosing-ellipsoid
    # https://stackoverflow.com/questions/1768197/bounding-ellipse/1768440#1768440
    x = []
    y = []
    z = []
    cell = sample.cells_[cellID]
    for vertexID in cell.vertices_:
        pos = sample.vertices_[vertexID].position_
        x.append(pos[0])
        y.append(pos[1])
        z.append(pos[2])
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    P = np.vstack((x,y,z))
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
    center = c.flatten()
    l = 1.0/np.sqrt(UQV[1][2])
    orientation = UQV[2][2,:]
    return (center, l, orientation)
