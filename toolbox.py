import numpy as np
import scipy
import matplotlib.pyplot as plt
import pandas as pd
from pyvtk import *
import scipy.stats as stats
import multiprocessing
import time as timer
import copy
import os



def mapmaker(dictionary):
    result={}
    for i, key in enumerate(dictionary):
        result[key]=i
    return result

def mk_coordinates_and_edges_dict(file):
    with open(file,"r") as f:
        lines=f.readlines()
        POINTSFlag=False#flag for coordinates
        LINESFlag=False#flag for edges
        TENSIONSFlag=False#flag for tensions
        coordinates=[]
        edges=[]
        tensions=[]
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
                if line.split()[0]=="LOOKUP_TABLE":
                    TENSIONSFlag=True
                    continue
                if POINTSFlag:
                    coordinates.append([float(line.split()[0]),float(line.split()[1]),float(line.split()[2])])
                if LINESFlag:
                    edges.append([int(line.split()[1]),int(line.split()[2])])

                if TENSIONSFlag:
                    tensions.append(float(line.split()[0]))

    coordinates_dict={i: np.array(coordinates[i]) for i in range(len(coordinates))}
    #print("Time to make coordinates dict: ",time.time()-start)
    edges_dict={i: fiber_edge(edges[i],tensions[i]) for i in range(len(edges))}
    #print("Time to make coordinates and edges dict: ",time.time()-start)
    return coordinates_dict,edges_dict

# Given an initial fiber network object, a final network object and an origin (to be center of mass of spheroid),
# this updates the displacements_dict, which is a dictionary of the form
# {radius (from r_bins): [displacements of nodes in that radius]}

def get_radius_and_avg_displacement_and_err_and_samplesize_arrays(dir_list,initial_time=1000,final_time=20000,r_bins=np.linspace(3.5,32,41)):

    displacements_dict={i:[] for i in r_bins}
    for dir in dir_list:
        # A dictionary of nodes for the initial network where each connected node is assigned to a radial bin
        nodes_dict={i:[] for i in r_bins}

        initial_network=fiber_network(config_dir=dir,time=initial_time)
        final_network=fiber_network(config_dir=dir,time=final_time)
        origin=initial_network.origin_
        #populate nodes_dict
        for node in initial_network.connected_nodes_:
            radial_distance = np.linalg.norm(initial_network.coordinates_dict_[node]-origin)
            if radial_distance>max(r_bins) or radial_distance<min(r_bins): continue
            bin = min(r_bins, key=lambda x: abs(x - radial_distance))
            nodes_dict[bin].append(node)


        
        for radius in r_bins:
            for node in nodes_dict[radius]:
                if len(nodes_dict[radius]) and node in final_network.connected_nodes_:
                    displacements_dict[radius].append(np.linalg.norm(initial_network.coordinates_dict_[node]-final_network.coordinates_dict_[node]))
        print("Processed ",dir)
    radius=[]
    displacement=[]
    samplesize=[]
    err=[]

    for key in displacements_dict.keys():
        if len(displacements_dict[key])>1:
            radius.append(key)
            displacement.append(np.mean(displacements_dict[key]))
            samplesize.append(len(displacements_dict[key]))
            err.append(stats.sem(displacements_dict[key]))

    return radius,displacement,err,samplesize

# Calculate the configuration tensor (in spherical coordinates) for a given fiber network (in the form of a dictionary of edges).
# The tensor is calculated for edges with tension above a certain cutoff, and within a spherical shell of radii between radial_distance_min and radial_distance_max.

def get_radius_and_omega_rr_and_err_and_samplesize_arrays(dir_list,time=25000,r_bins=np.linspace(3.5,32,41)):

    omega_dict={i:[] for i in r_bins}
    for dir in dir_list:
        # A dictionary of nodes for the initial network where each connected node is assigned to a radial bin
        radius_dict={i:[] for i in r_bins}

        network=fiber_network(config_dir=dir,time=time)
        network.calculate_edge_attributes()
        #populate nodes_dict
        for edgeID in network.edges_dict_:
            if network.edges_dict_[edgeID].tension_>network.tension_cutoff_:
                radial_distance = network.edges_dict_[edgeID].radial_distance_
                if radial_distance>max(r_bins) or radial_distance<min(r_bins):
                    continue
                bin = min(r_bins, key=lambda x: abs(x - radial_distance))
                radius_dict[bin].append(network.edges_dict_[edgeID])



        #print(radius_dict)
        for radius in radius_dict:
            if len(radius_dict[radius]):
                total_length=0
                omega_rr=0
                for edge in radius_dict[radius]:
                    total_length+=edge.length_
                    omega_rr+=edge.length_*edge.eta_spherical_[0]*edge.eta_spherical_[0]
                #print(omega_rr/total_length,radius)
                omega_dict[radius].append(omega_rr/total_length)
        print("Processed ",dir)
    #print(omega_dict)
    r=[]
    o_rr=[]
    samplesize=[]
    err=[]

    for key in omega_dict:
        if len(omega_dict[key])>1:
            r.append(key)
            o_rr.append(np.mean(omega_dict[key]))
            samplesize.append(len(omega_dict[key]))
            err.append(stats.sem(omega_dict[key]))

    return r,o_rr,err,samplesize

#Calculate density of nodes

def get_radius_and_densities_and_err_arrays(dir_list,time=25000,r_bins=np.linspace(3.5,32,41)):

    densities_dict={i:[] for i in r_bins}
    shell_width=r_bins[1]-r_bins[0]
    for dir in dir_list:
        # A dictionary of nodes for the network where each connected node is assigned to a radial bin
        radius_dict={i:0 for i in r_bins}

        network=fiber_network(config_dir=dir,time=time)
        network.calculate_edge_attributes()
        #populate nodes_dict
        for edgeID in network.edges_dict_:
            
            radial_distance = network.edges_dict_[edgeID].radial_distance_
            if radial_distance>max(r_bins) or radial_distance<min(r_bins):
                continue
            bin = min(r_bins, key=lambda x: abs(x - radial_distance))
            radius_dict[bin]+=1



        #print(radius_dict)
        for radius in radius_dict:
            if radius_dict[radius]>0:
                densities_dict[radius].append(radius_dict[radius]/((4/3)*np.pi*(radius+shell_width/2)**3-(4/3)*np.pi*(radius-shell_width/2)**3))

        print("Processed ",dir)
    #print(omega_dict)
    r=[]
    densities=[]
    err=[]

    for key in densities_dict:
        if len(densities_dict[key])>1:
            r.append(key)
            densities.append(np.mean(densities_dict[key]))
            
            err.append(stats.sem(densities_dict[key]))

    return r,densities,err

def calculate_orientation_tensor_rr(network,r_bins=np.linspace(3.5,32,41)):
    omega=np.zeros((3,3))
    if not network.edges_initialized_:network.calculate_edge_attributes()

    tension_cutoff=network.tension_cutoff_
    total_length=0
    edge_counter=0
    #make a dictionary of *** high tension edges *** in shells of radii r_bins
    radius_dict={i:[] for i in r_bins}
    for edge in list(network.edges_dict_.values()):
        if edge.tension_>tension_cutoff:
            bin = min(r_bins, key=lambda x: abs(x - edge.radial_distance_))
            radius_dict[bin].append(edge)
    r=[]
    omega_rr=[]
    for radius in radius_dict:
        if len(radius_dict[radius]):
            r.append(radius)
            total_length=0
            omega=0
            for edge in radius_dict[radius]:
                total_length+=edge.length_
                omega+=edge.length_*edge.eta_spherical_[0]*edge.eta_spherical_[0]
            omega_rr.append(omega/total_length)

    
    return r,omega_rr
# a somewhat slow function that I need to read text files such as topo.txt, cellVolumes.txt, etc.
# for the spheroid class. I would suggest finding a way to use it only once per run.
def mk_timevals_dict(textfile):
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

def make_sample_dot_topo_from_fully_loaded_classes(output_dir,vertices,edges,polygons,cells):


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
            file.write("\n")
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

## BOOK-KEEPING CLASSES ###
## These classes are mainly for bookkeeping purposes (storing or analyzing output data). 
## Their use is not related to, and does not interefere with, classes of the same name in the main code...
##
##

class Vertex:
    def __init__(self, id):
        self.id_ = id
        self.og_id_ = None
        self.position_ = [0., 0., 0.]
        self.cells_ = []
        self.is_surface_=False
    def setPosition(self, position):
        self.position_ = list(position)
    def addCell(self, cell):
        if cell not in self.cells_:
            self.cells_.append(cell)
    def deleteCell(self, cell):
        self.cells_.remove(cell)

class Edge:
    def __init__(self, id):
        self.id_ = id
        self.og_id_ = None
        self.dumpOn_ = True
        self.vertices_ = []
        self.is_surface_=False
    def addVertex(self, vertex):
        self.vertices_.append(vertex)

class Polygon:
    def __init__(self, id):
        self.id_ = id
        self.og_id_=None
        self.dumpOn_ = True
        self.edges_ = []
        self.vertices_ = []
        self.type_=0
        self.is_surface_=False
        self.normal_=None
        self.center_=None
        self.area_=None
    def addEdge(self, edge):
        self.edges_.append(edge)
    def addVertex(self, vertex):
        self.vertices_.append(vertex)

class Cell:
    def __init__(self, id):
        self.id_ = id
        self.og_id_ = None
        self.polygons_ = []
        self.is_surface_=False
        self.surface_area_=None
        self.shape_index_=None

        self.volume_=None
        self.center_=None
        self.stress_tensor_=None
        self.type_=None

    def addPolygon(self, polygon):
        self.polygons_.append(polygon)
    def deletePolygon(self, polygon):
        self.polygons_.remove(polygon)

# An edge object corresponding to a single edge in the fiber network. To be initialized with 
# 1. a list of two node ids,
# 2. a tension, 
# 3. a dictionary of coordinates for the entire network,
# 4. Origin. Consider making this the COM of the spheroid - which may be different from [0,0,0]. 
#       This will be the origin used to calculate spherical coordinates and spherical unit vectors at self.node_0_. We will not be shifting the cartesian coordinates of the nodes.

# Usage: mk_edges_dict() stores edge information in the form of these object. Then the fiber_network class is equipped with the dictionary from mk_edges_dict().

class fiber_edge:
    def __init__(self,nodes,tension):
        self.nodes_=nodes
        self.tension_=tension
        #Calculate these values at the fiber network level
        self.coordinates_=None
        self.length_=None
        self.radial_distance_=None
        self.eta_cartesian_=None
        self.r_hat_=None
        self.theta_hat_=None
        self.phi_hat_=None
        self.eta_spherical_=None
        self.eta_cartesian_=None

#extract all the fiber network information from a single vtk file
class fiber_network:
    def __init__(self,config_dir="build/",time=20000,origin=np.array([0,0,0])):
        self.time_=time
        self.origin_=origin
        self.file_dir_=config_dir+"{:07d}".format(time)+".ECM.vtk"
        self.coordinates_dict_,self.edges_dict_=mk_coordinates_and_edges_dict(self.file_dir_)
        self.connected_nodes_=np.unique([self.edges_dict_[edgeID].nodes_ for edgeID in self.edges_dict_])
        self.edges_initialized_=False
        #self.lengths_=[self.edges_dict_[i].length_ for i in self.edges_dict_]
        #self.radial_distances_=[self.edges_dict_[i].radial_distance_ for i in self.edges_dict_]
        #self.max_radial_distance_=np.max(self.radial_distances_)
        #self.min_radial_distance_=np.min(self.radial_distances_)
        #self.edge_counter_=len(self.edges_dict_)
        self.tensions_=None
        self.tension_cutoff_=None
        #Orientation tensor for all high tension edges, should you be interested...
        #self.omega_spherical_,self.high_tension_edge_counter_=calculate_omega_spherical(edges_dict=self.edges_dict_,
        #                                                                                tension_cutoff=self.tension_cutoff_,
        #                                                                                radial_distance_max=self.max_radial_distance_,
        #                                                                                radial_distance_min=self.min_radial_distance_)
        #self.high_tension_shells_=pd.DataFrame(columns=["radius","edge counter","O_rr"])
        #shell_width=2
        #shell_radii=np.arange(self.min_radial_distance_, self.max_radial_distance_+shell_width, shell_width)[:-1]
        #for r in shell_radii:
            #print(r,r-shell_width/2,r+shell_width/2)
        #    omega_spherical,edge_counter=calculate_omega_spherical(self.edges_dict_,tension_cutoff=self.tension_cutoff_,radial_distance_min=r-shell_width/2,radial_distance_max=r+shell_width/2)
        #    if edge_counter>0:
                #print(omega_spherical,edge_counter)
        #        self.high_tension_shells_=self.high_tension_shells_.append({"radius":r,"edge counter":edge_counter,"O_rr":omega_spherical[0,0]},ignore_index=True)
        #
    def calculate_edge_attributes(self):
        self.tensions_=[]
        for edgeID in self.edges_dict_:
            self.edges_dict_[edgeID].nodes_=sorted(self.edges_dict_[edgeID].nodes_, key=lambda id: np.linalg.norm(self.coordinates_dict_[id]-self.origin_))
            node_0_adjusted_coordinates_=self.coordinates_dict_[self.edges_dict_[edgeID].nodes_[0]]-self.origin_        
            self.edges_dict_[edgeID].length_=np.linalg.norm(self.coordinates_dict_[self.edges_dict_[edgeID].nodes_[0]]
                                    -self.coordinates_dict_[self.edges_dict_[edgeID].nodes_[1]])
            
            self.edges_dict_[edgeID].radial_distance_=np.linalg.norm(node_0_adjusted_coordinates_)
            coordinates_1=self.coordinates_dict_[self.edges_dict_[edgeID].nodes_[1]]
            coordinates_0=self.coordinates_dict_[self.edges_dict_[edgeID].nodes_[0]]
            self.edges_dict_[edgeID].eta_cartesian_=(coordinates_1-coordinates_0)/self.edges_dict_[edgeID].length_

            theta=np.arctan2(node_0_adjusted_coordinates_[1],node_0_adjusted_coordinates_[0])
            phi=np.arccos(node_0_adjusted_coordinates_[2]/self.edges_dict_[edgeID].radial_distance_)
            self.edges_dict_[edgeID].r_hat_=np.array([np.sin(phi)*np.cos(theta),
                                                        np.sin(phi)*np.sin(theta),
                                                        np.cos(phi)])
            self.edges_dict_[edgeID].theta_hat_=np.array([-np.sin(theta),
                                    np.cos(theta),
                                    0])
            self.edges_dict_[edgeID].phi_hat_=np.array([np.cos(phi)*np.cos(theta),
                                                        np.cos(phi)*np.sin(theta),
                                                        -np.sin(phi)])
            self.edges_dict_[edgeID].eta_spherical_=[np.dot(self.edges_dict_[edgeID].eta_cartesian_,self.edges_dict_[edgeID].r_hat_),
                                                        np.dot(self.edges_dict_[edgeID].eta_cartesian_,self.edges_dict_[edgeID].theta_hat_),
                                                        np.dot(self.edges_dict_[edgeID].eta_cartesian_,self.edges_dict_[edgeID].phi_hat_)]
            self.tensions_.append(self.edges_dict_[edgeID].tension_)
        
        self.tension_cutoff_=np.std(self.tensions_)
        self.edges_initialized_=True

        return

class spheroid:
    def __init__(self,config_dir="build/",time=20000):
        self.time_=time
        self.config_dir_=config_dir
        self.vertices_=None
        self.edges_=None
        self.polygons_=None
        self.cells_=None
        self.center_=None
        self.shape_index_=None
        self.gamma_=None
        self.kv_=None
        self.timevals_topo_=mk_timevals_dict(config_dir+"topo.txt")
        self.timevals_cellShapeIndex_=mk_timevals_dict(config_dir+"cellShapeIndex.txt")
        self.timevals_cellVolume_=mk_timevals_dict(config_dir+"cellVolume.txt")
        self.timevals_cellCenter_=mk_timevals_dict(config_dir+"cellCenter.txt")

        self.loadconfig()



        
    
### loadconfig(): This function loads the currently existing topo.txt in build/ into the class variables...
### ... specifically, self.vertices_,self.edges_,self.polygons_,self.cells_,self.cellIDs_

    def loadconfig(self):
        time=self.time_
        stopwatch_start=timer.time()
        with open(self.config_dir_+"topo.txt", "r") as file:
            verticesFlag = False
            edgesFlag = False
            polygonsFlag = False
            cellsFlag = False
            start=self.timevals_topo_[time][0]
            end=self.timevals_topo_[time][1]
            lines=file.readlines()
            for i in range (start,end):
                line=lines[i]
                ##sometimes there is a blank last line
                if len(line) <= 1:
                    break

                lineSplit = line.split()
                if lineSplit[0] == "time":
                    self.vertices_ = {}
                    self.edges_ = {}
                    self.polygons_ = {}
                    self.cells_ = {}
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
                    continue
                
                if verticesFlag:
                    vertexID = int(lineSplit[0])
                    x = float(lineSplit[1])
                    y = float(lineSplit[2])
                    z = float(lineSplit[3])
                    vertex = Vertex(vertexID)
                    vertex.setPosition([x,y,z])
                    self.vertices_[vertexID] = vertex

                if edgesFlag:
                    edgeID = int(lineSplit[0])
                    edge = Edge(edgeID)

                    for i in range(1, len(lineSplit)):edge.addVertex(int(lineSplit[i]))
                    
                    self.edges_[edgeID]=edge
            
                if polygonsFlag:
                    polygonID = int(lineSplit[0])
                    polygon=Polygon(polygonID)
                    for i in range(1, len(lineSplit)):polygon.addEdge(int(lineSplit[i]))

                    self.polygons_[polygonID]=polygon
                
                if cellsFlag:
                    cellID = int(lineSplit[0])
                    cell=Cell(cellID)
                    for i in range(1, len(lineSplit)-1):cell.addPolygon(int(lineSplit[i]))
                    cell.type_=int(lineSplit[-1])

                    self.cells_[cellID]=cell
            print("time to load topo.txt: ", (timer.time()-stopwatch_start))
            stopwatch_start=timer.time()
            self.load_conf_file()
            
            self.load_surface_attributes()
            self.loadcellattributes()
            self.load_polygon_attributes()
            #self.dump_surface_vtk()
            self.load_spheroid_attributes()
            self.calculate_stress_tensor()

            print("time to other attributes: ",timer.time()-stopwatch_start)
        return
    def load_conf_file(self):
        with open(self.config_dir_+"conf","r") as file:
            lines=file.readlines()
            for line in lines:
                if (line.split()[0])=="s0": 
                    self.shape_index_=(float(line.split()[1]))#s0
                    self.gamma_=(float(line.split()[2]))#gamma
                if (line.split()[0])=="kv": 
                    self.kv_=(float(line.split()[1]))#kv
        return
    def loadcellattributes(self):
        time=self.time_
        with open (self.config_dir_+"cellShapeIndex.txt", "r") as file:
            start=self.timevals_cellShapeIndex_[time][0]
            end=self.timevals_cellShapeIndex_[time][1]
            lines=file.readlines()
            for i in range (start+1,end):
                line=lines[i]
                lineSplit = line.split()
                cellID=int(lineSplit[0])
                self.cells_[cellID].shape_index_=float(lineSplit[1])
        with open (self.config_dir_+"cellVolume.txt", "r") as file:
            start=self.timevals_cellVolume_[time][0]
            end=self.timevals_cellVolume_[time][1]
            lines=file.readlines()
            for i in range (start+1,end):
                line=lines[i]
                lineSplit = line.split()
                cellID=int(lineSplit[0])
                self.cells_[cellID].volume_=float(lineSplit[1])
        with open (self.config_dir_+"cellCenter.txt", "r") as file:
            start=self.timevals_cellCenter_[time][0]
            end=self.timevals_cellCenter_[time][1]
            lines=file.readlines()
            for i in range (start+1,end):
                line=lines[i]
                lineSplit = line.split()
                cellID=int(lineSplit[0])
                self.cells_[cellID].center_=[float(lineSplit[1]),float(lineSplit[2]),float(lineSplit[3])]
        self.calculate_cell_surface_areas()
        return
    
    # Here I use information from cellVolume.txt and cellShapeIndex.txt: A=s*V^2/3
    # instead of calculating areas from triangular polygon patches

    def calculate_cell_surface_areas(self):
        for cellID in self.cells_:
            if self.cells_[cellID].type_:
                self.cells_[cellID].surface_area_=self.cells_[cellID].volume_**(2/3)*self.cells_[cellID].shape_index_
        return

    '''
    ## Strictly speaking this is not necessary. Just use polygon center and add up triangle areas...
    def load_polygon_vertices(self):
        for polygonID in self.polygons_:
            if self.polygons_[polygonID].type_:
                #print(polygonID,len(sample.polygons_[polygonID].edges_))
                init_edge=self.polygons_[polygonID].edges_[0]
                vertices=[self.edges_[init_edge].vertices_[0]]
                #print(polygonID,len(self.polygons_[polygonID].edges_),vertices)
                for edgeID in self.polygons_[polygonID].edges_:
                    #(self.edges_[edgeID].vertices_)

                    while len(vertices)<len(self.polygons_[polygonID].edges_):
                        for edgeID in self.polygons_[polygonID].edges_:
                            if vertices[-1] in self.edges_[edgeID].vertices_:
                                next_vertex = self.edges_[edgeID].vertices_[0] if self.edges_[edgeID].vertices_[0] != vertices[-1] else self.edges_[edgeID].vertices_[1]
                                if next_vertex in vertices: pass
                                else: 
                                    vertices.append(next_vertex)
                #print(vertices)
                self.polygons_[polygonID].vertices_=vertices
        return
        '''
    # As defined in Okuda et al 
    def calculate_polygon_centers(self):
        for polygonID in self.polygons_:
            total_length=0
            self.polygons_[polygonID].center_=[0,0,0]
            for edgeID in self.polygons_[polygonID].edges_:
                node0=self.vertices_[self.edges_[edgeID].vertices_[0]].position_
                node1=self.vertices_[self.edges_[edgeID].vertices_[1]].position_
                length=np.linalg.norm(np.subtract(node0,node1))
                edge_center=np.add(node0,node1)/2
                self.polygons_[polygonID].center_=np.add(self.polygons_[polygonID].center_,edge_center*length)
                total_length+=length
            self.polygons_[polygonID].center_=self.polygons_[polygonID].center_/total_length
        return

    def calculate_polygon_areas(self):
        for polygonID in self.polygons_:

            self.polygons_[polygonID].area_=0
            for edgeID in self.polygons_[polygonID].edges_:
                v_i=np.subtract(self.vertices_[self.edges_[edgeID].vertices_[0]].position_,self.polygons_[polygonID].center_)
                v_j=np.subtract(self.vertices_[self.edges_[edgeID].vertices_[1]].position_,self.polygons_[polygonID].center_)
                self.polygons_[polygonID].area_+=0.5*np.linalg.norm(np.cross(v_i,v_j))
        return
    

                            
    def load_surface_attributes(self):
        for id in self.polygons_:
            cell_types=[]
    

            for cellID in self.cells_:
                if id in self.cells_[cellID].polygons_: cell_types.append(self.cells_[cellID].type_)


            if cell_types==[0,1] or cell_types==[1,0]:self.polygons_[id].is_surface_=True

        for id in self.cells_:
            if self.cells_[id].type_:
                for polygonID in self.cells_[id].polygons_:
                    if self.polygons_[polygonID].is_surface_:self.cells_[id].is_surface_=True
        for id in self.polygons_:
            if self.polygons_[id].is_surface_==True:
                for edgeID in self.polygons_[id].edges_:
                    self.edges_[edgeID].is_surface_=True
                    for vertexID in self.edges_[edgeID].vertices_:
                        self.vertices_[vertexID].is_surface_=True
    def load_polygon_attributes(self):
        for cellID in self.cells_:
            if self.cells_[cellID].type_:
                for polygonID in self.cells_[cellID].polygons_:
                    self.polygons_[polygonID].type_=1
        self.calculate_polygon_centers()
        self.calculate_polygon_areas()
        return

    def load_spheroid_attributes(self):
        self.center_=[]
        for cellID in self.cells_:
            if self.cells_[cellID].type_==1:
                self.center_.append(self.cells_[cellID].center_)
        self.center_=np.mean(self.center_,axis=0)
    def calculate_stress_tensor(self):
        for cellID in self.cells_:
            if self.cells_[cellID].type_:
                volume=self.cells_[cellID].volume_
                pressure=2*self.kv_*(volume-1)
                self.cells_[cellID].stress_tensor_=[[pressure,0,0],[0,pressure,0],[0,0,pressure]]
                for polygonID in self.cells_[cellID].polygons_:
                    if self.polygons_[polygonID].is_surface_:
                        tension=2*(self.cells_[cellID].surface_area_-self.shape_index_)+self.gamma_
                        #print("this is a surface polygon")
                    else:
                        neighboring_cellID=None
                        #find neighboring cell, if it exists...
                        for check_cellID in self.cells_:
                            if polygonID in self.cells_[check_cellID].polygons_:
                                neighboring_cellID=check_cellID

                        tension=2*(self.cells_[cellID].surface_area_-self.shape_index_)+2*(self.cells_[neighboring_cellID].surface_area_-self.shape_index_)
                    polygon_center=self.polygons_[polygonID].center_
                    for edgeID in self.polygons_[polygonID].edges_:
                        v0=np.subtract(self.vertices_[self.edges_[edgeID].vertices_[0]].position_,polygon_center)
                        v1=np.subtract(self.vertices_[self.edges_[edgeID].vertices_[1]].position_,polygon_center)
                        unit_vector=np.cross(v0,v1)/np.linalg.norm(np.cross(v0,v1))
                        #orient the unit vector so that it has a positive projection onto the vector from the cell center to the polygon center
                        vector_to_polygon_from_cell_center=np.subtract(polygon_center,self.cells_[cellID].center_)
                        #print(np.sign(np.dot(vector_to_polygon_from_cell_center,unit_vector)))
                        unit_vector*=np.sign(np.dot(vector_to_polygon_from_cell_center,unit_vector))
                        for i in range(3):
                            for j in range(3):
                                self.cells_[cellID].stress_tensor_[i][j]+=tension*unit_vector[i]*unit_vector[j]/(2*volume)
        return
    def dump_surface_vtk(self):

        tmp_cells={}
        tmp_polygons={}
        tmp_edges={}
        tmp_vertices={}

        for cellID in self.cells_:
            if self.cells_[cellID].is_surface_:
                tmp_cells[cellID]=self.cells_[cellID]
                
                for polygonID in self.cells_[cellID].polygons_:
                    if polygonID in tmp_polygons: pass
                    else: 
                        tmp_polygons[polygonID]=self.polygons_[polygonID]
                        for edgeID in self.polygons_[polygonID].edges_:
                            if edgeID in tmp_edges: pass
                            else:
                                tmp_edges[edgeID]=self.edges_[edgeID]
                                
                                for vertexID in self.edges_[edgeID].vertices_:
                                    if vertexID in tmp_vertices:pass
                                    else: 
                                        tmp_vertices[vertexID]=self.vertices_[vertexID]


        v_map=mapmaker(tmp_vertices)
        e_map=mapmaker(tmp_edges)
        p_map=mapmaker(tmp_polygons)



        new_vertices={}
        for i, key in enumerate(tmp_vertices):
            vertex=Vertex(i)
            vertex.setPosition(tmp_vertices[key].position_)
            new_vertices[i]=vertex

        new_edges={}
        for i, key in enumerate(tmp_edges):
            edge=Edge(i)
            for vertexID in tmp_edges[key].vertices_:
                edge.addVertex(v_map[vertexID])
            new_edges[i]=edge

        new_polygons={}
        for i, key in enumerate(tmp_polygons):
            polygon=Polygon(i)
            for edgeID in tmp_polygons[key].edges_:
                polygon.addEdge(e_map[edgeID])
            new_polygons[i]=polygon

        new_cells={}
        for i, key in enumerate(tmp_cells):
            cell=Cell(i)
            for polygonID in tmp_cells[key].polygons_:
                cell.addPolygon(p_map[polygonID])
            new_cells[i]=cell

        Points_=[]
        Polygons_=[]
        cellscalars=[]
        for key in new_vertices:
            Points_.append(new_vertices[key].position_)

        for i, cellID in enumerate(new_cells):
            
            for polygonID in new_cells[cellID].polygons_:
                cellscalars.append(i)
                tmp_polygon=[]
                for edgeID in new_polygons[polygonID].edges_:
                    for vertexID in new_edges[edgeID].vertices_:
                        tmp_polygon.append(vertexID)
                Polygons_.append(tmp_polygon)

        structure = PolyData(points=Points_,polygons=Polygons_)
        celldata = CellData(\
            Scalars(cellscalars,
                    name='cell_scalars'))
        vtk = VtkData(structure,celldata)
        vtk.tofile(self.config_dir_+"{:07d}.surface".format(self.time_),'ascii')        
        return

class periodic_box:
    def __init__(self,config_dir="build/",time=20000):
        self.time_=time
        self.config_dir_=config_dir
        self.vertices_=None
        self.edges_=None
        self.polygons_=None
        self.cells_=None
        self.center_=None
        self.shape_index_=None
        self.kv_=None
        self.timevals_topo_=mk_timevals_dict(config_dir+"topo.txt")
        self.timevals_cellShapeIndex_=mk_timevals_dict(config_dir+"cellShapeIndex.txt")
        self.timevals_cellVolume_=mk_timevals_dict(config_dir+"cellVolume.txt")
        self.timevals_cellCenter_=mk_timevals_dict(config_dir+"cellCenter.txt")

        self.loadconfig()



        
    
### loadconfig(): This function loads the currently existing topo.txt in build/ into the class variables...
### ... specifically, self.vertices_,self.edges_,self.polygons_,self.cells_,self.cellIDs_

    def loadconfig(self):
        time=self.time_
        stopwatch_start=timer.time()
        with open(self.config_dir_+"topo.txt", "r") as file:
            verticesFlag = False
            edgesFlag = False
            polygonsFlag = False
            cellsFlag = False
            start=self.timevals_topo_[time][0]
            end=self.timevals_topo_[time][1]
            lines=file.readlines()
            for i in range (start,end):
                line=lines[i]
                ##sometimes there is a blank last line
                if len(line) <= 1:
                    break

                lineSplit = line.split()
                if lineSplit[0] == "time":
                    self.vertices_ = {}
                    self.edges_ = {}
                    self.polygons_ = {}
                    self.cells_ = {}
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
                    continue
                
                if verticesFlag:
                    vertexID = int(lineSplit[0])
                    x = float(lineSplit[1])
                    y = float(lineSplit[2])
                    z = float(lineSplit[3])
                    vertex = Vertex(vertexID)
                    vertex.setPosition([x,y,z])
                    self.vertices_[vertexID] = vertex

                if edgesFlag:
                    edgeID = int(lineSplit[0])
                    edge = Edge(edgeID)

                    for i in range(1, len(lineSplit)):edge.addVertex(int(lineSplit[i]))
                    
                    self.edges_[edgeID]=edge
            
                if polygonsFlag:
                    polygonID = int(lineSplit[0])
                    polygon=Polygon(polygonID)
                    for i in range(1, len(lineSplit)):polygon.addEdge(int(lineSplit[i]))

                    self.polygons_[polygonID]=polygon
                
                if cellsFlag:
                    cellID = int(lineSplit[0])
                    cell=Cell(cellID)
                    for i in range(1, len(lineSplit)):cell.addPolygon(int(lineSplit[i]))

                    self.cells_[cellID]=cell
            print("time to load topo.txt: ", (timer.time()-stopwatch_start))
            stopwatch_start=timer.time()
            self.load_conf_file()
            self.loadcellattributes()
            #self.load_polygon_attributes()
            #self.dump_surface_vtk()



            print("time to other attributes: ",timer.time()-stopwatch_start)
        return
    def load_conf_file(self):
        with open(self.config_dir_+"conf","r") as file:
            lines=file.readlines()
            for line in lines:
                if (line.split()[0])=="s0": 
                    self.shape_index_=(float(line.split()[1]))#s0
                if (line.split()[0])=="kv": 
                    self.kv_=(float(line.split()[1]))#kv
        return
    def loadcellattributes(self):
        time=self.time_
        with open (self.config_dir_+"cellShapeIndex.txt", "r") as file:
            start=self.timevals_cellShapeIndex_[time][0]
            end=self.timevals_cellShapeIndex_[time][1]
            lines=file.readlines()
            for i in range (start+1,end):
                line=lines[i]
                lineSplit = line.split()
                cellID=int(lineSplit[0])
                self.cells_[cellID].shape_index_=float(lineSplit[1])
        with open (self.config_dir_+"cellVolume.txt", "r") as file:
            start=self.timevals_cellVolume_[time][0]
            end=self.timevals_cellVolume_[time][1]
            lines=file.readlines()
            for i in range (start+1,end):
                line=lines[i]
                lineSplit = line.split()
                cellID=int(lineSplit[0])
                self.cells_[cellID].volume_=float(lineSplit[1])
        with open (self.config_dir_+"cellCenter.txt", "r") as file:
            start=self.timevals_cellCenter_[time][0]
            end=self.timevals_cellCenter_[time][1]
            lines=file.readlines()
            for i in range (start+1,end):
                line=lines[i]
                lineSplit = line.split()
                cellID=int(lineSplit[0])
                self.cells_[cellID].center_=[float(lineSplit[1]),float(lineSplit[2]),float(lineSplit[3])]
        self.calculate_cell_surface_areas()
        return
    
    # Here I use information from cellVolume.txt and cellShapeIndex.txt: A=s*V^2/3
    # instead of calculating areas from triangular polygon patches

    def calculate_cell_surface_areas(self):
        for cellID in self.cells_:

            self.cells_[cellID].surface_area_=self.cells_[cellID].volume_**(2/3)*self.cells_[cellID].shape_index_
        return


    # As defined in Okuda et al 
    def calculate_polygon_centers(self):
        for polygonID in self.polygons_:
            total_length=0
            self.polygons_[polygonID].center_=[0,0,0]
            for edgeID in self.polygons_[polygonID].edges_:
                node0=self.vertices_[self.edges_[edgeID].vertices_[0]].position_
                node1=self.vertices_[self.edges_[edgeID].vertices_[1]].position_
                length=np.linalg.norm(np.subtract(node0,node1))
                edge_center=np.add(node0,node1)/2
                self.polygons_[polygonID].center_=np.add(self.polygons_[polygonID].center_,edge_center*length)
                total_length+=length
            self.polygons_[polygonID].center_=self.polygons_[polygonID].center_/total_length
        return

    def calculate_polygon_areas(self):
        for polygonID in self.polygons_:

            self.polygons_[polygonID].area_=0
            for edgeID in self.polygons_[polygonID].edges_:
                v_i=np.subtract(self.vertices_[self.edges_[edgeID].vertices_[0]].position_,self.polygons_[polygonID].center_)
                v_j=np.subtract(self.vertices_[self.edges_[edgeID].vertices_[1]].position_,self.polygons_[polygonID].center_)
                self.polygons_[polygonID].area_+=0.5*np.linalg.norm(np.cross(v_i,v_j))
        return
    

                            

    def load_polygon_attributes(self):
        self.calculate_polygon_centers()
        self.calculate_polygon_areas()
        return

    def dump_surface_vtk(self):

        tmp_cells={}
        tmp_polygons={}
        tmp_edges={}
        tmp_vertices={}

        for cellID in self.cells_:
            if self.cells_[cellID].is_surface_:
                tmp_cells[cellID]=self.cells_[cellID]
                
                for polygonID in self.cells_[cellID].polygons_:
                    if polygonID in tmp_polygons: pass
                    else: 
                        tmp_polygons[polygonID]=self.polygons_[polygonID]
                        for edgeID in self.polygons_[polygonID].edges_:
                            if edgeID in tmp_edges: pass
                            else:
                                tmp_edges[edgeID]=self.edges_[edgeID]
                                
                                for vertexID in self.edges_[edgeID].vertices_:
                                    if vertexID in tmp_vertices:pass
                                    else: 
                                        tmp_vertices[vertexID]=self.vertices_[vertexID]


        v_map=mapmaker(tmp_vertices)
        e_map=mapmaker(tmp_edges)
        p_map=mapmaker(tmp_polygons)



        new_vertices={}
        for i, key in enumerate(tmp_vertices):
            vertex=Vertex(i)
            vertex.setPosition(tmp_vertices[key].position_)
            new_vertices[i]=vertex

        new_edges={}
        for i, key in enumerate(tmp_edges):
            edge=Edge(i)
            for vertexID in tmp_edges[key].vertices_:
                edge.addVertex(v_map[vertexID])
            new_edges[i]=edge

        new_polygons={}
        for i, key in enumerate(tmp_polygons):
            polygon=Polygon(i)
            for edgeID in tmp_polygons[key].edges_:
                polygon.addEdge(e_map[edgeID])
            new_polygons[i]=polygon

        new_cells={}
        for i, key in enumerate(tmp_cells):
            cell=Cell(i)
            for polygonID in tmp_cells[key].polygons_:
                cell.addPolygon(p_map[polygonID])
            new_cells[i]=cell

        Points_=[]
        Polygons_=[]
        cellscalars=[]
        for key in new_vertices:
            Points_.append(new_vertices[key].position_)

        for i, cellID in enumerate(new_cells):
            
            for polygonID in new_cells[cellID].polygons_:
                cellscalars.append(i)
                tmp_polygon=[]
                for edgeID in new_polygons[polygonID].edges_:
                    for vertexID in new_edges[edgeID].vertices_:
                        tmp_polygon.append(vertexID)
                Polygons_.append(tmp_polygon)

        structure = PolyData(points=Points_,polygons=Polygons_)
        celldata = CellData(\
            Scalars(cellscalars,
                    name='cell_scalars'))
        vtk = VtkData(structure,celldata)
        vtk.tofile(self.config_dir_+"{:07d}.surface".format(self.time_),'ascii')        
        return
    
##
##
### END OF BOOK-KEEPING CLASSES ###

