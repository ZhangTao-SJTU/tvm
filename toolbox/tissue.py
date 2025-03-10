import copy
import numpy as np
import os
from toolbox import functions, topology
# Class representing the tissue sample. 
# For now, this can handle either spheroids or 3-torus periodic tissue.

class Sample:
    def __init__(self):
        self.tissueType_ = None
        self.time_ = None
        self.config_dir_ = None
        self.file_ = None
        self.vertices_:dict[int,topology.Vertex] = {}
        self.edges_:dict[int,topology.Edge] = {}
        self.polygons_:dict[int,topology.Polygon] = {}
        self.cells_:dict[int,topology.Cell] = {}
        self.cell_neighbors_ = {}
        self.s0_ = None
        self.gamma_ = None
        self.kv_ = None
        self.boxSize_ = None
        self.sample_total_volume_ = None

    @classmethod
    def from_topology(cls,vertices,edges,polygons,cells):
        sample = cls()
        sample.vertices_ = vertices
        sample.edges_ = edges
        sample.polygons_ = polygons
        sample.cells_ = cells
        return sample
    
    # Given an unordered set of edges that make up a polygon,
    # We first load the node IDs for the edges into an list of the form
    # [[edge0node0,edge0node1],[edge1node0,edge1node1],...]

    # This function take such an input array, then returns an array of node ids 
    # [node0,node1..] arranged in either clockwise or anticlockwise order 
    # around the polygon.

    # Note: this function assumes that the polygon is simply connected. 
    # It still returns a result if the polygon is not, but it will print a warning.
    @staticmethod
    def resolve_polygon_edge_connectivity(array:list[list[int]]) -> list[int]:
        for sublist in array:
            if not len(sublist) == 2:
                raise ValueError("Input list must be of the form [[int,int],[int,int],...]")

        test = copy.deepcopy(array)
        result = [test[0][0],test[0][1]]
        findme = result[-1]
        test.pop(0)

        while len(test):
            flag=False
            for i,item in enumerate(test):
                if findme in item:
                    flag=True
                    if item[0]==findme: result.append(item[1])
                    else: result.append(item[0])
                    findme=result[-1]
                    test.pop(i)
                    break
            if not flag:
                print("Warning: (while using functions.arrange_polygon())")
                print("Polygon is not simply connected! Take a look:")
                print(array)
                result.append(test[0][0])
                result.append(test[0][1])
                findme = result[-1]
                test.pop(0)
        del result[-1]
        return result
    def set_file(self,file):
        if not os.path.isfile(self.config_dir_ + file):
            raise ValueError("Error in tissueSample.Sample: file must exist in self._config_dir")
        self.file_ = self.config_dir_ + file
   
    def load_config(self):
        if self.file_ is None:
            print("self.file_ not set. If you are loading info from topo.txt,")
            print("use self.load_config_from_topo()")
            raise ValueError("Error using tissueSample.Sample.load_config(): self.file_ not set.")
        # Load the topology. When loading cell polygons, note that there are
        # virtual cells (type=0) and real cells (type=1).
        self.vertices_ = {}
        self.edges_ = {}
        self.polygons_ = {}
        self.cells_ = {}
        with open(self.file_, "r") as file:
            verticesFlag = False
            edgesFlag = False
            polygonsFlag = False
            cellsFlag = False
            for line in file.readlines():
                #skip blank lines
                if not len(line.split()):
                    continue
                lineSplit = line.split()                
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
                    vertex = topology.Vertex(vertexID)
                    vertex.setPosition([x,y,z])
                    self.vertices_[vertexID] = vertex
                if edgesFlag:
                    edgeID = int(lineSplit[0])
                    edge = topology.Edge(edgeID)
                    for i in range(1, len(lineSplit)):
                        edge.addVertex(int(lineSplit[i]))
                    self.edges_[edgeID] = edge
                if polygonsFlag:
                    polygonID = int(lineSplit[0])
                    polygon = topology.Polygon(polygonID)
                    for i in range(1, len(lineSplit)):
                        polygon.addEdge(int(lineSplit[i]))
                    self.polygons_[polygonID] = polygon
                if cellsFlag:
                    cellID = int(lineSplit[0])
                    cell = topology.Cell(cellID)
                    # For spheroids, under cellsFlag in topo.txt, 
                    # the last element in the list of polygons for each cellID 
                    # is the cell type.
                    # 0 for virtual cells, 1 for real cells.
                    if self.tissueType_ == "spheroid":
                        for i in range(1, len(lineSplit)-1):
                            cell.addPolygon(int(lineSplit[i]))
                        cell.type_ = int(lineSplit[-1])

                    elif self.tissueType_ == "periodic":
                        for i in range(1, len(lineSplit)):
                            cell.addPolygon(int(lineSplit[i]))
                    self.cells_[cellID] = cell
    
    def load_conf_file(self):
        if not os.path.isfile(self.config_dir_ + "conf"):
            raise ValueError("conf file not found in config_dir")
        with open(self.config_dir_ + "conf","r") as file:
            lines = file.readlines()
            for line in lines:
                if (line.split()[0]) == "s0": 
                    self.s0_ = (float(line.split()[1]))#s0
                    if self.tissueType_ == "spheroid":
                        self.gamma_=(float(line.split()[2]))#gamma
                    for cellID,cell in self.cells_.items():
                        cell.v0_ = 1
                        cell.s0_ = copy.deepcopy(self.s0_)
                if (line.split()[0]) == "kv": 
                    self.kv_=(float(line.split()[1]))#kv
                if (self.tissueType_ == "periodic") and (line.split()[0]) == "box":
                    self.boxSize_ = float(line.split()[1])
    
    def load_cell_vertices(self):
        # Note: for spheroids we only care about real cells (type 1)
        # For periodic tissue, we care about all cells.
        for cellID, cell in self.cells_.items():
            if self.tissueType_ == "spheroid" and not cell.type_:
                continue
            cell.vertices_ = []
            for polygonID in cell.polygons_:
                for edgeID in self.polygons_[polygonID].edges_:
                    for vertexID in self.edges_[edgeID].vertices_:
                        cell.vertices_.append(vertexID)
            cell.vertices_ = np.unique(cell.vertices_)
    
    # Calculate cell centers from vertex positions.
    def calculate_cell_centers(self):
        for cellID, cell in self.cells_.items():
            if self.tissueType_ == "spheroid" and not cell.type_:
                continue
            if self.tissueType_ == "periodic" and cell.crossBoundary_:
                continue
            cell.center_ = np.zeros(3)
            for vertexID in cell.vertices_:
                cell.center_ = np.add(cell.center_, self.vertices_[vertexID].position_)
            cell.center_ = np.divide(cell.center_, len(cell.vertices_))
               
    # Load cell attributes from {self.time_}.cellInfo.txt
    # Create this file if it does not exist.
    def load_cell_attributes(self):
        # Load cell vertices. 
        self.load_cell_vertices()
        # Create {time}.cellInfo.txt if it does not exist
        if not os.path.isfile(self.config_dir_ 
                              + "{:07d}".format(self.time_) 
                              + ".cellInfo.txt"):
            functions.writeTimeCellInfo(self.config_dir_, self.time_)
        with open (self.config_dir_ 
                   + "{:07d}".format(self.time_) 
                   + ".cellInfo.txt", "r") as file:
            for i, line in enumerate(file.readlines()):
                # skip the header line
                if i == 0: continue
                # skip blank lines
                if not len(line.split()): continue
                lineSplit = line.split()
                cellID = int(lineSplit[0])
                self.cells_[cellID].center_ = [float(lineSplit[1]),
                                               float(lineSplit[2]),
                                               float(lineSplit[3])]
                self.cells_[cellID].volume_ = float(lineSplit[4])
                self.cells_[cellID].shape_index_ = float(lineSplit[5])
        # Calculating cell surface area from cell volume and shape index:
        for cellID, cell in self.cells_.items():
            if self.tissueType_ == "spheroid" and not cell.type_:
                continue
            cell.surface_area_ = (cell.shape_index_ * pow(cell.volume_,2/3))

    # As defined in Okuda et al 
    def calculate_polygon_centers_and_perimeters(self):
        for polygonID, polygon in self.polygons_.items():
            total_length = 0
            polygon.center_ = [0,0,0]
            for edgeID in polygon.edges_:
                length = np.linalg.norm(
                    np.subtract(
                        self.vertices_[self.edges_[edgeID].vertices_[0]].position_,
                        self.vertices_[self.edges_[edgeID].vertices_[1]].position_))
                edge_center = np.multiply(
                    1/2,np.add(
                        self.vertices_[self.edges_[edgeID].vertices_[0]].position_,
                        self.vertices_[self.edges_[edgeID].vertices_[1]].position_))
                polygon.center_ = np.add(self.polygons_[polygonID].center_,
                                         np.multiply(edge_center,length))
                total_length += length
            polygon.center_ = np.multiply(polygon.center_,
                                          1/total_length)
            polygon.perimeter_ = total_length

    # Equip polygon.center_ with the COM center coordinates (average of vertex positions)    
    def calculate_COM_polygon_centers(self):
        for polygonID, polygon in self.polygons_.items():
            if self.tissueType_ == "periodic" and polygon.crossBoundary_:
                continue
            if not len(polygon.vertices_):
                continue
            polygon.center_ = np.zeros(3)
            for vertexID in polygon.vertices_:
                polygon.center_ = np.add(polygon.center_, self.vertices_[vertexID].position_)
            polygon.center_ = np.divide(polygon.center_, len(polygon.vertices_))


    # Calculate polygon areas by breaking up into triangular patches.
    # Note that this requires that we first calculate polygon centers.
    def calculate_polygon_areas(self):
        for polygonID, polygon in self.polygons_.items():
            if self.tissueType_ == "periodic" and polygon.crossBoundary_:
                continue
            polygon.area_=0
            for edgeID in polygon.edges_:
                v_i = np.subtract(
                    self.vertices_[self.edges_[edgeID].vertices_[0]].position_,
                    polygon.center_)
                v_j = np.subtract(
                    self.vertices_[self.edges_[edgeID].vertices_[1]].position_,
                    polygon.center_)
                polygon.area_ += 0.5 * np.linalg.norm(np.cross(v_i, v_j))

    def calculate_cell_surface_areas(self):
        for cellID, cell in self.cells_.items():
            if self.tissueType_ == "spheroid" and not cell.type_:
                continue
            if self.tissueType_ == "periodic" and cell.crossBoundary_:
                continue
            cell.surface_area_ = 0
            for polygonID in cell.polygons_:
                cell.surface_area_ += self.polygons_[polygonID].area_
        return
    
    def calculate_cell_shape_indices(self):
        for cellID, cell in self.cells_.items():
            if self.tissueType_ == "spheroid" and not cell.type_:
                continue
            if self.tissueType_ == "periodic" and cell.crossBoundary_:
                continue
            cell.shape_index_ = cell.surface_area_ / pow(cell.volume_,2/3)

    def calculate_cell_volumes(self):
        self.arrange_polygon_vertices()
        for cellID, cell in self.cells_.items():
            if self.tissueType_ == "spheroid" and not cell.type_:
                continue
            if self.tissueType_ == "periodic" and cell.crossBoundary_:
                continue
            cell.volume_ = 0
            for polygonID in cell.polygons_:
                vertices_this_polygon = self.polygons_[polygonID].vertices_
                for i, vertexID in enumerate(vertices_this_polygon):
                    v0 = self.vertices_[vertexID].position_
                    v1ID = vertices_this_polygon[(i+1)%len(vertices_this_polygon)]
                    v1 = self.vertices_[v1ID].position_
                    v2 = self.polygons_[polygonID].center_
                    cell.volume_ += (1/6)*abs(
                        np.dot(
                            np.cross(
                                np.subtract(v0,cell.center_),
                                np.subtract(v1,cell.center_)),
                                np.subtract(v2,cell.center_)))

    # the static method resolve_polygon_edge_connectivity() 
    # orders the polygon vertices for a single polygon.

    # arrange_polygon_vertices() uses this method to order the vertices
    # for all relevant polygons in the sample.
    def arrange_polygon_vertices(self):
        for cellID, cell in self.cells_.items():
            if self.tissueType_ == "spheroid" and not cell.type_:
                continue
            # if self.tissueType_ == "periodic" and cell.crossBoundary_: 
            #     continue
            for polygonID in cell.polygons_:
                if self.tissueType_ == "periodic" and self.polygons_[polygonID].crossBoundary_:
                    continue
                tmp_vertices=[]
                for edgeID in self.polygons_[polygonID].edges_:
                    tmp_vertices.append(self.edges_[edgeID].vertices_)
                self.polygons_[polygonID].vertices_ = self.resolve_polygon_edge_connectivity(tmp_vertices)
    
    def evaluate_cell_neighbors(self):
        self.cell_neighbors_ = {}
        for cellID,cell in self.cells_.items():
            if self.tissueType_ == "spheroid" and not cell.type_:
                continue
            self.cell_neighbors_[cellID] = []
        for i, cellID in enumerate(self.cell_neighbors_):
            for polygonID in self.cells_[cellID].polygons_:
                for j in range(i+1, len(self.cell_neighbors_)):
                    test_cellID = list(self.cell_neighbors_.keys())[j]
                    if not polygonID in self.cells_[test_cellID].polygons_:
                        continue
                    self.cell_neighbors_[cellID].append(test_cellID)
                    self.cell_neighbors_[test_cellID].append(cellID)
                    # if polygonID in spheroid.cells_[list(cellID_to_neighbors.keys())[j]].polygons_:
                    #     cellID_to_neighbors[cellID].append(list(cellID_to_neighbors.keys())[j])
                    #     cellID_to_neighbors[list(cellID_to_neighbors.keys())[j]].append(cellID)
    
    def extract_cell(self,cellID):
        cell = copy.deepcopy(self.cells_[cellID])
        vertices = {}
        edges = {}
        polygons = {}
        cells = {cellID:cell}
        for polygonID in cell.polygons_:
            polygons[polygonID]=copy.deepcopy(self.polygons_[polygonID])
            for edgeID in polygons[polygonID].edges_:
                edges[edgeID]=copy.deepcopy(self.edges_[edgeID])
                for vertexID in edges[edgeID].vertices_:
                    vertices[vertexID]=copy.deepcopy(self.vertices_[vertexID])
        for polygonID in cell.polygons_:
            polygon = polygons[polygonID]
            sum_of_cross_product = np.array([0,0,0])
            for i in range(len(polygon.vertices_)):
                v1 = polygon.vertices_[i]
                v2 = polygon.vertices_[(i+1)%len(polygon.vertices_)]
                # print("     ", v1, v2)
                # print("before",vertices[v1].position_)
                vector1 = np.subtract(vertices[v1].position_, cell.center_)
                vector2 = np.subtract(vertices[v2].position_, cell.center_)
                # print("after",vertices[v1].position_)
                cross_product = np.cross(vector1, vector2)
                sum_of_cross_product = np.add(sum_of_cross_product, cross_product)

            sign = np.sign(np.dot(sum_of_cross_product,np.subtract(polygon.center_, cell.center_)))
            if sign==-1:
                polygon.vertices_=polygon.vertices_[::-1]

        for vertexID,vertex in vertices.items():
            for polygonID,polygon in polygons.items():
                if vertexID in polygon.vertices_:
                    vertex.polygons_.append(polygonID)
        single_cell = Sample.from_topology(vertices,edges,polygons,cells)
        single_cell.config_dir_ = self.config_dir_
        single_cell.time_ = self.time_
        single_cell.s0_ = cell.s0_
        single_cell.v0_ = cell.v0_
        single_cell.kv_ = self.kv_
        return single_cell

    def write_cell_collection_vtk(self,cells_array,filename):
        vertices = []
        polygons = []  
        for cellID in cells_array:
            cell = self.cells_[cellID]
            if cell.crossBoundary_:
                continue
            for polygonID in cell.polygons_:
                polygon = self.polygons_[polygonID]
                if not polygon.vertices_:
                    continue
                if polygon.crossBoundary_:
                    continue
                polygons.append(polygonID)
                for vertexID in polygon.vertices_:
                    vertices.append(vertexID)
        vertices = np.unique(vertices)
        polygons = np.unique(polygons)
        total_polygons = len(polygons)
        total_polygon_data_points = 0
        for polygonID in polygons:
            total_polygon_data_points += len(self.polygons_[polygonID].vertices_) + 1
        v_map = functions.mapmaker(vertices)
        with open("{}{}".format(self.config_dir_,filename),"w") as f:
            f.write("# vtk DataFile Version 2.0\n")
            f.write("polydata\n")
            f.write("ASCII\n")
            f.write("DATASET POLYDATA\n")
            f.write("POINTS {} double\n".format(len(vertices)))
            for vertexID in vertices:
                for i in range(3):
                    f.write("{} ".format(self.vertices_[vertexID].position_[i]))
                f.write("\n")
            f.write("POLYGONS {} {}\n".format(total_polygons,total_polygon_data_points))
            for polygonID in polygons:
                polygon = self.polygons_[polygonID]
                f.write("{} ".format(len(polygon.vertices_)))
                for vertexID in polygon.vertices_:
                    f.write("{} ".format(v_map[vertexID]))
                f.write("\n")

# class SingleCell:
#     def __init__(self,vertices,edges,polygons,cell):
#         self.config_dir_ = None
#         self.time_ = None
#         self.id_ = cell.id_
#         self.vertices_ = vertices
#         self.edges_ = edges
#         self.polygons_ = polygons
#         self.cells_ = {cell.id_:cell}
#         self.center_ = cell.center_
#         self.s0_ = None
#         self.kv_ = None
#         self.gamma_ = None

#     def update_center(self):
#         center = []
#         for _,vertex in self.vertices_.items():
#             center.append(vertex.position_)
#         center = np.mean(center,axis = 0)
#         self.center_ = center
#         self.cells_[self.id_].center_ = center
#     def dump_vtk(self,filename):
#         v_map = functions.mapmaker(self.vertices_)
#         # e_map = functions.mapmaker(self.edges_)
#         # p_map = functions.mapmaker(self.polygons_)
#         totalPolygonDataPoints = 0
#         for polygonID,polygon in self.polygons_.items():
#             totalPolygonDataPoints += (len(polygon.vertices_) + 1)
#         # with open(self.config_dir_
#                 #   + "{:07d}.cell_{}.vtk".format(self.time_,self.id_),'w') as file:
#         with open(filename, "w") as file: 
#             file.write("# vtk DataFile Version 2.0\n")
#             file.write("polydata\n")
#             file.write("ASCII\n")
#             file.write("DATASET POLYDATA\n")
#             file.write("POINTS {} double\n".format(len(self.vertices_)))
            
#             for vertexID, vertex in self.vertices_.items():
#                 file.write("{:12.5e} {:12.5e} {:12.5e}\n".format(
#                     vertex.position_[0],
#                     vertex.position_[1],
#                     vertex.position_[2]))
#             file.write("POLYGONS {} {}\n".format(
#                 len(self.polygons_),totalPolygonDataPoints))
        
#             for polygonID, polygon in self.polygons_.items():
#                 file.write("{:<7d}".format(len(polygon.edges_)))
#                 for vID in polygon.vertices_: 
#                     file.write("{:<7d}".format(v_map[vID]))
#                 file.write("\n")
#             file.close()
