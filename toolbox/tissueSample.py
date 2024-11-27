import copy
import numpy as np
import os
from toolbox import functions, topology

# Class representing the tissue sample. 
# For now, this can handle either spheroids or 3-torus periodic tissue.

class Sample:
    # def __init__(self,
    #              configDir:str = "samples/",
    #              simulationTime = 20000,
    #              tissueType:str = "spheroid"):
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

        # calculate from cellVolume.txt
        self.sample_total_volume_ = None
        # If the sample is a spheroid, these are relevant quantities
        # if self.tissueType_ == "spheroid":
        self.sample_center_ = None
        self.sample_surface_area_ = None
        # self.load_config()
        # self.load_conf_file()
        # self.load_cell_attributes()
        # self.calculate_cell_surface_areas()
        # if self.tissueType_ == "periodic":
        #     self.loadCrossBoundaryAttributes()
        # self.calculate_polygon_centers_and_perimeters()
        # self.calculate_polygon_areas()
        # self.arrange_polygon_vertices()
        # if self.tissueType_ == "spheroid":
        #     self.identify_surface_polygons_and_cells()
        #     self.load_spheroid_attributes()

        return

    @classmethod
    def from_topology(cls,vertices,edges,polygons,cells):
        sample = cls()
        sample.vertices_ = vertices
        sample.edges_ = edges
        sample.polygons_ = polygons
        sample.cells_ = cells
        # sample.center_ = cell.center_
        sample.s0_ = None
        sample.kv_ = None
        return sample
    
    @classmethod
    # def periodic_tissue(cls,configDir:str = "samples/",simulationTime = 20000):
    def periodic_tissue(cls, config_dir, input_filename = "sample.topo"):
        # Validate the input: 
        if not config_dir.endswith("/"):
            print("config_dir must end with a '/'. Attempting to fix this...")
            config_dir = config_dir + "/"
        if not os.path.isdir(config_dir):
            raise ValueError("config_dir must be a valid directory")
        if not os.path.isfile(config_dir+input_filename):
            raise ValueError("input_filename must exist and be in the config_dir")
        
        sample = cls()
        sample.tissueType_ = "periodic"
        sample.time_ = 0
        sample.config_dir_ = config_dir

        sample.load_config()
        sample.load_conf_file()
        sample.load_cell_attributes()
        # sample.calculate_cell_surface_areas()
        sample.loadCrossBoundaryAttributes()
        sample.calculate_polygon_centers_and_perimeters()
        sample.calculate_polygon_areas()
        sample.arrange_polygon_vertices()
        return sample
    
    # @classmethod 
    # def periodic_tissue(cls,configDir:str = "samples/",simulationTime = 20000):
    #     # Validate the input: 
    #     if not configDir.endswith("/"):
    #         raise ValueError("config_dir must end with a '/'")
        
    #     if not os.path.isdir(configDir):
    #         raise ValueError("config_dir must be a valid directory")
        
    #     cls.tissueType_ = "periodic"
    #     cls.time_ = simulationTime
    #     cls.config_dir_ = configDir

    #     Sample.load_config(cls)
    #     Sample.load_conf_file(cls)
    #     Sample.load_cell_attributes(cls)
    #     Sample.calculate_cell_surface_areas(cls)
    #     Sample.loadCrossBoundaryAttributes(cls)
    #     Sample.calculate_polygon_centers_and_perimeters(cls)
    #     Sample.calculate_polygon_areas(cls)
    #     Sample.arrange_polygon_vertices(cls)
    #     return cls
    # loadconfig(): given self.time_, this function first checks if
    # {time}.topo.txt exists in self.config_dir_. If not, it creates this file -
    # that is, it mines the topology at this time from topo.txt
    
    # Similarly, it checks and ensures the existence of {time}.cellInfo.txt
    # in self.config_dir_.
    
    # Then it loads this data into the dictionaries
    # self.vertices_,self.edges_,self.polygons_,self.cells_,self.cellIDs_

    def load_config(self):
        if self.file_ is None:
            if os.path.isfile(self.config_dir_ + "topo.txt"):
                self.file_ = "{}{:07}.topo.txt".format(self.config_dir_,self.time_)
        if not os.path.isfile(self.file_):
            if os.path.isfile(self.config_dir_ + "topo.txt"):    
                functions.make_time_topo(self.config_dir_, self.time_)
            else:
                raise ValueError("topo.txt not found in config_dir.")

        # Load the topology. When loading cell polygons, note that there are
        # virtual cells (type=0) and real cells (type=1).
        with open(self.file_, "r") as file:
            verticesFlag = False
            edgesFlag = False
            polygonsFlag = False
            cellsFlag = False

            for line in file.readlines():
                #skip blank lines
                if not len(line.split()): continue

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
                    # In other words, if you do a boolean cast on cell.type_,
                    # False is a virtual cell, True is a real cell.

                    if self.tissueType_ == "spheroid":
                        for i in range(1, len(lineSplit)-1):
                            cell.addPolygon(int(lineSplit[i]))
                        cell.type_ = int(lineSplit[-1])

                    elif self.tissueType_ == "periodic":
                        for i in range(1, len(lineSplit)):
                            cell.addPolygon(int(lineSplit[i]))

                    self.cells_[cellID] = cell

        return
    
    def load_conf_file(self):
        with open(self.config_dir_ + "conf","r") as file:
            lines = file.readlines()
            for line in lines:
                if (line.split()[0]) == "s0": 
                    self.s0_ = (float(line.split()[1]))#s0
                    if self.tissueType_ == "spheroid":
                        self.gamma_=(float(line.split()[2]))#gamma
                if (line.split()[0]) == "kv": 
                    self.kv_=(float(line.split()[1]))#kv
                if (self.tissueType_ == "periodic") and (line.split()[0]) == "box":
                    self.boxSize_ = float(line.split()[1])
        return
    
    # Load cell attributes from {self.time_}.cellInfo.txt
    # Create this file if it does not exist.
    
    def load_cell_attributes(self):
        # Load cell vertices. Note: we only care about real cells (type 1)
        for cellID, cell in self.cells_.items():
            if bool(cell.type_):
                # Cell.vertices_ is a list of vertex ids that make up the cell.
                cell.vertices_ = []
                for polygonID in cell.polygons_:
                    for edgeID in self.polygons_[polygonID].edges_:
                        for vertexID in self.edges_[edgeID].vertices_:
                            cell.vertices_.append(vertexID)
                cell.vertices_ = np.unique(cell.vertices_)

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

        return
    
    # Here I use information originally from 
    # cellVolume.txt and cellShapeIndex.txt: 
    # A=s*V^2/3
    # instead of calculating areas from triangular polygon patches

    def calculate_cell_surface_areas(self):
        for cellID, cell in self.cells_.items():

            if bool(cell.type_):
                cell.surface_area_ = (cell.shape_index_
                                      * pow(cell.volume_,2/3))
        return

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
        return

    # Calculate polygon areas by breaking up into triangular patches.
    # Note that this requires that we first calculate polygon centers.
    def calculate_polygon_areas(self):
        for polygonID, polygon in self.polygons_.items():
            polygon.area_=0
            for edgeID in polygon.edges_:
                v_i = np.subtract(
                    self.vertices_[self.edges_[edgeID].vertices_[0]].position_,
                    polygon.center_)
                v_j = np.subtract(
                    self.vertices_[self.edges_[edgeID].vertices_[1]].position_,
                    polygon.center_)
                polygon.area_ += 0.5 * np.linalg.norm(np.cross(v_i, v_j))
        return
    
    def calculate_cell_centers(self):
        for cellID, cell in self.cells_.items():
            # if not cell.type_:
            #     continue
            cell.center_ = np.zeros(3)
            for vertexID in cell.vertices_:
                cell.center_ = np.add(cell.center_, self.vertices_[vertexID].position_)
            cell.center_ = np.divide(cell.center_, len(cell.vertices_))
        return
    
    def calculate_cell_surface_areas(self):
        for cellID, cell in self.cells_.items():
            # if not cell.type_:
            #     continue
            cell.surface_area_ = 0
            for polygonID in cell.polygons_:
                cell.surface_area_ += self.polygons_[polygonID].area_
        return
    
    def calculate_cell_volumes(self):
        self.arrange_polygon_vertices()
        for cellID, cell in self.cells_.items():
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
        return

    # identify polygon and cells that are on the surface of the spheroid                     
    def identify_surface_polygons_and_cells(self):
        for polygonID, polygon in self.polygons_.items():
            cell_types = []

            for cellID, cell in self.cells_.items():
                if polygonID in cell.polygons_:
                    cell_types.append(cell.type_)

            if set(cell_types) == {0,1}:
                polygon.is_surface_ = True

        for cellID, cell in self.cells_.items():
            if bool(cell.type_):
                for polygonID in cell.polygons_:
                    if self.polygons_[polygonID].is_surface_:
                        cell.is_surface_ = True
        return


    # this function both loads and arranges the vertices of the polygons. Note that
    # the vertices may either be arranged clockwise or counterclockwise.
    # dont resolve this here, because each polygon is shared between two cells,
    # and for either cell the anticlockwise arrangement of
    # vertices on this polygon is opposed.

    def arrange_polygon_vertices(self):
        for cellID, cell in self.cells_.items():
            if bool(cell.type_) or self.tissueType_ == "periodic":
                for polygonID in cell.polygons_:
                    tmp_vertices=[]
                    for edgeID in self.polygons_[polygonID].edges_:
                        tmp_vertices.append(self.edges_[edgeID].vertices_)
                    self.polygons_[polygonID].vertices_ = functions.arrange_polygon(tmp_vertices)
                    
        return
    
    def load_spheroid_attributes(self):
        self.sample_center_ = []
        self.sample_surface_area_ = 0
        self.sample_total_volume_ = 0
        for cellID ,cell in self.cells_.items():
            if bool(cell.type_):
                self.sample_center_.append(cell.center_)
                self.sample_total_volume_ += cell.volume_
        self.sample_center_ = np.mean(self.sample_center_, axis=0)
        for polygonID, polygon in self.polygons_.items():
            if polygon.is_surface_:
                self.sample_surface_area_ += polygon.area_
        return
    
    def extract_cell(self,cellID):
        vertices={}
        edges={}
        polygons={}
        cell=copy.deepcopy(self.cells_[cellID])
        cells = {cellID:cell}
        for polygonID in cell.polygons_:
            polygons[polygonID]=copy.deepcopy(self.polygons_[polygonID])
            for edgeID in polygons[polygonID].edges_:
                edges[edgeID]=copy.deepcopy(self.edges_[edgeID])
                for vertexID in edges[edgeID].vertices_:
                    vertices[vertexID]=copy.deepcopy(self.vertices_[vertexID])

        for polygonID in cell.polygons_:
            # polygon=copy.deepcopy(polygons[polygonID])
            polygon=polygons[polygonID]

            # print("polygonID", polygonID,polygon.vertices_)

            sum_of_cross_product=np.array([0,0,0])
            for i in range(len(polygon.vertices_)):
                v1=polygon.vertices_[i]
                v2=polygon.vertices_[(i+1)%len(polygon.vertices_)]
                # print("     ", v1, v2)
                # print("before",vertices[v1].position_)
                vector1=np.subtract(vertices[v1].position_, cell.center_)
                vector2=np.subtract(vertices[v2].position_, cell.center_)
                # print("after",vertices[v1].position_)
                cross_product=np.cross(vector1, vector2)
                sum_of_cross_product=np.add(sum_of_cross_product, cross_product)

            sign=np.sign(np.dot(sum_of_cross_product,np.subtract(polygon.center_, cell.center_)))
            if sign==-1:
                polygon.vertices_=polygon.vertices_[::-1]

        for vertexID,vertex in vertices.items():
            for polygonID,polygon in polygons.items():
                if vertexID in polygon.vertices_:
                    vertex.polygons_.append(polygonID)
        
        single_cell = Sample.from_topology(vertices,edges,polygons,cells)
        single_cell.config_dir_ = self.config_dir_
        single_cell.time_ = self.time_
        single_cell.s0_ = self.s0_
        single_cell.kv_ = self.kv_
        return single_cell
    
    # This function dumps the surface layer of the spheroid to a vtk file.
    def dump_surface_vtk(self):
        if not self.tissueType_ == "spheroid":
            raise ValueError("This function is only for spheroids")

        # Step 1: Extract the surface cells, polygons, edges and vertices.
        tmp_cells:dict[int,topology.Cell] = {}
        tmp_polygons:dict[int,topology.Polygon] = {}
        tmp_edges:dict[int,topology.Edge] = {}
        tmp_vertices:dict[int,topology.Vertex] = {}

        for cellID,cell in self.cells_.items():
            if cell.is_surface_:
                tmp_cells[cellID] = cell
                for polygonID in self.cells_[cellID].polygons_:
                    if polygonID in tmp_polygons: pass
                    else: 
                        tmp_polygons[polygonID] = self.polygons_[polygonID]
                        for edgeID in self.polygons_[polygonID].edges_:
                            if edgeID in tmp_edges: pass
                            else:
                                tmp_edges[edgeID]=self.edges_[edgeID]
                                
                                for vertexID in self.edges_[edgeID].vertices_:
                                    if vertexID in tmp_vertices:pass
                                    else: 
                                        tmp_vertices[vertexID] = self.vertices_[vertexID]

        # At this point, we have the surface cells, polygons, edges and vertices
        # in the dictionaries tmp_cells, tmp_polygons, tmp_edges and tmp_vertices
        # respectively.
        
        # However, the vertex IDs in tmp_vertices are not contiguous, so we will
        # remap these IDs when listing the polygon vertices in the vtk file.
        
        v_map = functions.mapmaker(tmp_vertices)
        
        # This variable is the number  of data values in the 
        # vtk file associated with polygons.
        
        # In each POLYGONS line, the first value is the number of vertices in the
        # polygon, and the subsequent values are the vertex numbers (in the order
        # they are listed in the POINTS section of the vtk file).

        totalPolygonDataPoints = 0
        for polygonID,polygon in tmp_polygons.items():
            totalPolygonDataPoints += (len(polygon.vertices_)
                                       + 1)
        # Step 2: Write the vtk file
        with open(self.config_dir_
                  + "{:07d}.surface.vtk".format(self.time_),'w') as file:
            
            file.write("# vtk DataFile Version 2.0\n")
            file.write("polydata\n")
            file.write("ASCII\n")
            file.write("DATASET POLYDATA\n")
            file.write("POINTS {} double\n".format(len(tmp_vertices)))
            
            for vertexID, vertex in tmp_vertices.items():
                file.write("{:12.5e} {:12.5e} {:12.5e}\n".format(
                    vertex.position_[0],
                    vertex.position_[1],
                    vertex.position_[2]))
            file.write("POLYGONS {} {}\n".format(len(tmp_polygons),
                                                 totalPolygonDataPoints))
        
    
            for polygonID, polygon in tmp_polygons.items():
                file.write("{:<7d}".format(len(polygon.edges_)))
                for vID in polygon.vertices_: 
                    file.write("{:<7d}".format(v_map[vID]))
                file.write("\n")

            file.write("CELL_DATA {}\n".format(len(tmp_polygons)))
            file.write("SCALARS shape_index double\n")
            file.write("LOOKUP_TABLE default\n")
            for polygonID, polygon in tmp_polygons.items():
                for cellID, cell in tmp_cells.items():
                    if polygonID in cell.polygons_:
                        file.write("{:<12.6f}\n".format(cell.shape_index_))
                        break
            # You can include more scalars here if you want.
            file.close()

        return

    def loadCrossBoundaryAttributes(self):
        for edgeID, edge in self.edges_.items():
            edge.length_ = np.linalg.norm(
                np.subtract(
                    self.vertices_[edge.vertices_[0]].position_,
                    self.vertices_[edge.vertices_[1]].position_))
            if edge.length_ > self.boxSize_/2:
                edge.crossBoundary_ = True
        for polygonID,polygon in self.polygons_.items():
            for edgeID in polygon.edges_:
                if self.edges_[edgeID].crossBoundary_:
                    polygon.crossBoundary_ = True
                    break
        for cellID,cell in self.cells_.items():
            for polygonID in cell.polygons_:
                if self.polygons_[polygonID].crossBoundary_:
                    cell.crossBoundary_ = True
                    break
        return

    def dump_vtk(self,filename):
        v_map = functions.mapmaker(self.vertices_)
        # e_map = functions.mapmaker(self.edges_)
        # p_map = functions.mapmaker(self.polygons_)
        totalPolygonDataPoints = 0
        for polygonID,polygon in self.polygons_.items():
            totalPolygonDataPoints += (len(polygon.vertices_) + 1)
        # with open(self.config_dir_
                #   + "{:07d}.cell_{}.vtk".format(self.time_,self.id_),'w') as file:
        with open(filename, "w") as file: 
            file.write("# vtk DataFile Version 2.0\n")
            file.write("polydata\n")
            file.write("ASCII\n")
            file.write("DATASET POLYDATA\n")
            file.write("POINTS {} double\n".format(len(self.vertices_)))
            
            for vertexID, vertex in self.vertices_.items():
                file.write("{:12.5e} {:12.5e} {:12.5e}\n".format(
                    vertex.position_[0],
                    vertex.position_[1],
                    vertex.position_[2]))
            file.write("POLYGONS {} {}\n".format(len(self.polygons_),
                                                 totalPolygonDataPoints))
        
    
            for polygonID, polygon in self.polygons_.items():
                file.write("{:<7d}".format(len(polygon.edges_)))
                for vID in polygon.vertices_: 
                    file.write("{:<7d}".format(v_map[vID]))
                file.write("\n")

            # file.write("CELL_DATA {}\n".format(len(tmp_polygons)))
            # file.write("SCALARS shape_index double\n")
            # file.write("LOOKUP_TABLE default\n")
            # for polygonID, polygon in tmp_polygons.items():
            #     for cellID, cell in tmp_cells.items():
            #         if polygonID in cell.polygons_:
            #             file.write("{:<12.6f}\n".format(cell.shape_index_))
            #             break
            # You can include more scalars here if you want.
            file.close()
            
    #returns[del_nu(c^s)]j. a 3x3 matrix. the rows (first index) are s, which are components of polygon center.
    # def okuda_derivative_tensor(self, polygonID, vertexID):
    #     # self.arrange_cell_polygons(cellID)
    #     tensor = np.zeros((3,3))
    #     current_vertex_index = self.polygons_[polygonID].vertices_.index(vertexID)
    #     previous_vertex_index = current_vertex_index - 1
    #     next_vertex_index = (current_vertex_index + 1) % len(self.polygons_[polygonID].vertices_)
    #     current_vertex_id = self.polygons_[polygonID].vertices_[current_vertex_index]
    #     previous_vertex_id = self.polygons_[polygonID].vertices_[previous_vertex_index]
    #     next_vertex_id = self.polygons_[polygonID].vertices_[next_vertex_index]
    #     # print(previous_vertex_id, current_vertex_id, next_vertex_id)
    #     current_vertex = self.vertices_[current_vertex_id].position_
    #     previous_vertex = self.vertices_[previous_vertex_id].position_
    #     next_vertex = self.vertices_[next_vertex_id].position_
    #     current_edge = np.subtract(next_vertex,current_vertex)
    #     previous_edge = np.subtract(current_vertex,previous_vertex)

    #     # (1/2P)(||l_v||+||l_nu-1||)delta(sj)
    #     for i in range(3):
    #         tensor[i][i] += (np.linalg.norm(current_edge) + np.linalg.norm(previous_edge))/(2 * self.polygons_[polygonID].perimeter_)

    #     # (1/2P)(l_nu-1^j/||l_nu-1||(r_nu-1+r_v)s-(l_v^j/||l_v||(r_v+r_nu+1)^s)
    #     for s in range(3):
    #         for j in range(3):
    #             term = 0
    #             term += (previous_edge[j]/np.linalg.norm(previous_edge)) * (previous_vertex[s] + current_vertex[s])
    #             term -= (current_edge[j]/np.linalg.norm(current_edge)) * (current_vertex[s] + next_vertex[s])
    #             term /= (2 * self.polygons_[polygonID].perimeter_)
    #             tensor[s][j] += term

    #     # -1/Pc^s(l_nu-1^j/||l_nu-1||-l_v^j/||l_v||)
    #     for s in range(3):
    #         for j in range(3):
    #             term = 0
    #             term -= (previous_edge[j]/np.linalg.norm(previous_edge))
    #             term += (current_edge[j]/np.linalg.norm(current_edge))
    #             term *= (self.polygons_[polygonID].center_[s] / self.polygons_[polygonID].perimeter_)
    #             tensor[s][j] += term

    #     return tensor

class SingleCell:
    def __init__(self,vertices,edges,polygons,cell):
        self.config_dir_ = None
        self.time_ = None
        self.id_ = cell.id_
        self.vertices_ = vertices
        self.edges_ = edges
        self.polygons_ = polygons
        self.cells_ = {cell.id_:cell}
        self.center_ = cell.center_
        self.s0_ = None
        self.kv_ = None
        self.gamma_ = None

    def update_center(self):
        center = []
        for _,vertex in self.vertices_.items():
            center.append(vertex.position_)
        center = np.mean(center,axis = 0)
        self.center_ = center
        self.cells_[self.id_].center_ = center
    def dump_vtk(self,filename):
        v_map = functions.mapmaker(self.vertices_)
        # e_map = functions.mapmaker(self.edges_)
        # p_map = functions.mapmaker(self.polygons_)
        totalPolygonDataPoints = 0
        for polygonID,polygon in self.polygons_.items():
            totalPolygonDataPoints += (len(polygon.vertices_) + 1)
        # with open(self.config_dir_
                #   + "{:07d}.cell_{}.vtk".format(self.time_,self.id_),'w') as file:
        with open(filename, "w") as file: 
            file.write("# vtk DataFile Version 2.0\n")
            file.write("polydata\n")
            file.write("ASCII\n")
            file.write("DATASET POLYDATA\n")
            file.write("POINTS {} double\n".format(len(self.vertices_)))
            
            for vertexID, vertex in self.vertices_.items():
                file.write("{:12.5e} {:12.5e} {:12.5e}\n".format(
                    vertex.position_[0],
                    vertex.position_[1],
                    vertex.position_[2]))
            file.write("POLYGONS {} {}\n".format(len(self.polygons_),
                                                 totalPolygonDataPoints))
        
    
            for polygonID, polygon in self.polygons_.items():
                file.write("{:<7d}".format(len(polygon.edges_)))
                for vID in polygon.vertices_: 
                    file.write("{:<7d}".format(v_map[vID]))
                file.write("\n")

            # file.write("CELL_DATA {}\n".format(len(tmp_polygons)))
            # file.write("SCALARS shape_index double\n")
            # file.write("LOOKUP_TABLE default\n")
            # for polygonID, polygon in tmp_polygons.items():
            #     for cellID, cell in tmp_cells.items():
            #         if polygonID in cell.polygons_:
            #             file.write("{:<12.6f}\n".format(cell.shape_index_))
            #             break
            # You can include more scalars here if you want.
            file.close()

