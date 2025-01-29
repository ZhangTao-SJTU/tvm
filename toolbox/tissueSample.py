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

    @classmethod
    def from_topology(cls,vertices,edges,polygons,cells):
        sample = cls()
        sample.vertices_ = vertices
        sample.edges_ = edges
        sample.polygons_ = polygons
        sample.cells_ = cells
        sample.s0_ = None
        sample.kv_ = None
        return sample
    
    @classmethod
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
        sample.file_ = config_dir + input_filename
        sample.load_config()
        sample.load_conf_file()
        sample.load_cross_boundary_attributes()
        sample.load_cell_vertices()
        sample.arrange_polygon_vertices()
        sample.calculate_cell_centers()
        sample.calculate_periodic_sample_center()
        return sample
    
    # loadconfig(): given self.time_, this function first checks if
    # {time}.topo.txt exists in self.config_dir_. If not, it creates this file -
    # that is, it mines the topology at this time from topo.txt
    
    # Similarly, it checks and ensures the existence of {time}.cellInfo.txt
    # in self.config_dir_.
    
    # Then it loads this data into the dictionaries
    # self.vertices_,self.edges_,self.polygons_,self.cells_,self.cellIDs_
    def load_config_from_topo(self):
        self.file_ = "{}{:07}.topo.txt".format(self.config_dir_,self.time_)
        if not os.path.isfile(self.file_):
            functions.make_time_topo(self.config_dir_, self.time_)
        self.load_config()
    
    def load_config(self):
        if self.file_ is None:
            print("self.file_ not set. If you are loading info from topo.txt,")
            print("use self.load_config_from_topo()")
            raise ValueError("self.file_ not set.")
        # Load the topology. When loading cell polygons, note that there are
        # virtual cells (type=0) and real cells (type=1).
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
            if cell.type_ or self.tissueType_ == "periodic":
                cell.vertices_ = []
                for polygonID in cell.polygons_:
                    for edgeID in self.polygons_[polygonID].edges_:
                        for vertexID in self.edges_[edgeID].vertices_:
                            cell.vertices_.append(vertexID)
                cell.vertices_ = np.unique(cell.vertices_)

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
    
    # Here I use information originally from 
    # cellVolume.txt and cellShapeIndex.txt: 
    # A=s*V^2/3
    # instead of calculating areas from triangular polygon patches
    def calculate_cell_surface_areas(self):
        for cellID, cell in self.cells_.items():
            if bool(cell.type_):
                cell.surface_area_ = (
                    cell.shape_index_
                    * pow(cell.volume_,2/3))

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
    
    def calculate_cell_shape_indices(self):
        for cellID, cell in self.cells_.items():
            if self.tissueType_ == "spheroid" and not cell.type_:
                continue
            if self.tissueType_ == "periodic" and cell.crossBoundary_:
                continue
            cell.shape_index_ = cell.surface_area_ / pow(cell.volume_,2/3)
    
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
    
    def calculate_periodic_sample_center(self):
        center = []
        for cellID,cell in self.cells_.items():
            if cell.crossBoundary_:
                continue
            center.append(cell.center_)
        self.sample_center_=np.mean(center, axis = 0)

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
                self.polygons_[polygonID].vertices_ = functions.arrange_polygon(tmp_vertices)
    
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
    
    # Packaging the call of some functions from tissueSample.Sample
    # This loads cell.volume_, cell.shape_index_ etc
    def load_periodic_tissue_cell_properties(self):
        self.calculate_COM_polygon_centers()
        self.calculate_cell_volumes()
        self.calculate_polygon_areas()
        self.calculate_cell_surface_areas()
        self.calculate_cell_shape_indices()

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

    def load_cross_boundary_attributes(self):
        if not self.tissueType_ == "periodic":
            print("This function, load_cross_boundary_attributes()")
            print("is only supposed to be for periodic tissue")
            return
        
        for edgeID, edge in self.edges_.items():
            edge.length_ = np.linalg.norm(
                np.subtract(
                    self.vertices_[edge.vertices_[0]].position_,
                    self.vertices_[edge.vertices_[1]].position_))
            if edge.length_ > self.boxSize_/2:
                edge.crossBoundary_ = True
            # for vertexID in edge.vertices_:
            #     for coordinate in self.vertices_[vertexID].position_:
            #         if coordinate < 0 or coordinate > self.boxSize_:
            #             edge.crossBoundary_ = True
                
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

    def write_periodic_vtk(self,filename, use_scalar = False):
        if not self.tissueType_ == "periodic":
            print("this function is for periodic tissue only")
            return
        vertices = []
        polygons = []
        for cellID,cell in self.cells_.items():
            if cell.crossBoundary_:
                continue
            for polygonID in cell.polygons_:
                polygons.append(polygonID)
        polygons = np.unique(polygons)

        total_polygons = len(polygons)
        total_polygon_data_points = 0
        # for polygonID,polygon in self.polygons_.items():
        for polygonID in polygons:
            polygon = self.polygons_[polygonID]
            # if polygon.crossBoundary_:
            #     continue
            # total_polygons += 1
            total_polygon_data_points += len(polygon.vertices_) + 1
            for vertex in polygon.vertices_:
                vertices.append(vertex)
        vertices = np.unique(vertices)
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
            # for polygonID,polygon in self.polygons_.items():
            #     if polygon.crossBoundary_:
            #         continue
            for polygonID in polygons:
                polygon = self.polygons_[polygonID]
                f.write("{} ".format(len(polygon.vertices_)))
                for vertexID in polygon.vertices_:
                    f.write("{} ".format(v_map[vertexID]))
                f.write("\n")
            
            if not use_scalar:
                return
            
            f.write("CELL_DATA {}\n".format(total_polygons))
            f.write("SCALARS scalar_1 double\n")
            f.write("LOOKUP_TABLE default\n")
            # for polygonID,polygon in self.polygons_.items():
            #     if polygon.crossBoundary_:
            #         continue
            for polygonID in polygons:
                polygon = self.polygons_[polygonID]
                f.write("{:12.6f}\n".format(polygon.vtk_scalar_))

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
    
    def dump_vtk(self,filename):
        v_map = functions.mapmaker(self.vertices_)
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
            file.write("POLYGONS {} {}\n".format(
                len(self.polygons_),totalPolygonDataPoints))
        
            for polygonID, polygon in self.polygons_.items():
                file.write("{:<7d}".format(len(polygon.edges_)))
                for vID in polygon.vertices_: 
                    file.write("{:<7d}".format(v_map[vID]))
                file.write("\n")
            file.close()
