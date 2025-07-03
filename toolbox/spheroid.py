from toolbox.tissue import Sample
from toolbox import functions
from toolbox import topology
import numpy as np
import os

class Spheroid(Sample):
    def __init__(self):
        super().__init__()
        self.tissueType_ = "spheroid"
        self.spheroid_center_ = None
        self.spheroid_surface_area_ = None
        self.spheroid_volume_ = None
        
    @classmethod
    def from_config(cls, config_dir, simulation_time):
        sample = cls()
        sample.time_ = simulation_time
        sample.config_dir_ = config_dir
        sample.file_ = "{}{:07}.topo.txt".format(sample.config_dir_,sample.time_)
        # Validate the input:
        if not config_dir.endswith("/"):
            print("config_dir must end with a '/'. Attempting to fix this...")
            config_dir = config_dir + "/"
        if not os.path.isdir(config_dir):
            raise ValueError("sample.config_dir_ = {} must be a valid directory".format(config_dir))
        sample.load_config_from_topo()
        sample.load_conf_file()
        sample.load_cell_attributes()
        sample.calculate_polygon_centers_and_perimeters()
        sample.arrange_polygon_vertices()
        sample.calculate_polygon_areas()
        sample.identify_surface_polygons_and_cells()
        sample.load_spheroid_attributes()
        return sample

    # load_config_from_topo(): given self.time_, this function first checks if
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
    def identify_surface_polygons_and_cells(self):
        for polygonID, polygon in self.polygons_.items():
            cell_types = []
            for _, cell in self.cells_.items():
                if polygonID in cell.polygons_:
                    cell_types.append(cell.type_)
            if set(cell_types) == {0,1}:
                polygon.is_surface_ = True
        for _, cell in self.cells_.items():
            if not cell.type_:
                continue
            if cell.is_surface_:
                continue
            for polygonID in cell.polygons_:
                if self.polygons_[polygonID].is_surface_:
                    cell.is_surface_ = True
                    break
    def load_spheroid_attributes(self):
        self.spheroid_center_ = []
        self.spheroid_surface_area_ = 0
        self.spheroid_volume_ = 0
        for _,cell in self.cells_.items():
            if not cell.type_:
                continue
            self.spheroid_center_.append(cell.center_)
            self.spheroid_volume_ += cell.volume_
        self.spheroid_center_ = np.mean(self.spheroid_center_, axis=0)
        for _, polygon in self.polygons_.items():
            if polygon.is_surface_:
                self.spheroid_surface_area_ += polygon.area_

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
            if not cell.is_surface_:
                continue
            tmp_cells[cellID] = cell
            for polygonID in self.cells_[cellID].polygons_:
                if polygonID in tmp_polygons: 
                    continue

            tmp_polygons[polygonID] = self.polygons_[polygonID]
            for edgeID in self.polygons_[polygonID].edges_:
                if edgeID in tmp_edges:
                    continue
                tmp_edges[edgeID]=self.edges_[edgeID]
                    
                for vertexID in self.edges_[edgeID].vertices_:
                    if vertexID in tmp_vertices:
                        continue
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