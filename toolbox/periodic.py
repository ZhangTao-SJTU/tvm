from toolbox.tissue import Sample
from toolbox import functions
import numpy as np
import os

class PeriodicTissue(Sample):
    def __init__(self):
        super().__init__()
        self.tissueType_ = "periodic"

    @classmethod
    def from_config(cls, config_dir, input_filename = "sample.topo"):
        # Validate the input: 
        if not config_dir.endswith("/"):
            print("config_dir must end with a '/'. Attempting to fix this...")
            config_dir = config_dir + "/"
        if not os.path.isdir(config_dir):
            raise ValueError("config_dir must be a valid directory")
        sample = cls()
        sample.time_ = 0
        sample.config_dir_ = config_dir
        sample.load_periodic_tissue_from_file(input_filename)
        return sample
    
    # Packaging the loading of the periodic tissue from file
    def load_periodic_tissue_from_file(self,filename):
        self.set_file(filename)
        self.load_config()
        self.load_conf_file()
        self.load_cell_vertices()
        self.load_cross_boundary_attributes()
        self.arrange_polygon_vertices()
        self.calculate_cell_centers()
        self.calculate_periodic_sample_center()
        self.calculate_COM_polygon_centers()
        self.calculate_cell_volumes()
        self.calculate_polygon_areas()
        self.calculate_cell_surface_areas()
        self.calculate_cell_shape_indices()
        self.calculate_boundary_cell_attributes()
    def load_cross_boundary_attributes(self):
        for edgeID, edge in self.edges_.items():
            edge.length_ = np.linalg.norm(
                np.subtract(
                    self.vertices_[edge.vertices_[0]].position_,
                    self.vertices_[edge.vertices_[1]].position_))
            if edge.length_ > self.boxSize_/2:
                edge.crossBoundary_ = True
                
        for _,polygon in self.polygons_.items():
            for edgeID in polygon.edges_:
                if self.edges_[edgeID].crossBoundary_:
                    polygon.crossBoundary_ = True
                    break
        for _,cell in self.cells_.items():
            for polygonID in cell.polygons_:
                if self.polygons_[polygonID].crossBoundary_:
                    cell.crossBoundary_ = True
                    break

    def calculate_boundary_cell_attributes(self):
        boundary_cells = []
        for cellID,cell in self.cells_.items():
            if cell.crossBoundary_:
                boundary_cells.append(cellID)

        for num,testCellID in enumerate(boundary_cells):
            test_b_cell = self.extract_cell(testCellID)
            for cellID,cell in test_b_cell.cells_.items():
                cell.crossBoundary_ = False
            for polygonID, polygon in test_b_cell.polygons_.items():
                polygon.crossBoundary_ = False
            for edgeID, edge in test_b_cell.edges_.items():
                edge.crossBoundary_ = False

            axes_to_flip = []  # Flip all axes
            for edgeID,edge in test_b_cell.edges_.items():
                v0 = test_b_cell.vertices_[edge.vertices_[0]].position_
                v1 = test_b_cell.vertices_[edge.vertices_[1]].position_
                edge_vector = np.subtract(v1, v0)
                for i in range(3):
                    if abs(edge_vector[i]) > self.boxSize_/2 :
                        axes_to_flip.append(i)
            axes_to_flip = list(set(axes_to_flip))  # Remove duplicates
            # print("Axes to flip:", axes_to_flip)
            for vertexID, vertex in test_b_cell.vertices_.items():
                for i in axes_to_flip:
                    if vertex.position_[i] < self.boxSize_ / 2:
                        continue
                    vertex.position_[i] -= self.boxSize_

            # for vertexID, vertex in test_b_cell.vertices_.items():
            #     for i,coordinate in enumerate(vertex.position_):
            #         vertex.position_[i] -= boxSize*np.floor(coordinate / boxSize)
            test_b_cell.arrange_polygon_vertices()
            test_b_cell.calculate_COM_polygon_centers()
            test_b_cell.calculate_cell_centers()
            test_b_cell.calculate_cell_volumes()
            test_b_cell.calculate_polygon_areas()
            test_b_cell.calculate_cell_surface_areas()
            test_b_cell.calculate_cell_shape_indices()
            # test_b_cell.write_cell_collection_vtk([testCellID],"{}.cell.vtk".format(num))
            cell = test_b_cell.cells_[testCellID]
            self.cells_[testCellID].volume_ = cell.volume_
            self.cells_[testCellID].surface_area_ = cell.surface_area_
            self.cells_[testCellID].shape_index_ = cell.shape_index_
            # print(cell.id_,cell.volume_,cell.surface_area_,cell.shape_index_)                   
    def write_periodic_vtk(self,filename, use_scalar = False):
        vertices = []
        polygons = []
        for _,cell in self.cells_.items():
            if cell.crossBoundary_:
                continue
            for polygonID in cell.polygons_:
                polygons.append(polygonID)
        polygons = np.unique(polygons)
        total_polygons = len(polygons)
        total_polygon_data_points = 0
        for polygonID in polygons:
            polygon = self.polygons_[polygonID]
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
            for polygonID in polygons:
                polygon = self.polygons_[polygonID]
                f.write("{:12.6f}\n".format(polygon.vtk_scalar_))

    def calculate_periodic_sample_center(self):
        center = []
        for cellID,cell in self.cells_.items():
            if cell.crossBoundary_:
                continue
            center.append(cell.center_)
        self.sample_center_=np.mean(center, axis = 0)