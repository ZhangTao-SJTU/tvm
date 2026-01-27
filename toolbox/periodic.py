from toolbox.tissue import Sample
from toolbox import functions
import numpy as np
import os

class PeriodicTissue(Sample):
    def __init__(self):
        super().__init__()
        self.tissueType_ = "periodic"
        self.periodic_sample_center_ = None

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
        sample.set_config_dir(config_dir)
        sample.set_file(input_filename)
        sample.load_periodic_tissue_from_file()
        return sample
    
    # Packaging the loading of the periodic tissue from file
    def load_periodic_tissue_from_file(self):
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
    
    def calculate_periodic_sample_center(self):
        self.periodic_sample_center_ = np.mean([vertex.position_ for _,vertex in self.vertices_.items()], axis = 0)

    def write_periodic_vtk(self,filename, use_scalar = False):
        cells_array = [cellID for cellID,cell in self.cells_.items() if not cell.crossBoundary_]
        self.write_cell_collection_vtk(cells_array, filename, use_scalar)
        