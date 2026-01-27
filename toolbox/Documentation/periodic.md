Class Hierarchy:

Parent: Sample

Class structure:

1. periodic.PeriodicTissue:

    Contains periodic tissue specific add ons to the parent class. Initialize:
        self.tissueType_ = "periodic"
        self.periodic_sample_center_ = None


@classmethod
from_config(config_dir,input_filename):
    sample.time_ = 0
    sample.set_config_dir(config_dir)
    sample.set_file(input_filename)
    sample.load_periodic_tissue_from_file()

Routines:
    
load_periodic_tissue_from_file():
    (a) loads the periodic tissue, 
    (b) calculates cell attributes for ALL vertices: namely
        (i) cell.vertices_
        (ii) cell.center_
        (iii) cell.volume_
        (iv) cell.surface_

calculate_periodic_sample_center():
    average of positions of vertices

write_periodic_vtk():
    shortcut use of write_cell_collection_vtk from parent class.