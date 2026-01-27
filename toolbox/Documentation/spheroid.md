Class Hierarchy:

Parent: Sample

Class structure:

1. spheroid.Spheroid:
    Contains spheroid specific add ons to the parent class. Initialize:
        self.tissueType_ = "spheroid"
        self.spheroid_center_ = None
        self.spheroid_surface_area_ = None
        self.spheroid_volume_ = None

@classmethod
from_dir_and_time(config_dir, simulation_time):
    ...

@classmethod
from_config(config_dir,input_filename):
    sample.time_ = 0
    sample.set_config_dir(config_dir)
    sample.set_file(input_filename)
    sample.load_spheroid_from_file(input_filename)

Routines:
    
load_spheroid_from_file():
    (a) loads the spheroid, 
    (b) calculates cell attributes for spheroid cells, namely
        (i) cell.vertices_
        (ii) cell.center_
        (iii) cell.volume_
        (iv) cell.surface_

load_config_from_topo():
    ...

identify_surface_polygons_and_cells():
    surface polygons are shared between cells of types {0,1}. 
    surface cells have surface polygons

calculate_spheroid_attributes():
    calculates
        self.spheroid_center_
        self.spheroid_surface_area_
        self.spheroid_volume_

write_surface_vtk(filename)    
write_spheroid_vtk(filename)