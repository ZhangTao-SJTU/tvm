Class Hierarchy:

Parent: Sample
Children: PeriodicTissue, Spheroid

Class structure:

1. tissueSample.Sample:

Contains attributes and routines that are sensible for spheroids, periodic tissue and arbitrary collections of cells (which could even be a single cell...)

Nothing is initialized and no routines are run in init.

@classmethods:
from_topology: initialize a sample from dictionaries; may be useful for cell collection initializations.

Attributes:
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

Setters:
def set_config_dir(self,dir):
def set_file(self,dir):
Routines:
    load_config:
    load_conf_file:
    load_cell_vertices:
    calculate_cell_centers
    load_cell_attributes
    calculate_polygon_centers_and_perimeters
    calculate_COM_polygon_centers
    calculate_polygon_areas
    calculate_cell_surface_areas
    calculate_cell_shape_indices
    calculate_cell_volumes
    arrange_polygon_vertices
    evaluate_cell_neighbors
    extract_cell
    write_cell_collection_vtk

1. periodic.PeriodicTissue(tissueSample.Sample):

Attributes:
    self.tissueType = "periodic"
@classmethod
from_config(config_dir, input_filename = "sample.topo"):

Routines:
load_periodic_tissue_from_file
load_cross_boundary_attributes
write_periodic_vtk
calculate_periodic_sample_center


