Class Hierarchy:

Parent: Sample
Children: PeriodicTissue, Spheroid

Class structure:

1. tissueSample.Sample:

Contains attributes and routines that are sensible for spheroids, periodic tissue and arbitrary collections of cells (which could even be a single cell...)

Nothing is initialized and no routines are run in init.

Attributes and their initializations:
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
        self.v0_ = None
        self.gamma_ = None
        self.kv_ = None
        self.boxSize_ = None
        self.sample_total_volume_ = None

@classmethod
from_topology: initialize a sample from dictionaries; may be useful for cell collection initializations.


Setters:
set_config_dir(self,dir):
set_file(self,filename):
    Note: filename must exist in self.config_dir_
set_origin:

Routines:
    
    load_config():
        loads config from self.file_ (format like "sample.topo)
        Currently can handle:
        (a) periodic tissue
        (b) spheroid (type 0 cells are virtual)
    
    load_conf_file():
        Currently can handle:
        (a) periodic tissue
        (b) spheroid (s0,gamma)

    load_cell_vertices():
        once the configuration is loaded (e.g with load_config()), it evaluates cell.vertices_ by working backwards (cell->polygons->edges-.vertices)

    calculate_cell_centers():
        average of positions of cell vertices.
    
    load_cell_attributes():
        First, calls load_cell_vertices()
        Then, from {time}.cellInfo.txt, loads the following cell attributes
        (a) center
        (b) volume
        (c) shape index
        (d) surface area (calcuate from vol,shape index)
    
    load_cross_boundary_attributes():
        Once 
        (a) the topology is loaded (e.g with load_config())
        (b) boxSize_ is loaded (e.g with load_conf_file())
        
        This function determines
            edge.crossBoundary_
            polygon.crossBoundary_
            cell.crossBoundary_
        by evaluating all edge lengths and comparing to boxsize/2
    calculate_boundary_cell_attributes():
        *** USING THE COM POLYGON CENTER ***
        First it finds relevant boundary cells (for spheroid, it ignores cell type 0 or empty cells)
        Each boundary cell is extracted and vertices are transformed so the regular algorithms can be used.
        Then the cell's
        (a) Volume
        (b) Surface Area
        (c) Shape Index
        are calculated; and stored back in the main sample.

    calculate_polygon_centers_and_perimeters()
        These quantities are calculated as prescribed in Okuda
    
    calculate_COM_polygon_centers():
        COM polygon centers

    calculate_polygon_areas:
        calculate polygon areas, provided polygon vertices exist. Note: may be doing extra work for spheroid, calculating empty cell polygon areas

    calculate_cell_surface_areas():
        Once all polygon areas are available (e.g with calculate_polygon_areas), adds up the total surface area

    calculate_cell_shape_indices():
        once cell.volume_ and cell.surface_area_ are available...

    calculate_cell_volumes():
        needs cell_center,and polygon vertices in clockwise/anticlockwise order (although in theory there is a way to calculate without that)

    arrange_polygon_vertices():
        uses the static method resoluve_polygon_edge_connectivity(); arranges for non cross boundary (and for spheroid, type 1) cells
    
    evaluate_cell_neighbors():
        Evaluates the dictionary self.cell_neighbors_

    extract_cell():
        Returns a sample with only one cell. Useful in calculating stuff for cross boundary cells (and much more...)
    
    write_cell_collection_vtk(cells_array,filename,use_scalar = False):
        *** cannot handle crossBoundary cells yet ***
        write vtk to self.config_dir_/filename
        If use_scalar, uses polygon.vtk_scalar_




