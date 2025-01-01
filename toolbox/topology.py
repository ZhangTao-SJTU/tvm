class Vertex:
    def __init__(self, id):
        self.id_ = id
        self.og_id_ = None
        self.position_ = [0., 0., 0.]
        self.cells_:list[int] = []
        self.polygons_=[]
        self.volume_derivative_ = None
        self.surface_area_derivative_ = None
        self.boundary_area_derivative_ = None
        self.volume_force_ = None
        self.surface_force_ = None
        self.is_driven_ = False
        self.force_ = None
        self.is_surface_ = False
        self.is_daughter_ = False
        return
    def setPosition(self, position):
        self.position_ = list(position)
    def addCell(self, cell):
        if cell not in self.cells_:
            self.cells_.append(cell)

class Edge:
    def __init__(self, id):
        self.id_ = id
        self.og_id_ = None
        self.dumpOn_ = True
        self.vertices_ = []
        self.length_ = None
        self.crossBoundary_ = False
        self.is_driven_ = False
        self.is_surface_ = False
        self.is_mother_ = False
        self.is_daughter_ = False
        self.is_in_dividing_polygon_ = False
        self.mother_id_ = None
        self.mother_polygon_id_ = None
        self.divide_flag_ = False
        self.intersection_vertex_ = None
        return
    def addVertex(self, vertex):
        self.vertices_.append(vertex)

class Polygon:
    def __init__(self, id):
        self.id_ = id
        self.og_id_ = None
        self.dumpOn_ = True
        self.edges_ = []
        self.vertices_ = []
        self.type_ = 0
        self.crossBoundary_ = False
        self.is_surface_ = False
        self.normal_ = None
        self.center_ = None
        self.perimeter_ = None
        self.area_ = None
        self.is_mother_ = False
        self.is_daughter_ = False
        self.is_dividing_polygon_ = False
        self.intersection_edge_ = None
        self.mother_id_ = None
        self.divide_flag_ = False
        self.daughter_vertices_ = None
        self.daughter_edges_ = None
        self.scalar_ = None
        self.is_fixed_ = False
        return
    def addEdge(self, edge):
        self.edges_.append(edge)
    def addVertex(self, vertex):
        self.vertices_.append(vertex)
    
class Cell:
    def __init__(self, id):
        self.id_ = id
        self.og_id_ = None
        self.polygons_ = []
        self.vertices_ = None
        self.crossBoundary_ = False
        self.is_surface_ = False
        self.surface_area_ = None
        self.s0_ = None
        self.boundary_area_ = None
        self.shape_index_ = None
        self.shape_index_change_ = None
        self.volume_ = None
        self.v0_ = None
        self.center_ = None
        self.inertia_tensor_ = None
        self.stress_tensor_ = None
        self.type_:int = 0
        self.is_mother_ = False
        self.is_daughter_ = False
        self.is_in_chain_ = False
        self.mother_id_ = None
        self.is_fixed_ = False
        return
    def addPolygon(self, polygon):
        self.polygons_.append(polygon)
    def deletePolygon(self, polygon):
        self.polygons_.remove(polygon)

# An edge object corresponding to a single edge in the fiber network. To be initialized with 
# 1. a list of two node ids,
# 2. a strain, 
# 3. a dictionary of coordinates for the entire network,
# 4. Origin. Consider making this the COM of the spheroid - which may be different from [0,0,0]. 
#       This will be the origin used to calculate spherical coordinates and spherical unit vectors at self.node_0_. We will not be shifting the cartesian coordinates of the nodes.

# Usage: mk_edges_dict() stores edge information in the form of these object. Then the fiber_network class is equipped with the dictionary from mk_edges_dict().

class fiberEdge:
    def __init__(self,nodes,strain):
        self.nodes_ = nodes
        self.tension_ = None
        self.coordinates_ = None
        self.length_ = None
        self.strain_ = strain
        self.radialDistance_ = None
        self.rHat_ = None
        self.thetaHat_ = None
        self.phiHat_ = None
        self.etaSpherical_ = None
        self.etaCartesian_ = None
