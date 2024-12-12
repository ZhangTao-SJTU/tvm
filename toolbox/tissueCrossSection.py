import numpy as np
import pyvtk
from toolbox import topology,tissueSample, functions

def makeSampleCrossSection(sample:tissueSample.Sample,
                           normal:np.ndarray = np.array([1,0,0]),
                           scalar:str = "shape_index"):


    #####################################################################
    # Step 1: given the normal, the cross section plane passes through 
    # the center of the spheroid. We first identify the cells that 
    # are intersected by the cross section plane. 
    
    # We will call these cells the "mother cells": (cell.is_mother_ = True)

    # We accomplish this by checking the signs of the dot products
    # of the vectors from the spheroid center to the vertices of the cell
    # with the normal vector.

    # If any of the cell vertices lie on the cross section plane,
    # we raise an error.
    #####################################################################

    for cellID,cell in sample.cells_.items():
        if sample.tissueType_ == "spheroid" and not cell.type_:
            continue
        if sample.tissueType_ == "periodic" and cell.crossBoundary_:
            continue
        signs = []
        for vertexID in cell.vertices_:
            vector_from_spheroid_center = np.subtract(sample.vertices_[vertexID].position_,
                                                        sample.sample_center_)
            sign = np.sign(np.dot(normal,vector_from_spheroid_center))
            signs.append(sign)
            if int(sign) == 0:
                raise ValueError("Error: a vertex of cell {} on cross section".format(cellID))
        if len(set(signs)) == 1: pass
        else: cell.is_mother_ = True

    #####################################################################
    # Step 2: Given mother cells, we identify the polygons that are intersected 
    # by the cross section plane. We will call these polygons the mother polygons:
    # polygon.is_mother_ = True
    # and append the polygonID to the list mother_polygons
    #####################################################################
            
    mother_polygons = []
    for cellID,cell in sample.cells_.items():
        if cell.is_mother_:
            for polygonID in cell.polygons_:
                signs = []
                for vertexID in sample.polygons_[polygonID].vertices_:
                    vector_from_spheroid_center = np.subtract(sample.vertices_[vertexID].position_,
                                                              sample.sample_center_)
                    sign = np.sign(np.dot(normal,vector_from_spheroid_center))
                    signs.append(sign)
                    # redundant check: this was already raised in step 1
                    # if int(sign)==0:
                    #     print("Error: vertex on cross section")

                if len(set(signs)) == 1: pass
                else: 
                    sample.polygons_[polygonID].is_mother_=True
                    mother_polygons.append(polygonID)
    mother_polygons = np.unique(mother_polygons)

    #####################################################################
    # Step 3: We identify mother vertices as vertices that belong to mother polygons
    # this is useful because we will break up the mother polygons 
    # into triangular polygon sections. That will involve creating new vertices at the
    # center of the mother polygons.
    #####################################################################

    mother_vertices = []
    for polygonID,polygon in sample.polygons_.items():
        if polygon.is_mother_:
            for vertexID in polygon.vertices_:
                mother_vertices.append(vertexID)
    mother_vertices = np.unique(mother_vertices)
    
    #####################################################################
    # Step 4: Given the mother vertices and mother polygons
    # we create new exclusive dictionaries of these
    # These dictionaries are called new_vertices and new_polygons respectively
    #####################################################################

    v_map = functions.mapmaker(mother_vertices)
    # print("v_map: ",v_map)
    new_vertices = {}
    for vertexID in mother_vertices:
        new_vertices[v_map[vertexID]] = topology.Vertex(v_map[vertexID])
        new_vertices[v_map[vertexID]].setPosition(sample.vertices_[vertexID].position_)
        new_vertices[v_map[vertexID]].og_id_ = vertexID

    p_map = functions.mapmaker(mother_polygons)
    new_polygons = {}
    for polygonID in mother_polygons:
        new_polygons[p_map[polygonID]] = topology.Polygon(p_map[polygonID])
        new_polygons[p_map[polygonID]].og_id_ = polygonID
        new_polygons[p_map[polygonID]].center_ = sample.polygons_[polygonID].center_
        for vertexID in sample.polygons_[polygonID].vertices_:
            new_polygons[p_map[polygonID]].addVertex(v_map[vertexID])

    #####################################################################
    # Step 5: We break up the mother polygons into triangular polygons.
    # This involves creating new vertices at the polygon center, and we will append those to the 
    # new_vertices dictionary. We will also create new triangular polygons and append those to the
    # triangular_polygons dictionary.
    #####################################################################
            
    triangular_polygons = {}

    for polygonID,polygon in new_polygons.items():
        new_vertexID = len(new_vertices)
        new_vertex = topology.Vertex(new_vertexID)
        new_vertex.setPosition(polygon.center_)
        new_vertices[new_vertexID] = new_vertex
        for i,vertexID in enumerate(polygon.vertices_):
            new_triangular_polygonID = len(triangular_polygons)
            new_triangular_polygon = topology.Polygon(new_triangular_polygonID)
            new_triangular_polygon.og_id_ = polygon.og_id_
            next_vertexID = polygon.vertices_[(i+1)%len(polygon.vertices_)]
            new_triangular_polygon.addVertex(vertexID)
            new_triangular_polygon.addVertex(new_vertexID)
            new_triangular_polygon.addVertex(next_vertexID)
            signs = []
            for vertexID in new_triangular_polygon.vertices_:
                vector_from_spheroid_center = np.subtract(new_vertices[vertexID].position_,
                                                          sample.sample_center_)
                sign = np.sign(np.dot(normal,vector_from_spheroid_center))
                signs.append(sign)
                if int(sign) == 0:
                    raise ValueError(
                        "Error: the polygon center for polygon {} on cross section".format(polygonID))
                if len(set(signs)) == 1: pass
                else: 
                    triangular_polygons[new_triangular_polygonID] = new_triangular_polygon
    
    #####################################################################
    # Step 6: Create new edges for the triangular polygons
    #####################################################################
    edge_pairs = []

    for polygonID,polygon in triangular_polygons.items():
        edge_pairs.append([polygon.vertices_[0],polygon.vertices_[1]])
        edge_pairs.append([polygon.vertices_[1],polygon.vertices_[2]])
        edge_pairs.append([polygon.vertices_[2],polygon.vertices_[0]])

    unique_pairs_set = {tuple(sorted(pair)) for pair in edge_pairs}

    # Convert the set back to a list of lists

    edge_pairs = [list(pair) for pair in unique_pairs_set]
        
    new_edges = {}
    for pair in edge_pairs:
        new_edgeID = len(new_edges)
        new_edge = topology.Edge(new_edgeID)
        new_edge.addVertex(pair[0])
        new_edge.addVertex(pair[1])
        new_edges[new_edgeID] = new_edge

    for polygonID,polygon in triangular_polygons.items():
        pair1 = [polygon.vertices_[0],polygon.vertices_[1]]
        pair2 = [polygon.vertices_[1],polygon.vertices_[2]]
        pair3 = [polygon.vertices_[2],polygon.vertices_[0]]
        counter = 0
        for edgeID,edge in new_edges.items():
            if {tuple(sorted(edge.vertices_))} == {tuple(sorted(pair1))}: 
                polygon.addEdge(edgeID)
                counter += 1
            if {tuple(sorted(edge.vertices_))} == {tuple(sorted(pair2))}: 
                polygon.addEdge(edgeID)
                counter += 1
            if {tuple(sorted(edge.vertices_))} == {tuple(sorted(pair3))}: 
                polygon.addEdge(edgeID)
                counter += 1
            if counter == 3: break

    #####################################################################
    # Step 7: Create intersection vertices, edges and polygons:
            
    # Intersection vertices: vertices that lie on the cross section plane
    # We find the intersection of the edges of the triangular polygons
    # with the cross section plane. We create new vertices at these intersection points.
        
    # Intersection edges: edges that connect the intersection vertices
    # we connect the intersection vertices on the same triangular polygon 
    # by an intersection edge
            
    # Intersection polygons: polygons that are formed by the intersection edges
    # we connect the intersection edges on mother polygons of the same mother cell 
    # to form an intersection polygon. 
    #####################################################################
    intersection_vertices = {}
    for edgeID, edge in new_edges.items():
        r_0 = new_vertices[edge.vertices_[0]].position_
        r_1 = new_vertices[edge.vertices_[1]].position_
        dividing_factor = np.dot(normal,r_1) - np.dot(normal,r_0)
        if dividing_factor == 0:
            raise ValueError("Error: edge {} parallel to cross section plane".format(edgeID))
                                                
        t = ((np.dot(normal,sample.sample_center_) - np.dot(normal,r_0))
             / (np.dot(normal,r_1) - np.dot(normal,r_0)))
        if t <=1 and t>=0:
            id = len(intersection_vertices)
            position = np.add(r_0,np.multiply(t,np.subtract(r_1,r_0)))
            edge.intersection_vertex_ = id
            vertex_object = topology.Vertex(id)
            vertex_object.setPosition(position)
            intersection_vertices[id] = vertex_object

    intersection_edges = {}
    for polygonID,polygon in triangular_polygons.items():
        
        id = len(intersection_edges)
        polygon.intersection_edge_=id
        intersection_edges[id]=topology.Edge(id)

        for edgeID in polygon.edges_:
            if new_edges[edgeID].intersection_vertex_ != None:
                intersection_edges[id].addVertex(new_edges[edgeID].intersection_vertex_)

    new_cells = {}
    for cellID, cell in sample.cells_.items():
        if cell.is_mother_:
            new_cellID = len(new_cells)
            new_cell = topology.Cell(new_cellID)
            if scalar == "shape_index":
                new_cell.shape_index_change_ = cell.shape_index_change_
            elif scalar == "principal_radial_stress":
                new_cell.principal_radial_stress_ = cell.principal_radial_stress_
            new_cell.polygons_ = []
            for polygonID in cell.polygons_:
                for daughter_polygonID, polygon in triangular_polygons.items():
                    if polygon.og_id_ == polygonID:
                        new_cell.addPolygon(daughter_polygonID)
            new_cells[new_cellID] = new_cell

    intersection_polygons = {}

    intersection_polygon_scalars = {}

    for cellID, cell in new_cells.items():
        id = len(intersection_polygons)
        intersection_polygons[id] = topology.Polygon(id)
        if scalar == "shape_index":
            intersection_polygon_scalars[id]=cell.shape_index_change_
        elif scalar == "principal_radial_stress":
            intersection_polygon_scalars[id] = np.linalg.norm(cell.principal_radial_stress_)
        #intersection_polygon_scalars[id]=np.linalg.norm(cell.principal_radial_stress_)
        # intersection_polygon_scalars[id]=cell.shape_index_
    
        for polygonID in cell.polygons_:
            intersection_polygons[id].addEdge(triangular_polygons[polygonID].intersection_edge_)

    for polygonID,polygon in intersection_polygons.items():
        polygon.edges_ = np.unique(polygon.edges_)
        vertices = []
        for edgeID in polygon.edges_:
            vertices.append(intersection_edges[edgeID].vertices_)
        
        polygon.vertices_ = functions.arrange_polygon(vertices)

    Points_ = []
    Polygons_ = []
    cellscalars = []
    for vertexID, vertex in new_vertices.items():
        Points_.append(vertex.position_)
    for polygonID,polygon in triangular_polygons.items():
        Polygons_.append([vertexID for vertexID in polygon.vertices_])
        cellscalars.append(0)

    structure = pyvtk.PolyData(points=Points_,polygons=Polygons_)
    celldata = pyvtk.CellData(\
        pyvtk.Scalars(cellscalars,
                name = 'cell_scalars'))
    vtk = pyvtk.VtkData(structure,celldata)
    vtk.tofile(sample.config_dir_ + "{:07d}.triangles.vtk".format(sample.time_),'ascii')

    Points_=[]
    Polygons_=[]
    cellscalars=[]
    for vertexID, vertex in intersection_vertices.items():
        Points_.append(vertex.position_)
    for polygonID,polygon in intersection_polygons.items():
        Polygons_.append([vertexID for vertexID in polygon.vertices_])
        cellscalars.append(intersection_polygon_scalars[polygonID])

    structure = pyvtk.PolyData(points = Points_,polygons=Polygons_)
    celldata = pyvtk.CellData(\
        pyvtk.Scalars(cellscalars,
                name='cell_scalars'))
    vtk = pyvtk.VtkData(structure,celldata)
    vtk.tofile(sample.config_dir_+"{:07d}.crossSection.vtk".format(sample.time_),'ascii')
    
    return
