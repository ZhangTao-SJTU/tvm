import random
import numpy as np
from toolbox import tissueSample
from toolbox import functions
from toolbox import momentOfInertia
from toolbox import topology

def dumpVtk(sample:tissueSample.Sample):
    vertexMap = functions.mapmaker(sample.vertices_)
    with open(sample.config_dir_
                + "{:07d}.modifiedSample.vtk".format(sample.time_),'w') as file:
        
        file.write("# vtk DataFile Version 2.0\n")
        file.write("polydata\n")
        file.write("ASCII\n")
        file.write("DATASET POLYDATA\n")
        file.write("POINTS {} double\n".format(len(sample.vertices_)))
        for vertexID, vertex in sample.vertices_.items():
            file.write("{:12.5e} {:12.5e} {:12.5e}\n".format(
                vertex.position_[0],
                vertex.position_[1],
                vertex.position_[2]))
        totalPolygonDataPoints = 0
        totalPolygons = 0
        for polygonID,polygon in sample.polygons_.items():
            if not polygon.crossBoundary_:
                totalPolygons += 1
                totalPolygonDataPoints += len(polygon.vertices_) + 1

           
        file.write("POLYGONS {} {}\n".format(totalPolygons,
                                                totalPolygonDataPoints))
    

        for polygonID, polygon in sample.polygons_.items():
            if not polygon.crossBoundary_: 
                file.write("{:<7d}".format(len(polygon.vertices_)))
                for vID in polygon.vertices_: 
                    file.write("{:<7d}".format(vertexMap[vID]))
                file.write("\n")

        file.write("CELL_DATA {}\n".format(totalPolygons))
        file.write("SCALARS shapeIndex double\n")
        file.write("LOOKUP_TABLE default\n")
        for polygonID, polygon in sample.polygons_.items():
            if not polygon.crossBoundary_:
                for cellID, cell in sample.cells_.items():
                    if polygonID in cell.polygons_:
                        file.write("{:<12.6f}\n".format(cell.shape_index_))
                        break
    return
    
def dumpCellVtk(sample:tissueSample.Sample, cellID:int):
    cell = sample.cells_[cellID]
    tmp_vertices = {}
    tmp_polygons = {polygonID:sample.polygons_[polygonID] for polygonID in cell.polygons_}
    for polygonID, polygon in tmp_polygons.items():
        for vertexID in polygon.vertices_:
            tmp_vertices[vertexID] = sample.vertices_[vertexID]
    vertexMap = functions.mapmaker(tmp_vertices)
    with open(
        sample.config_dir_
        + "{:07d}.cell{:03d}.vtk".format(sample.time_, cell.id_),'w') as file:
        
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
        totalPolygonDataPoints = 0
        totalPolygons = 0
        for polygonID,polygon in tmp_polygons.items():
            totalPolygons += 1
            totalPolygonDataPoints += len(polygon.vertices_) + 1
        file.write("POLYGONS {} {}\n".format(totalPolygons,
                                                totalPolygonDataPoints))
        for polygonID, polygon in tmp_polygons.items():
            file.write("{:<7d}".format(len(polygon.vertices_)))
            for vID in polygon.vertices_: 
                file.write("{:<7d}".format(vertexMap[vID]))
            file.write("\n")
    return

def calculateElongationAxis(sample:tissueSample.Sample, cellID:int):
    inertia_tensor = momentOfInertia.calculate_moment_of_inertia_tensor(sample,cellID)
    eigenvalues, eigenvectors = np.linalg.eig(inertia_tensor)
    N = eigenvectors[:, np.argmin(eigenvalues)]
    return N


def evaluatePostDivisionTopology(sample:tissueSample.Sample, cellID:int):
    cell = sample.cells_[cellID]
    newVertices = {}
    newEdges = {}
    newPolygons = {}
    newCells = {}
    motherEdgeIDToDaughterEdgeIDs = {}
    motherPolygonIDToDaughterPolygonIDs = {}

    N = calculateElongationAxis(sample,cell.id_)

    # Identify mother edges and mother polygons. Create intersection vertices on mother edges.
    for polygonID in cell.polygons_:
        for edgeID in sample.polygons_[polygonID].edges_:
            edge = sample.edges_[edgeID]
            v0 = sample.vertices_[edge.vertices_[0]].position_
            v1 = sample.vertices_[edge.vertices_[1]].position_
            t = np.dot(np.subtract(cell.center_,v0),N)/np.dot(np.subtract(v1,v0),N)
            if t>0 and t<1:
                sample.polygons_[polygonID].is_mother_ = True
                if not polygonID in motherPolygonIDToDaughterPolygonIDs:
                    motherPolygonIDToDaughterPolygonIDs[polygonID] = []
                if edgeID not in motherEdgeIDToDaughterEdgeIDs:
                    motherEdgeIDToDaughterEdgeIDs[edgeID] = []
                    new_vertex = topology.Vertex(id=np.max(list(sample.vertices_.keys()))+1+len(newVertices))
                    new_vertex.setPosition(np.add(v0,np.multiply(t,np.subtract(v1,v0))))
                    new_vertex.is_daughter_ = True
                    newVertices[new_vertex.id_] = new_vertex
                    sample.edges_[edgeID].is_mother_ = True
                    sample.edges_[edgeID].intersection_vertex_ = new_vertex.id_


    # Create new dividing edges for the dividing polygon. Add to newEdges.
    # Note:There are more daugher edges that are created in the next step.            
    for polygonID in motherPolygonIDToDaughterPolygonIDs:
        newDividingEdge = topology.Edge(id=np.max(list(sample.edges_.keys()))+1+len(newEdges))
        newDividingEdge.is_in_dividing_polygon_ = True
        newDividingEdge.mother_polygon_id_ = polygonID
        vertices = []
        for edgeID in sample.polygons_[polygonID].edges_:
            if sample.edges_[edgeID].is_mother_:
                vertices.append(sample.edges_[edgeID].intersection_vertex_)
        for vertexID in vertices:
            newDividingEdge.addVertex(vertexID)
        newEdges[newDividingEdge.id_] = newDividingEdge

    # Create daugher edges for the mother edges. Each mother edge is divided into two daughter edges.
    # Add these to newEdges.
    for edgeID in motherEdgeIDToDaughterEdgeIDs:
        motherEdgeIDToDaughterEdgeIDs[edgeID] = [np.max(list(sample.edges_.keys()))+1+len(newEdges),
                                                np.max(list(sample.edges_.keys()))+2+len(newEdges)]
        newDaughterEdge1 = topology.Edge(id = motherEdgeIDToDaughterEdgeIDs[edgeID][0])
        newDaughterEdge2 = topology.Edge(id = motherEdgeIDToDaughterEdgeIDs[edgeID][1])
        newDaughterEdge1.is_daughter_ = True
        newDaughterEdge2.is_daughter_ = True
        newDaughterEdge1.mother_id_ = edgeID
        newDaughterEdge2.mother_id_ = edgeID
        vertices1 = [sample.edges_[edgeID].vertices_[0],sample.edges_[edgeID].intersection_vertex_]
        vertices2 = [sample.edges_[edgeID].vertices_[1],sample.edges_[edgeID].intersection_vertex_]
        for vertexID in vertices1:
            newDaughterEdge1.addVertex(vertexID)
        for vertexID in vertices2:
            newDaughterEdge2.addVertex(vertexID)
        newEdges[newDaughterEdge1.id_] = newDaughterEdge1
        newEdges[newDaughterEdge2.id_] = newDaughterEdge2

    # Create new dividing polygon. Add the dividing edges. to this polygon.
    # Add to newPolygons.
    dividingPolygon = topology.Polygon(id=np.max(list(sample.polygons_.keys()))+1+len(newPolygons))
    dividingPolygon.is_dividing_polygon_ = True
    for edgeID,edge in newEdges.items():
        if edge.is_in_dividing_polygon_:
            dividingPolygon.addEdge(edgeID)     
    newPolygons[dividingPolygon.id_] = dividingPolygon               

    # For each mother polygon there should be two daughter polygons. Each daughter polygon
    # has the corresponding dividing edge, two daughter edges and the other edges of the mother polygon.
    for polygonID in motherPolygonIDToDaughterPolygonIDs:
        motherPolygonIDToDaughterPolygonIDs[polygonID] = [np.max(list(sample.polygons_.keys()))+1+len(newPolygons),
                                                        np.max(list(sample.polygons_.keys()))+2+len(newPolygons)]
        newDaughterPolygon1 = topology.Polygon(id = motherPolygonIDToDaughterPolygonIDs[polygonID][0])
        newDaughterPolygon2 = topology.Polygon(id = motherPolygonIDToDaughterPolygonIDs[polygonID][1])
        newDaughterPolygon1.is_daughter_ = True
        newDaughterPolygon2.is_daughter_ = True
        newDaughterPolygon1.mother_id_ = polygonID
        newDaughterPolygon2.mother_id_ = polygonID
        # Add the dividing edges to both daughter polygons.
        for edgeID,edge in newEdges.items():
            if edge.is_in_dividing_polygon_ and edge.mother_polygon_id_ == polygonID:
                newDaughterPolygon1.addEdge(edgeID)
                newDaughterPolygon2.addEdge(edgeID)
                break
        # Add the daughter edges to the corresponding daughter polygons.
        # newDaughterPolygon1 should be "above" the dividing polygon plane.
        for edgeID,edge in newEdges.items():
            if edge.is_daughter_ and edge.mother_id_ in sample.polygons_[polygonID].edges_:
                if np.dot(
                    np.subtract(
                        sample.vertices_[edge.vertices_[0]].position_,cell.center_),
                        N) > 0:
                    newDaughterPolygon1.addEdge(edgeID)
                elif np.dot(
                    np.subtract(
                        sample.vertices_[edge.vertices_[0]].position_,cell.center_),
                        N) < 0:
                    newDaughterPolygon2.addEdge(edgeID)
                else: 
                    print("Error: The non intersecting vertex of the daughter edge is on the dividing plane.")
        
        # add the other edges from the mother polygon to the corresponding daughter polygons.
        for edgeID in sample.polygons_[polygonID].edges_:
            if not sample.edges_[edgeID].is_mother_:
                if np.dot(
                    np.subtract(
                        sample.vertices_[sample.edges_[edgeID].vertices_[0]].position_,cell.center_),
                        N) > 0:
                    newDaughterPolygon1.addEdge(edgeID)
                elif np.dot(
                    np.subtract(
                        sample.vertices_[sample.edges_[edgeID].vertices_[0]].position_,cell.center_),
                        N) < 0:
                    newDaughterPolygon2.addEdge(edgeID)
                else: 
                    print("Error: One of the edge vertices in the mother polygon is on the dividing plane.")
        newPolygons[newDaughterPolygon1.id_] = newDaughterPolygon1
        newPolygons[newDaughterPolygon2.id_] = newDaughterPolygon2
    # Arrange polygon vertices for each polygon in newPolygons.
    for polygonID,polygon in newPolygons.items():
        tmp_vertices=[]
        for edgeID in polygon.edges_:
            if edgeID in newEdges:
                tmp_vertices.append(newEdges[edgeID].vertices_)
            if edgeID in sample.edges_:
                tmp_vertices.append(sample.edges_[edgeID].vertices_)
        polygon.vertices_ = functions.arrange_polygon(tmp_vertices)

    # two new daughter cells.
    newDaughterCell1 = topology.Cell(id=np.max(list(sample.cells_.keys()))+1)
    newDaughterCell2 = topology.Cell(id=np.max(list(sample.cells_.keys()))+2)

    for polygonID,polygon in newPolygons.items():
        # Add the dividing polygon to both daughter cells.
        if polygon.is_dividing_polygon_:
            newDaughterCell1.addPolygon(polygonID)
            newDaughterCell2.addPolygon(polygonID)
        # Add the daughter polygons to the corresponding daughter cells.    
        signSet = set()
        for vertexID in polygon.vertices_:
            if vertexID in newVertices:
                if newVertices[vertexID].is_daughter_:
                    continue
                signSet.add(np.sign(np.dot(np.subtract(newVertices[vertexID].position_,cell.center_),N)))
            elif vertexID in sample.vertices_:
                if sample.vertices_[vertexID].is_daughter_:
                    continue
                signSet.add(np.sign(np.dot(np.subtract(sample.vertices_[vertexID].position_,cell.center_),N)))
        if len(signSet) > 1:
            print("Error: The polygon has vertices on both sides of the dividing plane.")
        elif len(signSet) == 1:
            if signSet.pop() > 0:
                newDaughterCell1.addPolygon(polygonID)
            else:
                newDaughterCell2.addPolygon(polygonID)
    # Add the remaining polygons from the mother cell to the corresponding daughter cells.
    for polygonID in cell.polygons_:
        polygon = sample.polygons_[polygonID]
        if not polygon.is_mother_:
            signSet = set()
            for vertexID in polygon.vertices_:
                if vertexID in newVertices:
                    if newVertices[vertexID].is_daughter_:
                        continue
                    signSet.add(np.sign(np.dot(np.subtract(newVertices[vertexID].position_,cell.center_),N)))
                elif vertexID in sample.vertices_:
                    if sample.vertices_[vertexID].is_daughter_:
                        continue
                    signSet.add(np.sign(np.dot(np.subtract(sample.vertices_[vertexID].position_,cell.center_),N)))
            if len(signSet) > 1:
                print("Error: The polygon has vertices on both sides of the dividing plane.")
            elif len(signSet) == 1:
                if signSet.pop() > 0:
                    newDaughterCell1.addPolygon(polygonID)
                else:
                    newDaughterCell2.addPolygon(polygonID)

    newCells[newDaughterCell1.id_] = newDaughterCell1
    newCells[newDaughterCell2.id_] = newDaughterCell2

    # update the topology in samples

    sample.vertices_.update(newVertices)
    sample.edges_.update(newEdges)
    for edgeID in motherEdgeIDToDaughterEdgeIDs:
        del sample.edges_[edgeID]
    sample.polygons_.update(newPolygons)
    for polygonID in motherPolygonIDToDaughterPolygonIDs:
        del sample.polygons_[polygonID]
    for polygonID,polygon in sample.polygons_.items():
        rearrangePolygonVerticesFlag = False
        for edgeID in polygon.edges_:
            if edgeID in motherEdgeIDToDaughterEdgeIDs:
                rearrangePolygonVerticesFlag = True
                polygon.edges_.remove(edgeID)
                for daughterEdgeID in motherEdgeIDToDaughterEdgeIDs[edgeID]:
                    polygon.addEdge(daughterEdgeID)
        if rearrangePolygonVerticesFlag:
            tmp_vertices=[]
            for edgeID in polygon.edges_:
                if edgeID in sample.edges_:
                    tmp_vertices.append(sample.edges_[edgeID].vertices_)
            polygon.vertices_ = functions.arrange_polygon(tmp_vertices)
    sample.cells_.update(newCells)
    del sample.cells_[cell.id_]
    for cellID,cell in sample.cells_.items():
        for polygonID in cell.polygons_:
            if polygonID in motherPolygonIDToDaughterPolygonIDs:
                cell.polygons_.remove(polygonID)
                for daughterPolygonID in motherPolygonIDToDaughterPolygonIDs[polygonID]:
                    cell.addPolygon(daughterPolygonID)
    print("Daughter cell IDs: ", newDaughterCell1.id_, newDaughterCell2.id_)
    dumpCellVtk(sample, newDaughterCell1.id_)
    dumpCellVtk(sample, newDaughterCell2.id_)

    return sample

def dumpSample(sample):
    with open("sample.topo", "w") as file:
        file.write("vertices {:d}\n".format(len(sample.vertices_)))
        for key, vertex in sample.vertices_.items():
            id = vertex.id_
            x = vertex.position_[0]
            y = vertex.position_[1]
            z = vertex.position_[2]
            file.write("{:6d} {:12.5e} {:12.5e} {:12.5e}\n".format(id, x, y, z))
        file.write("edges {:d}\n".format(len(sample.edges_)))

        for key, edge in sample.edges_.items():
            file.write("{:d}".format(edge.id_))
            for vertexID in edge.vertices_:
                file.write(" {:6d}".format(vertexID))
            file.write("\n")
        
        file.write("polygons {:d}\n".format(len(sample.polygons_)))
        for key,polygon in sample.polygons_.items():
            file.write("{:d}".format(polygon.id_))
            for edgeID in polygon.edges_:
                file.write(" {:6d}".format(edgeID))
            file.write("\n")

        file.write("cells {:d}\n".format(len(sample.cells_)))
        for key,cell in sample.cells_.items():
            file.write("{:d}".format(cell.id_))
            for polygonID in cell.polygons_:
                file.write(" {:6d}".format(polygonID))
            file.write("\n")
    return

def main():
    sample = tissueSample.Sample(configDir = "samples/", simulationTime = 500, tissueType = "periodic")
    crossBoundary = True
    while crossBoundary:
        cellID = random.choice(list(sample.cells_.keys()))
        crossBoundary = sample.cells_[cellID].crossBoundary_
    sample.cells_[cellID].is_mother_ = True
    print("Mother cell ID: ", cellID)
    dumpCellVtk(sample, cellID)
    sample = evaluatePostDivisionTopology(sample, cellID)
    dumpSample(sample)
    return

if __name__ == "__main__":
    main()