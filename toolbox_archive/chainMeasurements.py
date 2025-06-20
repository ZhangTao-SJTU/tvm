from toolbox import tissue
from toolbox import shape
from toolbox import stress
from toolbox import functions
import os
import numpy as np


class Chain:
    def __init__(self,id:int,sample:tissue.Sample, cellIDs:list[int]):
        self.id_ = id
        self.cellIDs_ = cellIDs
        self.polygonIDs_ = []
        for cellID in cellIDs:
            self.polygonIDs_.extend(sample.cells_[cellID].polygons_)
        self.polygonIDs_ = list(set(self.polygonIDs_))
        self.volume_ = None
        self.surface_area_ = None
        self.calculate_volume(sample)
        self.calculate_surface_area(sample)
        self.shape_index_ = self.surface_area_/self.volume_**(2/3)
        return
    

    def calculate_volume(self, sample:tissue.Sample):
        self.volume_ = 0
        for cellID in self.cellIDs_:
            self.volume_ += sample.cells_[cellID].volume_
        return
    
    def calculate_surface_area(self, sample:tissue.Sample):
        # first identify the boundary polygons
        # then add up their areas.
        self.surface_area_ = 0
        surface_polygons = []
        for cellID in self.cellIDs_:
            intersection = set(sample.cells_[cellID].polygons_).intersection(surface_polygons)
            for polygonID in sample.cells_[cellID].polygons_:
                if polygonID not in intersection:
                    surface_polygons.append(polygonID)
        for polygonID in surface_polygons:
            self.surface_area_ += sample.polygons_[polygonID].area_
        return
    
    def dump_vtk(self, sample:tissue.Sample):
        output_dir = sample.config_dir_+"chains/"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        chain_vertices = []
        chain_polygons = []

        for cellID in self.cellIDs_:
            cell = sample.cells_[cellID]
            for polygonID in cell.polygons_:
                chain_vertices.extend(sample.polygons_[polygonID].vertices_)
                chain_polygons.append(polygonID)

        chain_vertices = list(set(chain_vertices))
        chain_polygons = list(set(chain_polygons))
    
        # Construct vMap, this will be useful in writing the VTK file.
        vMap = functions.mapmaker(chain_vertices)
        nPolygons = len(chain_polygons)
        nVertices = 0
        for polygonID in chain_polygons:
            nVertices += len(sample.polygons_[polygonID].vertices_)
        
        with open(output_dir + "{:07d}.{}.chain.vtk".format(sample.time_, self.id_),"w") as file:
            file.write("# vtk DataFile Version 2.0\n")
            file.write("polydata\n")
            file.write("ASCII\n")
            file.write("DATASET POLYDATA\n")
            file.write("POINTS {} double\n".format(len(chain_vertices)))    
            for vertexID in chain_vertices:
                file.write("{} {} {}".format(
                    sample.vertices_[vertexID].position_[0],
                    sample.vertices_[vertexID].position_[1],
                    sample.vertices_[vertexID].position_[2]))
                file.write("\n")
            file.write("\nPOLYGONS {} {}\n".format(
                nPolygons, nPolygons + nVertices))
            for polygonID in chain_polygons:
                file.write("{} ".format(len(sample.polygons_[polygonID].vertices_)))
                for vertexID in sample.polygons_[polygonID].vertices_:
                    file.write("{} ".format(vMap[vertexID]))
                file.write("\n")    
            # file.write( "\nCELL_DATA {}\n".format(nPolygons))
            # file.write( "SCALARS type double\n")
            # file.write( "LOOKUP_TABLE default\n")
            # for polygonID in chain_polygons:
            #     file.write("{}\n".format(abs(sample.polygons_[polygonID].scalar_)))
        return

# mark cells in sample (via side effect)
def mark_chain_cells(sample:tissue.Sample):
    for cellID,cell in sample.cells_.items():
        if cell.type_:
            shp = shape.calculate_shape_tensor(sample,cellID)
            str = stress.calculate_stress_tensor(sample,cellID)
            _, shape_egvecs = np.linalg.eigh(shp)
            stress_egvals, stress_egvecs = np.linalg.eigh(str)
            proj = np.dot(stress_egvecs[-1],shape_egvecs[0])
            if abs(proj) > 0.8:
                cell.is_in_chain_ = True
            # max_shear = stress_egvals[-1]-stress_egvals[0]
            # if max_shear>0.02:
            #     cell.is_in_chain_ = True
    return

def evaluate_chains_dict(sample:tissue.Sample):
    # mark cells that will be part of a chain
    mark_chain_cells(sample)
    # Create a roster for cells that are part of some chain
    cell_assignment = {}
    for cellID,cell in sample.cells_.items():
        if cell.type_ and cell.is_in_chain_:
            cell_assignment[cellID] = False
    chains = {}
    for i,currentID in enumerate(list(cell_assignment.keys())):
        if cell_assignment[currentID]:
            continue
        cell_assignment[currentID] = True
        chain = []
        chain.append(currentID)
        # check all the following cells for potential chain membership.
        # if they are already part of a chain skip them.
        # otherwise they are a potential candidate.
        for j in range(i+1,len(cell_assignment)):
            candidateID = list(cell_assignment.keys())[j]
            if cell_assignment[candidateID]:
                continue
            candidateCell = sample.cells_[candidateID]
            # Now check if the candidate cell is a neighbor of any of the cells in the chain.
            # If it is, add it to the chain.
            # If it is not, continue to the next candidate.
            for assignedID in chain:
                assignedCell = sample.cells_[assignedID]

                if len(set(assignedCell.polygons_).intersection(candidateCell.polygons_)):
                    chain.append(candidateID)
                    cell_assignment[candidateID] = True
                    break
        # if this chain has more than one cell, include it in the "chains" dictionary.
        if len(chain)>1:
            chains[len(chains)] = chain

    # remove chains that are only one cell long
    # chains = {k:v for k,v in chains.items() if len(v)>1}

    # Transform the dictionary so that the entry is a chain object
    # (rather than a list of chain cell IDS) 
    for id, chain_list in chains.items():
        chains[id] = Chain(id=id,sample=sample,cellIDs=chain_list)
    return chains
