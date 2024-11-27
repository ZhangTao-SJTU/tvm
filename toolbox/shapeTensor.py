from toolbox import tissueSample
import numpy as np

def calculate_shape_tensor(sample:tissueSample.Sample, cellID:int):
    cell = sample.cells_[cellID]
    shapeTensor = np.zeros((3,3))
    for vertexID in cell.vertices_:
        r_prime = np.subtract(sample.vertices_[vertexID].position_,cell.center_)
        shapeTensor = np.add(
            shapeTensor,
            np.outer(r_prime,r_prime))
    # for polygonID in cell.polygons_:
    #     r_prime = np.subtract(sample.polygons_[polygonID].center_,cell.center_)
    #     shapeTensor = np.add(
    #         shapeTensor,
    #         np.outer(r_prime,r_prime))
    # shapeTensor /= len(cell.vertices_)+len(cell.polygons_)
    return shapeTensor