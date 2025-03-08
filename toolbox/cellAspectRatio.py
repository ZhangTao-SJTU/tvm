import numpy as np
from toolbox import momentOfInertia
from toolbox import tissue

def calculate_aspect_ratio(sample:tissue.Sample, cellID:int):
    inertiaTensor = momentOfInertia.calculate_moment_of_inertia_tensor(sample, cellID)
    eigenvalues = np.linalg.eigvalsh(inertiaTensor)
    # np.linalg.eigvalsh returns eigenvalues in ascending order.
    return np.sqrt(eigenvalues[-1]/eigenvalues[0])

def calculate_shape_tensor(sample:tissue.Sample, cellID:int):
    cell = sample.cells_[cellID]
    shapeTensor = np.zeros((3,3))
    for vertexID in cell.vertices_:
        r_prime = np.subtract(sample.vertices_[vertexID].position_,cell.center_)
        shapeTensor = np.add(
            shapeTensor,
            np.outer(r_prime,r_prime))
    shapeTensor /= len(cell.vertices_)
    return shapeTensor
# Shape tensor as seen in Nestor-Bergmann et al. 2021
def calculate_aspect_ratio_from_shape_tensor(sample:tissue.Sample, cellID:int):
    shapeTensor = calculate_shape_tensor(sample, cellID)
    eigenvalues = np.linalg.eigvalsh(shapeTensor)
    return np.sqrt(eigenvalues[-1]/eigenvalues[0])