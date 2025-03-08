import numpy as np
from toolbox import tissue
# Checked with Ligesh on 07/30/2024
# Found and corrected factor errors in a',b',c' calculations
def calculate_tetrahedral_inertia_tensor(vertices:list[list[float]]):
    # Formula from
    # https://docsdrive.com/pdfs/sciencepublications/jmssp/2005/8-11.pdf

    # Calculates the moment of inertia tensor about origin
    # of a tetrahedron with vertices at the origin,
    # and at the positions in the "vertices" argument. 
    # 
    # Hence, the "vertices" argument in the function should be of the form
    # vertices = [[x0,y0,z0],[x1,y1,z1],[x2,y2,z2]]

    determinant = np.linalg.det(vertices)
    x0, y0, z0 = vertices[0][0], vertices[0][1], vertices[0][2]
    x1, y1, z1 = vertices[1][0], vertices[1][1], vertices[1][2]
    x2, y2, z2 = vertices[2][0], vertices[2][1], vertices[2][2]

    a = determinant * (
        y0**2 + y1**2 + y2**2
        + z0**2 + z1**2 + z2**2
        + y0*y1 + y0*y2 + y1*y2
        + z0*z1+ z0*z2 + z1*z2) / 60
    b = determinant * (
        x0**2 + x1**2 + x2**2
        + z0**2 + z1**2 + z2**2
        + x0*x1 + x0*x2 + x1*x2
        + z0*z1 + z0*z2 + z1*z2) / 60
    c = determinant * (
        x0**2 + x1**2 + x2**2
        + y0**2 + y1**2 + y2**2
        + x0*x1 + x0*x2 + x1*x2
        + y0*y1 + y0*y2 + y1*y2) / 60
    a_prime = determinant * (
        2*y0*z0 + y0*z1 + y0*z2
        + y1*z0 + 2*y1*z1 + y1*z2
        + y2*z0 + y2*z1 + 2*y2*z2) / 120
    b_prime = determinant * (
        2*x0*z0 + x0*z1 + x0*z2
        + x1*z0 + 2*x1*z1 + x1*z2
        + x2*z0 + x2*z1 + 2*x2*z2) / 120
    c_prime = determinant * (
        2*x0*y0 + x0*y1 + x0*y2
        + x1*y0 + 2*x1*y1 + x1*y2
        + x2*y0 + x2*y1 + 2*x2*y2) / 120

    tetrahedral_inertia_tensor = np.array([
        [a, - b_prime, - c_prime],
        [- b_prime, b, - a_prime],
        [- c_prime, - a_prime, c]])
    return tetrahedral_inertia_tensor

def calculate_moment_of_inertia_tensor(sample:tissue.Sample, cellID:int):
    inertiaTensor = np.zeros((3,3))
    cell = sample.cells_[cellID]
    for polygonID in cell.polygons_:
        polygon_center = np.subtract(sample.polygons_[polygonID].center_,cell.center_)
        for edgeID in sample.polygons_[polygonID].edges_:
            id1 = sample.edges_[edgeID].vertices_[0]
            id2 = sample.edges_[edgeID].vertices_[1]
            vertex1 = np.subtract(sample.vertices_[id1].position_,cell.center_)
            vertex2 = np.subtract(sample.vertices_[id2].position_,cell.center_)
            tetrahedral_vertices = [polygon_center,vertex1,vertex2]
            if np.sign(np.linalg.det(tetrahedral_vertices)) == -1:
                tetrahedral_vertices[1] , tetrahedral_vertices[2] = tetrahedral_vertices[2] , tetrahedral_vertices[1]
            inertiaTensor = np.add(
                inertiaTensor,
                calculate_tetrahedral_inertia_tensor(tetrahedral_vertices))
    return inertiaTensor
