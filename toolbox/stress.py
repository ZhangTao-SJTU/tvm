import numpy as np

def calculate_max_shear_stress(sample, cellID:int):
    stress = calculate_stress_tensor(sample,cellID)
    egvals = np.linalg.eigvalsh(stress)
    max_shear = 0.5 * abs(egvals[-1] - egvals[0])
    return max_shear

def calculate_stress_tensor(sample, cellID:int):
    cell = sample.cells_[cellID]
    if cell.crossBoundary_:
        return calculate_stress_tensor_boundary_cell(sample, cellID)
    tensor_term = np.zeros((3,3))
    for polygonID in cell.polygons_:
        polygon = sample.polygons_[polygonID]
        # cell_center_to_polygon = np.subtract(polygon.center_,cell.center_)
        for edgeID in polygon.edges_:
            edge = sample.edges_[edgeID]
            v1 = np.subtract(
                sample.vertices_[edge.vertices_[0]].position_,
                polygon.center_)
            v2 = np.subtract(
                sample.vertices_[edge.vertices_[1]].position_,
                polygon.center_)
            a_poly = 0.5 * np.cross(v1, v2)
            ## No need to resolve orientation here; the tensor term is quadratic in a_poly
            tensor_term = np.add(
                tensor_term,
                np.outer(a_poly, a_poly) / np.linalg.norm(a_poly))
    factor = 2 * (cell.surface_area_ - cell.s0_) / cell.volume_
    stress_tensor = -2 * sample.kv_ * (cell.volume_ - cell.v0_) * np.identity(3)
    stress_tensor = np.add(stress_tensor, -1 * factor * cell.surface_area_ * np.identity(3))
    stress_tensor = np.add(stress_tensor, factor * tensor_term)
    return stress_tensor

def calculate_stress_tensor_boundary_cell(sample, cellID:int):
    test_b_cell = sample.extract_cell(cellID)
    for _,cell in test_b_cell.cells_.items():
        cell.crossBoundary_ = False
    for _, polygon in test_b_cell.polygons_.items():
        polygon.crossBoundary_ = False
    for _, edge in test_b_cell.edges_.items():
        edge.crossBoundary_ = False
    axes_to_flip = []  # Flip all axes
    for _,edge in test_b_cell.edges_.items():
        v0 = test_b_cell.vertices_[edge.vertices_[0]].position_
        v1 = test_b_cell.vertices_[edge.vertices_[1]].position_
        edge_vector = np.subtract(v1, v0)
        for i in range(3):
            if abs(edge_vector[i]) > sample.boxSize_/2 :
                axes_to_flip.append(i)
    axes_to_flip = list(set(axes_to_flip))  # Remove duplicates
    for _, vertex in test_b_cell.vertices_.items():
        for i in axes_to_flip:
            if vertex.position_[i] < sample.boxSize_ / 2:
                continue
            vertex.position_[i] -= sample.boxSize_
    test_b_cell.arrange_polygon_vertices()
    test_b_cell.calculate_COM_polygon_centers()
    test_b_cell.calculate_cell_centers()
    test_b_cell.calculate_cell_volumes()
    test_b_cell.calculate_polygon_areas()
    test_b_cell.calculate_cell_surface_areas()
    test_b_cell.calculate_cell_shape_indices()
    return calculate_stress_tensor(test_b_cell, cellID)
