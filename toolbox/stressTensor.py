import numpy as np
import math
from toolbox import tissueSample

# corrected second term in calculate_surface_term with Ligesh 7/30/2024

# Assuming that the polygon vertices are arranged in cyclic order, 
# this function ensures that they are in anticlockwise order
# (with respect to the vector from cell center to polygon center)
# Polygon vertices in sample are rearranged as a side effect.
# This function returns nothing.
def rearrange_polygon_vertices_for_cell(
        sample:tissueSample.Sample,
        cellID:int) -> None:
    cell = sample.cells_[cellID]
    for polygonID in cell.polygons_:
        polygon = sample.polygons_[polygonID]
        sum_of_cross_products = np.zeros(3)
        for i in range(len(polygon.vertices_)):
            v1 = polygon.vertices_[i]
            v2 = polygon.vertices_[(i + 1) % len(polygon.vertices_)]
            vector1 = np.subtract(
                sample.vertices_[v1].position_,
                cell.center_)
            vector2 = np.subtract(
                sample.vertices_[v2].position_,
                cell.center_)
            cross_product = np.cross(vector1, vector2)
            sum_of_cross_products = np.add(
                sum_of_cross_products,
                cross_product)
        sign = np.sign(
            np.dot(
                sum_of_cross_products,
                np.subtract(polygon.center_, cell.center_)))
        if sign == -1:
            # print("Reversing order polygon vertices for polygon {}".format(polygonID))
            polygon.vertices_ = polygon.vertices_[::-1]
    return

def calculate_volume_term(sample:tissueSample.Sample, cellID:int):
    cell = sample.cells_[cellID]
    # Initialize the volume term with identity matrix
    tensor_product = np.identity(3)
    # Add contributions from each polygon
    for polygonID in cell.polygons_:
        polygon = sample.polygons_[polygonID]
        P = polygon.perimeter_
        # For this polygon, evaluate the sum of r_i x r_{i+1}
        sum = np.zeros(3)
        for i in range(len(polygon.vertices_)):
            nowID = polygon.vertices_[i]
            nextID = polygon.vertices_[(i + 1) % len(polygon.vertices_)]
            nowVertex = sample.vertices_[nowID].position_
            nextVertex = sample.vertices_[nextID].position_
            sum = np.add(sum,np.cross(nowVertex,nextVertex))
        # For this polygon, evaluate the addition to the tensor product
        for i in range(len(polygon.vertices_)):
            nowID = polygon.vertices_[i]
            nextID = polygon.vertices_[(i + 1) % len(polygon.vertices_)]
            nowVertex = sample.vertices_[nowID].position_
            nextVertex = sample.vertices_[nextID].position_
            edgeVector = np.subtract(nextVertex,nowVertex)
            edgeCenter = np.add(nowVertex,nextVertex) / 2
            diffCenter = np.subtract(edgeCenter,polygon.center_)
            l = np.linalg.norm(edgeVector)
            scalarMultiple = np.dot(diffCenter, sum) / (P * l * 6)
            tensor_product = np.add(
                tensor_product,
                np.outer(edgeVector,edgeVector) * scalarMultiple)
            
    return 2 * sample.kv_ * (cell.volume_ - 1) * tensor_product

def calculate_surface_area_term(sample:tissueSample.Sample, cellID:int):
    cell = sample.cells_[cellID]
    # Initialize the surface term S*I
    tensor_product = cell.surface_area_ * np.identity(3)
    # Subtract the contributions from each polygon in the second term

    for polygonID in cell.polygons_:
        polygon = sample.polygons_[polygonID]
        for i in range(len(polygon.vertices_)):
            nowID = polygon.vertices_[i]
            nextID = polygon.vertices_[(i + 1) % len(polygon.vertices_)]
            nowVertex = sample.vertices_[nowID].position_
            nextVertex = sample.vertices_[nextID].position_
            areaVector = np.cross(
                np.subtract(nowVertex,polygon.center_),
                np.subtract(nextVertex,polygon.center_))/ 2
            area = np.linalg.norm(areaVector)

            tensor_product = np.subtract(
                tensor_product,
                np.outer(areaVector,areaVector) / area)

    # Subtract the contributions from each polygon in the third (i.e final) term
    for polygonID in cell.polygons_:
        polygon = sample.polygons_[polygonID]
        P = polygon.perimeter_
        # For this polygon, evaluate the sum of l_i x a_i / |a_i|
        sum = np.zeros(3)
        for i in range(len(polygon.vertices_)):
            nowID = polygon.vertices_[i]
            nextID = polygon.vertices_[(i + 1) % len(polygon.vertices_)]
            nowVertex = sample.vertices_[nowID].position_
            nextVertex = sample.vertices_[nextID].position_
            edgeVector = np.subtract(nextVertex, nowVertex)
            # No need to divide by 2 here, as we will use the unit vector)
            areaVector = np.cross(
                np.subtract(nowVertex,polygon.center_),
                np.subtract(nextVertex,polygon.center_))
            sum = np.add(
                sum,
                np.cross(edgeVector, areaVector) / np.linalg.norm(areaVector))
        # Now add the contribution from this polygon to the tensor product
        for i in range(len(polygon.vertices_)):
            nowID = polygon.vertices_[i]
            nextID = polygon.vertices_[(i + 1) % len(polygon.vertices_)]
            nowVertex = sample.vertices_[nowID].position_
            nextVertex = sample.vertices_[nextID].position_
            edgeVector = np.subtract(nextVertex, nowVertex)
            edgeCenter = np.add(nowVertex, nextVertex) / 2
            diffCenter = np.subtract(edgeCenter, polygon.center_)
            l = np.linalg.norm(edgeVector)
            scalarMultiple = np.dot(diffCenter, sum) / (P * l * 2)
            tensor_product = np.subtract(
                tensor_product,
                np.outer(edgeVector, edgeVector) * scalarMultiple)
    return 2 * (cell.surface_area_ - sample.s0_) * tensor_product

def calculate_gamma_term(sample:tissueSample.Sample, cellID:int):
    cell = sample.cells_[cellID]
    if not cell.is_surface_:
        return np.zeros((3,3))
    # Initialize the surface term S*I
    tensor_product = cell.surface_area_ * np.identity(3)
    # Subtract the contributions from each polygon in the second term
    for polygonID in cell.polygons_:
        polygon = sample.polygons_[polygonID]
        if polygon.is_surface_:
            for i in range(len(polygon.vertices_)):
                nowID = polygon.vertices_[i]
                nextID = polygon.vertices_[(i + 1) % len(polygon.vertices_)]
                nowVertex = sample.vertices_[nowID].position_
                nextVertex = sample.vertices_[nextID].position_
                areaVector = np.cross(
                    np.subtract(nowVertex,polygon.center_),
                    np.subtract(nextVertex,polygon.center_))/ 2
                area = np.linalg.norm(areaVector)
                tensor_product = np.subtract(
                    tensor_product,
                    np.outer(areaVector,areaVector) / area)

    # Subtract the contributions from each polygon in the third (i.e final) term
    for polygonID in cell.polygons_:
        polygon = sample.polygons_[polygonID]
        if polygon.is_surface_:
            P = polygon.perimeter_
            # For this polygon, evaluate the sum of l_i x a_i / |a_i|
            sum = np.zeros(3)
            for i in range(len(polygon.vertices_)):
                nowID = polygon.vertices_[i]
                nextID = polygon.vertices_[(i + 1) % len(polygon.vertices_)]
                nowVertex = sample.vertices_[nowID].position_
                nextVertex = sample.vertices_[nextID].position_
                edgeVector = np.subtract(nextVertex, nowVertex)
                areaVector = np.cross(
                    np.subtract(nowVertex,polygon.center_),
                    np.subtract(nextVertex,polygon.center_))
                sum = np.add(
                    sum,
                    np.cross(edgeVector, areaVector) / np.linalg.norm(areaVector))
            # Now add the contribution from this polygon to the tensor product
            for i in range(len(polygon.vertices_)):
                nowID = polygon.vertices_[i]
                nextID = polygon.vertices_[(i + 1) % len(polygon.vertices_)]
                nowVertex = sample.vertices_[nowID].position_
                nextVertex = sample.vertices_[nextID].position_
                edgeVector = np.subtract(nextVertex, nowVertex)
                edgeCenter = np.add(nowVertex, nextVertex) / 2
                diffCenter = np.subtract(edgeCenter, polygon.center_)
                l = np.linalg.norm(edgeVector)
                scalarMultiple = np.dot(diffCenter, sum) / (P * l * 2)
                tensor_product = np.subtract(
                    tensor_product,
                    np.outer(edgeVector, edgeVector) * scalarMultiple)
    return sample.gamma_ * tensor_product

def calculate_stress_tensor(sample:tissueSample.Sample, cellID:int):
    rearrange_polygon_vertices_for_cell(sample, cellID)
    stressTensor = np.add(
        calculate_volume_term(sample, cellID),
        calculate_surface_area_term(sample, cellID))
    
    # if sample.cells_[cellID].is_surface_:
    #     stressTensor = np.add(
    #         stressTensor,
    #         calculate_gamma_term(sample, cellID))
    # sample.cells_[cellID].stress_tensor_ = stressTensor
    return stressTensor

def calculate_stress_invariants(stress_tensor):
    egvals = np.linalg.eigvalsh(stress_tensor)
    hydrostatic = (egvals[0] + egvals[1] + egvals[2]) / 3
    max_shear = 0.5 * abs(egvals[-1] - egvals[0])
    # please check!
    von_mises = math.sqrt(
        0.5 * (
            (egvals[1] - egvals[0])**2
            + (egvals[2] - egvals[1])**2
            + (egvals[0] - egvals[2])**2))
    return {
        "hydrostatic" : hydrostatic,
        "max_shear" : max_shear,
        "von_mises" : von_mises}

def calculate_stress_tensor_COM_center(sample:tissueSample.Sample, cellID:int):
    rearrange_polygon_vertices_for_cell(sample, cellID)
    stress_tensor = np.zeros((3,3))
    cell = sample.cells_[cellID]
    if not cell.volume_:
        raise ValueError("Cell volume is zero or not calculated")
    # Add volume term:
    v_term = 2*sample.kv_*(cell.volume_-cell.v0_)*np.identity(3)
    
    # Calculate surface area term:

    s_term = cell.surface_area_ * np.identity(3)
    # Subtract the contributions from each polygon in the second term
    total_area = 0
    for polygonID in cell.polygons_:
        polygon = sample.polygons_[polygonID]
        for i in range(len(polygon.vertices_)):
            nowID = polygon.vertices_[i]
            nextID = polygon.vertices_[(i + 1) % len(polygon.vertices_)]
            nowVertex = sample.vertices_[nowID].position_
            nextVertex = sample.vertices_[nextID].position_
            areaVector = np.cross(
                np.subtract(nowVertex,polygon.center_),
                np.subtract(nextVertex,polygon.center_))/ 2
            area = np.linalg.norm(areaVector)
            total_area += area
            s_term = np.subtract(
                s_term,
                np.outer(areaVector,areaVector) / area)
    
    s_term = (2/cell.volume_)*(cell.surface_area_-cell.s0_)*s_term
    if not (cell.surface_area_ - total_area < 1e-4):
        raise ValueError("sanity check failed cell area: {} total area from triangles: {}".format(cell.surface_area_, total_area))

    stress_tensor = (-1) * (v_term + s_term)
    return stress_tensor

