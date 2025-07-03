from toolbox import tissueSample
from toolbox import stressTensor
from toolbox import cellAspectRatio
import vtk
import numpy as np
import math
import copy

def create_vtk_arrow(
        origin = np.array([0,0,0]),
        direction = np.array([1,0,0]),
        filename = "arrow.vtk"):
    
    x = [1,0,0]
    angle = math.degrees(math.acos(np.dot(x,direction)))
    axis = np.cross(x,direction)

    # This is an arrow with base at [0,0,0] and tip at [1,0,0]
    arrow_source = vtk.vtkArrowSource()

    # Create a transform to 
    # (a) rotate the arrow towards the direction, and 
    # (b) translate the base to the center
    transform = vtk.vtkTransform()
    transform.Translate(origin)
    transform.RotateWXYZ(angle, axis[0], axis[1], axis[2])  
    # Create a transform filter to apply the transform
    transform_filter = vtk.vtkTransformPolyDataFilter()
    transform_filter.SetTransform(transform)
    transform_filter.SetInputConnection(arrow_source.GetOutputPort())
    transform_filter.Update()

    # Write the output to a VTK file
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName("{}".format(filename))
    writer.SetInputConnection(transform_filter.GetOutputPort())
    writer.Write()
    return

def create_eigenvector_arrows(origin, matrix, filename_base):
    _,egvecs = np.linalg.eigh(matrix)
    for i,vec in enumerate(egvecs):
        create_vtk_arrow(
            origin = origin,
            direction = vec, 
            filename = "{}{}.vtk".format(filename_base,i))
    