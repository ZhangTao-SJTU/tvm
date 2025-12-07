from toolbox.periodic import PeriodicTissue
from toolbox.training import Training
from toolbox.overlap import calculate_Q
from toolbox import stress
import numpy as np
import os
import random
import pandas as pd

class CellDriver(Training):
    def __init__(self):
        super().__init__()
        self._target_cell = None
        self._target_cell_neighbors = None
        self._direction = None
        self._target_position = None
        self._distance = None
        # self._clamping_max_iters = 10
        # self._clamping_s0_lower_limit = 4.6
        # self._clamping_correction_factor = 0.1
        # self._clamping_FIRE_only = True
        # self._clamping_tolerance = 1e-15
        
    @classmethod
    def periodic_tissue(cls,tissue):
        inst = super().periodic_tissue(tissue)
        inst._modified_cells = list(inst._config.cells_.keys())
        return inst
    def set_direction(self,direction):
        norm = np.linalg.norm(direction)
        if norm == 0:
            raise ValueError("Direction vector cannot be zero")
        self._direction = direction / norm
    def set_target_cell(self,target_cell):
        self._target_cell = target_cell
        # For checking the above functionality with vtk:
        for polygonID,polygon in self._config.polygons_.items():
            polygon.vtk_scalar_ = 0
        cell = self._config.cells_[self._target_cell]
        for polygonID in cell.polygons_:
            polygon = self._config.polygons_[polygonID]
            polygon.vtk_scalar_ = 1
        self._config.write_periodic_vtk(filename = "target_cells.vtk", use_scalar=True)
        self._config.write_cell_collection_vtk(cells_array=[self._target_cell], filename="target_cell_only.vtk")
    def set_central_target_cell(self):
        # Pick the cell closest to the center of the box
        box_center = np.array([self._config.boxSize_/2,self._config.boxSize_/2,self._config.boxSize_/2])
        min_distance = float('inf')
        central_cell_id = None
        for cellID,cell in self._config.cells_.items():
            if cell.crossBoundary_:
                continue
            distance = np.linalg.norm(np.subtract(cell.center_,box_center))
            if distance < min_distance:
                min_distance = distance
                central_cell_id = cellID
        self.set_target_cell(central_cell_id)
    def set_target_position(self, position):
        self._target_position = position
    def evaluate_target_cell_neighbors(self, n_layers = 1):
        self._config.evaluate_cell_neighbors()
        self._target_cell_neighbors = self._config.cell_neighbors_[self._target_cell]
        if n_layers == 1:
            return
        for _ in range(1, n_layers):
            new_neighbors = set()
            for neighbor_cell_id in self._target_cell_neighbors:
                neighbor_neighbors = self._config.cell_neighbors_[neighbor_cell_id]
                for nn_id in neighbor_neighbors:
                    if nn_id != self._target_cell and nn_id not in self._target_cell_neighbors:
                        new_neighbors.add(nn_id)
            self._target_cell_neighbors.extend(new_neighbors)
    def set_s0_gradient(self):
        for cellID,cell in self._config.cells_.items():
            if cell.center_ is None:
                continue
            if cellID == self._target_cell:
                continue
            vector_from_target = np.subtract(cell.center_, self._config.cells_[self._target_cell].center_)
            distance_along_direction = vector_from_target.dot(self._direction)
            if distance_along_direction > 0 and distance_along_direction < 1:
                cell.s0_ = 5.6
            if distance_along_direction>1:
                cell.s0_ = 5.6
            if distance_along_direction < 0 and distance_along_direction > -1:
                cell.s0_ = 4.8
            if distance_along_direction < -1:
                cell.s0_ = 4.8
        self.write_cell_parameters()
    def set_target_cell_neighbors_s0(self):
        #   set the s0 of the target cell neighbors to:
        #   5.3 if distance > 0
        #   4.9 if distance < 0
        for neighbor_cell_id in self._target_cell_neighbors:
            neighbor_cell = self._config.cells_[neighbor_cell_id]
            vector_from_target = np.subtract(neighbor_cell.center_, self._config.cells_[self._target_cell].center_)
            neighbor_cell.s0_ = 5.6
            # if vector_from_target.dot(self._direction) > 0:
            #     neighbor_cell.s0_ = 5.6
            # else:
            #     neighbor_cell.s0_ = 4.8
        self.write_cell_parameters()
        for cellID in self._target_cell_neighbors:
            cell = self._config.cells_[cellID]
            cell.vtk_scalar_ = cell.s0_
            for polygonID in cell.polygons_:
                polygon = self._config.polygons_[polygonID]
                polygon.vtk_scalar_ = cell.s0_

        self._config.write_cell_collection_vtk(cells_array=self._target_cell_neighbors, filename="target_cell_neighbors_s0.vtk", use_scalar=True)
        

    def evaluate_distance(self):
        return     

    def initialize(self):
        if os.path.isfile("{}cellParameters.input".format(self._dir)):
            os.system("cp {}cellParameters.input {}cellParameters.init.input".format(self._dir,self._dir))
        else:
            self.write_cell_parameters("cellParameters.init.input")

        os.system("cp {}minimized.txt {}init_config.txt".format(self._dir,self._dir))
        self.set_initial_config(PeriodicTissue.from_config(self._dir,"init_config.txt".format(self._dir)))
        self._initial_config.evaluate_cell_neighbors()

    def single_iteration(self):
        self.set_s0_gradient()
        self.minimize_config()
        for _,cell in self._config.cells_.items():
            cell.vtk_scalar_ = cell.s0_
            for polygonID in cell.polygons_:
                polygon = self._config.polygons_[polygonID]
                polygon.vtk_scalar_ = cell.vtk_scalar_
        self._config.write_periodic_vtk("{:07d}.bulk.vtk".format(self._iter_counter),use_scalar=True)
        self._config.write_cell_collection_vtk(cells_array=[self._target_cell], filename="{:07d}.target.vtk".format(self._iter_counter), use_scalar=True)
        os.system("cp {}minimized.txt {}{:07d}.bulk.txt".format(self._dir,self._dir,self._iter_counter))
        os.system("cp {}cellParameters.input {}{:07d}.cellParameters.txt".format(self._dir,self._dir,self._iter_counter))
        self._iter_counter += 1
