from toolbox.periodic import PeriodicTissue
from toolbox.training import Training
from toolbox.overlap import calculate_Q
from toolbox import stress
import numpy as np
import os
import random
import pandas as pd

class cell_driver(Training):
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
    def set_target_cell(self,target_cell):
        self._target_cell = target_cell
        # For checking the above functionality with vtk:
        for polygonID,polygon in self._config.polygons_.items():
            polygon.vtk_scalar_ = 0
        for cellID in self._target_cell_to_stress:
            cell = self._config.cells_[cellID]
            for polygonID in cell.polygons_:
                polygon = self._config.polygons_[polygonID]
                polygon.vtk_scalar_ = 1
        self._config.write_periodic_vtk(filename = "target_cells.vtk", use_scalar=True)
    def set_target_position(self, position):
        self._target_position = position
    def set_target_cell_neighbors(self):
        self._config.evaluate_cell_neighbors()
        self._target_cell_neighbors = self._config.cell_neighbors_[self._target_cell]

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
        np.savetxt("{}initial_cost.txt".format(self._dir), [self.evaluate_cost()], fmt='%.2e')

        initial_stresses = {}
        self.calculate_max_shear_stresses()
        for cellID in self._target_cell_to_stress:
            cell = self._config.cells_[cellID]
            initial_stresses[cellID] = cell.max_shear_stress_
        df = pd.DataFrame(list(initial_stresses.items()), columns=['CellID', 'Current'])        
        # df = pd.DataFrame(self._target_cell_to_stress.items(), columns=['cellID', 'target_stress'])
        df.to_csv("{}initial_stress.csv".format(self._dir), index=False)
        # print("Initial Cost: {:.2e}".format(cost))
        self._cost_values = []
        self._q_values = []


