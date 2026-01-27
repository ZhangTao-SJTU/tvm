from toolbox.periodic import PeriodicTissue
from toolbox.spheroid import Spheroid
from toolbox.training import Training
from toolbox.overlap import calculate_Q
from toolbox import stress
import numpy as np
import os
import random
import pandas as pd

class Gradient(Training):
    def __init__(self):
        super().__init__()
        self._s0_min = 4.8
        self._s0_max = 5.2
        self._direction = np.array([0.0,0.0,1.0])

    @classmethod
    def spheroid(cls,tissue):
        inst = super().spheroid(tissue)
        inst._modified_cells = []
        for cellID,cell in inst._config.cells_.items():
            if cell.type_:
                inst._modified_cells.append(cellID)
        return inst
    
    def set_s0_range(self,s0_min,s0_max):
        self._s0_min = s0_min
        self._s0_max = s0_max
    
    def assign_linear_gradient_s0(self):
        m = (self._s0_max-self._s0_min)/2
        c = (self._s0_max+self._s0_min)/2
        for cellID,cell in self._config.cells_.items():
            if not cell.type_:
                continue
            r_hat = np.subtract(cell.center_,self._config.spheroid_center_)
            r_hat/=np.linalg.norm(r_hat)
            cell.s0_ = m*np.dot(r_hat,self._direction) + c
        self.write_cell_parameters()
        