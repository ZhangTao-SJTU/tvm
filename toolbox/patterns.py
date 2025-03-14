from toolbox.training import Training
from toolbox import stressTensor
import numpy as np
import os
import random

class Patterns(Training):
    def __init__(self):
        super().__init__()
        self._target_cells = None
        self._target_stress = None
    @classmethod
    def periodic_tissue(cls,tissue):
        inst = super().periodic_tissue(tissue)
        inst._modified_cells = list(inst._config.cells_.keys())
        return inst
    def set_target_stress(self,stress):
        self._target_stress = stress
    def set_target_cells(self,target_cells):
        self._target_cells = target_cells
        # For checking the above functionality with vtk:
        for polygonID,polygon in self._config.polygons_.items():
            polygon.vtk_scalar_ = 0
        for cellID in self._target_cells:
            cell = self._config.cells_[cellID]
            for polygonID in cell.polygons_:
                polygon = self._config.polygons_[polygonID]
                polygon.vtk_scalar_ = 1
        self._config.write_periodic_vtk(filename = "target_cells.vtk", use_scalar=True)
    def set_target_cells_spheroid(self, spheroid_radius = 0.5):
        # spheroid_radius_min = 2
        # spheroid_radius_max = 3.5
        self._target_cells = []
        for cellID,cell in self._config.cells_.items():
            if cell.crossBoundary_:
                continue
            r = np.subtract(cell.center_,self._config.sample_center_)
            r = np.linalg.norm(r)
            # if r < spheroid_radius_max and r > spheroid_radius_min:
            if r<spheroid_radius:
                self._target_cells.append(cellID)
        # For checking the above functionality with vtk:
        for polygonID,polygon in self._config.polygons_.items():
            polygon.vtk_scalar_ = 0
        for cellID in self._target_cells:
            cell = self._config.cells_[cellID]
            for polygonID in cell.polygons_:
                polygon = self._config.polygons_[polygonID]
                polygon.vtk_scalar_ = 1
        self._config.write_periodic_vtk(filename = "target_cells.vtk", use_scalar=True)
    def set_central_target_cells(self, n_cells = 1):
        r_lim = 0.9
        self._target_cells = []
        for cellID,cell in self._config.cells_.items():
            if cell.crossBoundary_:
                continue
            r = np.subtract(cell.center_,self._config.sample_center_)
            r = np.linalg.norm(r)
            if r<r_lim:
                self._target_cells.append(cellID)
            if len(self._target_cells) == n_cells:
                break
        # For checking the above functionality with vtk:
        for polygonID,polygon in self._config.polygons_.items():
            polygon.vtk_scalar_ = 0
        for cellID in self._target_cells:
            cell = self._config.cells_[cellID]
            for polygonID in cell.polygons_:
                polygon = self._config.polygons_[polygonID]
                polygon.vtk_scalar_ = 1
        self._config.write_periodic_vtk(filename = "target_cells.vtk", use_scalar=True)
    
    def set_random_target_cells(self, n_cells = 1):
        self._target_cells = []
        while len(self._target_cells)<n_cells:
            cellID = random.choice(list(self._config.cells_.keys()))
            cell = self._config.cells_[cellID]
            if cell.crossBoundary_:
                continue
            if cellID in self._target_cells:
                continue
            self._target_cells.append(cellID)
            if len(self._target_cells) == n_cells:
                break
        # For checking the above functionality with vtk:
        for polygonID,polygon in self._config.polygons_.items():
            polygon.vtk_scalar_ = 0
        for cellID in self._target_cells:
            cell = self._config.cells_[cellID]
            for polygonID in cell.polygons_:
                polygon = self._config.polygons_[polygonID]
                polygon.vtk_scalar_ = 1
        self._config.write_periodic_vtk(filename = "target_cells.vtk", use_scalar=True)
        
    def calculate_max_shear_stresses(self):
        for cellID in self._target_cells:
            cell = self._config.cells_[cellID]
            stress = stressTensor.calculate_stress_tensor_COM_center(self._config,cellID)
            egvals = np.linalg.eigvalsh(stress)
            max_shear = 0.5 * abs(egvals[-1] - egvals[0])
            cell.max_shear_stress_ = max_shear

    def evaluate_cost(self):
        multiplier = 1
        self._cost = 0
        # self._config.load_periodic_tissue_cell_properties()
        self.calculate_max_shear_stresses()
        for cellID in self._target_cells:
            cell = self._config.cells_[cellID]
            self._cost += multiplier * (cell.max_shear_stress_ - self._target_stress) ** 2

    def single_iteration(self,clamp_tol = None):
        if not len(self._target_cells):
            print("No target cells, iteration terminated.")
            return
        free_state_areas = {}
        # Starting with a minimized config with loaded cell properties...
        print("\n\n=====================================\n\n")  
        print("Starting iteration: {:d}".format(self._iter_counter))
        print("\n\n=====================================\n\n")

        print("Number of modifiable cells: {:d}".format(len(self._modified_cells)))
        print("Target cells: ", self._target_cells)

        print("\n\n-------------------------------------")
        print("Step 1: Evaluate and store the current (free state) areas of hidden (non-target) cells")
        # The stored cell areas will be used to calculate learning DOF changes
        for cellID,cell in self._config.cells_.items():
            #TO DO: fix stress calculation for boundary cells
            if cell.crossBoundary_:
                continue
            if cellID in self._target_cells:
                continue
            free_state_areas[cellID] = cell.surface_area_

        print("\n\n-------------------------------------")
        print("Step 2: CLAMPING")
        self.clamp_target_cells(clamping_tolerance = clamp_tol)

        print("\n\n-------------------------------------")
        print("Step 3: Use the clamped state areas to update all learning degrees of freedom.")

        for cellID, cell in self._config.cells_.items():
            if cell.crossBoundary_:
                continue
            if not cellID in free_state_areas:
                continue
            del_area = cell.surface_area_ - free_state_areas[cellID]
            s0_change = self._learning_rate * del_area
            self._config.cells_[cellID].s0_ -= s0_change
        
        print("\n\n-------------------------------------")
        print("Step 4: UNCLAMP target cells; write cell parameters")
        print("-------------------------------------\n\n")
        for cellID in self._target_cells:
            cell = self._config.cells_[cellID]
            cell.s0_ = self._config.s0_
        self.write_cell_parameters()
        self.minimize_config(FIRE_only = True)
        # Step 5: Logging, etc
        self.write_configuration(filename = "{}.bulk.txt".format(self._iter_counter))
        for cellID,cell in self._config.cells_.items():
            if cell.crossBoundary_:
                continue
            for polygonID in cell.polygons_:
                self._config.polygons_[polygonID].vtk_scalar_ = cell.s0_
        self._config.write_periodic_vtk(filename = "{}.bulk.vtk".format(self._iter_counter),use_scalar=True)
        os.system("cp {}cellParameters.input {}cellParameters.{}.input".format(
            self._dir,self._dir,self._iter_counter))
    

    def clamp_target_cells(self,clamping_tolerance = 1e-3, max_iters = 20):
        correction_factor = 10
        for iter in range(max_iters):
            print("Clamping iteration {}".format(iter))
            self.calculate_max_shear_stresses()
            needs_clamping = []
            for cellID in self._target_cells:
                cell = self._config.cells_[cellID]
                current_stress = cell.max_shear_stress_
                print("Cell: {}, Max shear stress: {}, s0: {}, surface area:{}".format(cellID,cell.max_shear_stress_, cell.s0_, cell.surface_area_))
                if abs(current_stress - self._target_stress) > clamping_tolerance:
                    needs_clamping.append(cellID)
            if not len(needs_clamping):
                print("Clamping successful to tolerance")
                return

            for cellID in needs_clamping:
                cell = self._config.cells_[cellID]
                current_stress = cell.max_shear_stress_
                # The target stress should be an overshoot/undershoot of self._target_stress
                target_stress = self._target_stress * (1 + correction_factor * (self._target_stress - current_stress))
                print("Cell {}, Target stress {} ".format(cellID, target_stress))
                # Desired change in stress:
                del_stress = self._lambda * (target_stress - current_stress)
                upper_bound = 5.6
                lower_bound = 4.8
                stress_change = 0
                # Binary search for s0 that produces the right stress change
                while abs(stress_change - del_stress) > self._tolerance:
                    cell.s0_ = (upper_bound + lower_bound) / 2
                    stress = stressTensor.calculate_stress_tensor_COM_center(self._config,cellID)
                    egvals = np.linalg.eigvalsh(stress)
                    max_shear = 0.5 * abs(egvals[-1] - egvals[0])
                    stress_change = max_shear - current_stress
                    if stress_change < del_stress:
                        upper_bound = cell.s0_
                    else:
                        lower_bound = cell.s0_
                cell.s0_ = (upper_bound + lower_bound) / 2
            self.write_cell_parameters()
            self.minimize_config(FIRE_only = True)
            
    def run(self):
        self._cost_values = []
        self.evaluate_cost()
        self._cost_values.append(self._cost)
        print("Initial Cost: {:.2e}".format(self._cost))
        while self._cost > self._tolerance:
            if self._cost<1e-6:
                self.single_iteration(clamp_tol = 1e-5)
            else:
                self.single_iteration(clamp_tol = 1e-3)
            self.evaluate_cost()
            self._cost_values.append(self._cost)
            print("Iteration: {:d}, Cost: {:.2e}".format(self._iter_counter,self._cost))
            self._iter_counter += 1

        print("Optimization finished at iteration: {:d}".format(self._iter_counter-1))
        print("Final cost: {:.2e}".format(self._cost))

