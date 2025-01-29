from toolbox.training import Training
from toolbox import stressTensor
import numpy as np
import os

class Patterns(Training):
    def __init__(self):
        super().__init__()
        self._target_cells = None
        self._target_stress = None
        self._cost_values = None
    @classmethod
    def from_config(cls, config_dir, input_filename):
        inst = super().from_config(config_dir,input_filename)
        inst._config.load_periodic_tissue_cell_properties()
        return inst
    def set_target_stress(self,stress):
        self._target_stress = stress
    
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
    
    def calculate_max_shear_stresses(self):
        for cellID, cell in self._config.cells_.items():
            if cell.crossBoundary_:
                continue
            stress = stressTensor.calculate_stress_tensor_COM_center(self._config,cellID)
            egvals = np.linalg.eigvalsh(stress)
            max_shear = 0.5 * abs(egvals[-1] - egvals[0])
            cell.max_shear_stress_ = max_shear

    def evaluate_cost(self):
        multiplier = 1
        self._cost = 0
        self._config.load_periodic_tissue_cell_properties()
        self.calculate_max_shear_stresses()
        for cellID in self._target_cells:
            cell = self._config.cells_[cellID]
            self._cost += multiplier * (cell.max_shear_stress_ - self._target_stress) ** 2

    def single_iteration(self):
        if not len(self._target_cells):
            print("No target cells, iteration terminated.")
            return
        self._modified_cells = []
        free_state_areas = {}
        clamped_state_areas = {}
        # Starting with a minimized config with loaded cell properties...
        print("\n\n=====================================\n\n")  
        print("Starting iteration: {:d}".format(self._iter_counter))
        print("\n\n=====================================\n\n")  
        
        print("Step 0: Identify all non-boundary cells as self._modified_cells. Load cell parameters.")
        # Note: Actually all cells are in fact modifiable; we are just ignoring any updates on boundary cells
        # Hopefully by making target cells closer to box center, their required updates are ignorable.
        for cellID, cell in self._config.cells_.items():
            if cell.crossBoundary_:
                continue
            self._modified_cells.append(cellID)

        print("Number of modifiable cells: {:d}".format(len(self._modified_cells)))
        print("Target cells: ", self._target_cells)

        # if cellParameters.input exists, load the parameters at this point
        if os.path.isfile("{}cellParameters.input".format(self._dir)):
            with open("{}cellParameters.input".format(self._dir),"r") as f:
                lines = f.readlines()
                for line in lines:
                    if not len(line.split()):
                        continue
                    if not len(line.split()) == 4:
                        print("Error in {}cellParameters.input".format(self._dir))
                        return
                    tmp_id = int(line.split()[0])
                    tmp_v0 = float(line.split()[1])
                    tmp_s0 = float(line.split()[2])
                    tmp_is_fixed = bool(int(line.split()[3]))
                    if not tmp_id in self._modified_cells:
                        continue
                    cell = self._config.cells_[tmp_id]
                    cell.v0_ = tmp_v0
                    cell.s0_ = tmp_s0
                    cell.is_fixed_ = tmp_is_fixed
                    if cell.is_fixed_:
                        for polygonID in cell.polygons_:
                            self._config.polygons_[polygonID].is_fixed_ = True
        print("\n\n-------------------------------------")
        print("Step 1: Evaluate and store the current (free state) stresses of all cells")
        # (Again, meaning all non-boundary cells)
        # The stored cell areas will be used to calculate learning DOF changes
        
        self.calculate_max_shear_stresses()
        for cellID,cell in self._config.cells_.items():
            if cell.crossBoundary_:
                continue
            free_state_areas[cellID] = cell.surface_area_

        print("\n\n-------------------------------------")
        print("Step 2: CLAMPING and minimizing")
        # Evaluate the stress difference for target cells
        # Solve for clamping s0 value to reproduce this stress difference (upto tolerence)
        # Conclude clamping by minimizing configuration with these s0 values

        for cellID in self._target_cells:
            cell = self._config.cells_[cellID]
            free_state_stress = cell.max_shear_stress_
            target_del_stress = self._lambda * (self._target_stress - cell.max_shear_stress_)

            upper_bound = 5.6
            lower_bound = 4.8

            stress_change = 0
            while abs(stress_change - target_del_stress) > self._tolerance:
                cell.s0_ = (upper_bound + lower_bound) / 2
                stress = stressTensor.calculate_stress_tensor_COM_center(self._config,cellID)
                egvals = np.linalg.eigvalsh(stress)
                max_shear = 0.5 * abs(egvals[-1] - egvals[0])
                stress_change = max_shear - free_state_stress
                # print(cell.s0_, stress_change)
                if stress_change < target_del_stress:
                    upper_bound = cell.s0_
                else:
                    lower_bound = cell.s0_
            cell.s0_ = (upper_bound + lower_bound) / 2
        self.write_cell_parameters()
        self.minimize_config(FIRE_only=True)
        self._config.load_periodic_tissue_cell_properties()
        
        print("\n\n-------------------------------------")
        print("Step 3: Find all clamped state areas. Hence update all learning degrees of freedom.")
        
        for cellID,cell in self._config.cells_.items():
            if cell.crossBoundary_:
                continue
            clamped_state_areas[cellID] = cell.surface_area_
        # note, here we are changing target cell s0 as well. However this will be 
        # corrected in the following unclamping step
        for cellID, clamped_area in clamped_state_areas.items():
            if not cellID in free_state_areas:
                continue
            del_area = clamped_area - free_state_areas[cellID]
            s0_change = self._learning_rate * del_area
            self._config.cells_[cellID].s0_ += s0_change
        
        print("\n\n-------------------------------------")
        print("Step 4: UNCLAMP target cells, minimize configuration")
        print("-------------------------------------\n\n")
        for cellID in self._target_cells:
            cell = self._config.cells_[cellID]
            cell.s0_ = self._config.s0_
        self.write_cell_parameters()
        # self.minimize_config()
        # self._config.load_periodic_tissue_cell_properties()
        self.calculate_max_shear_stresses()
        # Step 5: Logging, etc
        self.write_configuration(filename = "{}.bulk.txt".format(self._iter_counter))
        for cellID,cell in self._config.cells_.items():
            if cell.crossBoundary_:
                continue
            for polygonID in cell.polygons_:
                self._config.polygons_[polygonID].vtk_scalar_ = cell.max_shear_stress_
        self._config.write_periodic_vtk(filename = "{}.bulk.vtk".format(self._iter_counter),use_scalar=True)
        os.system("cp {}cellParameters.input {}cellParameters.{}.input".format(
            self._dir,self._dir,self._iter_counter))
        
    
    def run(self):
        self.minimize_config()
        self._cost_values = []
        self.evaluate_cost()
        self._cost_values.append(self._cost)
        print("Initial Cost: {:.2e}".format(self._cost))
        while self._cost > self._tolerance:
            self.single_iteration()
            self.evaluate_cost()
            self._cost_values.append(self._cost)
            print("Iteration: {:d}, Cost: {:.2e}".format(self._iter_counter,self._cost))
            self._iter_counter += 1

        print("Optimization finished at iteration: {:d}".format(self._iter_counter-1))
        print("Final cost: {:.2e}".format(self._cost))


        # for i in range(3):
        #     self.single_iteration()
        #     self.evaluate_cost()
        #     print("Iteration: {:d}, Cost: {:.4f}".format(self._iter_counter,self._cost))

