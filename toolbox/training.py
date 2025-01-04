from toolbox import tissueSample
from toolbox import functions
import copy
import os
import random
import numpy as np
from toolbox.minimization import FIREminimization

class Training(FIREminimization):
    def __init__(self):
        super().__init__()
        #pulse drive parameters
        self._iter_counter = 0

    @classmethod
    def from_config(cls, config_dir, input_filename):
        inst = super().from_config(config_dir,input_filename)
        return inst
    
    def set_iter_counter(self,iter_counter):
        self._iter_counter = iter_counter
    
    # Writes the file cellParameters.input,
    # with current state of cells from self._modified_cells
    # In other words; parameters to be written from self._config.cells_
    def write_cell_parameters(self):
        if not self._modified_cells:
            print("No modified cells")
            return
        
        # write the modified cell IDs to a file
        with open("{}cellParameters.input".format(self._dir),"w") as f:
            for cellID in self._modified_cells:
                cell = self._config.cells_[cellID]
                f.write("{:d} {} {} {}\n".format(
                    cell.id_,
                    cell.v0_,
                    cell.s0_,
                    int(cell.is_fixed_)))

    # pick a cell; cellID appended to self._modified_cells
    # Criteria:
    # 1. No cross boundary cells
    # 2. No moidified cellls should share polygons
    # 3. Should not be too close to the boundary (for visualization).
    def pick_random_modified_cell(self, v0 = 1, s0 = 5.2, is_fixed = False):
        sample = self._config
        cell_found = False
        while (not cell_found):
            test_cell_id = random.choice(list(sample.cells_.keys()))
            cell = sample.cells_[test_cell_id]
            
            if cell.crossBoundary_:
                continue
            for polygonID in cell.polygons_:
                if sample.polygons_[polygonID].crossBoundary_:
                    continue
            if self._modified_cells:
                if cell.id_ in self._modified_cells:
                    continue
                shared_polygons = False
                for cellID in self._modified_cells:
                    for polygonID in sample.cells_[cellID].polygons_:
                        if polygonID in cell.polygons_:
                            shared_polygons = True
                            break
                    if shared_polygons:
                        break
                if shared_polygons:
                    continue
            boundary = False
            for coordinate in cell.center_:
                if coordinate < 1:
                    boundary = True
                    break
                if coordinate > sample.boxSize_ - 1:
                    boundary = True
                    break
            if boundary:
                continue
            # print(cell.center_)
            self._modified_cells.append(cell.id_)
            cell.v0_ = v0
            cell.s0_ = s0
            cell.is_fixed_ = is_fixed
            if cell.is_fixed_:
                for polygonID in cell.polygons_:
                    self._config.polygons_[polygonID].is_fixed_ = True
            cell_found = True

    # Edit global parameters of the configuration               
    def edit_conf(self,**kwargs):
        init_time = None
        final_time = None
        euler_time = None
        dump_vtk = None
        log = None
        s0 = None
        Lth = None
        temp = None
        kv = None
        box_l = None
        box_periodic = None

        with open("{}conf".format(self._dir),"r") as f:
            lines = f.readlines()
            for line in lines:
                if not len(line.split()):
                    continue
                if line.startswith("time"):
                    init_time = line.split()[1]
                    final_time = line.split()[2]
                    euler_time = line.split()[3]
                if line.startswith("dump vtk"):
                    dump_vtk = line.split()[2]
                if line.startswith("log"):
                    log = line.split()[1]
                if line.startswith("s0"):
                    s0 = line.split()[1]
                if line.startswith("Lth"):
                    Lth = line.split()[1]
                if line.startswith("T"):
                    temp = line.split()[1]
                if line.startswith("kv"):
                    kv = line.split()[1]
                if line.startswith("box"):
                    box_l = line.split()[1]
                    box_periodic = line.split()[4]

        for key in kwargs:
            if key == "s0":
                s0 = kwargs[key]
            if key == "kv":
                kv = kwargs[key]

        with open("{}conf".format(self._dir),"w") as f:
            f.write("time {} {} {}\n".format(init_time,final_time,euler_time))
            f.write("dump vtk {}\n".format(dump_vtk))
            f.write("log {}\n".format(log))
            f.write("s0 {}\n".format(s0))
            f.write("Lth {}\n".format(Lth))
            f.write("T {}\n".format(temp))
            f.write("kv {}\n".format(kv))
            f.write("box {} {} {} {} {} {}\n".format(
                box_l,box_l,box_l,box_periodic,box_periodic,box_periodic))            
class pulseDrive(Training):
    def __init__(self):
        super().__init__()

    @classmethod
    def from_config(cls, config_dir, input_filename):
        inst = super().from_config(config_dir,input_filename)
        return inst
    
    # A single iteration, consisting of:
    # 1. writing current state of (driven) cell parameters
    # 2. Minimizing configuration (subject to current state of driven cells)
    # 3. Writing resulting configuration, bulk vtk and input cells vtk
    #   (labelled by iteration number)
    # 4. Incrementing the iteration counter (this is mainly helpful for resuming runs)
    def single_iteration(self):
        self.write_cell_parameters()
        self.minimize_config()
        self.write_configuration(filename = "{}.bulk.txt".format(self._iter_counter))
        self._config.write_periodic_vtk(filename = "{}.bulk.vtk".format(self._iter_counter))
        self._config.write_cell_collection_vtk(
            cells_array = self._modified_cells,
            filename = "{}.input.vtk".format(self._iter_counter))
        self._iter_counter += 1
    
    # Pulse the s0 of driven cells between min_s0 and max_s0
    def run_s0_pulsing(self, min_s0 = 5, max_s0 = 5.6, iterations = 10):
        for iter in range(iterations):
            for cellID in self._modified_cells:
                cell = self._config.cells_[cellID]
                if iter%2:
                    cell.s0_ = max_s0
                else:
                    cell.s0_ = min_s0
            self.single_iteration()

class stressPatterns(Training):
    def __init__(self):
        super().__init__()

    @classmethod
    def from_config(cls, config_dir, input_filename):
        inst = super().from_config(config_dir,input_filename)
        return inst
    


