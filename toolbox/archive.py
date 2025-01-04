from toolbox import tissueSample
from toolbox import functions
import copy
import os
import random
import numpy as np
from toolbox.minimization import FIREminimization

class Training:
    def __init__(self):
        self._config = None
        self._dir = None
        self._iter_counter = 0
        self._fixed_cells = None
        self._direction = None
        self._stepsize = None
        self._expansion_factor = None
        self._separation = None
        self._minimum_separation = None
        self._maximum_separation = None

    @classmethod
    def from_config(cls, config_dir, input_filename):
        sample = cls()
        sample._config = tissueSample.Sample.periodic_tissue(
            config_dir = config_dir,
            input_filename = input_filename)
        sample._dir = config_dir
        sample._maximum_separation = sample._config.boxSize_/2
        return sample
    def set_stepsize(self,stepsize):
        self._stepsize = stepsize
    def set_expansion_factor(self, expansion_factor):
        self._expansion_factor = expansion_factor
    def set_minimum_separation(self,minimum_separation):
        self._minimum_separation = minimum_separation
    def set_maximum_separation(self,maximum_separation):
        self._maximum_separation = maximum_separation
    def set_iter_counter(self,iter_counter):
        self._iter_counter = iter_counter

    # This function write the current state of self._config to
    # self._dir/sample.topo
    # The topology is then minimized using the tvm program.
    # Finally, the minimized topology is read back into self._config
    def minimize_config(self):
        self.write_configuration("sample.topo")
        # tvm produces a new minimized.txt in self._dir
        # Any file of the same name must be therefore first removed.
        # Otherwise, tvm will append to the existing file.
        if os.path.isfile("{}minimized.txt".format(self._dir)):
            os.remove("{}minimized.txt".format(self._dir))
        os.system("cd {} && ../build/tvm".format(self._dir))
        self._config = tissueSample.Sample.periodic_tissue(
            config_dir = self._dir,
            input_filename = "minimized.txt")
        self.load_fixed_cells()

    def write_configuration(self,filename = "sample.topo"):
        sample = self._config
        with open("{}{}".format(self._dir,filename), "w") as file:
            file.write("vertices {:d}\n".format(len(sample.vertices_)))
            for key,vertex in sample.vertices_.items():
                id = vertex.id_
                x = vertex.position_[0]
                y = vertex.position_[1]
                z = vertex.position_[2]
                file.write("{:6d} {:.14f} {:.14f} {:.14f}\n".format(id, x, y, z))
            file.write("edges {:d}\n".format(len(sample.edges_)))
            for key,edge in sample.edges_.items():
                file.write("{:d}".format(edge.id_))
                for vertexID in edge.vertices_:
                    file.write(" {:6d}".format(vertexID))
                file.write("\n")
            file.write("polygons {:d}\n".format(len(sample.polygons_)))
            for key, polygon in sample.polygons_.items():
                file.write("{:d}".format(polygon.id_))
                for edgeID in polygon.edges_:
                    file.write(" {:6d}".format(edgeID))
                file.write("\n")
            file.write("cells {:d}\n".format(len(sample.cells_)))
            for key, cell in sample.cells_.items():
                file.write("{:d}".format(cell.id_))
                for polygonID in cell.polygons_:
                    file.write(" {:6d}".format(polygonID))
                file.write("\n")
    
    def write_periodic_vtk(self,filename):
        sample = self._config
        vertices = []
        total_polygons = 0
        total_polygon_data_points = 0
        for polygonID,polygon in sample.polygons_.items():
            if polygon.crossBoundary_:
                continue
            total_polygons += 1
            total_polygon_data_points += len(polygon.vertices_) + 1
            for vertex in polygon.vertices_:
                vertices.append(vertex)
        vertices = list(set(vertices))
        v_map = functions.mapmaker(vertices)
        with open("{}{}".format(self._dir,filename),"w") as f:
            f.write("# vtk DataFile Version 2.0\n")
            f.write("polydata\n")
            f.write("ASCII\n")
            f.write("DATASET POLYDATA\n")
            f.write("POINTS {} double\n".format(len(vertices)))
            for vertexID in vertices:
                for i in range(3):
                    f.write("{} ".format(sample.vertices_[vertexID].position_[i]))
                f.write("\n")
            f.write("POLYGONS {} {}\n".format(total_polygons,total_polygon_data_points))
            for polygonID,polygon in sample.polygons_.items():
                if polygon.crossBoundary_:
                    continue
                f.write("{} ".format(len(polygon.vertices_)))
                for vertexID in polygon.vertices_:
                    f.write("{} ".format(v_map[vertexID]))
                f.write("\n")
            f.write("CELL_DATA {}\n".format(total_polygons))
            f.write("SCALARS is_fixed double\n")
            f.write("LOOKUP_TABLE default\n")
            for polygonID,polygon in sample.polygons_.items():
                if polygon.crossBoundary_:
                    continue
                if polygon.is_fixed_:
                    f.write("1\n")
                else:
                    f.write("0\n")
            # f.write("SCALARS shape_index double\n")
            # f.write("LOOKUP_TABLE default\n")
            # for polygonID,polygon in sample.polygons_.items():
            #     if polygon.crossBoundary_:
            #         continue
            #     for cellID,cell in sample.cells_.items():
            #         # if cell.crossBoundary_:
            #         #     continue
            #         if polygonID in cell.polygons_:
            #             if cell.shape_index_ is None:
            #                 f.write("0\n")
            #                 break
            #             else:
            #                 f.write("{:.14f}\n".format(cell.shape_index_))
            #                 break
            # f.write("SCALARS volume double\n")
            # f.write("LOOKUP_TABLE default\n")
            # for polygonID,polygon in sample.polygons_.items():
            #     if polygon.crossBoundary_:
            #         continue
            #     for cellID,cell in sample.cells_.items():
            #         # if cell.crossBoundary_:
            #         #     continue
            #         if polygonID in cell.polygons_:
            #             if cell.shape_index_ is None:
            #                 f.write("0\n")
            #                 break
            #             else:
            #                 f.write("{:.14f}\n".format(cell.volume_))
            #                 break
    
    def pick_fixed_cell(self):
        self._fixed_cells = []
        for cellID,cell in self._config.cells_.items():
            # cell = sample.cells_[random.choice(list(sample.cells_.keys()))]
            if cell.crossBoundary_:
                continue
            
            good_cell = True
            for i in range(3):
                if cell.center_[i] < 0.4*self._config.boxSize_:
                    good_cell = False
                if cell.center_[i] > 0.6*self._config.boxSize_:
                    good_cell = False
            if good_cell:
                self._fixed_cells.append(cell.id_)
                break        
        # write the fixed cell IDs to a file
        with open("{}fixed.topo".format(self._dir),"w") as f:
            for cellID in self._fixed_cells:
                f.write("{:d}\n".format(cellID))
        self.load_fixed_cells()

    # pick two fixed cells at random and write fixed.topo
    def pick_fixed_pair(self,n_fixed_cells = 2):
        sample = self._config
        self._fixed_cells = []
        while len(self._fixed_cells) < n_fixed_cells:
            cell = sample.cells_[random.choice(list(sample.cells_.keys()))]
            if cell.crossBoundary_:
                continue
            if cell.is_fixed_:
                continue
            self._fixed_cells.append(cell.id_)
            if len(self._fixed_cells) == 2:
                self.calculate_separation()
                if (self._separation < self._maximum_separation 
                    or self._separation > np.sqrt(3)*self._maximum_separation):
                    self._fixed_cells = []

        # write the fixed cell IDs to a file
        with open("{}fixed.topo".format(self._dir),"w") as f:
            for cellID in self._fixed_cells:
                f.write("{:d}\n".format(cellID))
        # fix the cells from self._fixed_cells in self._config
        self.load_fixed_cells()

    def shrink_cell(self,cellID,factor,inwards = True):
        cell = self._config.cells_[cellID]
        for vertexID in cell.vertices_:
            vertex = self._config.vertices_[vertexID]
            direction = vertex.position_ - cell.center_
            direction /= np.linalg.norm(direction)
            if inwards:
                vertex.position_ -= factor * direction
            else:
                vertex.position_ += factor * direction
    
    # Calculate the unit vector between two fixed cells
    def calculate_direction(self):
        if self._fixed_cells is None:
            print("fixed_cells is not set")
            return
        if not (len(self._fixed_cells) == 2):
            print("There are not exactly two fixed cells")
            return
        cell_1 = self._config.cells_[self._fixed_cells[0]]
        cell_2 = self._config.cells_[self._fixed_cells[1]]
        self._direction = cell_2.center_ - cell_1.center_
        self._direction /= np.linalg.norm(self._direction)

    # Radial separation between two vectors
    def calculate_separation(self):
        if self._fixed_cells is None:
            print("nothing to calculate. Need to initialize self._fixed_cells")
            return
        cell_1 = self._config.cells_[self._fixed_cells[0]]
        cell_2 = self._config.cells_[self._fixed_cells[1]]
        self._separation = np.linalg.norm(cell_2.center_ - cell_1.center_)
    
    def shift_cell(self,cellID,vector):
        cell = self._config.cells_[cellID]
        if not cell.is_fixed_:
            print("WARNING: shift_cells() should only be used on the fixed cells")
        for vertexID in cell.vertices_:
            self._config.vertices_[vertexID].position_ += vector

    # This function loads the fixed cell IDs into self._config
    # If the fixed cell IDs have not already been loaded, 
    # (i.e if self._fixed_cells is None)
    # the fixed.topo file and loads the fixed cell IDs into self._fixed_cells
    def load_fixed_cells(self):
        if self._fixed_cells is None:
            if not os.path.isfile("{}fixed.topo".format(self._dir)):
                print("fixed.topo does not exist")
                return
            self._fixed_cells = []
            with open("{}fixed.topo".format(self._dir),"r") as f:
                for line in f:
                    self._fixed_cells.append(int(line))

        #fix the cells from self._fixed_cells in self._config
        for cellID in self._fixed_cells:
            cell = self._config.cells_[cellID]
            cell.is_fixed_ = True
            for polygonID in cell.polygons_:
                self._config.polygons_[polygonID].is_fixed_ = True
        # self.set_direction()
    
    def linear_displace_fixed_pair(self, stepsize, inwards = True):
        if self._fixed_cells is None:
            print("nothing to drive. Need to initialize self._fixed_cells")
            return
        self.calculate_direction()
        if inwards:
            self.shift_cell(self._fixed_cells[0], stepsize * self._direction)
            self.shift_cell(self._fixed_cells[1], -1 * stepsize * self._direction)
        else:
            self.shift_cell(self._fixed_cells[0], -1 * stepsize * self._direction)
            self.shift_cell(self._fixed_cells[1], stepsize * self._direction)

    # Packages the linear displacement of fixed pair, 
    # followed by minimization, writing and iteration increment
    def pair_drive_iteration(self, stepsize, inwards = True):
        self.linear_displace_fixed_pair(stepsize = stepsize, inwards = inwards)
        self.minimize_config()
        self.calculate_separation()
        print("Current separation: {}".format(self._separation))
        self.write_configuration("{:07d}.sample.topo".format(self._iter_counter))
        self.write_periodic_vtk("{:07d}.sample.vtk".format(self._iter_counter))
        self._iter_counter += 1
        
    def pair_drive(self,n_oscillations = 1):
        self.write_configuration("initial.topo")
        self.write_periodic_vtk("initial.vtk")
        self.load_fixed_cells()
        self.calculate_separation()
        
        # Initialize to the maximum separation. This is the 0th iteration.
        if self._separation > self._maximum_separation:
            stepsize = (self._separation - self._maximum_separation)/2
            self.pair_drive_iteration(stepsize = stepsize, inwards = True)
        for _ in range(n_oscillations):
        #Step 1: Drive configuration inwards
            while (self._separation >= self._minimum_separation + 2 * self._stepsize):
                self.pair_drive_iteration(stepsize = self._stepsize, inwards = True)

            # Step 1 Continued: Final step to reach the minimum separation
            if self._separation > self._minimum_separation:
                stepsize = (self._separation - self._minimum_separation)/2
                self.pair_drive_iteration(stepsize = stepsize, inwards = True)

            #Step 2: Drive configuration outwards
            while (self._separation + 2 * self._stepsize <= self._maximum_separation):
                self.pair_drive_iteration(stepsize = self._stepsize, inwards = False)

            # Step 2 Continued: Final step to reach the maximum separation
            if self._separation < self._maximum_separation:
                stepsize = (self._maximum_separation - self._separation)/2
                self.pair_drive_iteration(stepsize = stepsize, inwards = False)

    def pulse_drive(self, n_iterations):
        cellID = self._fixed_cells[0]
        self.minimize_config()
        self.write_configuration("{:07d}.sample.topo".format(self._iter_counter))
        self.write_periodic_vtk("{:07d}.sample.vtk".format(self._iter_counter))
        self._iter_counter += 1
        # # Stage 1: Shrink Cell
        # for i in range(n_iterations):
        #     self.shrink_cell(cellID = cellID, factor = self._expansion_factor, inwards = True)
        #     self.minimize_config()
        #     self.write_configuration("{:07d}.sample.topo".format(self._iter_counter))
        #     self.write_periodic_vtk("{:07d}.sample.vtk".format(self._iter_counter))
        #     self._iter_counter += 1
        # # Stage 2: Expand Cell
        # for i in range(n_iterations):
        #     self.shrink_cell(cellID = cellID, factor = self._expansion_factor, inwards = False)
        #     self.minimize_config()
        #     self.write_configuration("{:07d}.sample.topo".format(self._iter_counter))
        #     self.write_periodic_vtk("{:07d}.sample.vtk".format(self._iter_counter))
        #     self._iter_counter += 1

        
        for _ in range(n_iterations):
            # One cycle of shrink and expand
            for cellID in self._fixed_cells:
                self.shrink_cell(cellID = cellID, factor = self._expansion_factor, inwards = True)
            self.minimize_config()
            self.write_configuration("{:07d}.sample.topo".format(self._iter_counter))
            self.write_periodic_vtk("{:07d}.sample.vtk".format(self._iter_counter))
            self._iter_counter += 1
            for cellID in self._fixed_cells:
                self.shrink_cell(cellID = cellID, factor = self._expansion_factor, inwards = False)
            self.minimize_config()
            self.write_configuration("{:07d}.sample.topo".format(self._iter_counter))
            self.write_periodic_vtk("{:07d}.sample.vtk".format(self._iter_counter))
            self._iter_counter += 1

