from toolbox.periodic import PeriodicTissue
from toolbox.spheroid import Spheroid
import os
import copy
class FIREminimization:
    def __init__(self):
        self._config = None
        self._dir = None
        self._cpp_executable_dir = None

    @classmethod
    def from_sample(cls, tissue):
        sample = cls()
        sample._config = tissue
        sample._dir = tissue.config_dir_
        return sample    

    def set_cpp_executable_dir(self, cpp_executable_dir:str):
        self._cpp_executable_dir = cpp_executable_dir
    def set_config(self, tissue:PeriodicTissue):
        self._config = copy.deepcopy(tissue)
    def get_dir(self):
        return self._dir
    def minimize_config(self, FIRE_only = False):
        self.write_configuration("sample.topo")
        # tvm produces a new minimized.txt in self._dir
        # Any file of the same name must be therefore first removed or renamed.
        # Otherwise, tvm will append to the existing file.
        if os.path.isfile("{}minimized.txt".format(self._dir)):
            os.system("mv {}minimized.txt {}minimized_old.txt".format(self._dir,self._dir))
        # Run tvm executable
        if FIRE_only:
            os.system("cd {} && {}tvm FIRE_only".format(self._dir,self._cpp_executable_dir))
        else:
            os.system("cd {} && {}tvm".format(self._dir,self._cpp_executable_dir))
        if not os.path.isfile("{}minimized.txt".format(self._dir)):
            raise ValueError("FIREminimization failed: minimized.txt not found")
        self._config.set_file("minimized.txt")
        if self._config.tissueType_ == "spheroid":
            self._config.load_spheroid_from_file()
        elif self._config.tissueType_ == "periodic":
            self._config.load_periodic_tissue_from_file()
        self.load_cell_parameters()

    def write_configuration(self,filename = "sample.topo"):
        sample = self._config
        with open("{}{}".format(self._dir,filename), "w") as file:
            file.write("vertices {:d}\n".format(len(sample.vertices_)))
            for _,vertex in sample.vertices_.items():
                id = vertex.id_
                x = vertex.position_[0]
                y = vertex.position_[1]
                z = vertex.position_[2]
                file.write("{:6d} {:.14f} {:.14f} {:.14f}\n".format(id, x, y, z))
            file.write("edges {:d}\n".format(len(sample.edges_)))
            for _,edge in sample.edges_.items():
                file.write("{:d}".format(edge.id_))
                for vertexID in edge.vertices_:
                    file.write(" {:6d}".format(vertexID))
                file.write("\n")
            file.write("polygons {:d}\n".format(len(sample.polygons_)))
            for _, polygon in sample.polygons_.items():
                file.write("{:d}".format(polygon.id_))
                for edgeID in polygon.edges_:
                    file.write(" {:6d}".format(edgeID))
                file.write("\n")
            file.write("cells {:d}\n".format(len(sample.cells_)))
            for _, cell in sample.cells_.items():
                file.write("{:d}".format(cell.id_))
                for polygonID in cell.polygons_:
                    file.write(" {:6d}".format(polygonID))
                if self._config.tissueType_ == "spheroid":
                    file.write(" {:d}".format(cell.type_))
                file.write("\n")

    def load_cell_parameters(self, filename = "cellParameters.input"):
        if not os.path.isfile("{}{}".format(self._dir,filename)):
            print("{} does not exist".format(filename))
            return
        with open("{}{}".format(self._dir,filename),"r") as f:
            lines = f.readlines()
            for line in lines:
                if not len(line.split()):
                    continue
                if not len(line.split()) == 4:
                    print("Error in {}{}".format(self._dir,filename))
                    return
                tmp_id = int(line.split()[0])
                tmp_v0 = float(line.split()[1])
                tmp_s0 = float(line.split()[2])
                tmp_is_fixed = bool(int(line.split()[3]))
                cell = self._config.cells_[tmp_id]
                cell.v0_ = tmp_v0
                cell.s0_ = tmp_s0
                cell.is_fixed_ = tmp_is_fixed
                if cell.is_fixed_:
                    for polygonID in cell.polygons_:
                        self._config.polygons_[polygonID].is_fixed_ = True
