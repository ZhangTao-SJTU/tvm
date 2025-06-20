import random
from toolbox.minimization import FIREminimization
class Training(FIREminimization):
    def __init__(self):
        super().__init__()
        self._iter_counter = 0
        self._learning_rate = 10
        self._lambda = 1.1
        self._tolerance = 1e-8
        self._cost = None
        self._cost_values = None
    @classmethod
    def periodic_tissue(cls, tissue):
        inst = super().periodic_tissue(tissue)
        return inst
    
    def set_iter_counter(self,iter_counter):
        self._iter_counter = iter_counter
    def set_lambda(self,lam):
        self._lambda = lam
    def set_learning_rate(self,learning_rate):
        self._learning_rate = learning_rate
    def set_tolerance(self,tolerance):
        self._tolerance = tolerance

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

    def pick_random_modified_cell(self, v0 = 1, s0 = 5.2, is_fixed = False):
        sample = self._config
        cell_found = False
        while (not cell_found):
            test_cell_id = random.choice(list(sample.cells_.keys()))
            cell = sample.cells_[test_cell_id]
            if cell.crossBoundary_:
                continue
            cross_boundary_polygons = False
            for polygonID in cell.polygons_:
                if sample.polygons_[polygonID].crossBoundary_:
                    cross_boundary_polygons = True
                    break
            if cross_boundary_polygons:
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
            self._modified_cells.append(cell.id_)
            cell.v0_ = v0
            cell.s0_ = s0
            cell.is_fixed_ = is_fixed
            if cell.is_fixed_:
                for polygonID in cell.polygons_:
                    self._config.polygons_[polygonID].is_fixed_ = True
            cell_found = True

    # Edit global parameters of the configuration by changing the conf file.             
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
            elif key == "kv":
                kv = kwargs[key]
            else:
                raise ValueError("Invalid key, or not supported yet.")
            
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
