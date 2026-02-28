import random
import copy
from toolbox.minimization import FIREminimization

class Training(FIREminimization):
    def __init__(self):
        super().__init__()
        self._iter_counter = 0
        self._tolerance = 1e-8
        self._cost = None
        self._cost_values = None
        self._q_values = None
        self._initial_config = None
    @classmethod
    def from_sample(cls, tissue):
        inst = super().from_sample(tissue)
        return inst
    def set_iter_counter(self,iter_counter):
        self._iter_counter = iter_counter
    def set_tolerance(self,tolerance):
        self._tolerance = tolerance
    def set_initial_config(self,config):
        self._initial_config = config
        self._initial_config.evaluate_cell_neighbors()
    def get_iter_counter(self):
        return self._iter_counter
    def get_last_overlap(self):
        return self._q_values[-1]
    def get_tolerance(self):
        return self._tolerance
    
    def write_cell_parameters(self,filename = "cellParameters.input"):
        with open("{}{}".format(self._dir,filename),"w") as f:
            for cellID, cell in self._config.cells_.items():
                if self._config.tissueType_ == "spheroid" and not cell.type_:
                    continue
                f.write("{:d} {} {} {}\n".format(
                    cellID,
                    cell.v0_,
                    cell.s0_,
                    int(cell.is_fixed_)))

    # Edit global parameters of the configuration by changing the conf file.             
    def edit_conf(self,**kwargs):
        if self._config.tissueType_ == "spheroid":
            self.edit_conf_spheroid(**kwargs)
            return
        # Initialize variables
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
        # Load existing conf file
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

        if "s0" in kwargs:
            s0 = kwargs["s0"]
            print("s0 edited:", s0)
        if "kv" in kwargs:
            kv = kwargs["kv"]
            print("kv edited:", kv)

        if "final_time" in kwargs:
            final_time = kwargs["final_time"]
            print("final_time edited:", final_time)
        if "log" in kwargs:
            log = kwargs["log"]
            print("log edited:", log)
            
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

    def edit_conf_spheroid(self,**kwargs):
        # Initialize variables
        init_time = None
        final_time = None
        euler_time = None
        dump_vtk = None
        log = None
        s0 = None
        gamma = None
        Lth = None
        temp = None
        kv = None
        box_l = None
        box_periodic = None
        fire_equilibrium_tolerance = None
        pull_max = 3
        pull_period = 10
        fire_iter_max = None
        # Load existing conf file
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
                    gamma = line.split()[2]
                if line.startswith("Lth"):
                    Lth = line.split()[1]
                if line.startswith("T"):
                    temp = line.split()[1]
                if line.startswith("kv"):
                    kv = line.split()[1]
                if line.startswith("box"):
                    box_l = line.split()[1]
                    box_periodic = line.split()[4]
                if line.startswith("pull"):
                    fire_dt_max = line.split()[1]
                    pull_max = line.split()[2]
                    pull_period = line.split()[3]
                    fire_iter_max = line.split()[4]

        if "s0" in kwargs:
            s0 = kwargs["s0"]
            print("s0 edited:", s0)
        if "kv" in kwargs:
            kv = kwargs["kv"]
            print("kv edited:", kv)
        if "final_time" in kwargs:
            final_time = kwargs["final_time"]
            print("final_time edited:", final_time)
        if "log" in kwargs:
            log = kwargs["log"]
            print("log edited:", log)
        if "gamma" in kwargs:
            gamma = kwargs["gamma"]
            print("gamma edited:", gamma)
        if "box_l" in kwargs:
            box_l = kwargs["box_l"]
            print("box_l edited:", box_l)
        if "fire_equilibrium_tolerance" in kwargs:
            fire_equilibrium_tolerance = kwargs["fire_equilibrium_tolerance"]
            print("fire_equilibrium_tolerance edited:", fire_dt_max)
        if "fire_iter_max" in kwargs:
            fire_iter_max = kwargs["fire_iter_max"]
            print("fire_iter_max edited:", fire_iter_max)          
        with open("{}conf".format(self._dir),"w") as f:
            f.write("time {} {} {}\n".format(init_time,final_time,euler_time))
            f.write("dump vtk {}\n".format(dump_vtk))
            f.write("log {}\n".format(log))
            f.write("s0 {} {}\n".format(s0,gamma))
            f.write("Lth {}\n".format(Lth))
            f.write("T {}\n".format(temp))
            f.write("kv {}\n".format(kv))
            f.write("box {} {} {} {} {} {}\n".format(
                box_l,box_l,box_l,box_periodic,box_periodic,box_periodic))            
            f.write("pull {} {} {} {}\n".format(
                fire_equilibrium_tolerance,pull_max,pull_period,fire_iter_max))