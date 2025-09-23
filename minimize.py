from toolbox.periodic import PeriodicTissue
from toolbox.training import Training
import os
for i in range(0,10):
    dir = "init_homogeneous/{:03d}/".format(i)
    file = "sample.topo"
    tissue = PeriodicTissue.from_config(dir,file)
    minimization_instance = Training.periodic_tissue(tissue)
    minimization_instance.set_cpp_executable_dir("/home/shabeeb/Projects/tvm-fire/build/")
    minimization_instance.edit_conf(log = 500, final_time = 50000)
    minimization_instance.minimize_config()
    minimization_instance.edit_conf(log = 100)
    os.system("rm {}FE.txt".format(dir))