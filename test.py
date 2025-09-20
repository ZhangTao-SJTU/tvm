from toolbox.periodic import PeriodicTissue
from toolbox.minimization import FIREminimization
from toolbox.training import Training
from toolbox.pattern_trainer import train_random_cells
from toolbox.periodic import PeriodicTissue
from toolbox import stress
import numpy as np

dir = "tests/6_5.0_4_cells/"
sample = PeriodicTissue.from_config(dir,"sample.topo")

min_instance = Training.periodic_tissue(sample)
min_instance.set_cpp_executable_dir("/Users/shabeebameen/Projects/tvm-fire/build/")
min_instance.minimize_config()
sample =  PeriodicTissue.from_config(dir,"minimized.txt")
stresses = [stress.calculate_max_shear_stress(sample,cellID) for cellID in sample.cells_]
train_random_cells(dir,average_cells_only = True, target_stress = np.mean(stresses),n_cells = 4,cpp_executable_dir = "/Users/shabeebameen/Projects/tvm-fire/build/",max_iters = 2000)