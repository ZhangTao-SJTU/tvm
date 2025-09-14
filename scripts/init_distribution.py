# calculate_stresses
code = """
from toolbox import stress
from toolbox.periodic import PeriodicTissue
from toolbox.training import Training
import numpy as np
stresses = []
main_dir = "init_homogeneous/"
for i in range(50):
    dir = main_dir+"7_{}/".format(i)
    print("Processing {}".format(dir))
    sample = PeriodicTissue.from_config(dir,"minimized.txt")
    # training_instance = Training.periodic_tissue(sample)
    for cellID,cell in sample.cells_.items():
        stresses.append(stress.calculate_max_shear_stress(sample,cellID))
np.savetxt(main_dir+"stresses.txt",stresses)
"""
exec(code)


code = """
import os
for i in range(50,100):
    os.system("rm -rf init_homogeneous/7_{}/".format(i))
"""
exec(code)