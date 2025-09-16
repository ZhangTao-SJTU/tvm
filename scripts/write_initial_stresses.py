from toolbox import stress
from toolbox.periodic import PeriodicTissue
from toolbox.training import Training
import numpy as np
import glob
import os

def calculate_stresses(dir,file):
    stresses = []
    sample = PeriodicTissue.from_config(dir,file)
    for cellID in sample.cells_:
        stresses.append(stress.calculate_max_shear_stress(sample,cellID))
    return stresses

def write_stresses(output_dir,dir_list, **kwargs):
    file = "minimized.txt"
    if "file" in kwargs:
        file = kwargs["file"]
    stresses = []
    for dir in dir_list:
        print(dir)
        stresses.extend(calculate_stresses(dir,file))
    np.savetxt("{}stresses.txt".format(output_dir),stresses)
    
def main():
    for s0 in [4.8,4.9,5.0,5.1,5.2,5.3]:
        output_dir = "init_homogeneous/7_{:.1f}/".format(s0)
        dir_list = []
        # for dir in sorted(glob.glob(output_dir+"*")):
        #     dir += "/"
        for i in range(100):
            dir = output_dir + "{:03d}/".format(i)
            if not os.path.isfile(dir+"minimized.txt"):
                continue
            dir_list.append(dir)
        if not len(dir_list):
            continue
        write_stresses(output_dir,dir_list)

if __name__ == "__main__":
    main()
