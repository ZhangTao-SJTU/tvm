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

# write stresses of all cells, useful for histogram etc...
def write_stresses(dir_list, filename, **kwargs):
    file = "minimized.txt"
    if "input_file" in kwargs:
        file = kwargs["input_file"]
    stresses = []
    for dir in dir_list:
        print(dir)
        stresses.extend(calculate_stresses(dir,file))
    np.savetxt(filename,stresses)
    
def main():
    # for s0 in [4.8,4.9,5.0,5.1,5.2,5.3]:
    for length in [4,5]:
        output_dir = "init/init_homogeneous_{}/".format(length)
        dir_list = []
        # for dir in sorted(glob.glob(output_dir+"*")):
        #     dir += "/"
        for i in range(100):
            dir = output_dir + "{:03d}/".format(i)
            if not os.path.isfile(dir+"minimized.txt"):
                raise ValueError("no input file in {}!".format(dir))
            dir_list.append(dir)
        write_stresses(output_dir,dir_list)

if __name__ == "__main__":
    main()

# from toolbox.periodic import PeriodicTissue
# from toolbox import stress
# from toolbox.training import Training
# import os
# import numpy as np
# stresses = []
# for i in range(100):
#     dir = "init/init_homogeneous_4/{:03d}/".format(i)
#     # os.system("scp shabeebameen@macbookpro.lan:/Users/shabeebameen/Projects/tvm-fire/{}sample.topo {}".format(dir,dir))
#     file = "minimized.txt"
# #     if not os.path.isfile (dir+file):
# #         continue
#     tissue = PeriodicTissue.from_config(dir,file)
#     stresses.extend([stress.calculate_max_shear_stress(tissue,cellID) for cellID in tissue.cells_])
# print(np.mean(stresses), np.std(stresses))
# np.savetxt("init_homogeneous/stresses.txt", stresses)