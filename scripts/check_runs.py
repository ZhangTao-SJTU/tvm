import os
import numpy as np
experiments = ["1_cell_increase/", "1_cell_decrease/","2_cells/","4_cells/"]
s0_vals = [4.8,4.9,5.0,5.1,5.2,5.3]
for exp in experiments:
    for s0 in s0_vals:
        s0_dir = "{}7_{:.1f}/".format(exp,s0)
        for dir in [s0_dir +"{:03d}/".format(i) for i in range(100)]:
            if not os.path.isdir(dir):
                continue
            if not os.path.isfile("{}costs.txt".format(dir)):
                continue
            with open("{}costs.txt".format(dir)) as f:
                lines = f.readlines()
            

            print(dir, len(lines),lines[-1])