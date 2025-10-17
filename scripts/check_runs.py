import os
import numpy as np
experiments = ["2_cells_mean_l_4/","2_cells_mean_l_5/","2_cells_mean_l_6/","4_cells_mean_l_4/","4_cells_mean_l_5/","4_cells_mean_l_6/"]
for exp in experiments:
    for i in range(100):
        dir = "/home/mameen/{}{:03d}/".format(exp,i)
        if not os.path.isfile(dir + "initial_cost.txt"):
            print(dir, "no initial cost file")
            # continue
        if not os.path.isfile(dir + "costs.txt"):
            print(dir, "run started but no costs file")
            continue
        if os.path.isfile(dir + "error.txt"):
            print(dir, "error file exists")
        with open(dir + "costs.txt", "r") as f:
            lines = f.readlines()
            print(dir, len(lines), lines[-1])
                

        # with open(dir + "costs.txt", "r") as f:
        #     lines = f.readlines()
        #     print(dir, len(lines), lines[-1])

        #     print(len(lines), lines[-1])
        # for pattern in ["patternA/", "patternB/"]:
        #     if os.path.isfile(dir + pattern + "distances.txt"):
        #         with open(dir + pattern + "distances.txt", "r") as f:
        #             lines = f.readlines()
        #             print("   ", pattern, len(lines), lines[-1])
        #     else:
        #         print("   ", pattern, "No file")