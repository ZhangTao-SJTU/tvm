import pandas as pd
import numpy as np
import math
import os
import sys

def write_SD_s0_in_dir_spheroid(dir):
    if not os.path.isfile(dir + "costs.txt"):
        print("costs.txt not found in {}".format(dir))
        return
    costs = np.genfromtxt(dir + "costs.txt")

    if costs.shape == ():
        return
    if costs.shape == (0,):
        return
    if math.isnan(costs[-1]):
        return
    iter_to_SD_s0 = {}
    for i in range(50000):
        file = dir + f"files/{i:07d}.cellParameters.input"
        if not os.path.isfile(file):
            continue
        s0_vals = []
        df = pd.read_csv(file,sep = " ", header=None)
        spheroid_cells = np.loadtxt(dir+"spheroid_cells.txt", dtype=int)
        for _, row in df.iterrows():
            s0_vals.append(row[2])
        iter_to_SD_s0[i] = np.std(s0_vals)
    for i in range(50000):
        file = dir + f"{i:07d}.cellParameters.input"
        if not os.path.isfile(file):
            continue
        s0_vals = []
        df = pd.read_csv(file,sep = " ", header=None)
        spheroid_cells = np.loadtxt(dir+"spheroid_cells.txt", dtype=int)
        for _, row in df.iterrows():
            if row[0] in spheroid_cells:
                s0_vals.append(row[2])
        iter_to_SD_s0[i] = np.std(s0_vals)
    pd.DataFrame({"iter": list(iter_to_SD_s0.keys()), "SD_s0": list(iter_to_SD_s0.values())}).to_csv(f"{dir}SD_s0.csv", index=False)

def write_SD_s0_in_dir(dir):
    if not os.path.isfile(dir + "costs.txt"):
        print("costs.txt not found in {}".format(dir))
        return
    costs = np.genfromtxt(dir + "costs.txt")

    if costs.shape == ():
        return
    if costs.shape == (0,):
        return
    if math.isnan(costs[-1]):
        return
    iter_to_SD_s0 = {}
    for i in range(50000):
        file = dir + f"files/{i:07d}.cellParameters.input"
        if not os.path.isfile(file):
            continue
        s0_vals = pd.read_csv(file,sep = " ", header=None)[2].to_list()
        iter_to_SD_s0[i] = np.std(s0_vals)
    for i in range(50000):
        file = dir + f"{i:07d}.cellParameters.input"
        if not os.path.isfile(file):
            continue
        s0_vals = pd.read_csv(file,sep = " ", header=None)[2].to_list()
        iter_to_SD_s0[i] = np.std(s0_vals)
    pd.DataFrame({"iter": list(iter_to_SD_s0.keys()), "SD_s0": list(iter_to_SD_s0.values())}).to_csv(f"{dir}SD_s0.csv", index=False)


# def main():
#     l = int(sys.argv[1])
#     n = int(sys.argv[2])

#     for i in range(100):
#         dir = f"kv_10_l_{l}_n_sp_{n:03d}/{i:03d}/"
#         write_SD_s0_in_dir(dir)
def main():
    l = int(sys.argv[1])
    n = int(sys.argv[2])

    for i in range(100):
        dir = f"kv_10_l_{l}_n_{n}/{i:03d}/"
        write_SD_s0_in_dir(dir)
if __name__ == "__main__":    
    main()