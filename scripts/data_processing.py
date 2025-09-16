from toolbox import stress
from toolbox.periodic import PeriodicTissue
from toolbox.minimization import FIREminimization
import os
import glob
import numpy as np
import pandas as pd
from scipy import stats

def write_patterns_histogram_data(experiments_list):
    for experiment in experiments_list:
        for pattern in ["PatternA","PatternB"]:
            init_stresses = []
            final_stresses = []
            init_s0 = []
            final_s0 = []
            for i in range(20):
                # dir = "/home/mameen/init_homogeneous/7_{}/".format(i)
                dir = "{}/run_{}/{}/".format(experiment,i,pattern)
                if not os.path.isfile(dir + "costs.txt"):
                    continue
                with open(dir + "costs.txt", "r") as f:
                    lines = f.readlines()
                    if len(lines) < 2:
                        continue
                    if float(lines[-1]) > 1e-5:
                        continue
                    print(dir)
                init_file = "init_config.txt"
                final_file = sorted(glob.glob(dir + "*.bulk.txt"))[-1].split("/")[-1]
                print("Final file: ", final_file)
                init_tissue = PeriodicTissue.from_config(dir, init_file)
                final_tissue = PeriodicTissue.from_config(dir, final_file)
                init = FIREminimization.periodic_tissue(init_tissue)
                final = FIREminimization.periodic_tissue(final_tissue)
                if os.path.isfile("{}init_cellParameters.txt".format(dir)):
                    init.load_cell_parameters("init_cellParameters.txt")
                final.load_cell_parameters()

                for cellID,cell in init_tissue.cells_.items():
                    cell.max_shear_stress_ = stress.calculate_max_shear_stress(init_tissue,cellID)
                    init_stresses.append(cell.max_shear_stress_)
                    init_s0.append(cell.s0_)
                for cellID,cell in final_tissue.cells_.items():
                    cell.max_shear_stress_ = stress.calculate_max_shear_stress(final_tissue,cellID)
                    final_stresses.append(cell.max_shear_stress_)
                    final_s0.append(cell.s0_)
            data_dict = {"Initial_Stress":init_stresses, "Final_Stress": final_stresses, "Initial_s0": init_s0, "Final_s0": final_s0}
            df = pd.DataFrame(data_dict)
            df.to_csv("{}/{}_histogram_data.csv".format(experiment,pattern))

def write_patterns_costs(experiment_list):
    # experiment = "sp_1_cell_increase"
    for experiment in experiment_list:
        for pattern in ["PatternA","PatternB"]:
            iteration_to_costs = {i:[]for i in range(1000)}
            for i in range(20):
                # dir = "/home/mameen/init_homogeneous/7_{}/".format(i)
                # dir = "/home/mameen/{}/run_{}/".format(experiment,i)
                dir = "{}/run_{}/{}/".format(experiment,i,pattern)
                # print(dir)
                if not os.path.isfile(dir + "costs.txt"):
                    continue
                with open(dir + "costs.txt", "r") as f:
                    lines = f.readlines()
                    if len(lines) < 2:
                        continue
                    if float(lines[-1]) > 5e-6:
                        continue
                print(dir)
                initial_stress = pd.read_csv("{}initial_stress.csv".format(dir))
                # print(initial_stress)
                init = initial_stress["Stress"].to_numpy()
                targets = pd.read_csv("{}0000.stresses.csv".format(dir))["Target"].to_numpy()
                iteration_to_costs[0].append(np.mean(abs((init - targets))/init))
                for iter in range(len(lines)):
                    current_stress_file = "{}{:04d}.stresses.csv".format(dir,iter)
                    current_stress = pd.read_csv(current_stress_file)
                    mean_stress = np.mean(abs((current_stress["Current"].to_numpy() - targets))/current_stress["Current"].to_numpy())
                    # print(mean_stress)
                    iteration_to_costs[iter+1].append(mean_stress)
            mn = []
            sm = []
            for i, array in iteration_to_costs.items():
                if len(array)<3: continue
                mn.append(np.mean(array))
                sm.append(stats.sem(array))
            avg_stress = {"mean":mn, "sem":sm}
            df = pd.DataFrame(avg_stress)
            df.to_csv("{}/{}_costs.csv".format(experiment,pattern))

def write_distances(experiment_list):
    # experiment = "sp_1_cell_increase"
    for experiment in experiment_list:

        iteration_to_distances = {i:[]for i in range(1000)}
        for i in [0,7,8,15,16,18]:
            dir = "{}/run_{}/".format(experiment,i)
            print()
            if not os.path.isfile(dir + "distances.txt"):
                continue
            with open(dir + "distances.txt") as f:
                lines = f.readlines()
                if len(lines)<2:
                    continue
        
            init_dist = np.loadtxt("{}initial_distances.txt".format(dir))
            iteration_to_distances[0].append(init_dist)
            # print(initial_stress)

            distances = np.loadtxt("{}distances.txt".format(dir))
            print(distances)

            for iter, dist in enumerate(distances):
                iteration_to_distances[iter+1].append(dist)
        mn = []
        sm = []
        for i, array in iteration_to_distances.items():
            if len(array)<3: continue
            mn.append(np.mean(array))
            sm.append(stats.sem(array))
        avg_stress = {"mean":mn, "sem":sm}
        df = pd.DataFrame(avg_stress)
        df.to_csv("{}/distances.csv".format(experiment))

# histogram for stresses and s0s:

def write_histogram_data(experiment):
    init_stresses = []
    final_stresses = []
    init_s0 = []
    final_s0 = []
    for i in range(100):
        # dir = "/home/mameen/init_homogeneous/7_{}/".format(i)
        dir = "{}{:03d}/".format(experiment,i)
        if not os.path.isfile(dir + "costs.txt"):
            continue
        if not os.path.isfile(dir + "0000.bulk.txt"):
            continue
        with open(dir + "costs.txt", "r") as f:
            lines = f.readlines()
            if len(lines) < 2:
                continue
            # if float(lines[-1]) > 1e-6:
            #     continue
            print(dir)
        init_file = "init_config.txt"
        final_file = sorted(glob.glob(dir + "*.bulk.txt"))[-1].split("/")[-1]
        print("Final file: ", final_file)
        init_tissue = PeriodicTissue.from_config(dir, init_file)
        final_tissue = PeriodicTissue.from_config(dir, final_file)
        init = FIREminimization.periodic_tissue(init_tissue)
        final = FIREminimization.periodic_tissue(final_tissue)
        if os.path.isfile("{}init_cellParameters.txt".format(dir)):
            init.load_cell_parameters("init_cellParameters.txt")
        final.load_cell_parameters()

        for cellID,cell in init_tissue.cells_.items():
            cell.max_shear_stress_ = stress.calculate_max_shear_stress(init_tissue,cellID)
            init_stresses.append(cell.max_shear_stress_)
            init_s0.append(cell.s0_)
        for cellID,cell in final_tissue.cells_.items():
            cell.max_shear_stress_ = stress.calculate_max_shear_stress(final_tissue,cellID)
            final_stresses.append(cell.max_shear_stress_)
            final_s0.append(cell.s0_)
    data_dict = {"Initial_Stress":init_stresses, "Final_Stress": final_stresses, "Initial_s0": init_s0, "Final_s0": final_s0}
    df = pd.DataFrame(data_dict)
    df.to_csv("{}/histogram_data.csv".format(experiment))

def write_average_error(experiment):
    iteration_to_costs = {i:[]for i in range(1000)}
    for i in range(100):
        dir = "{}{:03d}/".format(experiment,i)
        if not os.path.isfile(dir + "costs.txt"):
            continue
        with open(dir + "costs.txt", "r") as f:
            lines = f.readlines()
            if len(lines) < 2:
                continue
            if float(lines[-1]) > 1e-6:
                continue
        initial_stress = pd.read_csv("{}initial_stress.csv".format(dir))
        init = initial_stress["Stress"].to_numpy()
        targets = pd.read_csv("{}0000.stresses.csv".format(dir))["Target"].to_numpy()
        # Append initial (pretraining value).
        # Note that 0000.stresses.csv has current stresses AFTER the 0th iteration of training...
        iteration_to_costs[0].append(np.mean(abs((init - targets))/init))
        for iter in range(len(lines)):
            current_stress_file = "{}{:04d}.stresses.csv".format(dir,iter)
            current_stress = pd.read_csv(current_stress_file)
            mean_stress = np.mean(abs((current_stress["Current"].to_numpy() - targets))/current_stress["Current"].to_numpy())
            print(mean_stress)
            iteration_to_costs[iter+1].append(mean_stress)
    mn = []
    sd = []
    for i, array in iteration_to_costs.items():
        if len(array)<3: continue
        mn.append(np.mean(array))
        sd.append(np.std(array))
        # sd.append(stats.sem(array))
    avg_stress = {"mean":mn, "sem":sd}
    df = pd.DataFrame(avg_stress)
    df.to_csv("{}error.csv".format(experiment))

def main():
    s0_vals = [4.8,4.9,5.0,5.1,5.2,5.3]
    d_list = ["1_cell_increase/", "1_cell_decrease/", "2_cells/","4_cells/"]
    for d in d_list:
        for s0 in s0_vals:
            experiment = d+"7_{:.1f}/".format(s0)
            write_histogram_data(experiment)
            write_average_error(experiment) 
    # write_histogram_data(experiments_list)
    # write_costs(experiments_list)
    
if __name__ == "__main__":
    main()
