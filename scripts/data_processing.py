from toolbox import stress
from toolbox.periodic import PeriodicTissue
from toolbox.minimization import FIREminimization
from toolbox.overlap import calculate_Q
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
                    if float(lines[-1]) > 1e-6:
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

def write_histogram_data(dir_list,filename):
    # init_stresses = []
    final_stresses = []
    # init_s0 = []
    final_s0 = []
    for dir in dir_list:
        print(dir)
        # init_file = "init_config.txt"
        costs = np.loadtxt(dir + "costs.txt")
        final_iter = len(costs)-1
        # final_file = sorted(glob.glob(dir + "*.bulk.txt"))[-1].split("/")[-1]
        final_file = "{:04d}.bulk.txt".format(final_iter)
        print("Final file: ", final_file)
        # init_tissue = PeriodicTissue.from_config(dir, init_file)
        final_tissue = PeriodicTissue.from_config(dir, final_file)
        # init = FIREminimization.periodic_tissue(init_tissue)
        final = FIREminimization.periodic_tissue(final_tissue)
        # if os.path.isfile("{}init_cellParameters.txt".format(dir)):
            # init.load_cell_parameters("init_cellParameters.txt")
        final.load_cell_parameters("{:04d}.cellParameters.input".format(final_iter))

        # for cellID,cell in init_tissue.cells_.items():
        #     cell.max_shear_stress_ = stress.calculate_max_shear_stress(init_tissue,cellID)
        #     init_stresses.append(cell.max_shear_stress_)
        #     init_s0.append(cell.s0_)
        for cellID,cell in final_tissue.cells_.items():
            cell.max_shear_stress_ = stress.calculate_max_shear_stress(final_tissue,cellID)
            final_stresses.append(cell.max_shear_stress_)
            final_s0.append(cell.s0_)
    # data_dict = {"Initial_Stress":init_stresses, "Final_Stress": final_stresses, "Initial_s0": init_s0, "Final_s0": final_s0}
    data_dict = {"Final_Stress": final_stresses, "Final_s0": final_s0}

    df = pd.DataFrame(data_dict)
    df.to_csv(filename)
#write overlap between consecutive epochs. set value at epoch 0 to be 1.

'''
def write_Q2(dir,**kwargs):
    intervals = np.arange(0,1000,100)
    if "intervals" in kwargs:
        intervals = kwargs["intervals"]
    if not os.path.isdir(dir):
        raise ValueError("Dir doesnt exist!")
    if not os.path.isfile(dir + "costs.txt"):
        print("No costs.txt: maybe run hasn't started yet.")
        return
    costs = np.loadtxt(dir + "costs.txt")
    init_config = PeriodicTissue.from_config(dir, "init_config.txt") 
    init_config.evaluate_cell_neighbors()
    init_cell_neighbors_dict = init_config.cell_neighbors_
    iter_to_sample = {i:PeriodicTissue.from_config(dir, "{:04d}.bulk.txt".format(i)) for i in intervals}
    for _,sample in iter_to_sample.items():
        sample.evaluate_cell_neighbors()
    
    iter_to_overlap = {i:[] for i in intervals}
    # for i,iter in enumerate(intervals):
    #     current_neighbors_dict = iter_to_sample[iter].cell_neighbors_
    #     if i == 0:
    #         iter_to_overlap[iter] = calculate_Q(current_neighbors_dict,init_cell_neighbors_dict)
    #         continue
    #     previous_iter = intervals[i-1]
    #     current_neighbors_dict = iter_to_sample[iter].cell_neighbors_
    #     previous_neighbors_dict = iter_to_sample[previous_iter].cell_neighbors_
    #     iter_to_overlap[iter] = calculate_Q(current_neighbors_dict,previous_neighbors_dict)
    for iter in intervals:
        current_neighbors_dict = iter_to_sample[iter].cell_neighbors_
        iter_to_overlap[iter] = calculate_Q(current_neighbors_dict,init_cell_neighbors_dict)
   
   
    df = pd.DataFrame({"iteration":list(iter_to_overlap.keys()),"overlap":list(iter_to_overlap.values())})
    df.to_csv("{}overlaps.csv".format(dir))
'''

def write_Q2(dir):
    if not os.path.isdir(dir):
        raise ValueError("Dir doesnt exist!")
    if not os.path.isfile(dir + "costs.txt"):
        print("No costs.txt: maybe run hasn't started yet.")
        return
    costs = np.loadtxt(dir + "costs.txt")
    init_config = PeriodicTissue.from_config(dir, "init_config.txt") 
    init_config.evaluate_cell_neighbors()
    init_cell_neighbors_dict = init_config.cell_neighbors_
    samples = [PeriodicTissue.from_config(dir, "{:04d}.bulk.txt".format(i)) for i in range(len(costs))]
    overlaps = [1]
    for sample in samples:
        sample.evaluate_cell_neighbors()
        overlaps.append(calculate_Q(init_cell_neighbors_dict,sample.cell_neighbors_))
    np.savetxt("{}overlaps.txt".format(dir),overlaps)

# Write average error - the 0th element being the initial (pretraining) error
def write_average_error(dir):
    if not os.path.isdir(dir):
        raise ValueError("Dir doesnt exist!")
    errors = []
    if not os.path.isfile(dir + "costs.txt"):
        print("No costs.txt: maybe run hasn't started yet.")
        return
    costs = np.loadtxt(dir + "costs.txt")
    initial_stress = pd.read_csv("{}initial_stress.csv".format(dir))
    init = initial_stress["Current"].to_numpy()
    targets = pd.read_csv("{}0000.stresses.csv".format(dir))["Target"].to_numpy()
    # Append initial (pretraining value).
    # Note that 0000.stresses.csv has current stresses AFTER the 0th iteration of training...
    errors.append(np.mean(abs((init - targets))/init))
    for iter in range(len(costs)):
        current_stress_file = "{}{:04d}.stresses.csv".format(dir,iter)
        current_stress = pd.read_csv(current_stress_file)
        mean_stress = np.mean(abs((current_stress["Current"].to_numpy() - targets))/current_stress["Current"].to_numpy())
        errors.append(mean_stress)
    np.savetxt("{}errors.txt".format(dir), errors)

def write_mean_error_from_dirlist(dirlist,filename):
    errors_dict = {0:[]}
    for dir in dirlist:
        if not os.path.isdir(dir):
            print("Dir doesnt exist!")
            continue
        if not os.path.isfile(dir + "costs.txt"):
            print("No costs.txt: maybe this run hasn't started yet.")
            continue
        
        with open(dir + "costs.txt", "r") as f:
            lines = f.readlines()
            if len(lines) < 2:
                print("Costs file too short.")
                continue
            if float(lines[-1]) > 1e-12:
                print("Run hasn't converged.")
                continue
        costs = np.loadtxt(dir + "costs.txt")
        initial_stress = pd.read_csv("{}initial_stress.csv".format(dir))
        init = initial_stress["Current"].to_numpy()
        targets = pd.read_csv("{}0000.stresses.csv".format(dir))["Target"].to_numpy()
        # Append initial (pretraining value).
        # Note that 0000.stresses.csv has current stresses AFTER the 0th iteration of training...
        errors_dict[0].extend((abs((init - targets))/init))
        for iter in range(len(costs)):
            if iter not in errors_dict:
                errors_dict[iter] = []
            current_stress_file = "{}{:04d}.stresses.csv".format(dir,iter)
            current_stress = pd.read_csv(current_stress_file)
            errors_dict[iter].extend(abs((current_stress["Current"].to_numpy() - targets))/current_stress["Current"].to_numpy())
    for i,array in errors_dict.items():
        print(i,array)
    errors = {"mean":[], "sem":[]}
    for i, array in errors_dict.items():
        print(i,array,"\n\n\n")
        if len(array)<10:
            continue
        errors["mean"].append(np.mean(array))
        errors["sem"].append(stats.sem(array))
    df = pd.DataFrame(errors)
    df.to_csv(filename)
'''
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
            if float(lines[-1])>1e-4:
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
            iteration_to_costs[iter+1].append(mean_stress)
    mn = []
    sem = []
    for i, array in iteration_to_costs.items():
        if len(array)<3: continue
        mn.append(np.mean(array))
        sem.append(stats.sem(array))
        # sd.append(stats.sem(array))
    avg_stress = {"mean":mn, "sem":sem}
    df = pd.DataFrame(avg_stress)
    df.to_csv("{}error.csv".format(experiment))
'''
def download(experiment):
    os.system("scp mameen@smatter-login.syr.edu:/home/mameen/{}error.csv {}error.csv".format(experiment,experiment))

# def download_all():
#     s0_vals = [4.8,4.9,5.0,5.1,5.2,5.3]
#     d_list = ["1_cell_increase/", "1_cell_decrease/", "2_cells/","4_cells/"]
#     for d in d_list:
#         os.makedirs(d, exist_ok=True)
#         for s0 in s0_vals:
#             experiment = d+"7_{:.1f}/".format(s0)
#             os.makedirs(experiment,exist_ok=True)
#             download(experiment)

def download_all(d):
    os.makedirs("data/"+d, exist_ok=True)
    os.system("scp mameen@smatter-login.syr.edu:/home/mameen/{}stresses.txt data/{}".format(d,d))
    os.system("scp mameen@smatter-login.syr.edu:/home/mameen/{}histogram_data.csv data/{}".format(d,d))

    for i in range(100):
        os.makedirs("data/{}{:03d}/".format(d,i), exist_ok=True)
        os.system(("scp mameen@smatter-login.syr.edu:/home/mameen/{}{:03d}/errors.txt data/{}{:03d}/".format(d,i,d,i)))
        os.system(("scp mameen@smatter-login.syr.edu:/home/mameen/{}{:03d}/target data/{}{:03d}/".format(d,i,d,i)))

# def main():
#     s0_vals = [4.8,4.9,5.0,5.1,5.2,5.3]
#     d_list = ["1_cell_increase/", "1_cell_decrease/", "2_cells/","4_cells/"]
#     for d in d_list:
#         for s0 in s0_vals:
#             experiment = d+"7_{:.1f}/".format(s0)
#             write_average_error(experiment)
#             # write_histogram_data(experiment)

def main():
    for experiment in ["1_cell_decrease_2_sigma/","1_cell_increase_2_sigma/","2_cells_mean/","4_cells_mean/"]:
        experiment = "/home/mameen/{}".format(experiment)
        dirlist = []
        for dir in ["{}{:03d}/".format(experiment,i) for i in range(100)]:
            if not os.path.isfile(dir+"costs.txt"):
                continue
            with open(dir+"costs.txt","r") as f:
                lines = f.readlines()
                if len(lines)<2:
                    continue
                if float(lines[-1])>1e-10:
                    continue
            dirlist.append(dir)
        write_histogram_data(dirlist,"{}/histogram_data.csv".format(experiment))    

if __name__ == "__main__":
    main()
