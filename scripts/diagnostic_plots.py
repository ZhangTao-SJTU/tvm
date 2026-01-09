from toolbox import manuscriptPlots
from toolbox.periodic import PeriodicTissue
from toolbox.minimization import FIREminimization
from toolbox.stress import calculate_max_shear_stress
import glob
import numpy as np
import pandas as pd
import os

def error_to_epoch(run_dir,n_cells,out_dir, filename="error_vs_epoch.png"):
    os.makedirs(out_dir,exist_ok=True)
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(1e-7,5)
    plotter.set_xlim(1,1500)
    plotter.set_xticks([500*i for i in range(150)])
    plotter.set_yticks([0.2*i for i in range(1,100)])
    plotter.set_xlabel("Epochs")
    plotter.set_ylabel(r"$\sum|1-\sigma_{target}/\sigma_{current}|$")

    plotter.set_yScaled()
    plotter.initialize_figure()
    plotter.set_yLog()
    df = pd.read_csv("{}info.csv".format(run_dir))
    err = df["Error"].to_numpy()
    iters = [(i+1)/n_cells for i in range(len(err))]
    plotter.plot_xy(iters,err,color = "red",alpha = 1)
    plotter.save_fig("{}/{}".format(out_dir,filename))

def single_cell_error_to_iters(run_dir,out_dir,filename="single_cell_error_to_iters.png"):
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(1e-9,5)
    plotter.set_xlim(1,50000)
    plotter.set_xticks([25000*i for i in range(150)])
    plotter.set_yticks([0.2*i for i in range(1,100)])
    plotter.set_xlabel("Iterations")
    plotter.set_ylabel("Single Cell Error")
    plotter.set_yScaled()
    plotter.initialize_figure()
    # plotter.set_xLog()
    plotter.set_yLog()
    err = np.loadtxt("{}costs.txt".format(run_dir))
    iters = [i+1 for i in range(len(err))]
    plotter.plot_xy(iters,err)
    plotter.save_fig("{}{}".format(out_dir,filename))

def overlap_to_iters(run_dir,out_dir,filename="overlap_to_iters.png"):
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(0,1)
    plotter.set_xlim(1,50000)
    plotter.set_xticks([25000*i for i in range(10)])
    plotter.set_yticks([0.2*i for i in range(1,100)])
    plotter.set_xlabel("Iterations")
    plotter.set_ylabel(r"$Q_n$")

    plotter.set_yScaled()
    plotter.initialize_figure()
    # plotter.set_xLog()
    # plotter.set_yLog()
    err = np.loadtxt("{}q_values.txt".format(run_dir))
    iters = [i+1 for i in range(len(err))]
    plotter.plot_xy(iters,err)
    plotter.save_fig("{}{}".format(out_dir,filename))

def parameter_space_distance_to_epochs(run_dir,out_dir,filename="distance_to_epochs.png"):
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(0,20)
    plotter.set_xlim(0,1500)
    plotter.set_xticks([500*i for i in range(10)])
    plotter.set_yticks([10*i for i in range(1,10)])
    plotter.set_xlabel("Epochs")
    plotter.set_ylabel(r"$\sum (s_0^{hidden}-s_0^{(init)})^2$")

    # plotter.set_yScaled()
    plotter.initialize_figure()
    df = pd.read_csv("{}info.csv".format(run_dir))
    err = df["Distance"].to_numpy()
    iters = [(i+1)/4 for i in range(len(err))]
    plotter.plot_xy(iters,err,color = "red",alpha = 1)
    plotter.save_fig("{}{}".format(out_dir,filename))

def s0_histogram(run_dir,out_dir,filename="s0_histogram.png"):
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(0,1.5)
    plotter.set_xlim(2,8)
    plotter.set_xticks([i for i in range(3,8)])
    plotter.set_yticks([i for i in range(1,10)])
    plotter.set_xlabel(r"$s_0^{hidden}$")

    # plotter.set_yScaled()
    plotter.initialize_figure()
    df = pd.read_csv("{}cellParameters.input".format(run_dir), sep = " ", header=None)
    plotter.histogram_from_dataframe(df[2],bins =20)
    plotter.save_fig("{}{}".format(out_dir,filename))

def last_epoch_error(run_dir,n_cells,out_dir,filename="last_epoch.png"):
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(1e-8,1e0)
    plotter.set_xlim(1,4)
    plotter.set_xticks([i for i in range(1,8)])
    plotter.set_yticks([0.00001*i for i in range(0,5)])
    plotter.set_xlabel("Training Step")
    plotter.set_ylabel("Error")
    plotter.set_yScaled()
    plotter.initialize_figure()
    df = pd.read_csv("{}info.csv".format(run_dir), sep = ",")
    print(df.columns)
    #check if value is the same for the last four rows in "epochs" column. if not, remove one row:
    epochs = df["Epoch"].to_numpy()
    while not np.all(epochs[-1*n_cells:] == epochs[-1]):
        df = df.iloc[:-1]
        epochs = df["Epoch"].to_numpy()
    
    iterations = df["Iter"].to_numpy()[-1*n_cells:]
    errors = df["Error"].to_numpy()[-1*n_cells:]

    target_cells = pd.read_csv("{}initial_stress.csv".format(run_dir), sep = ",")["CellID"].to_numpy()
    target_stress = pd.read_csv("{}initial_stress.csv".format(run_dir), sep = ",")["Target"].to_numpy()[0]
    target_cell_to_error ={i:[] for i in target_cells}

    # costs = np.loadtxt("{}costs.txt".format(run_dir))
    # # print(iter, costs[iter[0]+1:iter[-1]+1])
    # iters_to_costs = {}
    colors = ["purple", "blue", "green", "orange"]
    # for i in range(len(iter)-1):
    #     iters_array = np.arange(iter[i]+1, iter[i+1]+1)
    #     costs_array = costs[iter[i]+1:iter[i+1]+1]
    #     print(iters_array, costs_array)
    #     plotter.plot_xy(iters_array, costs_array, label = "Cell {}".format(i+1), color = colors[i])
    plotter.plot_xy([i for i in range(1,5)], errors, label = "Mean Error", color = "red")
    for num,i in enumerate(target_cell_to_error):
        for iter in iterations:
            sample = PeriodicTissue.from_config(run_dir+"files/","{:07d}.bulk.txt".format(iter))
            min = FIREminimization.periodic_tissue(sample)
            min.load_cell_parameters("{:07d}.cellParameters.input".format(iter))
            target_cell_to_error[i].append(abs(calculate_max_shear_stress(min._config,i)-target_stress)/target_stress)
        print(target_cell_to_error[i])  
        plotter.plot_scatter([i for i in range(1,5)], target_cell_to_error[i], color=colors[num], label = "Cell {}".format(num+1), alpha = 0.5)
    plotter.ax.legend(loc='upper right')
    plotter.set_yLog()
    plotter.save_fig("{}{}".format(out_dir,filename))

def displacement_histogram(run_dir,out_dir,filename="displacement_histogram.png"):
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(0,5)
    plotter.set_xlim(0,3)
    plotter.set_xticks([i for i in range(4)])
    plotter.set_yticks([5*i for i in range(1,5)])
    plotter.set_xlabel("Cell Displacement")

    # plotter.set_yScaled()
    plotter.initialize_figure()
    init_sample = PeriodicTissue.from_config(run_dir,"init_config.txt")
    # final_sample = PeriodicTissue.from_config(run_dir,sorted(glob.glob("{}*.bulk.txt".format(run_dir)))[-1])
    costs = np.loadtxt("{}costs.txt".format(run_dir))
    final_iter = len(costs)-1
    final_sample = PeriodicTissue.from_config(run_dir,"{:07d}.bulk.txt".format(final_iter))
    init_coordinates = {}
    displacements = []
    for cellID,cell in init_sample.cells_.items():
        if cell.center_ is None:
            continue
        init_coordinates[cellID] = cell.center_
    for cellID,cell in final_sample.cells_.items():
        if cell.center_ is None:
            continue
        if cellID in init_coordinates:
            displacement = np.linalg.norm(np.subtract(cell.center_,init_coordinates[cellID]))
            displacements.append(displacement)
    plotter.histogram_from_array(displacements, bins=5, color='blue', alpha=0.7)
    plotter.save_fig("{}{}".format(out_dir,filename))

def plot_all(run_dir,out_dir,n_cells):
    # error_to_epoch(run_dir,n_cells,out_dir)
    # single_cell_error_to_iters(run_dir,out_dir)
    # overlap_to_iters(run_dir,out_dir)
    # parameter_space_distance_to_epochs(run_dir,out_dir)
    # s0_histogram(run_dir,out_dir)
    # last_epoch_error(run_dir,n_cells,out_dir)
    displacement_histogram(run_dir,out_dir)
