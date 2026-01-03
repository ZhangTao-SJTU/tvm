from toolbox import manuscriptPlots
import numpy as np
import pandas as pd
import os

def error_to_epoch(run_dir,n_cells,out_dir, filename="error_vs_epoch.png"):
    os.makedirs(out_dir,exist_ok=True)
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(1e-7,5)
    plotter.set_xlim(1,500)
    plotter.set_xticks([100*i for i in range(150)])
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
    plotter.set_xlim(1,15000)
    plotter.set_xticks([5000*i for i in range(150)])
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
    plotter.set_xlim(1,15000)
    plotter.set_xticks([5000*i for i in range(150)])
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
    plotter.set_ylim(0,7)
    plotter.set_xlim(0,500)
    plotter.set_xticks([200*i for i in range(150)])
    plotter.set_yticks([2*i for i in range(1,10)])
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

def plot_all(run_dir,out_dir,n_cells):
    error_to_epoch(run_dir,n_cells,out_dir)
    single_cell_error_to_iters(run_dir,out_dir)
    overlap_to_iters(run_dir,out_dir)
    parameter_space_distance_to_epochs(run_dir,out_dir)
    s0_histogram(run_dir,out_dir)
