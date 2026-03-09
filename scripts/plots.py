from turtle import title

from toolbox import manuscriptPlots
from toolbox.stress import calculate_max_shear_stress
import math
import numpy as np
import pandas as pd
from scipy.stats import gamma
from scipy.stats import sem
import os
from toolbox.spheroid import Spheroid
from toolbox.periodic import PeriodicTissue
from toolbox.training import Training
from concurrent.futures import ProcessPoolExecutor, as_completed


# from scripts.plots import *
tolerance = 1e-6
lwmap = {1:15, 2:10, 3:5, 4:5, 5:5, 6:5}
lw_sp_map = {40:15,20:10,10:5, 5:20}
n_cells_color_map = {1:"blue",2:"orange",3:"purple", 4:"red" ,5:"yellow",6:"green"}
l_to_marker = {4:'o', 5:"s", 6:"^"}
l_to_alpha = {4:1, 5: 0.7, 6:0.5}
l_to_color = {4:"black", 5: "purple", 6:"green"}
n_runs = 100

def create_inset_plot(main_plot_file,inset_plot_file,inset_position,filename):
    import matplotlib.pyplot as plt
    from matplotlib.image import imread
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes    
    img_main = imread(main_plot_file)
    img_inset = imread(inset_plot_file)

    fig = plt.figure(figsize=(15,15))

    ax_main = fig.add_axes([0, 0, 1, 1])
    ax_main.imshow(img_main)
    ax_main.axis("off")

    # Inset image (x, y, width, height in figure coords)
    ax_inset = fig.add_axes(inset_position)
    ax_inset.imshow(img_inset)
    ax_inset.axis("off")

    plt.savefig(filename, transparent=True)

def write_final_s0(dir_list,savefile):
    s0_vals = []
    for dir in dir_list:
        df = pd.read_csv(dir+"cellParameters.input", sep = " ", header=None)
        s0_vals += df[2].to_list()
    np.savetxt(savefile,s0_vals)

# def write_stresses(dirlist,savefile):
#     stress_vals = {i:None for i in range(len(dirlist))}
#     for i,dir in enumerate(dirlist):
#         if not os.path.isfile(dir+"minimized.txt"):
#             continue
#         print(dir)
#         sample = PeriodicTissue.from_config(dir,"minimized.txt")
#         tr = Training.from_sample(sample)
#         tr.load_cell_parameters()
#         stress_vals[i] = [calculate_max_shear_stress(tr._config,cellID) for cellID in tr._config.cells_]
#     np.savetxt(savefile,list(stress_vals.values()), fmt="%.7f")
def process_single(i, dir):
        if not os.path.isfile(dir + "minimized.txt"):
            return i, None
        sample = PeriodicTissue.from_config(dir, "minimized.txt")
        tr = Training.from_sample(sample)
        tr.load_cell_parameters()
        stresses = [calculate_max_shear_stress(tr._config, cellID) for cellID in tr._config.cells_]
        return i, stresses

def write_stresses(dirlist, savefile, max_workers=8):
    stress_vals = {i: None for i in range(len(dirlist))}
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_single, i, d): i for i, d in enumerate(dirlist)}
        for future in as_completed(futures):
            i, result = future.result()
            stress_vals[i] = result
    np.savetxt(savefile, [val for val in stress_vals.values() if val is not None], fmt="%.7f")

def find_complete_runs(dirlist):
    complete_dirs =[]
    for test_dir in dirlist:
        if not os.path.isfile(test_dir + "costs.txt"):
            continue
        costs = np.genfromtxt(test_dir + "costs.txt")
        if costs.shape == ():
            continue
        if costs.shape == (0,):
            continue
        if costs[-1]>tolerance:
            continue
        if math.isnan(costs[-1]):
            continue
        print(test_dir,costs[-1])
        complete_dirs.append(test_dir)
    print("yield:", len(complete_dirs),complete_dirs[0],"\n\n")
    return complete_dirs

def find_complete_runs_multiple_patterns(dirlist):
    complete_dirs = []
    for test_dir in dirlist:
        if not os.path.isfile(test_dir+"info.csv"):
            continue
        info = pd.read_csv(test_dir+"info.csv")
        if info["Error"].to_list()[-1]>tolerance:
            continue
        complete_dirs.append(test_dir)
    print("yield:", len(complete_dirs),complete_dirs[0],"\n\n")

    return complete_dirs

def single_tracks_to_iters( label_to_data,savefile,
                            input_filename = "costs.txt",
                            title = None,
                            ylim = [1e-9,5],
                            xlim =[1,3000],
                            yticks =[5*i for i in range(1,100)],
                            ylabel =r"$<|1-\sigma_{T}/\sigma|>$",
                            ylog = True):
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(ylim)
    plotter.set_xlim(xlim)
    plotter.set_yticks(yticks)
    plotter.set_yticks([0.2*i for i in range(1,100)])
    plotter.set_xlabel("Iterations")
    plotter.set_ylabel(ylabel)
    # plotter.set_title(r"$n_{(total)} =$"+ "{}".format(l**3))
    if title is not None:
        plotter.set_title(title)
    plotter.initialize_figure()
    plotter.ax.set_xscale("log")
    if ylog:
        plotter.ax.set_yscale("log")
    for label, data in label_to_data.items():
        dirlist = data["dirlist"]
        label_added = False
        for test_dir in dirlist:
            y_vals = np.genfromtxt(test_dir+input_filename)
            if not label_added:
                plotter.plot_xy([i+1 for i in range(len(y_vals))],y_vals,label=label,color = data.get("color","black"),alpha = data.get("alpha",0.6),linewidth = data.get("linewidth",10) )
                label_added = True
            else: 
                plotter.plot_xy([i+1 for i in range(len(y_vals))],y_vals,label="_none",color = data.get("color","black"),alpha = data.get("alpha",0.6),linewidth = data.get("linewidth",10) )
    plotter.ax.legend()
    plotter.save_fig(savefile)

def histogram(label_to_data,
              input_filename = "s0_vals.txt",
              title = None,
              xlim = [3,7],
              ylim = [0,5],
              xlabel = r"$s_0$",
              xticks = [i for i in range(3,8)],
              yticks = [i for i in range(1,10)]):
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(ylim)
    plotter.set_xlim(xlim)
    plotter.set_xticks(xticks)
    plotter.set_yticks(yticks)
    plotter.set_xlabel(xlabel)
    if title is not None:
        plotter.set_title(title)
    plotter.initialize_figure()
    for label, data in label_to_data.items():
        array = np.genfromtxt(data["dir"]+input_filename)
    # plotter.ax.vlines(x = vline_x, ymin = vline_min, ymax = vline_max, linestyle= "dashed",color = plotter.neutralColor, label = r"$s_0^{(init)}$")
        plotter.histogram_from_array(array, fit_type= 'kde',bins = data.get("bins",30), label = label,alpha = data.get("alpha",0.3),color = data.get("color","black"),edgecolor = data.get("color","black"),linewidth = data.get("linewidth",5))
    return plotter

def scatter_plot(label_to_data,
            title = None,
            ylim = [0.2,1.05],
            xlim = [10,10000],
            xticks = [0.2*i for i in range(10)],
            yticks =[.2*i for i in range(1,10)],
            xlabel = r"$i_{final}$",
            ylabel =r"$Q_2$",
            xlog=True,
            ylog = False):
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(ylim)
    plotter.set_xlim(xlim)
    plotter.set_xticks(xticks)
    plotter.set_yticks(yticks)


    plotter.set_ylabel(ylabel)
    plotter.set_xlabel(xlabel)
    if title is not None:
        plotter.set_title(title)
    if xlog:
        plotter.set_xLog()
    if ylog:
        plotter.set_yLog()
    plotter.initialize_figure()
    # plotter.ax.set_yscale("log")


    # for label,data in label_to_data.items():
    #     plotter.plot_scatter(
    #         data["x_array"],
    #         data["y_array"],
    #         marker=data.get("marker","o"),
    #         # markeredgecolor="black",
    #         #
    #         # markeredgewidth=2,
    #         s=data.get("s",2000),
    #         lw = 5,
    #         color= data.get("color","black"),
    #         alpha = data.get("alpha",0.4),
    #         label = label,
    #         edgecolors = (0, 0, 0, 1.0),
    #     )
    for label,data in label_to_data.items():
        plotter.plot_scatter(
            data["x_array"],
            data["y_array"],
            marker=data.get("marker","o"),
            s=data.get("s",2000),
            color= data.get("color","black"),
            alpha = data.get("alpha",0.4),
            label = label,
        )
        plotter.plot_scatter(
            data["x_array"],
            data["y_array"],
            marker=data.get("marker","o"),
            facecolors = "none",
            edgecolors = "black",
            # edgecolors = data.get("color","black"),
            s=data.get("s",2000),
            lw = 1,
            label = "_none",
            alpha = 0.5,
        )

    return plotter

def stacked_plot(**kwargs):
    plotter = manuscriptPlots.plot_shared_x_axis()
    plotter.set_xlim(kwargs.get("xlim", [10,10000]))
    plotter.set_ylim_top(kwargs.get("ylim_top", [1e-8,1e1]))
    plotter.set_ylim_bottom(kwargs.get("ylim_bottom", [0.5,1.1]))
    plotter.set_ylim_bottom_right(kwargs.get("ylim_bottom_right", [0,0.5]))
    plotter.set_yticks_top(kwargs.get("yticks_top", None))
    plotter.set_yticks_bottom(kwargs.get("yticks_bottom", [0.6,0.8,1]))
    plotter.set_xLog(kwargs.get("xlog", True))
    plotter.set_yLog_top(kwargs.get("ylog_top", True))
    plotter.set_yLog_bottom(kwargs.get("ylog_bottom", False))
    plotter.set_yLog_bottom_right(kwargs.get("ylog_bottom_right", False))
    plotter.set_xlabel(kwargs.get("xlabel", "Iterations"))
    plotter.set_ylabel_top(kwargs.get("ylabel_top", r"$<|1-\sigma_{T}/\sigma|>$"))
    plotter.set_ylabel_bottom(kwargs.get("ylabel_bottom", r"$Q_2$"))
    plotter.set_ylabel_bottom_right(kwargs.get("ylabel_bottom_right", r"$SD(s_0)$"))
    plotter.set_title(kwargs.get("title", None))
    plotter.initialize_sharedx_figure()
    return plotter

