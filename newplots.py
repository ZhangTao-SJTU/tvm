from toolbox import manuscriptPlots
import numpy as np
import pandas as pd
from scipy.stats import gamma
import os

lwmap = {1:15,2:10,3:5,4:5,5:5,6:5}
lw_sp_map = {40:15,20:10,10:5}
n_spheroid_color_map = {40:"blue", 20:"orange", 10:"red" }
n_cells_color_map = {1:"blue",2:"orange",3:"green", 4:"red" ,5:"purple",6:"brown"}

def error_to_iters_n_spheroid():
    l = 5
    savefile = "spheroid_graphs/spheroid_error_to_iters.png"
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(1e-8,5)
    # plotter.set_ylim(0.6,1)

    plotter.set_xlim(1,2000)
    plotter.set_xticks([20*i for i in range(150)])
    plotter.set_yticks([5*i for i in range(1,100)])
    plotter.set_yticks([0.2*i for i in range(1,100)])

    plotter.set_xlabel("Iterations")
    plotter.set_ylabel(r"$<|1-\sigma_{T}/\sigma|>$")
    # plotter.set_ylabel(r"$Q_n$")

    # plotter.set_title(r"$n_{(target)} =$"+ "{}".format(n_target))
    plotter.set_yScaled()

    plotter.initialize_figure()
    plotter.ax.set_xscale("log")
    plotter.ax.set_yscale("log")
    for n_spheroid in [10,20,40]:
        label_added = False
        for i in range(4):
            test_dir = "runs/kv_10_l_{}_n_hidden_{:03d}/{:03d}/".format(l,n_spheroid,i)
            if not os.path.isfile(test_dir + "costs.txt"):
                continue
            costs = np.loadtxt(test_dir + "costs.txt")
            if not label_added:
                plotter.plot_xy([i+1 for i in range(len(costs))],costs, label=r"$n_{sp}$"+  "={}".format(n_spheroid),color = n_spheroid_color_map[n_spheroid], linewidth=lw_sp_map[n_spheroid])
                label_added = True
            else:
                plotter.plot_xy([i+1 for i in range(len(costs))],costs,label = "_none",color = n_spheroid_color_map[n_spheroid], linewidth=lw_sp_map[n_spheroid])

    plotter.ax.legend(loc = "lower left")
    plotter.save_fig(savefile)

def error_to_iters_periodic_system_size():
    # n_cells_color_map = {1:"purple",2:"blue",3:"orange", 4:"yellow" ,6:"red" }
    for l in [4,5,6]:
        savefile = "spheroid_graphs/error_to_iters_periodic_l_{}.png".format(l)
        plotter = manuscriptPlots.plot()
        plotter.set_ylim(1e-8,5)
        # plotter.set_ylim(0.6,1)

        plotter.set_xlim(1,1500)
        plotter.set_xticks([20*i for i in range(150)])
        plotter.set_yticks([5*i for i in range(1,100)])
        plotter.set_yticks([0.2*i for i in range(1,100)])

        plotter.set_xlabel("Iterations")
        plotter.set_ylabel(r"$<|1-\sigma_{T}/\sigma|>$")
        # plotter.set_ylabel(r"$Q_n$")

        # plotter.set_title(r"$n_{(target)} =$"+ "{}".format(n_target))
        plotter.set_yScaled()

        plotter.initialize_figure()
        plotter.ax.set_xscale("log")
        plotter.ax.set_yscale("log")
        for n_cells in [1,2,3,4,5,6]:
            label_added = False
            for i in range(20):
                test_dir = "data/kv_10_l_{}_n_{}/{:03d}/".format(l,n_cells,i)
                if not os.path.isfile(test_dir + "costs.txt"):
                    continue
                costs = np.loadtxt(test_dir + "costs.txt")
                # print(costs)
                if not label_added:
                    plotter.plot_xy([i+1 for i in range(len(costs))],costs, label=r"$n_{T}$"+  "={}".format(n_cells),color = n_cells_color_map[n_cells],linewidth  = lwmap[n_cells])
                    label_added = True
                else:
                    plotter.plot_xy([i+1 for i in range(len(costs))],costs,label = "_none",color = n_cells_color_map[n_cells],linewidth  = lwmap[n_cells])

        plotter.ax.legend(loc = "lower left")
        plotter.save_fig(savefile,transparent=False)

def overlap_to_iters_periodic_system_size():
    for l in [4,5,6]:
        savefile = "spheroid_graphs/overlap_to_iters_periodic_l_{}.png".format(l)
        plotter = manuscriptPlots.plot()
        # plotter.set_ylim(1e-8,5)
        plotter.set_ylim(0.55,1.05)

        plotter.set_xlim(1,1500)
        plotter.set_xticks([20*i for i in range(150)])
        plotter.set_yticks([5*i for i in range(1,100)])
        plotter.set_yticks([0.2*i for i in range(1,100)])

        plotter.set_xlabel("Iterations")
        plotter.set_ylabel(r"$Q_2$")
        # plotter.set_ylabel(r"$Q_n$")

        # plotter.set_title(r"$n_{(target)} =$"+ "{}".format(n_target))
        plotter.set_yScaled()

        plotter.initialize_figure()
        plotter.ax.set_xscale("log")
        # plotter.ax.set_yscale("log")
        for n_cells in [1,2,3,4,5,6]:
            label_added = False
            for i in range(20):
                test_dir = "data/kv_10_l_{}_n_{}/{:03d}/".format(l,n_cells,i)
                if not os.path.isfile(test_dir + "costs.txt"):
                    continue
                costs = np.loadtxt(test_dir + "q_values.txt")
                # print(costs)
                if not label_added:
                    plotter.plot_xy([i+1 for i in range(len(costs))],costs, label=r"$n_{T}$"+  "={}".format(n_cells),color = n_cells_color_map[n_cells],linewidth = lwmap[n_cells])
                    label_added = True
                else:
                    plotter.plot_xy([i+1 for i in range(len(costs))],costs,label = "_none",color = n_cells_color_map[n_cells],linewidth  = lwmap[n_cells])

        plotter.ax.legend(loc = "lower left")
        plotter.save_fig(savefile,transparent=False)

def overlap_to_iters_n_spheroid():
    l = 5
    savefile = "spheroid_graphs/spheroid_overlap_to_iters.png"
    plotter = manuscriptPlots.plot()
    # plotter.set_ylim(1e-7,5)
    plotter.set_ylim(0.5,1.05)
    plotter.set_xlim(1,2000)
    plotter.set_xticks([20*i for i in range(150)])
    plotter.set_yticks([5*i for i in range(1,100)])
    plotter.set_yticks([0.2*i for i in range(1,100)])

    plotter.set_xlabel("Iterations")
    plotter.set_ylabel(r"$Q_2$")
    # plotter.set_ylabel(r"$Q_n$")

    # plotter.set_title(r"$n_{(target)} =$"+ "{}".format(n_target))
    # plotter.set_yScaled()

    plotter.initialize_figure()
    plotter.ax.set_xscale("log")
    # plotter.ax.set_yscale("log")
    for n_spheroid in [10,20,40]:
        label_added = False
        for i in range(4):
            test_dir = "runs/kv_10_l_{}_n_hidden_{:03d}/{:03d}/".format(l,n_spheroid,i)
            if not os.path.isfile(test_dir + "costs.txt"):
                continue
            costs = np.loadtxt(test_dir + "q_values.txt")
            if not label_added:
                plotter.plot_xy([i+1 for i in range(len(costs))],costs, label=r"$n_{sp}$"+  "={}".format(n_spheroid),color = n_spheroid_color_map[n_spheroid], linewidth=lw_sp_map[n_spheroid])
                label_added = True
            else:
                plotter.plot_xy([i+1 for i in range(len(costs))],costs,label = "_none",color = n_spheroid_color_map[n_spheroid], linewidth=lw_sp_map[n_spheroid])

    plotter.ax.legend(loc = "lower left")
    plotter.save_fig(savefile)


def s0_histogram():
    for l in [4,5,6]:
        savefile = "spheroid_graphs/s0_histogram_periodic_l_{}.png".format(l)
        plotter = manuscriptPlots.plot()
        # plotter.set_ylim(1e-8,5)
        plotter.set_ylim(0.55,1.05)

        plotter.set_xlim(1,1500)
        plotter.set_xticks([20*i for i in range(150)])
        plotter.set_yticks([5*i for i in range(1,100)])
        plotter.set_yticks([0.2*i for i in range(1,100)])

        plotter.set_xlabel("Iterations")
        plotter.set_ylabel(r"$<|1-\sigma_{T}/\sigma|>$")
        # plotter.set_ylabel(r"$Q_n$")

        # plotter.set_title(r"$n_{(target)} =$"+ "{}".format(n_target))
        plotter.set_yScaled()

        plotter.initialize_figure()
        plotter.ax.set_xscale("log")
        # plotter.ax.set_yscale("log")
        for n_cells in [1,2,4]:
            label_added = False
            for i in range(20):
                test_dir = "data/kv_10_l_{}_n_{}/{:03d}/".format(l,n_cells,i)
                if not os.path.isfile(test_dir + "costs.txt"):
                    continue
                costs = np.loadtxt(test_dir + "q_values.txt")
                # print(costs)
                if not label_added:
                    plotter.plot_xy([i+1 for i in range(len(costs))],costs, label=r"$n_{T}$"+  "={}".format(n_cells),color = n_cells_color_map[n_cells],linewidth = lwmap[n_cells])
                    label_added = True
                else:
                    plotter.plot_xy([i+1 for i in range(len(costs))],costs,label = "_none",color = n_cells_color_map[n_cells],linewidth  = lwmap[n_cells])

        plotter.ax.legend(loc = "lower left")
        plotter.save_fig(savefile,transparent=False)


def main():
    error_to_iters_n_spheroid()
    overlap_to_iters_n_spheroid()
    error_to_iters_periodic_system_size()
    overlap_to_iters_periodic_system_size()
    
if __name__ == "__main__":
    main()