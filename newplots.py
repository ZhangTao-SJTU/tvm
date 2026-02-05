from toolbox import manuscriptPlots
import numpy as np
import pandas as pd
from scipy.stats import gamma
import os


def error_to_iters_n_spheroid():
    n_spheroid_color_map = {5:"blue",10:"orange", 20:"brown", 40:"red" }
    l = 5
    savefile = "spheroid_graphs/spheroid_error_to_iters.png"
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(1e-7,5)
    # plotter.set_ylim(0.6,1)

    plotter.set_xlim(1,2000)
    plotter.set_xticks([20*i for i in range(150)])
    plotter.set_yticks([5*i for i in range(1,100)])
    plotter.set_yticks([0.2*i for i in range(1,100)])

    plotter.set_xlabel("Iterations")
    plotter.set_ylabel(r"$<|1-\sigma_{target}/\sigma|>$")
    # plotter.set_ylabel(r"$Q_n$")

    # plotter.set_title(r"$n_{(target)} =$"+ "{}".format(n_target))
    plotter.set_yScaled()

    plotter.initialize_figure()
    plotter.ax.set_xscale("log")
    plotter.ax.set_yscale("log")
    for n_spheroid in [5,10,20,40]:
        label_added = False
        for i in range(4):
            test_dir = "runs/kv_10_l_{}_n_hidden_{:03d}/{:03d}/".format(l,n_spheroid,i)
            if not os.path.isfile(test_dir + "costs.txt"):
                continue
            costs = np.loadtxt(test_dir + "costs.txt")
            if not label_added:
                plotter.plot_xy([i+1 for i in range(len(costs))],costs, label=r"$n_{total}$"+  "={}".format(n_spheroid),color = n_spheroid_color_map[n_spheroid])
                label_added = True
            else:
                plotter.plot_xy([i+1 for i in range(len(costs))],costs,label = "_none",color = n_spheroid_color_map[n_spheroid])

    plotter.ax.legend(loc = "lower left")
    plotter.save_fig(savefile)

def overlap_to_iters_n_spheroid():
    n_spheroid_color_map = {5:"blue",10:"orange", 20:"brown", 40:"red" }
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
    for n_spheroid in [5,10,20,40]:
        label_added = False
        for i in range(4):
            test_dir = "runs/kv_10_l_{}_n_hidden_{:03d}/{:03d}/".format(l,n_spheroid,i)
            if not os.path.isfile(test_dir + "costs.txt"):
                continue
            costs = np.loadtxt(test_dir + "q_values.txt")
            if not label_added:
                plotter.plot_xy([i+1 for i in range(len(costs))],costs, label=r"$n_{total}$"+  "={}".format(n_spheroid),color = n_spheroid_color_map[n_spheroid])
                label_added = True
            else:
                plotter.plot_xy([i+1 for i in range(len(costs))],costs,label = "_none",color = n_spheroid_color_map[n_spheroid])

    plotter.ax.legend(loc = "lower left")
    plotter.save_fig(savefile)


def s0_histograms():
    
def main():
    error_to_iters_n_spheroid()

    overlap_to_iters_n_spheroid()
if __name__ == "__main__":
    main()