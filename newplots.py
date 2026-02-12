from toolbox import manuscriptPlots
from toolbox.stress import calculate_max_shear_stress
import numpy as np
import pandas as pd
from scipy.stats import gamma
from scipy.stats import sem
import os
from toolbox.spheroid import Spheroid
from toolbox.periodic import PeriodicTissue
from toolbox.training import Training
tolerance = 1e-6
lwmap = {1:15,2:10,3:5,4:5,5:5,6:5}
lw_sp_map = {40:15,20:10,10:5, 5:20}
n_spheroid_color_map = {40:"blue", 20:"orange", 10:"red", 5: "green" }
n_cells_color_map = {1:"blue",2:"orange",3:"purple", 4:"red" ,5:"yellow",6:"green"}

def error_to_iters_spheroid(l_vals = [5,6], n_spheroid = [5,10,20,40]):
    for l in l_vals:
        savefile = "spheroid_graphs/error_to_iters_spheroid_l_{}.png".format(l)
        plotter = manuscriptPlots.plot()
        plotter.set_ylim(1e-8,5)
        plotter.set_xlim(1,5000)
        plotter.set_xticks([20*i for i in range(150)])
        plotter.set_yticks([5*i for i in range(1,100)])
        plotter.set_yticks([0.2*i for i in range(1,100)])
        plotter.set_xlabel("Iterations")
        plotter.set_ylabel(r"$<|1-\sigma_{T}/\sigma|>$")
        plotter.set_title(r"$n_{(total)} =$"+ "{}".format(l**3))
        plotter.set_yScaled()
        plotter.initialize_figure()
        plotter.ax.set_xscale("log")
        plotter.ax.set_yscale("log")
        for n in n_spheroid:
            label_added = False
            for i in range(100):
                test_dir = "data/kv_10_l_{}_n_sp_{:03d}/{:03d}/".format(l,n,i)
                if not os.path.isfile(test_dir + "costs.txt"):
                    continue
                costs = np.loadtxt(test_dir + "costs.txt")
                if costs.shape == ():
                    continue
                if costs[-1]>tolerance:
                    continue
                if not label_added:
                    plotter.plot_xy([i+1 for i in range(len(costs))],costs,label=r"$n_{sp}$"+  "={}".format(n),color = n_spheroid_color_map[n], linewidth=lw_sp_map[n],alpha = 0.5)
                    label_added = True
                else:
                    plotter.plot_xy([i+1 for i in range(len(costs))],costs,label = "_none",color = n_spheroid_color_map[n], linewidth=lw_sp_map[n],alpha = 0.5)

        plotter.ax.legend(loc = "lower left")
        plotter.save_fig(savefile)

def overlap_to_iters_spheroid(l_vals = [5,6], n_spheroid = [5,10,20,40]):
    for l in l_vals:
        savefile = "spheroid_graphs/spheroid_overlap_to_iters.png"
        plotter = manuscriptPlots.plot()
        # plotter.set_ylim(1e-7,5)
        plotter.set_ylim(0.5,1.05)
        plotter.set_xlim(1,5000)
        plotter.set_xticks([20*i for i in range(150)])
        plotter.set_yticks([5*i for i in range(1,100)])
        plotter.set_yticks([0.2*i for i in range(1,100)])
        plotter.set_xlabel("Iterations")
        plotter.set_ylabel(r"$Q_2$")
        plotter.set_title(r"$n_{(total)} =$"+ "{}".format(l**3))
        plotter.initialize_figure()
        plotter.ax.set_xscale("log")
        # plotter.ax.set_yscale("log")
        for n in n_spheroid:
            label_added = False
            for i in range(4):
                test_dir = "data/kv_10_l_{}_n_sp_{:03d}/{:03d}/".format(l,n,i)
                if not os.path.isfile(test_dir + "costs.txt"):
                    continue
                costs = np.loadtxt(test_dir + "costs.txt")
                if costs.shape == ():
                    continue
                if costs[-1]>tolerance:
                    continue
                q_values = np.loadtxt(test_dir + "q_values.txt")
                if q_values.shape == ():
                    continue
                if not label_added:
                    plotter.plot_xy([i+1 for i in range(len(q_values))],q_values, label=r"$n_{sp}$"+  "={}".format(n),color = n_spheroid_color_map[n], linewidth=lw_sp_map[n])
                    label_added = True
                else:
                    plotter.plot_xy([i+1 for i in range(len(q_values))],q_values,label = "_none",color = n_spheroid_color_map[n], linewidth=lw_sp_map[n])
        plotter.ax.legend(loc = "lower left")
        plotter.save_fig(savefile)

def error_to_iters_periodic(l_vals = [4,5,6], n_target_cells = [1,2,3,4,5,6]):
    n_runs = 100
    for l in l_vals:
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

        plotter.set_title(r"$n_{(total)} =$"+ "{}".format(l**3))
        plotter.set_yScaled()

        plotter.initialize_figure()
        plotter.ax.set_xscale("log")
        plotter.ax.set_yscale("log")
        for n_cells in n_target_cells:
            label_added = False
            for i in range(n_runs):
                test_dir = "data/kv_10_l_{}_n_{}/{:03d}/".format(l,n_cells,i)
                if not os.path.isfile(test_dir + "costs.txt"):
                    continue
                costs = np.loadtxt(test_dir + "costs.txt")
                if costs.shape == ():
                    continue
                if costs[-1]>tolerance:
                    continue
                if not label_added:
                    plotter.plot_xy([i+1 for i in range(len(costs))],costs, label=r"$n_{T}$"+  "={}".format(n_cells),color = n_cells_color_map[n_cells],linewidth  = lwmap[n_cells])
                    label_added = True
                else:
                    plotter.plot_xy([i+1 for i in range(len(costs))],costs,label = "_none",color = n_cells_color_map[n_cells],linewidth  = lwmap[n_cells])

        plotter.ax.legend(loc = "lower left")
        plotter.save_fig(savefile)

def overlap_to_iters_periodic(l_vals = [4,5,6], n_target_cells = [1,2,3,4,5,6]):
    n_runs = 100
    for l in l_vals:
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

        plotter.set_title(r"$n_{(total)} =$"+ "{}".format(l**3))
        plotter.set_yScaled()

        plotter.initialize_figure()
        plotter.ax.set_xscale("log")
        # plotter.ax.set_yscale("log")
        for n_cells in n_target_cells:
            label_added = False
            for i in range(n_runs):
                test_dir = "data/kv_10_l_{}_n_{}/{:03d}/".format(l,n_cells,i)
                if not os.path.isfile(test_dir + "costs.txt"):
                    continue
                costs = np.loadtxt(test_dir + "costs.txt")
                if costs.shape == ():
                    continue
                if costs[-1]>tolerance:
                    continue

                q_values = np.loadtxt(test_dir + "q_values.txt")
                if q_values.shape == ():
                    continue
                if not label_added:
                    plotter.plot_xy([i+1 for i in range(len(q_values))],q_values, label=r"$n_{T}$"+  "={}".format(n_cells),color = n_cells_color_map[n_cells],linewidth = lwmap[n_cells])
                    label_added = True
                else:
                    plotter.plot_xy([i+1 for i in range(len(q_values))],q_values,label = "_none",color = n_cells_color_map[n_cells],linewidth  = lwmap[n_cells])

        plotter.ax.legend(loc = "lower left")
        plotter.save_fig(savefile)

def final_iterations_to_n_cells_spheroid(l_vals = [5,6], n_spheroid = [5,10,20,40]):
    n_runs = 100
    for l in l_vals:
        savefile = "spheroid_graphs/final_iterations_to_n_cells_spheroid_l_{}.png".format(l)
        n_cell_to_iters = {t:[] for t in n_spheroid}

        plotter = manuscriptPlots.plot()
        plotter.set_ylim(1,5e4)
        plotter.set_xlim(min(n_spheroid)-1, max(n_spheroid)+1)
        plotter.set_xticks(n_spheroid)
        plotter.set_xlabel(r"$n_{spheroid}$")
        plotter.set_ylabel(r"$i_{final}$")
        plotter.set_title(r"$n_{(total)} =$"+ "{}".format(l**3))
        plotter.set_yScaled()

        plotter.initialize_figure()
        plotter.ax.set_yscale("log")

        for t,iters in n_cell_to_iters.items():
            for i in range(n_runs):
                test_dir = "data/kv_10_l_{}_n_sp_{:03d}/{:03d}/".format(l,t,i)
                if not os.path.isfile(test_dir + "costs.txt"):
                    continue
                costs = np.loadtxt(test_dir + "costs.txt")
                if costs.shape == ():
                    continue
                if costs[-1]>tolerance:
                    continue
                iters.append(len(costs))

        mean_vals = np.array([np.median(iter) for iter in n_cell_to_iters.values()])
        min_vals  = np.array([min(iter) for iter in n_cell_to_iters.values()])
        max_vals  = np.array([max(iter) for iter in n_cell_to_iters.values()])

        plotter.plot_errorbar(
            n_spheroid,
            mean_vals,
            [mean_vals-min_vals, max_vals-mean_vals],
            capsize=20,
            capthick=10,
            marker='d',
            markersize=30,
            color="black",
            linestyle="--"
        )

        plotter.save_fig(savefile)

def final_overlap_to_n_cells_spheroid(l_vals = [5,6], n_spheroid = [5,10,20,40]):
    n_runs = 100
    for l in l_vals:
        savefile = "spheroid_graphs/final_overlap_to_n_cells_spheroid_l_{}.png".format(l)
        n_cell_to_final_overlaps = {t:[] for t in n_spheroid}

        plotter = manuscriptPlots.plot()
        plotter.set_ylim(0.5,1.05)
        plotter.set_xlim(min(n_spheroid)-1, max(n_spheroid)+1)
        plotter.set_xticks(n_spheroid)
        plotter.set_xlabel(r"$n_{spheroid}$")
        plotter.set_ylabel(r"$Q_2$")
        plotter.set_title(r"$n_{(total)} =$"+ "{}".format(l**3))
        plotter.set_yScaled()

        plotter.initialize_figure()

        for t,iters in n_cell_to_final_overlaps.items():
            for i in range(n_runs):
                test_dir = "data/kv_10_l_{}_n_sp_{:03d}/{:03d}/".format(l,t,i)
                if not os.path.isfile(test_dir + "costs.txt"):
                    continue
                costs = np.loadtxt(test_dir + "costs.txt")
                if costs.shape == ():
                    continue
                if costs[-1]>tolerance:
                    continue
                if not os.path.isfile(test_dir + "q_values.txt"):
                    continue
                q_values = np.loadtxt(test_dir + "q_values.txt")
                if q_values.shape == ():
                    continue
                iters.append(q_values[-1])

        mean_vals = np.array([np.median(iter) for iter in n_cell_to_final_overlaps.values()])
        min_vals  = np.array([min(iter) for iter in n_cell_to_final_overlaps.values()])
        max_vals  = np.array([max(iter) for iter in n_cell_to_final_overlaps.values()])

        plotter.plot_errorbar(
            n_spheroid,
            mean_vals,
            [mean_vals-min_vals, max_vals-mean_vals],
            capsize=20,
            capthick=10,
            marker='d',
            markersize=30,
            color="black",
            linestyle="--"
        )

        plotter.ax.hlines(y=1, xmin=min(n_spheroid)-1,
                          xmax=max(n_spheroid)+1,
                          linestyle="dotted",
                          color=plotter.neutralColor,
                          alpha=0.2)

        plotter.save_fig(savefile)

def final_iterations_to_n_cells_periodic(l_vals = [4,5,6], n_target_cells = [1,2,3,4,5,6]):
    n_runs = 100
    for l in l_vals:
        savefile = "spheroid_graphs/final_iteration_to_n_cells_l_{}.png".format(l)
        n_cell_to_iters = {t:[] for t in n_target_cells}
        plotter = manuscriptPlots.plot()
        # plotter.set_ylim(1e-8,5)
        plotter.set_ylim(10,5e3)

        plotter.set_xlim(0.5,6.5)
        plotter.set_xticks([i for i in range(10)])
        plotter.set_yticks([5*i for i in range(1,100)])
        # plotter.set_yticks([0.2*i for i in range(1,100)])

        plotter.set_xlabel(r"$n_T$")
        plotter.set_ylabel(r"$i_{final}$")
        # plotter.set_ylabel(r"$Q_n$")

        plotter.set_title(r"$n_{(total)} =$"+ "{}".format(l**3))
        plotter.set_yScaled()

        plotter.initialize_figure()
        plotter.ax.set_yscale("log")
        for t,iters in n_cell_to_iters.items():
            for i in range(n_runs):
                test_dir = "data/kv_10_l_{}_n_{}/{:03d}/".format(l,t,i)
                if not os.path.isfile(test_dir + "costs.txt"):
                    continue
                costs = np.loadtxt(test_dir + "costs.txt")
                if costs.shape == ():
                    continue
                if costs[-1]>tolerance:
                    continue
                iters.append(len(costs))
        # plotter.plot_xy(
        #     np.array(n_target_cells),
        #     np.array([np.mean(list(iter)) for iter in n_cell_to_iters.values()])
        # )
        mean_vals = np.array([np.median([iter]) for _,iter in n_cell_to_iters.items()])
        min_vals = np.array([min(list(iter)) for iter in n_cell_to_iters.values()])
        max_vals = np.array([max(list(iter)) for iter in n_cell_to_iters.values()])
        plotter.plot_errorbar(
            n_target_cells,
            mean_vals,
            [mean_vals-min_vals,max_vals-mean_vals],
            capsize=20,
            capthick=10,
            marker = 'd',
            markersize = 30,
            color = "black",
            linestyle = "--")
        # plotter.plot_max_min_fill(
        #     n_target_cells,
        #     min_vals,
        #     max_vals)
        
        # plotter.plot_errorfill(
        #     np.array(n_target_cells),
        #     np.array([np.mean(list(iter)) for iter in n_cell_to_iters.values()]),
        #     np.array([sem(list(iter)) for iter in n_cell_to_iters.values()]))
        # print(np.array([np.mean(list(iter)) for iter in n_cell_to_iters.values()]))
        # print(np.array([sem(list(iter)) for iter in n_cell_to_iters.values()]))


        plotter.save_fig(savefile)

def final_overlap_to_n_cells_periodic(l_vals = [4,5,6], n_target_cells = [1,2,3,4,5,6]):
    n_runs = 100
    for l in l_vals:
        savefile = "spheroid_graphs/final_overlap_to_n_cells_l_{}.png".format(l)
        n_cell_to_final_overlaps = {t:[] for t in n_target_cells}
        plotter = manuscriptPlots.plot()
        # plotter.set_ylim(1e-8,5)
        plotter.set_ylim(0.5,1.05)

        plotter.set_xlim(0.5,6.5)
        plotter.set_xticks([i for i in range(10)])
        plotter.set_yticks([0.2*i for i in range(1,100)])
        # plotter.set_yticks([0.2*i for i in range(1,100)])

        plotter.set_xlabel(r"$n_T$")
        plotter.set_ylabel(r"$Q_2$")
        # plotter.set_ylabel(r"$Q_n$")

        plotter.set_title(r"$n_{(total)} =$"+ "{}".format(l**3))
        plotter.set_yScaled()

        plotter.initialize_figure()
        # plotter.ax.set_yscale("log")
        for t,iters in n_cell_to_final_overlaps.items():
            for i in range(n_runs):
                test_dir = "data/kv_10_l_{}_n_{}/{:03d}/".format(l,t,i)
                if not os.path.isfile(test_dir + "costs.txt"):
                    continue
                costs = np.loadtxt(test_dir + "costs.txt")
                if costs.shape == ():
                    continue
                if costs[-1]>tolerance:
                    continue
                q_values = np.loadtxt(test_dir + "q_values.txt")
                if q_values.shape == ():
                    continue
                iters.append(q_values[-1])
        # plotter.plot_xy(
        #     np.array(n_target_cells),
        #     np.array([np.mean(list(iter)) for iter in n_cell_to_iters.values()])
        # )
        mean_vals = np.array([np.median([iter]) for _,iter in n_cell_to_final_overlaps.items()])
        err = np.array([np.std([iter]) for _,iter in n_cell_to_final_overlaps.items()])
        min_vals = np.array([min(list(iter)) for iter in n_cell_to_final_overlaps.values()])
        max_vals = np.array([max(list(iter)) for iter in n_cell_to_final_overlaps.values()])
        plotter.plot_errorbar(
            n_target_cells,
            mean_vals,
            [mean_vals-min_vals,max_vals-mean_vals],
            # err,
            capsize=20,
            capthick=10,
            marker = 'd',
            markersize = 30,
            color = "black",
            linestyle = "--")
        plotter.ax.hlines(y=1,xmin =0,xmax=7,linestyle = "dotted",color = plotter.neutralColor, alpha = 0.2)
        # plotter.plot_max_min_fill(
        #     n_target_cells,
        #     min_vals,
        #     max_vals)
        
        # plotter.plot_errorfill(
        #     np.array(n_target_cells),
        #     np.array([np.mean(list(iter)) for iter in n_cell_to_iters.values()]),
        #     np.array([sem(list(iter)) for iter in n_cell_to_iters.values()]))
        # print(np.array([np.mean(list(iter)) for iter in n_cell_to_iters.values()]))
        # print(np.array([sem(list(iter)) for iter in n_cell_to_iters.values()]))
        plotter.save_fig(savefile)

def s0_histogram_spheroid(l_vals = [5,6],n_spheroid_cells = [5,10,20,40]):
    n_runs = 100

    for l in l_vals:
        savefile = "spheroid_graphs/final_s0_histogram_spheroid_l_{}.png".format(l)
        plotter = manuscriptPlots.plot()
        plotter.set_ylim(0,4)
        plotter.set_xlim(3,7)
        plotter.set_xticks([i for i in range(3,8)])
        plotter.set_yticks([i for i in range(1,10)])
        plotter.set_xlabel(r"$s_0^{hidden}$")
        plotter.set_title("l = {}, ".format(l)+ r"$n_T = 1$")

        plotter.initialize_figure()
        plotter.ax.vlines(x = 5, ymin = 0, ymax = 3, linestyle= "dashed",color = plotter.neutralColor, label = r"$s_0^{(init)}$")

        for t in n_spheroid_cells:
            s0_vals = []
            for i in range(n_runs):
                test_dir ="data/kv_10_l_{}_n_sp_{:03d}/{:03d}/".format(l,t,i)
                testfile = test_dir+"cellParameters.input"
                if not os.path.isfile(testfile):
                    continue
                if not os.path.isfile(test_dir + "costs.txt"):
                    continue
                costs = np.loadtxt(test_dir + "costs.txt")
                if costs.shape == ():
                    continue
                if costs[-1]>tolerance:
                    continue
                df = pd.read_csv(testfile, sep = " ", header=None)
                spheroid_cells = np.loadtxt("data/kv_10_l_{}_n_sp_{:03d}/{:03d}/spheroid_cells.txt".format(l,t,i))
                for i,row in df.iterrows():
                    if not row[0] in spheroid_cells:
                        continue
                    s0_vals.append(row[2])
            if not len(s0_vals):
                continue
            plotter.histogram_from_array(s0_vals, fit_type= 'kde',bins = 30,label = r"$n_{sp}=$"+ "{}".format(t),alpha = 0.2,color = n_spheroid_color_map[t])
        plotter.ax.legend()
        plotter.save_fig(savefile)

def s0_histogram_periodic(l_vals = [4,5,6],n_target_cells = [1,2,3,4,5,6]):
    n_runs = 100

    for l in l_vals:
        savefile = "spheroid_graphs/final_s0_histogram_l_{}.png".format(l)
        
        plotter = manuscriptPlots.plot()
        plotter.set_ylim(0,4)
        plotter.set_xlim(3,7)
        plotter.set_xticks([i for i in range(3,8)])
        plotter.set_yticks([i for i in range(1,10)])
        plotter.set_xlabel(r"$s_0^{hidden}$")
        plotter.set_title(r"$n_{(total)} = $"+ "{}".format(l**3))

        plotter.initialize_figure()
        plotter.ax.vlines(x = 5, ymin = 0, ymax = 3, linestyle= "dashed",color = plotter.neutralColor, label = r"$s_0^{(init)}$")

        for t in n_target_cells:
            s0_vals = []
            for i in range(n_runs):
                test_dir = "data/kv_10_l_{}_n_{}/{:03d}/".format(l,t,i)
                testfile = test_dir+"cellParameters.input"
                if not os.path.isfile(testfile):
                    continue
                if not os.path.isfile(test_dir+"costs.txt"):
                    continue
                costs = np.loadtxt(test_dir + "costs.txt")
                if costs.shape == ():
                    continue
                if costs[-1]>tolerance:
                    continue
                s0_vals += pd.read_csv(testfile, sep = " ", header=None)[2].tolist()
            plotter.histogram_from_array(s0_vals, fit_type= 'kde',bins = 30,label = r"$n_T = {}$".format(t),alpha = 0.2,color = n_cells_color_map[t])
        plotter.ax.legend()
        plotter.save_fig(savefile)

def stress_histogram_spheroid(l_vals = [5,6],n_spheroid_cells = [5,10,20,40]):
    n_runs = 100

    for l in l_vals:
        savefile = "spheroid_graphs/final_stress_histogram_spheroid_l_{}.png".format(l)
        
        plotter = manuscriptPlots.plot()
        plotter.set_ylim(0,4)
        plotter.set_xlim(0,1)
        plotter.set_xticks([i for i in range(3,8)])
        plotter.set_yticks([i for i in range(1,10)])
        plotter.set_xlabel(r"$s_0^{hidden}$")
        plotter.set_title("l = {}, ".format(l)+ r"$n_T = 1$")

        plotter.initialize_figure()
        initial_stresses = np.loadtxt("init/kv_10_l_{}/stresses.txt".format(l))
        plotter.ax.vlines(x = np.mean(initial_stresses), ymin = 0, ymax = 3, linestyle= "dashed",color = plotter.neutralColor, label = r"$\sigma_{T}$")

        for t in n_spheroid_cells:
            stress_vals = []
            for i in range(n_runs):
                test_dir ="data/kv_10_l_{}_n_sp_{:03d}/{:03d}/".format(l,t,i)
                testfile = test_dir+"minimized.txt"
                if not os.path.isfile(testfile):
                    continue
                if not os.path.isfile(test_dir + "costs.txt"):
                    continue
                costs = np.loadtxt(test_dir + "costs.txt")
                if costs.shape == ():
                    continue
                if costs[-1]>tolerance:
                    continue
                spheroid_cells = np.loadtxt("data/kv_10_l_{}_n_sp_{:03d}/{:03d}/spheroid_cells.txt".format(l,t,i))
                sample = PeriodicTissue.from_config(test_dir,"minimized.txt")
                tr = Training.from_sample(sample)
                tr.load_cell_parameters()
                for cellID in spheroid_cells:
                    stress_vals.append(calculate_max_shear_stress(tr._config,cellID))

            if not len(stress_vals):
                continue
            np.savetxt("data/kv_10_l_{}_n_sp_{:03d}/final_stresses.txt".format(l,t),stress_vals)
            plotter.histogram_from_array(stress_vals, fit_type= 'kde',bins = 30,label = r"$n_{sp}=$"+ "{}".format(t),alpha = 0.2,color = n_spheroid_color_map[t])
        plotter.ax.legend()
        plotter.save_fig(savefile)



def main():
    stress_histogram_spheroid()
    # error_to_iters_spheroid()
    # overlap_to_iters_spheroid()
    # error_to_iters_periodic(n_target_cells=[1,2,4])
    # overlap_to_iters_periodic(n_target_cells=[1,2,4])
    # final_iterations_to_n_cells_spheroid()
    # final_overlap_to_n_cells_spheroid()
    # final_iterations_to_n_cells_periodic()
    # final_overlap_to_n_cells_periodic()
    # s0_histogram_spheroid()
    # s0_histogram_periodic(n_target_cells=[1,2,4])
    
if __name__ == "__main__":
    main()