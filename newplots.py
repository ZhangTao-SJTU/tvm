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
l_to_marker = {4:'o', 5:"s", 6:"^"}
l_to_alpha = {4:1, 5: 0.7, 6:0.5}
l_to_color = {4:"black", 5: "purple", 6:"green"}
n_runs = 100

def find_complete_runs(dirlist):
    complete_dirs =[]
    for test_dir in dirlist:
        if not os.path.isfile(test_dir + "costs.txt"):
            continue
        costs = np.loadtxt(test_dir + "costs.txt")
        if costs.shape == ():
            continue
        if costs[-1]>tolerance:
            continue
        complete_dirs.append(test_dir)
    # print(complete_dirs)
    return complete_dirs

def single_tracks_to_iters( label_to_data,savefile,
                            input_filename = "costs.txt",
                            title = None,
                            ylim = [1e-8,5],
                            yticks =[5*i for i in range(1,100)],
                            ylabel =r"$<|1-\sigma_{T}/\sigma|>$",
                            ylog = True):
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(ylim[0],ylim[1])
    plotter.set_xlim(0.5,5000)
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
            y_vals = np.loadtxt(test_dir+input_filename)
            if not label_added:
                plotter.plot_xy([i+1 for i in range(len(y_vals))],y_vals,label=label,color = data.get("color","black"),alpha = data.get("alpha",0.2) )
                label_added = True
            else: 
                plotter.plot_xy([i+1 for i in range(len(y_vals))],y_vals,label="_none",color = data.get("color","black"),alpha = data.get("alpha",0.2) )
    plotter.ax.legend()
    plotter.save_fig(savefile)

def histogram(label_to_data,
              input_filename = "s0_vals.txt",
              title = None,
              xlim = [3,7],
              ylim = [0,4],
              xlabel = r"$s_0$",
              xticks = [i for i in range(3,8)],
              yticks = [i for i in range(1,10)]):
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(ylim[0],ylim[1])
    plotter.set_xlim(xlim[0],xlim[1])
    plotter.set_xticks(xticks)
    plotter.set_yticks(yticks)
    plotter.set_xlabel(xlabel)
    if title is not None:
        plotter.set_title(title)
    plotter.initialize_figure()
    for label, data in label_to_data.items():
        array = np.loadtxt(data["dir"]+input_filename)
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
    plotter.set_ylim(ylim[0],ylim[1])
    plotter.set_xlim(xlim[0],xlim[1])
    plotter.set_xticks(xticks)
    plotter.set_yticks(yticks)


    plotter.set_ylabel(ylabel)
    plotter.set_xlabel(xlabel)
    if title is not None:
        plotter.set_title(title)
    # plotter.set_yScaled()

    plotter.initialize_figure()
    # plotter.ax.set_yscale("log")
    if xlog:
        plotter.set_xLog()
    if ylog:
        plotter.set_yLog()

    for label,data in label_to_data.items():
        plotter.plot_scatter(
            data["x_array"],
            data["y_array"],
            marker=data.get("marker","o"),
            s=data.get("s",500),
            color= data.get("color","black"),
            alpha = data.get("alpha",0.2),
            label = label,
        )
    return plotter
    # plotter.ax.legend(ncol=2)
    # plotter.save_fig(savefile)

def error_to_iters_single_pattern():
    for l in [4,5,6]:
        for n in [2,4]:
            dirlist = find_complete_runs(["data/kv_10_l_{}_n_{}/{:03d}/".format(l,n,i) for i in range(100)])
            savefile = "Panel_1/error_to_iters_periodic_l_{}_n_{}.png".format(l,n)
            label_to_data = {"_none":{"dirlist":dirlist, "color":"black"}}
            title = r"$n_{total} = $"+"{}".format(l**3)
            single_tracks_to_iters(label_to_data=label_to_data,title=title,savefile=savefile)

def error_to_iters_spheroid():
    for l in [5,6]:
        savefile = "Panel_5/error_to_iters_spheroid_l_{}.png".format(l)
        label_to_data = {}
        title = r"$n_{total} = $"+"{}".format(l**3)

        for n in [5,10,40]:
            dirlist = find_complete_runs(["data/kv_10_l_{}_n_sp_{:03d}/{:03d}/".format(l,n,i) for i in range(100)])
        
            label_to_data[r"$n_{sp}=$"+"{}".format(n)]={"dirlist":dirlist, "color":n_spheroid_color_map[n]}
        single_tracks_to_iters(label_to_data=label_to_data,title=title,savefile=savefile)

def overlap_to_iters_single_pattern():
    for l in [4,5,6]:
        for n in [2,4]:
            dirlist = find_complete_runs(["data/kv_10_l_{}_n_{}/{:03d}/".format(l,n,i) for i in range(100)])
            savefile = "Panel_1/overlap_to_iters_periodic_l_{}_n_{}.png".format(l,n)
            label_to_data = {"_none":{"dirlist":dirlist, "color":"black", "bins":30}}
            title = r"$n_{total} = $"+"{}".format(l**3)
            single_tracks_to_iters(label_to_data=label_to_data,title=title,savefile=savefile,
                                   input_filename="q_values.txt",
                                   ylabel = r"$Q_2$",
                                   ylog=False,
                                   ylim=[0.6,1.02],
                                   yticks=[0.6,0.8,1])

def s0_histogram_periodic():
    input_filename = "final_s0.txt"
    for l in [4,5,6]:
        for n in [2,4]:
            dir = "data/kv_10_l_{}_n_{}/".format(l,n)
            savefile = "Panel_1/s0_histogram_periodic_l_{}_n_{}.png".format(l,n)
            label_to_data = {r"$s_0^{(trained)}$":{"dir":dir, "color":"red", "bins":35}}
            title = r"$n_{total} = $"+"{}".format(l**3)
            hist = histogram(   label_to_data=label_to_data,
                                input_filename=input_filename,
                                title=title,
                                xlim = [3.5,6.5],
                                ylim = [0,3.5])
            hist.ax.vlines(x = 5, ymin = 0, ymax = 2.6, linestyle= "dashed",color = "black", label = r"$s_0^{(init)}$",linewidth =15)
            hist.ax.legend()
            hist.save_fig(savefile)

def s0_histogram_periodic_combined():
    input_filename = "final_s0.txt"
    for l in [4,5,6]:
        label_to_data = {}
        savefile = "Panel_2/s0_histogram_periodic_l_{}.png".format(l)

        for n in [1,2,4]:
            dir = "data/kv_10_l_{}_n_{}/".format(l,n)
            label_to_data[r"$n_T=$"+"{}".format(n)] = {"dir":dir, "color":n_cells_color_map[n], "bins":30}
        title = r"$n_{total} = $"+"{}".format(l**3)
        hist = histogram(   label_to_data=label_to_data,
                            input_filename=input_filename,
                            title=title,
                            xlim = [3.5,6.5],
                            ylim = [0,3.5])
        hist.ax.vlines(x = 5, ymin = 0, ymax = 3.2, linestyle= "dashed",color = "black", label = r"$s_0^{(init)}$",linewidth =15)
        hist.ax.legend()
        hist.save_fig(savefile)


def stress_histogram_periodic():
    input_filename = "final_stresses.txt"
    for l in [4,5,6]:
        for n in [2,4]:
            dir = "data/kv_10_l_{}_n_{}/".format(l,n)
            savefile = "Panel_1/stress_histogram_periodic_l_{}_n_{}.png".format(l,n)
            # label_to_data = {r"$\sigma_T = <\sigma^{(init)}>$":{"dir":dir, "color":"red", "bins":20}}
            label_to_data = {r"$\sigma^{(trained)}$":{"dir":dir, "color":"red", "bins":20}}

            title = r"$n_{total} = $"+"{}".format(l**3)
            hist = histogram(   label_to_data=label_to_data,
                                input_filename=input_filename,
                                title = title,
                                xlabel = r"$\sigma$",
                                xlim = [0,0.8],
                                xticks=[0.2*i for i in range(5)],
                                ylim = [0,5])
            hist.ax.vlines(x = np.mean(np.loadtxt(dir+"final_stresses.txt")), ymin = 0, ymax = 4.5, linestyle= "dotted",color = "black", label = r"$\sigma_T = <\sigma^{(init)}>$",linewidth =15)
            hist.ax.legend()
            hist.save_fig(savefile)

def final_iteration_to_final_overlap_periodic():
    for l in [4,5,6]:  
        savefile = "Panel_2/final_iterations_to_final_overlap_l_{}.png".format(l)
        label_to_dir = {}
        # for n in [1,2,3,4,5,6]:
        for n in [1,2,4]:

            dirlist = find_complete_runs(["data/kv_10_l_{}_n_{}/{:03d}/".format(l,n,i) for i in range(100)])
            x_array = [len(np.loadtxt(dir+"costs.txt"))for dir in dirlist]
            y_array = [np.loadtxt(dir+"q_values.txt")[-1] for dir in dirlist]
            label_to_dir[r"$n_T=$"+"{}".format(n)]= {"x_array":x_array,"y_array":y_array,"marker":"o","color":n_cells_color_map[n]}
        plotter = scatter_plot(label_to_dir,
                               ylim=[0.4,1.05],
                                title = r"$n_{total}=$"+"{}".format(l**3),
)
        plotter.ax.legend()
        plotter.save_fig(savefile)

def SD_s0_scatter_periodic():
    for l in [4,5,6]:  
        savefile = "Panel_2/SD_s0_scatter_l_{}.png".format(l)
        label_to_dir = {}
        # for n in [1,2,3,4,5,6]:
        #     dirlist = find_complete_runs(["data/kv_10_l_{}_n_{}/{:03d}/".format(l,n,i) for i in range(100)])
        x_array = [1,2,3,4,5,6]
        y_array = [np.std(np.loadtxt("data/kv_10_l_{}_n_{}/final_s0.txt".format(l,n))) for n in x_array]
        label_to_dir= {"_none": {"x_array":x_array,"y_array":y_array,"marker":"o","color":"black"}}
        plotter = scatter_plot(label_to_dir,
                     xlim =[0.5,6.5],
                     ylim=[0,0.5],
                     title = r"$n_{total}=$"+"{}".format(l**3),
                     xticks= x_array,
                     xlog=False,
                     xlabel=r"$n_T$",
                     ylabel=r"$SD(s_0)$")
        plotter.save_fig(savefile)
        
def SD_s0_scatter_periodic_all_sizes():
    savefile = "Panel_2/SD_s0_scatter.png"
    label_to_dir = {}
    for l in [4,5,6]:  
        
        # for n in [1,2,3,4,5,6]:
        #     dirlist = find_complete_runs(["data/kv_10_l_{}_n_{}/{:03d}/".format(l,n,i) for i in range(100)])
        x_array = [1,2,3,4,5,6]
        y_array = [np.std(np.loadtxt("data/kv_10_l_{}_n_{}/final_s0.txt".format(l,n))) for n in x_array]
        label_to_dir[r"$n_{total} = $"+"{}".format(l**3)] = {"x_array":x_array,"y_array":y_array,"marker":l_to_marker[l],"color":"black","alpha": l_to_alpha[l], "color":l_to_color[l]}
    plotter = scatter_plot(label_to_dir,
                    xlim =[0.5,6.5],
                    ylim=[0,0.5],
                    # title = r"$n_{total}=$"+"{}".format(l**3),
                    xticks= x_array,
                    xlog=False,
                    xlabel=r"$n_T$",
                    ylabel=r"$SD(s_0)$")
    plotter.ax.legend()
    plotter.save_fig(savefile)
# def error_to_iters_spheroid(l_vals = [5,6], n_spheroid = [5,10,20,40]):
#     for l in l_vals:
#         savefile = "spheroid_graphs/error_to_iters_spheroid_l_{}.png".format(l)
#         plotter = manuscriptPlots.plot()
#         plotter.set_ylim(1e-8,5)
#         plotter.set_xlim(1,5000)
#         plotter.set_xticks([20*i for i in range(150)])
#         plotter.set_yticks([5*i for i in range(1,100)])
#         plotter.set_yticks([0.2*i for i in range(1,100)])
#         plotter.set_xlabel("Iterations")
#         plotter.set_ylabel(r"$<|1-\sigma_{T}/\sigma|>$")
#         plotter.set_title(r"$n_{(total)} =$"+ "{}".format(l**3))
#         plotter.set_yScaled()
#         plotter.initialize_figure()
#         plotter.ax.set_xscale("log")
#         plotter.ax.set_yscale("log")
#         for n in n_spheroid:
#             label_added = False
#             for i in range(100):
#                 test_dir = "data/kv_10_l_{}_n_sp_{:03d}/{:03d}/".format(l,n,i)
#                 if not os.path.isfile(test_dir + "costs.txt"):
#                     continue
#                 costs = np.loadtxt(test_dir + "costs.txt")
#                 if costs.shape == ():
#                     continue
#                 if costs[-1]>tolerance:
#                     continue
#                 if not label_added:
#                     plotter.plot_xy([i+1 for i in range(len(costs))],costs,label=r"$n_{sp}$"+  "={}".format(n),color = n_spheroid_color_map[n], linewidth=lw_sp_map[n],alpha = 0.5)
#                     label_added = True
#                 else:
#                     plotter.plot_xy([i+1 for i in range(len(costs))],costs,label = "_none",color = n_spheroid_color_map[n], linewidth=lw_sp_map[n],alpha = 0.5)

#         plotter.ax.legend(loc = "lower left")
#         plotter.save_fig(savefile)

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

def final_iteration_to_final_overlap(l_vals = [5,6], n_cells = [5,10,20,40],spheroid = True):
    n_runs = 100
    for l in l_vals:
        if spheroid:
            savefile = "spheroid_graphs/final_iterations_to_final_overlap_spheroid_l_{}.png".format(l)
        else:
            savefile = "spheroid_graphs/final_iterations_to_final_overlap_l_{}.png".format(l)
        plotter = manuscriptPlots.plot()
        plotter.set_ylim(0.5,1.1)
        plotter.set_xlim(10,10000)
        # plotter.set_xticks([500*i for i in range(10)])
        plotter.set_xticks([0.5,1])

        plotter.set_ylabel(r"$q_2$")
        plotter.set_xlabel(r"$i_{final}$")
        plotter.set_title(r"$n_{(total)} =$"+ "{}".format(l**3))
        # plotter.set_yScaled()

        plotter.initialize_figure()
        # plotter.ax.set_yscale("log")
        plotter.ax.set_xscale("log")


        for t in n_cells:
            i_final = []
            q_final = []
            for i in range(n_runs):
                if spheroid:
                    test_dir = "data/kv_10_l_{}_n_sp_{:03d}/{:03d}/".format(l,t,i)
                else:
                    test_dir = "data/kv_10_l_{}_n_{}/{:03d}/".format(l,t,i)
                if not os.path.isfile(test_dir + "costs.txt"):
                    continue
                costs = np.loadtxt(test_dir + "costs.txt")
                if costs.shape == ():
                    continue
                if costs[-1]>tolerance:
                    continue
                i_final.append(len(costs))
                q_final.append(np.loadtxt(test_dir + "q_values.txt")[-1])


            if spheroid:
                plotter.plot_scatter(
                    i_final,
                    q_final,
                    marker='d',
                    markersize=30,
                    color=n_spheroid_color_map[t],
                    label = r"$n_{sp}=$"+"{}".format(t),
                )
            else:
                plotter.plot_scatter(
                    i_final,
                    q_final,
                    marker='d',
                    s=200,
                    color=n_cells_color_map[t],
                    label = r"$n_T=$"+"{}".format(t),
                )
        plotter.ax.legend()
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

# def s0_histogram_periodic(l_vals = [4,5,6],n_target_cells = [1,2,3,4,5,6]):
#     n_runs = 100

#     for l in l_vals:
#         savefile = "spheroid_graphs/final_s0_histogram_l_{}.png".format(l)
        
#         plotter = manuscriptPlots.plot()
#         plotter.set_ylim(0,4)
#         plotter.set_xlim(3,7)
#         plotter.set_xticks([i for i in range(3,8)])
#         plotter.set_yticks([i for i in range(1,10)])
#         plotter.set_xlabel(r"$s_0^{hidden}$")
#         plotter.set_title(r"$n_{(total)} = $"+ "{}".format(l**3))

#         plotter.initialize_figure()
#         plotter.ax.vlines(x = 5, ymin = 0, ymax = 3, linestyle= "dashed",color = plotter.neutralColor, label = r"$s_0^{(init)}$")

#         for t in n_target_cells:
#             s0_vals = []
#             for i in range(n_runs):
#                 test_dir = "data/kv_10_l_{}_n_{}/{:03d}/".format(l,t,i)
#                 testfile = test_dir+"cellParameters.input"
#                 if not os.path.isfile(testfile):
#                     continue
#                 if not os.path.isfile(test_dir+"costs.txt"):
#                     continue
#                 costs = np.loadtxt(test_dir + "costs.txt")
#                 if costs.shape == ():
#                     continue
#                 if costs[-1]>tolerance:
#                     continue
#                 s0_vals += pd.read_csv(testfile, sep = " ", header=None)[2].tolist()
#             plotter.histogram_from_array(s0_vals, fit_type= 'kde',bins = 30,label = r"$n_T = {}$".format(t),alpha = 0.2,color = n_cells_color_map[t])
#         plotter.ax.legend()
#         plotter.save_fig(savefile)

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

def calculate_surface_area(sample, spheroid_cells, frozen_cells):
    frozen_polygons = {
        polygon_id
        for cell_id in frozen_cells
        for polygon_id in sample.cells_[cell_id].polygons_
    }
    spheroid_polygons = {
        polygon_id
        for cell_id in spheroid_cells
        for polygon_id in sample.cells_[cell_id].polygons_
    }
    
    common_polygons = frozen_polygons & spheroid_polygons
    try:
        return sum(sample.polygons_[polygon_id].area_ for polygon_id in common_polygons)
    except:
        return None
def write_init_and_final_surface_areas(test_dir):
    if not os.path.isfile(test_dir+"minimized.txt"):
        return
    spheroid_cells = np.loadtxt(test_dir+"spheroid_cells.txt")
    frozen_cells = np.loadtxt(test_dir+"frozen_cells.txt")
    areas = [calculate_surface_area(PeriodicTissue.from_config(test_dir,file),spheroid_cells,frozen_cells) for file in ["initial_config.txt","minimized.txt"]]
    if None in areas:
        return
    np.savetxt(test_dir+"areas.txt",areas)

def write_all_areas():
    for l in [6]:
        for n in [5,10,20,40]:
            for i in range(100):
                test_dir ="data/kv_10_l_{}_n_sp_{:03d}/{:03d}/".format(l,n,i)
                write_init_and_final_surface_areas(test_dir)

def write_stresses_spheroid(l_vals = [5,6],n_spheroid_cells = [5,10,20,40]):
    n_runs = 100
    for l in l_vals:
        for t in n_spheroid_cells:
            stress_vals = []
            for i in range(n_runs):
                test_dir ="data/kv_10_l_{}_n_sp_{:03d}/{:03d}/".format(l,t,i)
                if not os.path.isfile(test_dir+"minimized.txt"):
                    continue
                if not os.path.isfile(test_dir + "cellParameters.input"):
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
                stress_vals += [calculate_max_shear_stress(tr._config,cellID) for cellID in spheroid_cells]
            if not len(stress_vals):
                continue
            np.savetxt("data/kv_10_l_{}_n_sp_{:03d}/final_stresses.txt".format(l,t),stress_vals)

def write_stresses_periodic(l_vals = [4,5,6],n_target_cells = [1,2,3,4,5,6]):
    n_runs = 100
    for l in l_vals:
        for t in n_target_cells:
            stress_vals = []
            for i in range(n_runs):
                test_dir ="data/kv_10_l_{}_n_{}/{:03d}/".format(l,t,i)
                if not os.path.isfile(test_dir+"minimized.txt"):
                    continue
                if not os.path.isfile(test_dir + "cellParameters.input"):
                    continue
                if not os.path.isfile(test_dir + "costs.txt"):
                    continue
                costs = np.loadtxt(test_dir + "costs.txt")
                if costs.shape == ():
                    continue
                if costs[-1]>tolerance:
                    continue
                sample = PeriodicTissue.from_config(test_dir,"minimized.txt")
                tr = Training.from_sample(sample)
                tr.load_cell_parameters()
                stress_vals += [calculate_max_shear_stress(tr._config,cellID) for cellID in tr._config.cells_]
            if not len(stress_vals):
                continue
            np.savetxt("data/kv_10_l_{}_n_{}/final_stresses.txt".format(l,t),stress_vals)

def write_final_s0(l_vals = [5,6],n_cells = [5,10,20,40], spheroid = True):
    for l in l_vals:
        for n in n_cells:
            s0_vals = []
            if spheroid: 
                savefile = "data/kv_10_l_{}_n_sp_{:03d}/final_s0.txt".format(l,n)
            else:
                savefile = "data/kv_10_l_{}_n_{}/final_s0.txt".format(l,n)
            for i in range(n_runs):
                if spheroid:
                    test_dir ="data/kv_10_l_{}_n_sp_{:03d}/{:03d}/".format(l,n,i)
                else:
                    test_dir ="data/kv_10_l_{}_n_{}/{:03d}/".format(l,n,i)
                if not os.path.isfile(test_dir+"minimized.txt"):
                    continue
                if not os.path.isfile(test_dir + "cellParameters.input"):
                    continue
                if not os.path.isfile(test_dir + "costs.txt"):
                    continue
                costs = np.loadtxt(test_dir + "costs.txt")
                if costs.shape == ():
                    continue
                if costs[-1]>tolerance:
                    continue
                if spheroid:
                    spheroid_cells = np.loadtxt(test_dir+"spheroid_cells.txt".format(l,n,i))
                    df = pd.read_csv(test_dir+"cellParameters.input", sep = " ", header=None)
                    for i,row in df.iterrows():
                        if not row[0] in spheroid_cells:
                            continue
                        s0_vals.append(row[2])
                else:
                    df = pd.read_csv(test_dir+"cellParameters.input", sep = " ", header=None)
                    s0_vals += df[2].to_list()

            if not len(s0_vals):
                continue
            np.savetxt(savefile,s0_vals)



def area_change_to_stress_change(l_vals = [5,6],n_spheroid_cells = [5,10,20,40]):
    for l in l_vals:
        savefile = "spheroid_graphs/area_change_to_stress_change_{}.png".format(l)
        plotter = manuscriptPlots.plot()
        plotter.set_xlim(-1,1)
        plotter.set_ylim(-.1,.1)
        plotter.set_xticks([0.5*i for i in range(-2,3)])
        plotter.set_yticks([0.1*i for i in range(-2,3)])
        plotter.set_xlabel(r"$(\sigma_T - \sigma_{0})/\sigma_{0}$")
        plotter.set_ylabel(r"$(S_{T} - S_{0})/S_{0}$")
        plotter.set_title("l = {}, ".format(l)+ r"$n_T = 1$")

        plotter.initialize_figure()
        for n in n_spheroid_cells:
            del_A = []
            del_stress = []
            for i in range(100):
                test_dir ="data/kv_10_l_{}_n_sp_{:03d}/{:03d}/".format(l,n,i)
                if not os.path.isfile(test_dir + "areas.txt"):
                    continue
                if not os.path.isfile(test_dir + "0000000.stresses.csv"):
                    continue
                areas = np.loadtxt(test_dir + "areas.txt")
                df = pd.read_csv(test_dir + "0000000.stresses.csv", sep = ",")
                del_A.append(areas[1]/areas[0]-1)
                del_stress.append(df["Target"].to_list()[0]/df["Current"].to_list()[0] - 1)
            plotter.plot_scatter(del_stress,del_A, color = n_spheroid_color_map[n], label = r"$n_{sp}=$"+ "{}".format(n))
        plotter.ax.legend()
        plotter.save_fig(savefile)

# def SD_s0_to_n_cells(l_vals = [4,5,6], n_cells = [1,2,3,4,5,6], spheroid = False):
#     for l in l_vals:
#         if spheroid:
#             savefile = "spheroid_graphs/sd_s0_to_n_spheroid_l_{}.png".format(l)
#             n_cell_to_final_s0 = {n:np.loadtxt("data/kv_10_l_{}_n_sp_{:03d}/final_s0.txt".format(l,n)) for n in n_cells}
#         else:
#             savefile = "spheroid_graphs/sd_s0_to_n_cells_l_{}.png".format(l)
#             n_cell_to_final_s0 = {n:np.loadtxt("data/kv_10_l_{}_n_{}/final_s0.txt".format(l,n)) for n in n_cells}
#         plotter = manuscriptPlots.plot()
#         plotter.set_xticks(n_cells)
#         if spheroid:
#             plotter.set_yticks([0.5*i for i in range(1,10)])
#             plotter.set_xlim(3,42)
#             plotter.set_ylim(0.2,1.3)
#             plotter.set_xlabel(r"$n_{sp}$")
#         else:
#             plotter.set_yticks([0.5*i for i in range(1,10)])
#             plotter.set_xlim(0,7)
#             plotter.set_ylim(0.1,0.5)
#             plotter.set_xlabel(r"$n_T$")
#         plotter.set_ylabel(r"$SD(s_0)$")
#         # plotter.set_ylabel(r"$Q_n$")

#         plotter.set_title(r"$n_{(total)} =$"+ "{}".format(l**3))
#         plotter.set_yScaled()
#         plotter.initialize_figure()
#         print([np.std(n_cell_to_final_s0[n]) for n in n_cells])
#         plotter.plot_scatter(
#             n_cells,
#             [np.std(n_cell_to_final_s0[n]) for n in n_cells],
#             marker = 'd',
#             s = 500,
#             color = "black")
#         plotter.save_fig(savefile)

def main():
    # stress_histogram_spheroid()
    # write_all_areas()
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
    # area_change_to_stress_change()
    # write_stresses_periodic()
    # write_final_s0(l_vals = [5,6],n_cells = [5,10,20,40], spheroid=True)
    # write_final_s0(l_vals = [4,5,6],n_cells = [1,2,3,4,5,6], spheroid=False)
    # SD_s0_to_n_cells(l_vals = [4,5,6], n_cells = [1,2,3,4,5,6], spheroid = False)
    # SD_s0_to_n_cells(l_vals = [5,6], n_cells = [5,10,20,40], spheroid = True)
    # final_iteration_to_final_overlap(l_vals=[4,5,6],n_cells = [1,2,4,6],spheroid=False)
    
    # Panel 1
    # error_to_iters_single_pattern()
    # overlap_to_iters_single_pattern()
    # s0_histogram_periodic()
    # stress_histogram_periodic()

    #Panel 2
    final_iteration_to_final_overlap_periodic()
    # SD_s0_scatter_periodic()
    # s0_histogram_periodic_combined()
    SD_s0_scatter_periodic_all_sizes()

    #Panel 5
    # error_to_iters_spheroid()
if __name__ == "__main__":
    main()