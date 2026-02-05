from toolbox import manuscriptPlots
import numpy as np
import pandas as pd
from scipy.stats import gamma
import os

alpha_map ={4:0.2,5:0.5, 6:1}

def error_to_iter_single_tracks():
    for n_target in [4,2]:
        graph_output = "new_graphs/error_to_iters_single_tracks_n_target_{}.jpg".format(n_target)

        plotter = manuscriptPlots.plot()
        plotter.set_ylim(1e-6,5e0)
        # plotter.set_ylim(0.6,1)

        plotter.set_xlim(1,2500)
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
        for size in [6]:
            legend_added = False
            for i in range(100):
                dir = "data/{}_cells_mean_l_{}/{:03d}/".format(n_target,size,i)
                if not os.path.isfile(dir+"errors.txt"):
                    continue
                errors = np.loadtxt(dir+"errors.txt")
                if errors[-1]>1e-4:
                    continue
                # errors = np.loadtxt(dir+"q_values.txt")

                x_array = [i+1 for i in range(len(errors))]
                plotter.plot_xy(x_array, errors, color = "black", alpha = 0.4,label = "_none")
                if legend_added:
                    continue
                plotter.plot_xy(x_array, errors, color = "black", alpha = 0.4, label = r"$n_{(target)}=$"+r'${}$'.format(n_target))
                legend_added = True


        # plotter.plot_errorfill(x_array,df_decrease["mean"].to_numpy(),df_decrease["sem"].to_numpy(),color = "black", alpha = 0.4, label = r"$\sigma_{hidden}^{(initial)}$")

        plotter.ax.legend(loc='lower left')
        plotter.save_fig(graph_output)

def overlap_to_iter_single_tracks():
    for n_target in [4,2]:
        graph_output = "new_graphs/overlap_to_iters_single_tracks_n_target_{}.jpg".format(n_target)

        plotter = manuscriptPlots.plot()
        # plotter.set_ylim(5e-6,5e0)
        plotter.set_ylim(0.6,1)

        plotter.set_xlim(1,2500)
        plotter.set_xticks([20*i for i in range(150)])
        plotter.set_yticks([5*i for i in range(1,100)])
        plotter.set_yticks([0.2*i for i in range(1,100)])

        plotter.set_xlabel("Iterations")
        # plotter.set_ylabel(r"$(|1-\sigma_{target}/\sigma|)$")
        plotter.set_ylabel(r"$Q_n$")

        # plotter.set_title(r"$n_{(target)} =$"+ "{}".format(n_target))
        plotter.set_yScaled()

        plotter.initialize_figure()
        for size in [6]:
            legend_added = False
            for i in range(100):
                dir = "data/{}_cells_mean_l_{}/{:03d}/".format(n_target,size,i)
                if not os.path.isfile(dir+"errors.txt"):
                    continue
                errors = np.loadtxt(dir+"q_values.txt")
                # if errors[-1]>1e-4:
                #     continue
                # errors = np.loadtxt(dir+"q_values.txt")

                x_array = [i+1 for i in range(len(errors))]
                plotter.plot_xy(x_array, errors, color = "black", alpha = 0.4,label = "_none")
                if legend_added:
                    continue
                plotter.plot_xy(x_array, errors, color = "black", alpha = 0.4, label = r"$n_{(target)}=$"+r'${}$'.format(n_target))
                legend_added = True


        # plotter.plot_errorfill(x_array,df_decrease["mean"].to_numpy(),df_decrease["sem"].to_numpy(),color = "black", alpha = 0.4, label = r"$\sigma_{hidden}^{(initial)}$")
        plotter.ax.set_xscale("log")
        # plotter.ax.set_yscale("log")
        plotter.ax.legend(loc='lower left')
        plotter.save_fig(graph_output)
def error_to_iter_single_pattern():
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(0,1.08)
    plotter.set_xlim(1,3000)
    plotter.set_xticks([1000*i for i in range(150)])
    plotter.set_yticks([0.2*i for i in range(1,100)])
    plotter.set_xlabel("Epochs")
    plotter.set_ylabel(r"$Q_n$")

    plotter.initialize_figure()
    # legend_added = False
    # for i in range(100):
    #     dir = "data/4_cells_mean/{:03d}/".format(i)
    #     if not os.path.isfile(dir+"errors.txt"):
    #         continue
    #     q_values = np.loadtxt(dir+"q_values.txt")
    #     # print(dir, errors[-1])
    #     x_array = [i+1 for i in range(len(q_values))]
    #     plotter.plot_xy(x_array, q_values, color ="black", alpha = 0.1,label = "_none")
        
    #     if legend_added:
    #         continue
    #     plotter.plot_xy(x_array, q_values, color ="black", alpha = 0.1, label = r"$\sigma_{target}=\sigma_\mu$")
    #     legend_added = True

    colors = {4:"red",5:"blue",6:"yellow"}
    for size in [4,5,6]:
        legend_added = False
        for i in range(100):
            dir = "data/4_cells_mean_l_{}/{:03d}/".format(size,i)
            if not os.path.isfile(dir+"q_values.txt"):
                continue
            with open(dir+"q_values.txt","r") as f:
                lines = f.readlines()
                if len(lines)<2:
                    continue
            errors = np.loadtxt(dir+"q_values.txt")
            x_array = [i+1 for i in range(len(errors))]
            plotter.plot_xy(x_array, errors, color = colors[size], alpha = 0.4,label = "_none")
            if legend_added:
                continue
            plotter.plot_xy(x_array, errors, color = colors[size], alpha = 0.4, label = "L = {}".format(size))
            legend_added = True
    # plotter.plot_errorfill(x_array,df_decrease["mean"].to_numpy(),df_decrease["sem"].to_numpy(),color = "black", alpha = 0.4, label = r"$\sigma_{hidden}^{(initial)}$")
    # plotter.ax.set_xscale("log")
    # plotter.ax.set_yscale("log")
    plotter.ax.legend(loc='upper right')
    plotter.save_fig("graphs/4_cells_overlap_all_sizes.jpg")

def error_to_iter_fixed_n_target():
    n_target = 4
    graph_output = "new_graphs/error_to_iters_n_target_{}.jpg".format(n_target)
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(5e-6,5)
    plotter.set_xlim(1,2500)
    plotter.set_xticks([20*i for i in range(150)])
    plotter.set_yticks([0.2*i for i in range(1,100)])
    plotter.set_xlabel("Iterations")
    plotter.set_ylabel(r"$|1-\sigma_{target}/\sigma|$")
    plotter.set_title(r"$n_{(target)}=$" +"{}".format(n_target))
    plotter.set_yScaled()
    plotter.initialize_figure()
    for size in [4,5,6]:
        errors_dict = {}
        max_len = 0
        for i in range(100):
            dir = "data/{}_cells_mean_l_{}/{:03d}/".format(n_target,size,i)
            if not os.path.isfile(dir+"errors.txt"):
                continue
            # errors = np.loadtxt(dir+"q_values.txt")
            errors = np.loadtxt(dir+"errors.txt")
            if errors[-1]>1e-4:
                continue

            errors_dict[i] = errors
            if len(errors)>max_len:
                max_len = len(errors)

        iter_to_errors = {i: [] for i in range(max_len)}
        for _, errors_array in errors_dict.items():
            for iter_idx in range(len(errors_array)):
                iter_to_errors[iter_idx].append(errors_array[iter_idx])
        mean_errors = []
        max_errors = []
        min_errors = []
        for iter_idx in range(max_len):
            mean_errors.append(np.mean(iter_to_errors[iter_idx]))
            max_errors.append(np.max(iter_to_errors[iter_idx]))
            min_errors.append(np.min(iter_to_errors[iter_idx]))

        x_array = [i+1 for i in range(len(mean_errors))]
        plotter.plot_max_min_fill(x_array, max_errors, min_errors, color = "blue", alpha = alpha_map[size], label = r"$n_{total}$"+r"$={}^3$".format(size))

    plotter.ax.set_xscale("log")
    plotter.ax.set_yscale("log")
    plotter.ax.legend(loc='lower left')
    plotter.save_fig(graph_output)

def overlap_to_iter_fixed_n_target():
    n_target = 4
    graph_output = "new_graphs/overlap_to_iters_n_target_{}.jpg".format(n_target)
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(0,1)
    plotter.set_xlim(1,2500)
    plotter.set_xticks([20*i for i in range(150)])
    plotter.set_yticks([0.2*i for i in range(1,100)])
    plotter.set_xlabel("Iterations")
    plotter.set_ylabel(r"$Q_n$")
    plotter.set_title(r"$n_{(target)}=$" +"{}".format(n_target))
    plotter.set_yScaled()
    plotter.initialize_figure()
    for size in [4,5,6]:
        errors_dict = {}
        max_len = 0
        for i in range(100):
            dir = "data/{}_cells_mean_l_{}/{:03d}/".format(n_target,size,i)
            if not os.path.isfile(dir+"errors.txt"):
                continue
            errors = np.loadtxt(dir+"q_values.txt")
            errors_dict[i] = errors
            if len(errors)>max_len:
                max_len = len(errors)

        iter_to_errors = {i: [] for i in range(max_len)}
        for _, errors_array in errors_dict.items():
            for iter_idx in range(len(errors_array)):
                iter_to_errors[iter_idx].append(errors_array[iter_idx])
        mean_errors = []
        max_errors = []
        min_errors = []
        for iter_idx in range(max_len):
            mean_errors.append(np.mean(iter_to_errors[iter_idx]))
            max_errors.append(np.max(iter_to_errors[iter_idx]))
            min_errors.append(np.min(iter_to_errors[iter_idx]))

        x_array = [i+1 for i in range(len(mean_errors))]
        plotter.plot_max_min_fill(x_array, max_errors, min_errors, color = "blue", alpha = alpha_map[size], label = r"$n_{total}$"+r"$={}^3$".format(size))

    plotter.ax.set_xscale("log")
    # plotter.ax.set_yscale("log")
    plotter.ax.legend(loc='lower left')
    plotter.save_fig(graph_output)

def error_to_iter_fixed_n_total():
    colormap = {6:"red", 4:"blue", 2:"purple"}
    alpha = 0.5
    size = 6
    graph_output = "new_graphs/error_to_iters_n_total_{}.jpg".format(size)
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(1e-7,5)
    plotter.set_xlim(1,5000)
    plotter.set_xticks([20*i for i in range(150)])
    plotter.set_yticks([0.2*i for i in range(1,100)])
    plotter.set_xlabel("Iterations")
    plotter.set_ylabel(r"$|1-\sigma_{target}/\sigma|$")
    plotter.set_title(r"$n_{total}=$" +r"${}^3$".format(size))
    plotter.set_yScaled()
    plotter.initialize_figure()
    for n_target in [2,4,6]:
        if n_target == 4:
            alpha = 1
        else:
            alpha =0.5
        errors_dict = {}
        max_len = 0
        for i in range(100):
            dir = "data/{}_cells_mean_l_{}/{:03d}/".format(n_target,size,i)
            if not os.path.isfile(dir+"errors.txt"):
                continue
            # errors = np.loadtxt(dir+"q_values.txt")
            errors = np.loadtxt(dir+"errors.txt")
            if errors[-1]>1e-4:
                continue

            errors_dict[i] = errors
            if len(errors)>max_len:
                max_len = len(errors)

        iter_to_errors = {i: [] for i in range(max_len)}
        for _, errors_array in errors_dict.items():
            for iter_idx in range(len(errors_array)):
                iter_to_errors[iter_idx].append(errors_array[iter_idx])
        mean_errors = []
        max_errors = []
        min_errors = []
        for iter_idx in range(max_len):
            mean_errors.append(np.mean(iter_to_errors[iter_idx]))
            max_errors.append(np.max(iter_to_errors[iter_idx]))
            min_errors.append(np.min(iter_to_errors[iter_idx]))

        x_array = [i+1 for i in range(len(mean_errors))]
        plotter.plot_max_min_fill(x_array, max_errors, min_errors, color = colormap[n_target], alpha = alpha, label = r"$n_{(target)}$"+r"$={}$".format(n_target))

    plotter.ax.set_xscale("log")
    plotter.ax.set_yscale("log")
    plotter.ax.legend(loc='lower left')
    plotter.save_fig(graph_output)

def overlap_to_iter_fixed_n_total():
    colormap = {6:"red", 4:"blue", 2:"purple"}
    alpha = 0.5
    size = 6
    graph_output = "new_graphs/overlap_to_iters_n_total_{}.jpg".format(size)
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(0,1)
    plotter.set_xlim(1,5000)
    plotter.set_xticks([20*i for i in range(150)])
    plotter.set_yticks([0.2*i for i in range(1,100)])
    plotter.set_xlabel("Iterations")
    plotter.set_ylabel(r"$Q_n$")
    plotter.set_title(r"$n_{total}=$" +r"${}^3$".format(size))
    plotter.set_yScaled()
    plotter.initialize_figure()
    for n_target in [2,4,6]:
        if n_target == 4:
            alpha = 1
        else:
            alpha =0.5
        errors_dict = {}
        max_len = 0
        for i in range(100):
            dir = "data/{}_cells_mean_l_{}/{:03d}/".format(n_target,size,i)
            if not os.path.isfile(dir+"errors.txt"):
                continue
            errors = np.loadtxt(dir+"q_values.txt")
            # errors = np.loadtxt(dir+"errors.txt")
            # if errors[-1]>1e-4:
            #     continue

            errors_dict[i] = errors
            if len(errors)>max_len:
                max_len = len(errors)

        iter_to_errors = {i: [] for i in range(max_len)}
        for _, errors_array in errors_dict.items():
            for iter_idx in range(len(errors_array)):
                iter_to_errors[iter_idx].append(errors_array[iter_idx])
        mean_errors = []
        max_errors = []
        min_errors = []
        for iter_idx in range(max_len):
            mean_errors.append(np.mean(iter_to_errors[iter_idx]))
            max_errors.append(np.max(iter_to_errors[iter_idx]))
            min_errors.append(np.min(iter_to_errors[iter_idx]))

        x_array = [i+1 for i in range(len(mean_errors))]
        plotter.plot_max_min_fill(x_array, max_errors, min_errors, color = colormap[n_target], alpha = alpha, label = r"$n_{(target)}$"+r"$={}$".format(n_target))

    plotter.ax.set_xscale("log")
    # plotter.ax.set_yscale("log")
    plotter.ax.legend(loc='lower left')
    plotter.save_fig(graph_output)

def error_to_iter_multiple_patterns():
    for a,b in [(1,1),(1,2),(2,2)]:
        graph_output = "new_graphs/error_to_iters_{}_{}.jpg".format(a,b)
        plotter = manuscriptPlots.plot()
        plotter.set_ylim(1e-7,5)
        plotter.set_xlim(1,30000)
        plotter.set_xticks([20*i for i in range(150)])
        plotter.set_yticks([0.2*i for i in range(1,100)])
        plotter.set_xlabel("Iterations")
        plotter.set_ylabel(r"$(|\sigma_{target}-\sigma|)/\sigma$")
        plotter.set_yScaled()
        plotter.initialize_figure()
        colors = {2:"purple",3:"red",4:"blue",5:"yellow",6:"green"}
        alphas = {2:0.2,3:0.3,4:0.4,5:0.5,6:0.6}
        for size in [4,5,6]:
            errors_dict = {}
            max_len = 0
            for i in range(100):
                dir = "data/{}_{}_l_{}/{:03d}/".format(a,b,size,i)
                if not os.path.isfile(dir+"costs.txt"):
                    continue
                errors = np.loadtxt(dir+"costs.txt")
                if not len(errors):
                    continue
                if len(errors)<2:
                    continue
                if errors[-1]>1e-4:
                    continue

                errors_dict[i] = errors
                if len(errors)>max_len:
                    max_len = len(errors)

            iter_to_errors = {i: [] for i in range(max_len)}
            for _, errors_array in errors_dict.items():
                for iter_idx in range(len(errors_array)):
                    iter_to_errors[iter_idx].append(errors_array[iter_idx])
            mean_errors = []
            for iter_idx in range(max_len):
                mean_errors.append(np.mean(iter_to_errors[iter_idx]))
            x_array = [i+1 for i in range(len(mean_errors))]
            errors = mean_errors
            # for i, errors in errors_dict.items():
            max_errors = []
            min_errors = []
            for iter_idx in range(max_len):
                max_errors.append(np.max(iter_to_errors[iter_idx]))
                min_errors.append(np.min(iter_to_errors[iter_idx]))

            plotter.plot_max_min_fill(x_array, max_errors, min_errors, color = colors[size], alpha = alphas[size], label = r"$n_{cells}^{(total)}$"+"={}".format(size**3))

        plotter.ax.set_xscale("log")
        plotter.ax.set_yscale("log")
        plotter.ax.legend(title = r"$n_A = {}, n_B = {}$".format(a,b),title_fontsize = 48,loc='lower left')
        plotter.save_fig(graph_output)

def overlap_to_iter_multiple_patterns():
    for a,b in [(1,1),(1,2),(2,2)]:
        graph_output = "new_graphs/overlap_to_iters_{}_{}.jpg".format(a,b)
        plotter = manuscriptPlots.plot()
        plotter.set_ylim(0,1)
        plotter.set_xlim(1,30000)
        plotter.set_xticks([10000*i for i in range(10)])
        plotter.set_yticks([0.2*i for i in range(6)])
        plotter.set_xlabel("Iterations")
        plotter.set_ylabel(r"$Q_n$")
        plotter.set_yScaled()
        plotter.set_xScaled()
        plotter.initialize_figure()
        plotter.set_xLog()
        colors = {2:"purple",3:"red",4:"blue",5:"yellow",6:"green"}
        alphas = {2:0.2,3:0.3,4:0.4,5:0.5,6:0.6}
        for size in [4,5,6]:
            errors_dict = {}
            max_len = 0
            for i in range(100):
                dir = "data/{}_{}_l_{}/{:03d}/".format(a,b,size,i)
                if not os.path.isfile(dir+"q_values.txt"):
                    continue
                errors = np.loadtxt(dir+"q_values.txt")
                if not len(errors):
                    continue
                if len(errors)<2:
                    continue

                errors_dict[i] = errors
                if len(errors)>max_len:
                    max_len = len(errors)

            iter_to_errors = {i: [] for i in range(max_len)}
            for _, errors_array in errors_dict.items():
                for iter_idx in range(len(errors_array)):
                    iter_to_errors[iter_idx].append(errors_array[iter_idx])
            mean_errors = []
            for iter_idx in range(max_len):
                mean_errors.append(np.mean(iter_to_errors[iter_idx]))
            x_array = [i+1 for i in range(len(mean_errors))]
            errors = mean_errors
            # for i, errors in errors_dict.items():
            max_errors = []
            min_errors = []
            for iter_idx in range(max_len):
                max_errors.append(np.max(iter_to_errors[iter_idx]))
                min_errors.append(np.min(iter_to_errors[iter_idx]))

            plotter.plot_max_min_fill(x_array, max_errors, min_errors, color = colors[size], alpha = alphas[size], label = r"$n_{cells}^{(total)}$"+"={}".format(size**3))

        # plotter.ax.set_xscale("log")
        # plotter.ax.set_yscale("log")
        plotter.ax.legend(title = r"$n_A = {}, n_B = {}$".format(a,b),title_fontsize = 48,loc='lower left')
        plotter.save_fig(graph_output)

def error_to_epoch_multiple_patterns():
    for a,b in [(1,1),(1,2),(2,2)]:
        graph_output = "new_graphs/error_to_epochs_{}_{}.jpg".format(a,b)
        plotter = manuscriptPlots.plot()
        plotter.set_ylim(1e-9,5)
        plotter.set_xlim(0.5,2000)
        plotter.set_xticks([1000*i for i in range(6)])
        plotter.set_yticks([0.2*i for i in range(6)])
        plotter.set_xlabel("Epochs")
        plotter.set_ylabel(r"$|1-\sigma_{target}/\sigma|$")
        plotter.set_yScaled()
        plotter.set_xScaled()
        plotter.set_title(r"$n_A = {}, n_B = {}$".format(a,b))

        plotter.initialize_figure()
        plotter.set_xLog()
        plotter.set_yLog()
        for size in [4,5,6]:
            errors_dict = {}
            max_len = 0
            for i in range(100):
                dir = "data/{}_{}_l_{}/{:03d}/".format(a,b,size,i)
                if not os.path.isfile(dir+"info.csv"):
                    continue
                errors = pd.read_csv(dir+"info.csv")["Error"].to_numpy()

                if not len(errors):
                    continue
                if len(errors)<2:
                    continue

                errors_dict[i] = errors
                if len(errors)>max_len:
                    max_len = len(errors)

            iter_to_errors = {i: [] for i in range(max_len)}
            for _, errors_array in errors_dict.items():
                for iter_idx in range(len(errors_array)):
                    iter_to_errors[iter_idx].append(errors_array[iter_idx])
            mean_errors = []
            for iter_idx in range(max_len):
                mean_errors.append(np.mean(iter_to_errors[iter_idx]))
            x_array = [(i+1)/2 for i in range(len(mean_errors))]
            max_errors = []
            min_errors = []
            for iter_idx in range(max_len):
                max_errors.append(np.max(iter_to_errors[iter_idx]))
                min_errors.append(np.min(iter_to_errors[iter_idx]))

            plotter.plot_max_min_fill(x_array, max_errors, min_errors, color = "blue", alpha = alpha_map[size], label = r"$n_{total}$"+r"$={}^3$".format(size))
        plotter.ax.legend(loc='lower left')
        plotter.save_fig(graph_output)

def distance_to_epoch_multiple_patterns():
    for a,b in [(1,1),(1,2),(2,2)]:
        graph_output = "new_graphs/distance_to_epochs_{}_{}.jpg".format(a,b)
        plotter = manuscriptPlots.plot()
        plotter.set_ylim(1e-2,50)
        plotter.set_xlim(0.5,2500)
        plotter.set_xticks([1000*i for i in range(6)])
        plotter.set_yticks([0.2*i for i in range(6)])
        plotter.set_xlabel("Epochs")
        plotter.set_ylabel(r"$|\vec{s}_{0,hidden}-\vec{s}_{0,hidden}^{(init)}|$")
        plotter.set_yScaled()
        plotter.set_xScaled()
        plotter.set_title(r"$n_A = {}, n_B = {}$".format(a,b))

        plotter.initialize_figure()
        plotter.set_xLog()
        plotter.set_yLog()
        for size in [4,5,6]:
            errors_dict = {}
            max_len = 0
            for i in range(100):
                dir = "data/{}_{}_l_{}/{:03d}/".format(a,b,size,i)
                if not os.path.isfile(dir+"info.csv"):
                    continue
                errors = pd.read_csv(dir+"info.csv")["Distance"].to_numpy()

                if not len(errors):
                    continue
                if len(errors)<2:
                    continue

                errors_dict[i] = errors
                if len(errors)>max_len:
                    max_len = len(errors)

            iter_to_errors = {i: [] for i in range(max_len)}
            for _, errors_array in errors_dict.items():
                for iter_idx in range(len(errors_array)):
                    iter_to_errors[iter_idx].append(errors_array[iter_idx])
            mean_errors = []
            for iter_idx in range(max_len):
                mean_errors.append(np.mean(iter_to_errors[iter_idx]))
            x_array = [(i+1)/2 for i in range(len(mean_errors))]
            max_errors = []
            min_errors = []
            for iter_idx in range(max_len):
                max_errors.append(np.max(iter_to_errors[iter_idx]))
                min_errors.append(np.min(iter_to_errors[iter_idx]))

            plotter.plot_max_min_fill(x_array, max_errors, min_errors, color = "blue", alpha = alpha_map[size], label = r"$n_{total}$"+r"$={}^3$".format(size))
        plotter.ax.legend(loc='lower right')
        plotter.save_fig(graph_output)

def overlap_to_epoch_multiple_patterns():
    for a,b in [(1,1),(1,2),(2,2)]:
        graph_output = "new_graphs/overlap_to_epochs_{}_{}.jpg".format(a,b)
        plotter = manuscriptPlots.plot()
        plotter.set_ylim(0,1)
        plotter.set_xlim(0.5,2000)
        plotter.set_xticks([1000*i for i in range(6)])
        plotter.set_yticks([0.2*i for i in range(6)])
        plotter.set_xlabel("Epochs")
        plotter.set_ylabel(r"$Q_n$")
        plotter.set_yScaled()
        plotter.set_xScaled()
        plotter.set_title(r"$n_A = {}, n_B = {}$".format(a,b))
        plotter.initialize_figure()
        plotter.set_xLog()
        # plotter.set_yLog()
        for size in [4,5,6]:
            errors_dict = {}
            max_len = 0
            for i in range(100):
                dir = "data/{}_{}_l_{}/{:03d}/".format(a,b,size,i)
                if not os.path.isfile(dir+"info.csv"):
                    continue
                iters = pd.read_csv(dir+"info.csv")["Iter"].to_numpy()
                all_overlaps = np.loadtxt(dir+"q_values.txt")
                if not len(all_overlaps):
                    continue
                errors = [all_overlaps[i] for i in iters]
                if not len(errors):
                    continue
                if len(errors)<2:
                    continue

                errors_dict[i] = errors
                if len(errors)>max_len:
                    max_len = len(errors)

            iter_to_errors = {i: [] for i in range(max_len)}
            for _, errors_array in errors_dict.items():
                for iter_idx in range(len(errors_array)):
                    iter_to_errors[iter_idx].append(errors_array[iter_idx])
            mean_errors = []
            max_errors = []
            min_errors = []
            for iter_idx in range(max_len):
                mean_errors.append(np.mean(iter_to_errors[iter_idx]))
                max_errors.append(np.max(iter_to_errors[iter_idx]))
                min_errors.append(np.min(iter_to_errors[iter_idx]))
            x_array = [(i+1)/2 for i in range(len(mean_errors))]
            plotter.plot_max_min_fill(x_array, max_errors, min_errors, color = "blue", alpha = alpha_map[size], label = r"$n_{total}$"+r"$={}^3$".format(size))
        plotter.ax.legend(loc='lower right')
        plotter.save_fig(graph_output)

def overlap_to_epoch_single_pattern_cellwise():
    for n in [3,4,5,6]:
        graph_output = "new_graphs/overlap_to_epochs_single_pattern_cellwise_n_{}.jpg".format(n)
        plotter = manuscriptPlots.plot()
        plotter.set_ylim(0,1)
        plotter.set_xlim(0.5,2000)
        plotter.set_xticks([1000*i for i in range(6)])
        plotter.set_yticks([0.2*i for i in range(6)])
        plotter.set_xlabel("Epochs")
        plotter.set_ylabel(r"$Q_2$")
        plotter.set_yScaled()
        plotter.set_xScaled()
        plotter.set_title(r"$n = {}$".format(n))
        plotter.initialize_figure()
        plotter.set_xLog()
        # plotter.set_yLog()
        # for size in [4,5,6]:
        errors_dict = {}
        max_len = 0
        for i in range(10):
            dir = "data/single_pattern_cellwise_n_{}/{:03d}/".format(n,i)
            if not os.path.isfile(dir+"info.csv"):
                continue
            iters = pd.read_csv(dir+"info.csv")["Iter"].to_numpy()
            all_overlaps = np.loadtxt(dir+"q_values.txt")
            if not len(all_overlaps):
                continue
            errors = [all_overlaps[i] for i in iters]
            if not len(errors):
                continue
            if len(errors)<2:
                continue

            errors_dict[i] = errors
            if len(errors)>max_len:
                max_len = len(errors)

        iter_to_errors = {i: [] for i in range(max_len)}
        for _, errors_array in errors_dict.items():
            for iter_idx in range(len(errors_array)):
                iter_to_errors[iter_idx].append(errors_array[iter_idx])
        mean_errors = []
        max_errors = []
        min_errors = []
        for iter_idx in range(max_len):
            mean_errors.append(np.mean(iter_to_errors[iter_idx]))
            max_errors.append(np.max(iter_to_errors[iter_idx]))
            min_errors.append(np.min(iter_to_errors[iter_idx]))
        x_array = [(i+1)/2 for i in range(len(mean_errors))]
        plotter.plot_max_min_fill(x_array, max_errors, min_errors, color = "blue", alpha = 0.8, label = r"$n_{target}$"+r"$={}$".format(n))
        plotter.ax.legend(loc='lower right')
        plotter.save_fig(graph_output)

def multiple_patterns_single_track(input_file,output_graph):
    colors = ["blue","orange","red","purple" ]
def error_to_epoch_single_pattern_cellwise():
    for n in [3,4,5,6]:
        graph_output = "new_graphs/error_to_epochs_single_pattern_cellwise_n_{}.jpg".format(n)
        plotter = manuscriptPlots.plot()
        plotter.set_ylim(1e-7,1)
        plotter.set_xlim(0.5,2000)
        plotter.set_xticks([1000*i for i in range(6)])
        plotter.set_yticks([0.2*i for i in range(6)])
        plotter.set_xlabel("Epochs")
        plotter.set_ylabel(r"$Q_2$")
        plotter.set_yScaled()
        plotter.set_xScaled()
        plotter.set_title(r"$n = {}$".format(n))
        plotter.initialize_figure()
        plotter.set_xLog()
        plotter.set_yLog()
        errors_dict = {}
        max_len = 0
        for i in range(10):
            dir = "data/single_pattern_cellwise_n_{}/{:03d}/".format(n,i)
            if not os.path.isfile(dir+"info.csv"):
                continue
            iters = pd.read_csv(dir+"info.csv")["Iter"].to_numpy()
            all_overlaps = np.loadtxt(dir+"costs.txt")
            if not len(all_overlaps):
                continue
            errors = [all_overlaps[i] for i in iters]
            if not len(errors):
                continue
            if len(errors)<2:
                continue

            errors_dict[i] = errors
            if len(errors)>max_len:
                max_len = len(errors)

        iter_to_errors = {i: [] for i in range(max_len)}
        for _, errors_array in errors_dict.items():
            for iter_idx in range(len(errors_array)):
                iter_to_errors[iter_idx].append(errors_array[iter_idx])
        mean_errors = []
        max_errors = []
        min_errors = []
        for iter_idx in range(max_len):
            mean_errors.append(np.mean(iter_to_errors[iter_idx]))
            max_errors.append(np.max(iter_to_errors[iter_idx]))
            min_errors.append(np.min(iter_to_errors[iter_idx]))
        x_array = [(i+1)/n for i in range(len(mean_errors))]
        plotter.plot_max_min_fill(x_array, max_errors, min_errors, color = "blue", alpha = 0.8, label = r"$n_{target}$"+r"$={}$".format(n))
        plotter.ax.legend(loc='lower right')
        plotter.save_fig(graph_output)

def error_to_iter_cellwise_single_tracks():
    for n_target in [3,4,5,6]:
        graph_output = "new_graphs/error_to_iters_cellwise_single_tracks_n_target_{}.jpg".format(n_target)

        plotter = manuscriptPlots.plot()
        plotter.set_ylim(1e-7,5)
        # plotter.set_ylim(0.6,1)

        plotter.set_xlim(1,100000)
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

        legend_added = False
        for i in range(10):
            dir = "data/single_pattern_cellwise_n_{}/{:03d}/".format(n_target,i)
            if not os.path.isfile(dir+"costs.txt"):
                continue
            errors = np.loadtxt(dir+"costs.txt")
            # if not errors[-1]<1e-4:
            #     continue
            # errors = np.loadtxt(dir+"q_values.txt")

            x_array = [i+1 for i in range(len(errors))]
            plotter.plot_xy(x_array, errors, color = "black", alpha = 0.4,label = "_none")
            if legend_added:
                continue
            plotter.plot_xy(x_array, errors, color = "black", alpha = 0.4, label = r"$n_{(target)}=$"+r'${}$'.format(n_target))
            legend_added = True


        # plotter.plot_errorfill(x_array,df_decrease["mean"].to_numpy(),df_decrease["sem"].to_numpy(),color = "black", alpha = 0.4, label = r"$\sigma_{hidden}^{(initial)}$")

        plotter.ax.legend(loc='lower left')
        plotter.save_fig(graph_output)

def error_to_iter_cellwise_sample_track():
    n=5
    i = 5
    graph_output = "new_graphs/error_to_iters_cellwise_sample_track_n_{}_{}.jpg".format(n,i)

    plotter = manuscriptPlots.plot()
    plotter.set_ylim(1e-7,5)
    # plotter.set_ylim(0.6,1)

    plotter.set_xlim(1,100000)
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

    legend_added = False
    # for i in range(10):

    dir = "data/single_pattern_cellwise_n_{}/{:03d}/".format(n,i)

    errors = np.loadtxt(dir+"costs.txt")
    # if not errors[-1]<1e-4:
    #     continue
    # errors = np.loadtxt(dir+"q_values.txt")

    x_array = [i+1 for i in range(len(errors))]
    plotter.plot_xy(x_array, errors, color = "black", alpha = 0.4,label = "_none")

    plotter.plot_xy(x_array, errors, color = "black", alpha = 0.4, label = r"$n_{(target)}=$"+r'${}$'.format(n))


    # plotter.plot_errorfill(x_array,df_decrease["mean"].to_numpy(),df_decrease["sem"].to_numpy(),color = "black", alpha = 0.4, label = r"$\sigma_{hidden}^{(initial)}$")

    plotter.ax.legend(loc='lower left')
    plotter.save_fig(graph_output)

    
    

def main():
    os.makedirs("new_graphs/",exist_ok=True)
    error_to_epoch_single_pattern_cellwise()
    overlap_to_epoch_single_pattern_cellwise()
    error_to_iter_cellwise_sample_track()

if __name__ == "__main__":
    main()