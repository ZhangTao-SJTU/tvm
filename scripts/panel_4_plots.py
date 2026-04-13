from plots import *
from panel_3_plots_new import construct_multiple_pattern_foldername, label_writer, n_cell_to_patterns
from toolbox import manuscriptPlots
from toolbox.periodic import PeriodicTissue
from toolbox.training import Training
from toolbox.stress import calculate_max_shear_stress
import os
output_folder = "Panels/Panel_4/"
color_dict = {"Error": "#D72638", r"$SD(s_0)$": "#F49D37", r"$Q_2$": "#3F88C5"}
linestyle_dict = {"Error": "-", r"$SD(s_0)$": ":", r"$Q_2$": "--"}
marker_dict = {"Error": "D", r"$SD(s_0)$": "o", r"$Q_2$": "^"}
alpha = 0.4
markersize = 700
def single_run_stacked_plot(label_to_data,**kwargs):
    plotter = stacked_plot(**kwargs)
    for label, data in label_to_data.items():
        if data["loc"] == "top":
            plotter.plot_xy_top(data["x"], data["y"], color = color_dict[label], label=label,alpha = alpha)
            plotter.ax_top.scatter(data["x"], data["y"], color = color_dict[label], label=label, marker = marker_dict[label], s=markersize,alpha = alpha)
            plotter.ax_top.scatter(data["x"], data["y"], color = color_dict[label], marker = marker_dict[label], s=markersize,facecolors = "none",edgecolors = "black")
        elif data["loc"] == "bottom":
            plotter.plot_xy_bottom(data["x"], data["y"], color = color_dict[label], label=label, alpha=alpha)
            plotter.ax_bottom.scatter(data["x"], data["y"], color = color_dict[label], label=label, marker = marker_dict[label], s=markersize, alpha=alpha)
            plotter.ax_bottom.scatter(data["x"], data["y"],  color = color_dict[label], marker = marker_dict[label],s=markersize,facecolors = "none",edgecolors = "black")

        elif data["loc"] == "bottom right":
            plotter.plot_xy_bottom_right(data["x"], data["y"], color = color_dict[label], label=label, alpha=alpha)
            plotter.ax_bottom_right.scatter(data["x"], data["y"], color = color_dict[label], label=label,marker = marker_dict[label], s=markersize,alpha = alpha)
            plotter.ax_bottom_right.scatter(data["x"], data["y"], color = color_dict[label], marker = marker_dict[label], s=markersize,facecolors = "none",edgecolors = "black")

    return plotter



def stacked_multiple_pattern_plots():
    l = 4
    # title = r"$n_{total}=$"+"{}".format(l**3)
    for n in [2,3,4]:
        for pattern in n_cell_to_patterns[n]:
            output_subfolder = construct_multiple_pattern_foldername(header = f"{output_folder}l_{l}", pattern = pattern)
            os.makedirs(output_subfolder,exist_ok=True)
            dirlist = find_complete_runs_multiple_patterns([construct_multiple_pattern_foldername("data/new_kv_10_l_{}_p".format(l),pattern)+"{:03d}/".format(j) for j in range(100)])
            for dir in dirlist:
                df = pd.read_csv(dir+"info.csv")
                iters = df["Iter"].to_list()
                if not sorted(iters) == iters:
                    print("iters not sorted in {}".format(dir))
                    continue
                costs = df["Error"].to_list()
                q_values = df["Overlap"].to_list()
                s0_overlap = df["Distance"].to_list()
                # add initial values to iters, costs, q_values, s0_overlap:
                # first element of iters is 1, 
                # first element of costs is initial cost
                # , first element of q_values is 1
                # , first element of s0_overlap is 0.
                def initial_cost():
                    init_stresses = dir+"0000000.stresses.csv"
                    df = pd.read_csv(init_stresses)
                    target = df["Target"].to_numpy()
                    current = df["Current"].to_numpy()
                    error = np.sqrt(np.sum((current-target)**2))
                    return error

                iters = [1] + iters
                costs = [initial_cost()] + costs
                q_values = [1] + q_values
                s0_overlap = [0] + s0_overlap
                # 1. error and q values from info.csv
                savefile = output_subfolder+"{:03d}_error_overlap.png".format(int(os.path.basename(os.path.dirname(dir))))
                label_to_data = {}
                label_to_data["Error"] = {"x": iters, "y": costs, "loc": "top"}
                label_to_data[r"$Q_2$"] = {"x": iters, "y": q_values, "loc": "bottom"}
                label_to_data[r"$SD(s_0)$"] = {"x": iters, "y": s0_overlap, "loc": "bottom right"}
                title = f"{label_writer(pattern)}"
                plotter = single_run_stacked_plot(label_to_data, title=title, xlim = [1,20000], ylim_top = [1e-8,1e1], ylim_bottom = [0,1.1], ylim_bottom_right = [0,0.6],
                                                   yticks_bottom = [0, 0.5, 1])
                # plotter.ax_bottom.hlines(y=1, xmin = plotter.xlim[0], xmax = plotter.xlim[1], color="black", linestyle=":",alpha = 0.2)

                plotter.legend_top()
                plotter.transparent=False
                plotter.save_fig(savefile)
                
# this function writes cellErrors.csv in dir.
# dir must have info.csv and files/
def write_cell_errors(dir):
    info = pd.read_csv(dir+"info.csv")
    stresses = pd.read_csv(dir+"files/0000000.stresses.csv")
    cellIDs = stresses["CellID"].to_list()
    target = stresses["Target"].to_list()[0]
    iter_to_cell_errors_list = []
    for iter in info["Iter"].to_list():
        iter = int(iter)
        file = f"{iter:07d}.bulk.txt"
        sample = PeriodicTissue.from_config(dir + "files/",file)
        training = Training.from_sample(sample)
        training.load_cell_parameters(f"{iter:07d}.cellParameters.input")
        iter_to_cell_errors = {"Iteration": iter}
        for i in range(len(cellIDs)):
            error = abs(calculate_max_shear_stress(training._config, cellIDs[i]) - target)/target
            iter_to_cell_errors[f"{i}"] = error
        iter_to_cell_errors_list.append(iter_to_cell_errors)
    df = pd.DataFrame(iter_to_cell_errors_list)
    df.to_csv(dir+"cellErrors.csv", index=False)
def epoch_cell_errors(dir = "data/004/", xlim = [34,273], offset = 20,savefile = "Panels/Panel_4/initial_epochs.png", title = None):
    colors = [ "#473144",  "#DF9B6D","#af1b3f",]
    markers = ["o", "s", "D"]
    
    plotter = manuscriptPlots.plot()
    plotter.set_ylim([1e-9,10])
    # x_0 = 34
    # x_1 = 273
    # x_0 = 34
    # x_1 = 273
    plotter.set_xlim([xlim[0]-offset, xlim[1]+offset])
    # plotter.set_xticks(xticks)
    # plotter.set_yticks(yticks)


    plotter.set_xlabel("Iteration")
    plotter.set_ylabel(r"$|1-\sigma_{T}/\sigma|$")
    if title is not None:
        plotter.set_title(title)

    # plotter.set_xLog()
    plotter.set_yLog()
    plotter.initialize_figure()
    info = pd.read_csv(dir + "info.csv")
    errors = pd.read_csv(dir + "cellErrors.csv")

    plotter.plot_xy(info["Iter"], info["Error"], label = "<Error>", color="black", alpha = 1, linewidth=10)
    for i in range(3):
        plotter.plot_xy(info["Iter"], errors[str(i)],label = "_none",alpha = 0.6, color=colors[i],linestyle=":", linewidth=10)
        plotter.plot_scatter(info["Iter"], errors[str(i)], label=f"Cell {i+1}",alpha = 0.8, color=colors[i], s=3000, marker=markers[i])
        plotter.plot_scatter(info["Iter"], errors[str(i)], label=f"_Cell {i+1}",alpha = 0.8, color=colors[i], s=3000, marker=markers[i], facecolors = "none", edgecolors = "black", linewidth=2)

            # plotter.plot_xy(info["Iter"], errors["1"], color="red",alpha=0.4, label="Cell 1")
    # plotter.save_fig("test.png")
    plotter.ax.hlines(y = 1e-6, xmin = plotter.xlim[0], xmax = plotter.xlim[1], color="black", linestyle="--",alpha = 0.5)
    plotter.ax.legend(ncol = 2)
    plotter.save_fig(savefile)

def main():
    os.makedirs(output_folder,exist_ok=True)
    stacked_multiple_pattern_plots()
    # epoch_cell_errors(xlim = [34,273])
    # epoch_cell_errors(xlim = [4668,4692], offset = 4, savefile = "Panels/Panel_4/final_epochs.png")

if __name__ == "__main__":
    main()