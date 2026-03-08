from plots import *
from panel_3_plots_new import construct_multiple_pattern_foldername, label_writer, n_cell_to_patterns
from toolbox import manuscriptPlots
from toolbox.periodic import PeriodicTissue
from toolbox.training import Training
from toolbox.stress import calculate_max_shear_stress
import os

output_folder = "Panels/Panel_4/"
def single_run_stacked_plot(label_to_data,**kwargs):
    plotter = stacked_plot(**kwargs)
    for label, data in label_to_data.items():
        if data["loc"] == "top":
            plotter.plot_xy_top(data["x"], data["y"], color = data.get("color", None), label=label)
        elif data["loc"] == "bottom":
            plotter.plot_xy_bottom(data["x"], data["y"], color = data.get("color", None), label=label)
        # plotter.plot_xy_top(iters, costs, marker="o", color="blue", )
        # plotter.plot_xy_bottom(iters, q_values, marker="o", color="blue")
        # plotter.ax_bottom.hlines(y=1, xmin = plotter.xlim[0], xmax = plotter.xlim[1], color="black", linestyle="--",alpha = 0.5)
        # plotter.transparent = False
        # plotter.save_fig(savefile)
    return plotter

def stacked_error_Q2_single_pattern():
    for l in [4,5,6]:
        title = r"$n_{total}=$"+"{}".format(l**3)
        for n in [1,2,3,4,5,6]:
            output_subfolder = f"{output_folder}l_{l}_n_{n}/"
            os.makedirs(output_subfolder,exist_ok=True)
            dirlist = find_complete_runs(["data/new_kv_10_l_{}_n_{}/{:03d}/".format(l,n,i) for i in range(100)])
            for dir in dirlist:
                label_to_data = {}
                costs = np.genfromtxt(dir+"costs.txt")
                iters = np.arange(1,len(costs)+1)
                q_values = np.genfromtxt(dir+"q_values.txt")
                savefile = output_subfolder+"{:03d}_error_overlap.png".format(int(os.path.basename(os.path.dirname(dir))))
                label_to_data["Error"] = {"x": iters, "y": costs, "loc": "top"}
                label_to_data["Distance"] = {"x": iters, "y": q_values, "loc": "bottom"}
                plotter = single_run_stacked_plot(label_to_data, title=title)
                plotter.save_fig(savefile)

def stacked_multiple_pattern_plots():
    l = 4
    title = r"$n_{total}=$"+"{}".format(l**3)
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
                # 1. error and q values from info.csv
                savefile = output_subfolder+"{:03d}_error_overlap.png".format(int(os.path.basename(os.path.dirname(dir))))
                label_to_data = {}
                label_to_data["Error"] = {"x": iters, "y": costs, "loc": "top"}
                label_to_data["Distance"] = {"x": iters, "y": q_values, "loc": "bottom"}
                plotter = single_run_stacked_plot(label_to_data, title=title, xlim = [10,20000], ylim_top = [1e-8,1e1], ylim_bottom = [0.5,1.1], yticks_bottom = [0.6,0.8,1])
                plotter.ax_bottom.hlines(y=1, xmin = plotter.xlim[0], xmax = plotter.xlim[1], color="black", linestyle="--",alpha = 0.5)
                plotter.save_fig(savefile)
                
                # 2. error and s0 from info.csv
                savefile = output_subfolder+"{:03d}_error_s0.png".format(int(os.path.basename(os.path.dirname(dir))))
                label_to_data = {}
                label_to_data["Error"] = {"x": iters, "y": costs, "loc": "top"}
                label_to_data["Distance"] = {"x": iters, "y": s0_overlap, "loc": "bottom", "color": "red"}
                plotter = single_run_stacked_plot(label_to_data, 
                                                  title=title, 
                                                  xlim = [10,20000], 
                                                  ylim_top = [1e-8,1e1], 
                                                  ylim_bottom = [0,1], 
                                                  yticks_bottom = [0.2*i for i in range(5)],
                                                  ylabel_bottom = r"$SD(s_0)$")
                plotter.save_fig(savefile)

                # 3. Overlap and s0 from info.csv
                savefile = output_subfolder+"{:03d}_overlap_s0.png".format(int(os.path.basename(os.path.dirname(dir))))
                label_to_data = {}
                label_to_data["Overlap"] = {"x": iters, "y": q_values, "loc": "top"}
                label_to_data["Distance"] = {"x": iters, "y": s0_overlap, "loc": "bottom", "color": "red"}
                plotter = single_run_stacked_plot(label_to_data, 
                                                  title=title, 
                                                  xlim = [10,20000], 
                                                  ylim_top = [0.5,1.1], 
                                                  ylim_bottom = [0,1],
                                                yticks_top = [0.6,0.8,1], 
                                                yticks_bottom = [0.2*i for i in range(1,5)],
                                                ylog_top = False, 
                                                ylabel_top = r"$Q_2$",
                                                ylabel_bottom = r"$SD(s_0)$")
                plotter.ax_top.hlines(y=1, xmin = plotter.xlim[0], xmax = plotter.xlim[1], color="black", linestyle="--",alpha = 0.5)
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

def main():
    os.makedirs(output_folder,exist_ok=True)
    # write_cell_errors("data/004/")

    # stacked_error_Q2_single_pattern()
    stacked_multiple_pattern_plots()

if __name__ == "__main__":
    main()