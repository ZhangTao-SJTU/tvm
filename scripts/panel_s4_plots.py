from plots import *
from panel_3_plots_new import construct_multiple_pattern_foldername, label_writer, n_cell_to_patterns
from panel_4_plots import single_run_stacked_plot
from toolbox import manuscriptPlots
from toolbox.periodic import PeriodicTissue
from toolbox.training import Training
from toolbox.stress import calculate_max_shear_stress
import os
output_folder = "Panels/Panel_S4/"
# color_dict = {"Error": "#D72638", r"$SD(s_0)$": "#F49D37", r"$Q_2$": "#3F88C5"}
color_dict = {"Error": "red", r"$SD(s_0)$": "orange", r"$Q_2$": "blue"}

# colors = ["#5e029c", "#faa23e","#0e4008"]
colors = ['purple','green',"orange"]

# color_dict = {"Error": "#360657", r"$SD(s_0)$": "#F49D37", r"$Q_2$": "#23591c"}
# 
linestyle_dict = {"Error": "-", r"$SD(s_0)$": ":", r"$Q_2$": "--"}
marker_dict = {"Error": "D", r"$SD(s_0)$": "o", r"$Q_2$": "^"}
alpha = 0.4
markersize = 700


def stacked_spheroid_run_plots():
    for l in [5,6]:
        for n in [5,10,20,40]:
            success_folder = output_folder+f"l_{l}_n_{n:03d}_success/"
            failure_folder = output_folder+f"l_{l}_n_{n:03d}_failure/"
            os.makedirs(success_folder, exist_ok=True)
            os.makedirs(failure_folder, exist_ok=True)
            for i in range(100):
                dir = f"data/kv_10_l_{l}_n_sp_{n:03d}/{i:03d}/"
                # print("Checking: ", dir)
                if not os.path.isfile(dir+"costs.txt"):
                    continue
                if not os.path.isfile(dir+"q_values.txt"):
                    continue

                # print("Processing: ", dir)
                costs = np.loadtxt(dir+"costs.txt")
                if costs.shape == ():
                    continue
                q_values = np.loadtxt(dir+"q_values.txt")
                df = pd.read_csv(dir+"SD_s0.csv")

                # def calculate_s0_overlap():
                #     overlaps = []
                #     for i in range(len(costs)):
                #         df = pd.read_csv(dir+"cellPara".format(i))

                # 1. error and q values from info.csv
                if costs[-1] < 1e-6:
                    savefile = success_folder+f"{i:03d}.png"
                else:
                    savefile = failure_folder+f"{i:03d}.png"
                label_to_data = {}
                label_to_data["Error"] = {"x": np.arange(1, len(costs)+1), "y": costs, "loc": "top"}
                label_to_data[r"$Q_2$"] = {"x": np.arange(1, len(q_values)+1), "y": q_values, "loc": "bottom"}
                label_to_data[r"$SD(s_0)$"] = {"x": df["iter"], "y": df["SD_s0"], "loc": "bottom right"}
                plotter = single_run_stacked_plot(label_to_data, title=None, xlog = True, xlim = [1,20000], ylim_top = [1e-8,1e1], ylim_bottom = [0,1.1], ylim_bottom_right = [0,2],
                                                    yticks_bottom = [0, 0.5, 1])
                # plotter.ax_bottom.hlines(y=1, xmin = plotter.xlim[0], xmax = plotter.xlim[1], color="black", linestyle=":",alpha = 0.2)
                plotter.xLog = False
                plotter.legend_top()
                plotter.transparent=False

                plotter.save_fig(savefile)


def main():
    os.makedirs(output_folder,exist_ok=True)
    stacked_spheroid_run_plots()
    
if __name__ == "__main__":
    main()