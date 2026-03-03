from plots import *
from panel_3_plots_new import construct_multiple_pattern_foldername, label_writer, n_cell_to_patterns
from toolbox import manuscriptPlots
import os


def single_run_error_overlap(dirlist,
                             output_folder,
                             xlim = [1,2500],
                             title = None):
    # for l in [4,5,6]:
    #     for n in [1,2,3,4,5,6]:
    #         output_folder = "Panel_4/l_{}_n_{}/".format(l,n)
    #         os.makedirs(output_folder,exist_ok=True)
    #         dirlist = find_complete_runs(["data/kv_10_l_{}_n_{}/{:03d}/".format(l,n,i) for i in range(100)])
    for i, dir in enumerate(dirlist):
        if os.path.isfile(dir +"info.csv"):
            info = pd.read_csv(dir+"info.csv")
            iters = info["Iter"].to_list()
            if not sorted(iters) == iters:
                print("iters not sorted in {}".format(dir))
                continue
            costs = info["Error"].to_list()
            q_values = info["Overlap"].to_list()
        else: 
            costs = np.loadtxt(dir+"costs.txt")
            q_values = np.loadtxt(dir+"q_values.txt")
            iters = np.arange(1,len(costs)+1)
        savefile = output_folder+"{:03d}_error_overlap.png".format(int(os.path.basename(os.path.dirname(dir))))
        plotter = manuscriptPlots.plot_shared_x_axis()
        plotter.set_xlim(xlim)
        plotter.set_ylim_top([1e-8,1e1])
        plotter.set_ylim_bottom([0.5,1.1])
        plotter.set_yticks_bottom([0.6,0.8,1])
        plotter.set_xLog(True)
        plotter.set_yLog_top(True)
        plotter.set_xlabel("Iterations")
        plotter.set_ylabel_top(r"$<|1-\sigma_{T}/\sigma|>$")
        plotter.set_ylabel_bottom(r"$Q_2$")
        plotter.set_title(title)
        plotter.initialize_sharedx_figure()
        plotter.plot_xy_top(iters, costs, marker="o", color="blue")
        plotter.plot_xy_bottom(iters, q_values, marker="o", color="blue")
        plotter.ax_bottom.hlines(y=1, xmin = plotter.xlim[0], xmax = plotter.xlim[1], color="black", linestyle="--",alpha = 0.5)
        plotter.transparent = False
        plotter.save_fig(savefile)

def single_run_error_s0_overlap(dirlist,
                             output_folder,
                             xlim = [1,2500],
                             title = None):
    # for l in [4,5,6]:
    #     for n in [1,2,3,4,5,6]:
    #         output_folder = "Panel_4/l_{}_n_{}/".format(l,n)
    #         os.makedirs(output_folder,exist_ok=True)
    #         dirlist = find_complete_runs(["data/kv_10_l_{}_n_{}/{:03d}/".format(l,n,i) for i in range(100)])
    for i, dir in enumerate(dirlist):
        if os.path.isfile(dir +"info.csv"):
            info = pd.read_csv(dir+"info.csv")
            iters = info["Iter"].to_list()
            if not sorted(iters) == iters:
                print("iters not sorted in {}".format(dir))
                continue
            costs = info["Error"].to_list()
            q_values = info["Distance"].to_list()
        else: 
            costs = np.loadtxt(dir+"costs.txt")
            q_values = np.loadtxt(dir+"q_values.txt")
            iters = np.arange(1,len(costs)+1)
        savefile = output_folder+"{:03d}_error_distance.png".format(int(os.path.basename(os.path.dirname(dir))))
        plotter = manuscriptPlots.plot_shared_x_axis()
        plotter.set_xlim(xlim)
        plotter.set_ylim_top([1e-8,1e1])
        plotter.set_ylim_bottom([0,0.5])
        plotter.set_yticks_bottom([0.2,0.4])
        plotter.set_xLog(True)
        plotter.set_yLog_top(True)
        plotter.set_xlabel("Iterations")
        plotter.set_ylabel_top(r"$<|1-\sigma_{T}/\sigma|>$")
        plotter.set_ylabel_bottom(r"$SD(s_0)$")
        plotter.set_title(title)
        plotter.initialize_sharedx_figure()
        plotter.plot_xy_top(iters, costs, marker="o", color="blue")
        plotter.plot_xy_bottom(iters, q_values, marker="o", color="red")
        plotter.ax_bottom.hlines(y=1, xmin = plotter.xlim[0], xmax = plotter.xlim[1], color="black", linestyle="--",alpha = 0.5)
        plotter.transparent = False
        plotter.save_fig(savefile)

def main():
    # os.makedirs("Panel_4",exist_ok=True)
    # for l in [4,5,6]:
    #     title = r"$n_{total}=$"+"{}".format(l**3)
    #     for n in [1,2,3,4,5,6]:
    #         output_folder = "Panel_4/l_{}_n_{}/".format(l,n)
    #         os.makedirs(output_folder,exist_ok=True)
    #         dirlist = find_complete_runs(["data/kv_10_l_{}_n_{}/{:03d}/".format(l,n,i) for i in range(100)])
    #         single_run_error_overlap(dirlist, output_folder, title=title)
    
    l = 4
    title = r"$n_{total}=$"+"{}".format(l**3)
    for n in [2,3,4]:
        for pattern in n_cell_to_patterns[n]:
            output_folder = construct_multiple_pattern_foldername(header = "Panel_4/l_{}".format(l), pattern = pattern)
            os.makedirs(output_folder,exist_ok=True)
            dirlist = find_complete_runs_multiple_patterns([construct_multiple_pattern_foldername("data/new_kv_10_l_{}_p".format(l),pattern)+"{:03d}/".format(j) for j in range(100)])
            single_run_error_s0_overlap(dirlist, output_folder,
                                     xlim = [10,20000],
                                    title = title)
    return
if __name__ == "__main__":
    main()