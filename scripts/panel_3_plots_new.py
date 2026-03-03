from plots import *
import os

n_cell_to_patterns = {2:[[1,1]], 3:[[1,1,1],[1,2]], 4:[[1,1,1,1],[1,1,2],[1,3],[2,2]]}
markersize = 1000

patterns = [{"pattern":[1,1],"marker":"o"},
            {"pattern":[1,1,1],"marker":"x"},
            {"pattern":[1,2],"marker":"h"},
            {"pattern":[1,1,1,1],"marker":"d"},
            {"pattern":[1,3],"marker":"P"},
            {"pattern":[2,2],"marker":"X"},
            {"pattern":[1,1,2],"marker":"s"}]

markers = ["P","*","s","^"]
linestyles = ["dashed","dotted","dashdot","densely dashed","loosely dashed"]
alphas = [0.4,0.6,0.8,1.0,0.2]

def label_writer(pattern):
    label = ""
    for i,p in enumerate(pattern):
        label += str(p)
        if i <len(pattern)-1:
            label +=r"\rightarrow"
    return r"${}$".format(label)

def construct_multiple_pattern_foldername(header = "", pattern = [1,1]):
    dir = "{}".format(header)
    for p in pattern:
        dir+="_{}".format(p)
    dir +="/"
    return dir


def write_final_s0_multiple_patterns():
    header = "data/new_kv_10_l_4_p"
    for _, pattern_list in n_cell_to_patterns.items():
        for pattern in pattern_list:
            experiment_dir = construct_multiple_pattern_foldername(header=header,pattern=pattern)
            dir_list = find_complete_runs_multiple_patterns([experiment_dir+"{:03d}/".format(i) for i in range(100)])
            savefile = experiment_dir+ "final_s0.txt"
            print(experiment_dir,savefile)
            write_final_s0(dir_list,savefile)

def s0_histogram_periodic_combined():
    input_filename = "final_s0.txt"
    for l in [4]:
        for n in [2,3,4]:
            savefile = "Panel_3_new/s0_histogram_periodic_l_{}_n_{}.png".format(l,n)
            label_to_data = {}
            label_to_data[r"$n_T=$"+"{}".format(n)]={"dir":"data/new_kv_10_l_{}_n_{}/".format(l,n), "color":n_cells_color_map[n], "bins":30, "alpha":0.1}
            for pattern in n_cell_to_patterns[n]:
                dir = construct_multiple_pattern_foldername("data/new_kv_10_l_{}_p".format(l),pattern)
                label_to_data[label_writer(pattern)] = {"dir":dir, "color":n_cells_color_map[n], "bins":30}
            title = r"$n_{total} = $"+"{}".format(l**3)
            hist = histogram(   label_to_data=label_to_data,
                                input_filename=input_filename,
                                title=title,
                                xlim = [3,7],
                                ylim = [0,5])
            hist.ax.vlines(x = 5, ymin = 0, ymax = 3.2, linestyle= "dashed",color = "black", label = r"$s_0^{(init)}$",linewidth =15)
            hist.ax.legend()
            hist.save_fig(savefile)


def error_to_iters_single_pattern():
    for l in [4]:
        label_to_data = {}
        for n in [1,2,4]:
            savefile = "Panel_3_new/error_to_iters_periodic_l_{}_n_{}.png".format(l,n)
            title = r"$n_{total} = $"+"{}".format(l**3)
            dirlist = find_complete_runs(["data/new_kv_10_l_{}_n_{}/{:03d}/".format(l,n,i) for i in range(100)])
            label_to_data[r"$n_T=$"+"{}".format(n)]={"dirlist":dirlist, "color":n_cells_color_map[n]}
        single_tracks_to_iters(label_to_data=label_to_data,title=title,savefile=savefile,)

def final_iteration_to_final_overlap_periodic():
    l = 4
    for n in [2,3,4]: 
        label_to_dir = {}
        savefile = "Panel_3_new/final_iterations_to_final_overlap_l_{}_n_{}.png".format(l,n)
        title = r"$n_{total} = $"+"{}".format(l**3)
        dirlist = find_complete_runs(["data/kv_10_l_{}_n_{}/{:03d}/".format(l,n,i) for i in range(100)])
        x_array = [len(np.loadtxt(dir+"costs.txt"))for dir in dirlist]
        y_array = [np.loadtxt(dir+"q_values.txt")[-1] for dir in dirlist]
        label_to_dir[r"$n_T=$"+"{}".format(n)]= {"x_array":x_array,"y_array":y_array,"marker":"o","color":n_cells_color_map[n],"s":3000,"alpha":0.5}
        for i,pattern in enumerate(n_cell_to_patterns[n]):
            print(n,pattern)
            dirlist = find_complete_runs_multiple_patterns([construct_multiple_pattern_foldername("data/new_kv_10_l_{}_p".format(l),pattern)+"{:03d}/".format(j) for j in range(100)])
            x_array = [len(np.loadtxt(dir+"costs.txt"))for dir in dirlist]
            y_array = [np.loadtxt(dir+"q_values.txt")[-1] for dir in dirlist]
            label_to_dir[label_writer(pattern)] = {"x_array":x_array,"y_array":y_array,"marker":markers[i],"color":n_cells_color_map[n],"alpha":0.7}

        
        plotter = scatter_plot(label_to_dir,
                            xlim = [50,20000],
                            ylim=[0.2,1.05],
                                title = r"$n_{total}=$"+"{}".format(l**3),
)
        plotter.ax.legend()
        plotter.save_fig(savefile)

#hard coded points...
def SD_s0_scatter_all_patterns():
    l = 4
    n_cells =[2,3,4]
    header = "data/new_kv_10_l_{}_p".format(l)
    savefile = "Panel_3_new/SD_s0_scatter_patterns_l_{}.png".format(l)
    label_to_dir = {}
    for n in n_cells:
        #single pattern
        y_array = [np.std(np.loadtxt("data/kv_10_l_{}_n_{}/final_s0.txt".format(l,n)))]
        label_to_dir[r"$n_T=$"+"{}".format(n)]={"x_array":[n],"y_array":y_array,"marker":"o","color":n_cells_color_map[n], "s":4000}
        for i, pattern in enumerate(n_cell_to_patterns[n]):
            y_array = [np.std(np.loadtxt(construct_multiple_pattern_foldername(header = header,pattern = pattern)+"final_s0.txt"))]
            label_to_dir[label_writer(pattern)]={"x_array":[n],"y_array":y_array,"marker":markers[i],"color":n_cells_color_map[n], "s":3000}
    plotter = scatter_plot(label_to_dir,
                     xlim =[1.5,4.5],
                     ylim=[0.2,0.42],
                     title = r"$n_{total}=$"+"{}".format(l**3),
                     xticks= n_cells,
                     xlog=False,
                     xlabel=r"$n_T$",
                     ylabel=r"$SD(s_0)$")
    # plotter.ax.legend()
    plotter.save_fig(savefile)
        # for n in [1,2,3,4,5,6]:
        #     dirlist = find_complete_runs(["data/new_kv_10_l_{}_n_{}/{:03d}/".format(l,n,i) for i in range(100)]))


def main():
    os.makedirs("Panel_3_new",exist_ok=True)
    write_final_s0_multiple_patterns()
    final_iteration_to_final_overlap_periodic()
    SD_s0_scatter_all_patterns()
if __name__ == "__main__":
    main()