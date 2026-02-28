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

tolerance = 1e-5
lwmap = {1:15,2:10,3:5,4:5,5:5,6:5}
lw_sp_map = {40:15,20:10,10:5, 5:20}
n_spheroid_color_map = {40:"blue", 20:"orange", 10:"red", 5: "green" }
n_cells_color_map = {1:"blue",2:"orange",3:"purple", 4:"red" ,5:"yellow",6:"green"}
n_runs = 100
n_cell_to_patterns = {2:[[1,1]], 3:[[1,2],[1,1,1]], 4:[[1,1,1,1],[1,1,2],[1,3],[2,2]]}
markersize = 1000

patterns = [{"pattern":[1,1],"marker":"o"},
            {"pattern":[1,1,1],"marker":"x"},
            {"pattern":[1,2],"marker":"h"},
            {"pattern":[1,1,1,1],"marker":"d"},
            {"pattern":[1,3],"marker":"P"},
            {"pattern":[2,2],"marker":"X"},
            {"pattern":[1,1,2],"marker":"s"}]

markers = ["o","o","P", "o","P","*","s"]
def construct_multiple_pattern_foldername(header = "", pattern = [1,1]):
    dir = "{}".format(header)
    for p in pattern:
        dir+="_{}".format(p)
    dir +="/"
    return dir

def write_final_s0(dir_list,savefile):
    s0_vals = []
    for dir in dir_list:
        df = pd.read_csv(dir+"cellParameters.input", sep = " ", header=None)
        s0_vals += df[2].to_list()
    np.savetxt(savefile,s0_vals)

def write_final_s0_multiple_patterns():
    for n_cells, pattern_list in n_cell_to_patterns.items():
        for subpattern in pattern_list:
            dir_list = find_complete_runs_multiple_patterns(subpattern)
            savefile = "data/" + construct_multiple_pattern_foldername(subpattern = subpattern) + "final_s0.txt"
            write_final_s0(dir_list,savefile)

def find_complete_runs_multiple_patterns(subpattern):
    test_dir = "data/" + construct_multiple_pattern_foldername(subpattern = subpattern)
    complete_dirs = []
    for i in range(n_runs):
        run_dir = test_dir+"{:03d}/".format(i)
        if not os.path.isfile(run_dir+"info.csv"):
            continue
        info = pd.read_csv(run_dir+"info.csv")
        if info["Error"].to_list()[-1]>tolerance:
            continue
        complete_dirs.append(run_dir)
    print(subpattern, "has {} complete runs".format(len(complete_dirs)))
    return complete_dirs

def find_complete_runs_single_pattern_periodic(l=4,n=2):
    complete_dirs =[]
    for i in range(n_runs):
        test_dir ="data/kv_10_l_{}_n_{}/{:03d}/".format(l,n,i)
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
        complete_dirs.append(test_dir)
    return complete_dirs

def s0_histogram(n_cells = 3):
    l = 4
    savefile = "multiple_pattern_graphs/final_s0_histogram_l_{}_n_{}.png".format(l,n_cells)
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(0,4)
    plotter.set_xlim(3,7)
    plotter.set_xticks([i for i in range(3,8)])
    plotter.set_yticks([i for i in range(1,10)])
    plotter.set_xlabel(r"$s_0^{hidden}$")
    plotter.set_title(r"$n_{(total)} = $"+ "{}".format(l**3))
    plotter.initialize_figure()
    plotter.ax.vlines(x = 5, ymin = 0, ymax = 3, linestyle= "dashed",color = plotter.neutralColor, label = r"$s_0^{(init)}$")
    full_pattern_array =np.loadtxt("data/kv_10_l_{}_n_{}/final_s0.txt".format(l,n_cells))
    plotter.histogram_from_array(full_pattern_array, fit_type="kde",bins = 50,label = "Full Pattern",color = n_cells_color_map[n_cells],alpha = 0.2)
    for subpattern in n_cell_to_patterns[n_cells]:
        array = np.loadtxt("data/"+construct_multiple_pattern_foldername(subpattern=subpattern)+"final_s0.txt")
        label = "" 
        for i in subpattern: 
            label+="{} ".format(i)
        plotter.histogram_from_array(array,label = label ,fit_type = "kde", bins = 50, color = "black", alpha = 0.8)

    plotter.ax.legend()
    plotter.save_fig(savefile)

#hard coded points...
def SD_s0_to_n_cells():
    def label_writer(pattern):
        label = ""
        for i,p in enumerate(pattern):
            label += str(p)
            if i <len(pattern)-1:
                label +=r"\rightarrow"
        return label
    l = 4
    n_cells =[2,3,4]
    header = "data/kv_10_l_{}_".format(l)
    savefile = "multiple_pattern_graphs/sd_s0_to_subpatterns_l_{}.png".format(l)
    plotter = manuscriptPlots.plot()
    plotter.set_xticks(n_cells)
    plotter.set_yticks([0.5*i for i in range(1,10)])
    plotter.set_xlim(1.5,4.5)
    plotter.set_ylim(0,0.5)
    plotter.set_xlabel(r"$n_T$")
    plotter.set_ylabel(r"$SD(s_0)$")
    # plotter.set_ylabel(r"$Q_n$")

    plotter.set_title(r"$n_{(total)} =$"+ "{}".format(l**3))
    # plotter.set_yScaled()
    plotter.initialize_figure()
    for n in n_cells:
        plotter.plot_scatter(
            [n],
            [np.std(np.loadtxt(header + "n_{}/final_s0.txt".format(n)))],
            marker = 'X',
            s = markersize,
            color = n_cells_color_map[n],
            label = r"$n_T =$"+"{}".format(n))
    # n = 2
    # plotter.plot_scatter(
    #     [2],
    #     [np.std(np.loadtxt(header + "p_1_1/final_s0.txt"))],
    #     marker = 'D',
    #     s = markersize,
    #     color = n_cells_color_map[2],
    #     label = "1->1")
    for pattern_dict in patterns:
        pattern = pattern_dict["pattern"]
        marker = pattern_dict["marker"]
        plotter.plot_scatter(
            [sum(pattern)],
            [np.std(np.loadtxt(construct_multiple_pattern_foldername(header+"p",pattern)+"final_s0.txt"))],
            marker = marker,
            s = markersize,
            color = n_cells_color_map[sum(pattern)],
            label = r"${}$".format(label_writer(pattern)))


    plotter.ax.legend(ncols=2)
    plotter.save_fig(savefile)


def main():
    write_final_s0_multiple_patterns()
    # s0_histogram(4)
    # SD_s0_to_n_cells()
if __name__ == "__main__":
    main()