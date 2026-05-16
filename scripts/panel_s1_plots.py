from plots import *
import os
import numpy as np
l=6
title = None

def find_complete_runs_old(dirlist):
    complete_dirlist = []
    for dir in dirlist:
        if not (os.path.isfile(dir+"q_values.txt") and os.path.isfile(dir+"costs.txt") and os.path.isfile(dir+"minimized.txt") and os.path.isfile(dir+"cellParameters.input")):
            continue
        errors = np.loadtxt(dir+"costs.txt")
        if len(errors) < 2:
            continue
        if errors[-1] > 5e-5:
            continue
        complete_dirlist.append(dir)
    return complete_dirlist
tolerance = 5e-5
def write_final_s0_periodic():
    for n in [2,4]:
        experiment_dir = f"data/{n}_cells_mean_l_6/"
        dir_list = find_complete_runs_old([experiment_dir+"{:03d}/".format(i) for i in range(100)])
        savefile = experiment_dir+ "final_s0.txt"
        write_final_s0(dir_list,savefile)
        for inc_or_dec in ["increase","decrease"]:
            experiment_dir = f"data/1_cell_{inc_or_dec}_2_sigma/"
            dir_list = find_complete_runs_old([experiment_dir+"{:03d}/".format(i) for i in range(100)])
            savefile = experiment_dir+ "final_s0.txt"
            write_final_s0(dir_list,savefile)

def write_stresses_periodic():
    for n in [4]:
        experiment_dir = f"data/{n}_cells_mean_l_6/"
        dir_list = find_complete_runs_old([experiment_dir+"{:03d}/".format(i) for i in range(100)])
        savefile = experiment_dir+ "final_stresses.txt"
        write_stresses(dir_list,savefile)
    for inc_or_dec in ["increase","decrease"]:
        experiment_dir = f"data/1_cell_{inc_or_dec}_2_sigma/"
        dir_list = find_complete_runs_old([experiment_dir+"{:03d}/".format(i) for i in range(100)])
        savefile = experiment_dir+ "final_stresses.txt"
        write_stresses(dir_list,savefile)

def error_to_iters_single_pattern():
        for n in [2,4]:
            dirlist = find_complete_runs_old([f"data/{n}_cells_mean_l_6/{i:03d}/" for i in range(100)])
            print("Yield: {}".format(len(dirlist)), dirlist[0])
            savefile = "Panels/Panel_S1/error_to_iters_periodic_n_{}.png".format(n)
            label_to_data = {"_"+r"$N_T=$"+"{}".format(n):{"dirlist":dirlist, "color":"black",}}
            # title = r"$n_{total} = $"+"{}".format(l**3)
            single_tracks_to_iters(label_to_data=label_to_data,title=title,savefile=savefile)

def error_to_iters_single_pattern_inc_dec():        
    # title = r"$n_{total} = $"+"{}".format(l**3)
    dirlist = find_complete_runs_old([f"data/1_cell_increase_2_sigma/{i:03d}/" for i in range(100)])
    savefile = "Panels/Panel_S1/error_to_iters_periodic_inc_dec.png"
    label_to_data = {r"$\sigma_T=\bar{\sigma}_{(0)}+2 \sigma_{(SD)}$":{"dirlist":dirlist, "color":"red","alpha":0.5}}
    dirlist = find_complete_runs_old([f"data/1_cell_decrease_2_sigma/{i:03d}/" for i in range(100)])
    label_to_data[r"$\sigma_T=\bar{\sigma}_{(0)}- 2 \sigma_{(SD)}$"]={"dirlist":dirlist, "color":"blue","alpha":0.5}
    single_tracks_to_iters(label_to_data=label_to_data,title=title,savefile=savefile,xlim=[1,10000])


def overlap_to_iters_single_pattern():
    for n in [2,4]:
        dirlist = find_complete_runs_old([f"data/{n}_cells_mean_l_6/{i:03d}/" for i in range(100)])
        savefile = "Panels/Panel_S1/overlap_to_iters_periodic_n_{}.png".format(n)
        # title = r"$n_{total} = $"+"{}".format(l**3)
        label_to_data = {"_"+r"$N_T=$"+"{}".format(n):{"dirlist":dirlist, "color":"black",}}
        single_tracks_to_iters(label_to_data=label_to_data,title=title,savefile=savefile,
                                input_filename="q_values.txt",
                                ylabel = r"$Q_2$",
                                ylog=False,
                                ylim=[0.55,1.02],
                                yticks=[0.6,0.8,1]
                                )
            
def overlap_to_iters_single_pattern_inc_dec():
        # title = r"$n_{total} = $"+"{}".format(l**3)
        n = 1
        savefile = "Panels/Panel_S1/overlap_to_iters_periodic_inc_dec.png"
        dirlist = find_complete_runs([f"data/1_cell_increase_2_sigma/{i:03d}/" for i in range(100)])
        label_to_data = {r"$\sigma_T=\bar{\sigma}_{(0)}+2 \sigma_{(SD)}$":{"dirlist":dirlist, "color":"red", "alpha":0.5}}
        dirlist = find_complete_runs([f"data/1_cell_decrease_2_sigma/{i:03d}/" for i in range(100)])
        label_to_data[r"$\sigma_T=\bar{\sigma}_{(0)} - 2 \sigma_{(SD)}$"]={"dirlist":dirlist, "color":"blue","alpha":0.5}
        single_tracks_to_iters(label_to_data=label_to_data,title=title,savefile=savefile,
                                    input_filename="q_values.txt",
                                    ylabel = r"$Q_2$",
                                    ylog=False,
                                    xlim=[1,10000],
                                    ylim=[0.55,1.02],
                                    yticks=[0.6,0.8,1])
def s0_histogram_periodic():
    input_filename = "final_s0.txt"

    for n in [2,4]:
        dir = "data/{}_cells_mean_l_6/".format(n)
        savefile = "Panels/Panel_S1/s0_histogram_periodic_n_{}.png".format(n)
        # label_to_data = {r"$s_0^{(trained)}$":{"dir":dir, "color":"red", "bins":35}}
        label_to_data = {"Trained":{"dir":dir, "color":"red", "bins":35}}
        # title = r"$n_{total} = $"+"{}".format(l**3)
        hist = histogram(   label_to_data=label_to_data,
                            input_filename=input_filename,
                            title=title,
                            xlim = [3.5,6.5],
                            ylim = [0,4])
        hist.ax.vlines(x = 5, ymin = 0, ymax = 2.6, linestyle= "dashed",color = "black", label = r"$s_0^{(0)}$",linewidth =15)
        hist.ax.legend()
        hist.save_fig(savefile)

def s0_histogram_periodic_inc_dec():
    input_filename = "final_s0.txt"
    label_to_data = {}
    savefile = "Panels/Panel_S1/s0_histogram_periodic_inc_dec.png"
    label_to_data["Trained: "+r"$\sigma_T=\bar{\sigma}_{(0)}+2 \sigma_{(SD)}$"]={"dir":"data/1_cell_increase_2_sigma/", "color":"red", "bins":100,"linewidth":5}
    label_to_data["Trained: "+r"$\sigma_T=\bar{\sigma}_{(0)}-2 \sigma_{(SD)}$"]={"dir":"data/1_cell_decrease_2_sigma/", "color":"blue", "bins":50,"linewidth":5}
    # title = r"$n_{total} = $"+"{}".format(l**3)
    hist = histogram(   label_to_data=label_to_data,
                        input_filename=input_filename,
                        title=title,
                        xlim = [4,6],
                        ylim = [0,25],
                        xticks=[4,5,6],
                        yticks = [0,10,20])
    hist.ax.vlines(x = 5, ymin = 0, ymax = 10, linestyle= "dashed",color = "black", label = r"$s_0^{(0)}$",linewidth =15)
    hist.ax.legend()
    hist.save_fig(savefile)

def stress_histogram_periodic():
    input_filename = "final_stresses.txt"
    for n in [2,4]:
        label_to_data = {}
        dir = "data/{}_cells_mean_l_6/".format(n)
        savefile = "Panels/Panel_S1/stress_histogram_periodic_n_{}.png".format(n)
        # label_to_data[r"$\sigma^T$"] = {"dir": "data/kv_10_l_{}_n_{}/".format(l,n), "color":"red", "bins":20}
        label_to_data["Trained"] = {"dir": "data/{}_cells_mean_l_6/".format(n), "color":"red", "bins":20}

        # title = r"$n_{total} = $"+"{}".format(l**3)
        hist = histogram(   label_to_data=label_to_data,
                            input_filename=input_filename,
                            title = title,
                            xlabel = r"$\sigma$",
                            xlim = [0,0.8],
                            xticks=[0.2*i for i in range(5)],
                            ylim = [0,5])
        hist.ax.vlines(x = np.mean(np.loadtxt(dir+"final_stresses.txt")), ymin = 0, ymax = 4.5, linestyle= "dotted",color = "black", label = r"$\sigma_T = \bar{\sigma}_{(0)}$",linewidth =25)
        hist.ax.legend()
        hist.save_fig(savefile)

def stress_histogram_periodic_inc_dec():
    input_filename = "final_stresses.txt"
    for l in [6]:
        n=1
        label_to_data ={}
        savefile = "Panels/Panel_S1/stress_histogram_periodic_inc_dec.png"
        # label_to_data = {r"$\sigma_T = \bar{\sigma}_{(0)}>$":{"dir":dir, "color":"red", "bins":20}}
        label_to_data[r"$\sigma_{(0)}$"]={"data":np.loadtxt("init/init_homogeneous_6/stresses.txt".format(l,n)), "color":"black", "bins":20,"linewidth":5,"alpha":0.5}

        label_to_data["_"+r"$\sigma_T=\bar{\sigma}_{(0)}+2 \sigma_{(SD)}$"]={"dir":"data/1_cell_increase_2_sigma/".format(l,n), "color":"red", "bins":50,"linewidth":5}
        label_to_data["_"+r"$\sigma_T=\bar{\sigma}_{(0)}-2 \sigma_{(SD)}$"]={"dir":"data/1_cell_decrease_2_sigma/".format(l,n), "color":"blue", "bins":20,"linewidth":5}
        # label_to_data[r"$\sigma_T=\bar{\sigma}_{(0)}+2 \sigma_{(SD)}$"]={"dir":"data/kv_10_l_{}_n_{}_increase/".format(l,n), "color":"red", "bins":50,"linewidth":5}
        # label_to_data[r"$\sigma_T=\bar{\sigma}_{(0)}-2 \sigma_{(SD)}$"]={"dir":"data/kv_10_l_{}_n_{}_decrease/".format(l,n), "color":"blue", "bins":20,"linewidth":5}
        # title = r"$n_{total} = $"+"{}".format(l**3)
        hist = histogram(   label_to_data=label_to_data,
                            input_filename=input_filename,
                            title = title,
                            xlabel = r"$\sigma$",
                            xlim = [0,1.2],
                            xticks=[0.2*i for i in range(10)],
                            ylim = [0,20],
                            yticks=[0,10,20])

        init_stresses = np.loadtxt("/Users/shabeebameen/Projects/tvm-fire/init/init_homogeneous_6/stresses.txt")
        hist.ax.vlines(x = np.mean(init_stresses)+2*np.std(init_stresses), ymin = 0, ymax = 15, linestyle= "dotted",color = "red",linewidth =25, label = r"$\sigma_T = \bar{\sigma}_{(0)} + 2 \sigma_{(SD)}$")
        hist.ax.vlines(x = np.mean(init_stresses)-2*np.std(init_stresses), ymin = 0, ymax = 15, linestyle= "dotted",color = "blue",linewidth =25, label = r"$\sigma_T = \bar{\sigma}_{(0)} - 2 \sigma_{(SD)}$")
        # hist.ax.vlines(x = np.mean(init_stresses)+2*np.std(init_stresses), ymin = 0, ymax = 5, linestyle= "dotted",color = "red",linewidth =25)
        # hist.ax.vlines(x = np.mean(init_stresses)-2*np.std(init_stresses), ymin = 0, ymax = 5, linestyle= "dotted",color = "blue",linewidth =25)
        hist.ax.legend()
        hist.save_fig(savefile)

def main():
    os.makedirs("Panels/Panel_S1",exist_ok=True)
    # write_stresses_periodic()
    # write_stresses_inc_dec()
    write_final_s0_periodic()
    error_to_iters_single_pattern()
    error_to_iters_single_pattern_inc_dec()
    overlap_to_iters_single_pattern_inc_dec()
    overlap_to_iters_single_pattern()
    s0_histogram_periodic()
    s0_histogram_periodic_inc_dec()
    stress_histogram_periodic()
    stress_histogram_periodic_inc_dec()
    
if __name__ == "__main__":
    # for i in range(100):
    #     os.system("cp init/kv_10_l_6/000/conf data/kv_10_l_6_n_1_increase/{:03d}/".format(i))
    #     os.system("cp init/kv_10_l_6/000/conf data/kv_10_l_6_n_1_decrease/{:03d}/".format(i))
    main()