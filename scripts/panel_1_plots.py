from plots import *
import os

tolerance = 1e-6
def write_final_s0_periodic():
    for l in [4,5,6]:
        header = "data/kv_10_l_{}_".format(l)
        for n in [1,2,3,4,5,6]:
            experiment_dir = header+"n_{}/".format(n)
            dir_list = find_complete_runs([experiment_dir+"{:03d}/".format(i) for i in range(100)])
            savefile = experiment_dir+ "final_s0.txt"
            write_final_s0(dir_list,savefile)

            if n == 1:
                        
                experiment_dir = header+"n_{}_increase/".format(n)
                dir_list = find_complete_runs([experiment_dir+"{:03d}/".format(i) for i in range(100)])
                # print(dir_list)
                savefile = experiment_dir+ "final_s0.txt"
                write_final_s0(dir_list,savefile)
                experiment_dir = header+"n_{}_decrease/".format(n)
                dir_list = find_complete_runs([experiment_dir+"{:03d}/".format(i) for i in range(100)])
                savefile = experiment_dir+ "final_s0.txt"
                write_final_s0(dir_list,savefile)

def write_stresses_periodic():
    for l in [6]:
        header = "data/kv_10_l_{}_".format(l)
        # for n in [1,2,3,4,5,6]:
        for n in [1]:
            if n == 1:
                experiment_dir = header+"n_{}_increase/".format(n)
                dir_list = find_complete_runs([experiment_dir+"{:03d}/".format(i) for i in range(100)])
                savefile = experiment_dir+ "final_stresses.txt"
                write_stresses(dir_list,savefile)
                experiment_dir = header+"n_{}_decrease/".format(n)
                dir_list = find_complete_runs([experiment_dir+"{:03d}/".format(i) for i in range(100)])
                savefile = experiment_dir+ "final_stresses.txt"
                write_stresses(dir_list,savefile)

            # experiment_dir = header+"n_{}/".format(n)
            # dir_list = find_complete_runs([experiment_dir+"{:03d}/".format(i) for i in range(100)])
            # savefile = experiment_dir+ "final_stresses.txt"
            # write_stresses(dir_list,savefile)

def error_to_iters_single_pattern():
    for l in [4,5,6]:
        for n in [2,4]:
            dirlist = find_complete_runs(["data/kv_10_l_{}_n_{}/{:03d}/".format(l,n,i) for i in range(100)])
            savefile = "Panel_1/error_to_iters_periodic_l_{}_n_{}.png".format(l,n)
            label_to_data = {r"$n_T=$"+"{}".format(n):{"dirlist":dirlist, "color":"black",}}
            title = r"$n_{total} = $"+"{}".format(l**3)
            single_tracks_to_iters(label_to_data=label_to_data,title=title,savefile=savefile)

def error_to_iters_single_pattern_inc_dec():
    for l in [4,5,6]:
        n=1
        title = r"$n_{total} = $"+"{}".format(l**3)
        dirlist = find_complete_runs(["data/kv_10_l_{}_n_{}_increase/{:03d}/".format(l,n,i) for i in range(100)])
        savefile = "Panel_1/error_to_iters_periodic_l_{}_n_{}_inc_dec.png".format(l,n)
        label_to_data = {r"$\sigma_T=\bar{\sigma}_{(0)}+2 \sigma^{(SD)}$":{"dirlist":dirlist, "color":"red","alpha":0.5}}
        dirlist = find_complete_runs(["data/kv_10_l_{}_n_{}_decrease/{:03d}/".format(l,n,i) for i in range(100)])
        label_to_data[r"$\sigma_T=\bar{\sigma}_{(0)}- 2 \sigma^{(SD)}$"]={"dirlist":dirlist, "color":"blue","alpha":0.5}
        single_tracks_to_iters(label_to_data=label_to_data,title=title,savefile=savefile,xlim=[1,100000])


def overlap_to_iters_single_pattern():
    for l in [4,5,6]:
        for n in [2,4]:
            dirlist = find_complete_runs(["data/kv_10_l_{}_n_{}/{:03d}/".format(l,n,i) for i in range(100)])
            savefile = "Panel_1/overlap_to_iters_periodic_l_{}_n_{}.png".format(l,n)
            title = r"$n_{total} = $"+"{}".format(l**3)
            label_to_data = {r"$n_T=$"+"{}".format(n):{"dirlist":dirlist, "color":"black",}}

            single_tracks_to_iters(label_to_data=label_to_data,title=title,savefile=savefile,
                                   input_filename="q_values.txt",
                                   ylabel = r"$Q_2$",
                                   ylog=False,
                                   ylim=[0.55,1.02],
                                   yticks=[0.6,0.8,1]
                                   )
def overlap_to_iters_single_pattern_inc_dec():
    for l in [4,5,6]:
        title = r"$n_{total} = $"+"{}".format(l**3)
        n = 1
        savefile = "Panel_1/overlap_to_iters_periodic_l_{}_n_{}_inc_dec.png".format(l,n)
        dirlist = find_complete_runs(["data/kv_10_l_{}_n_{}_increase/{:03d}/".format(l,n,i) for i in range(100)])
        label_to_data = {r"$\sigma_T=\bar{\sigma}_{(0)}+2 \sigma^{(SD)}$":{"dirlist":dirlist, "color":"red", "alpha":0.5}}
        dirlist = find_complete_runs(["data/kv_10_l_{}_n_{}_decrease/{:03d}/".format(l,n,i) for i in range(100)])
        label_to_data[r"$\sigma_T=\bar{\sigma}_{(0)} - 2 \sigma^{(SD)}$"]={"dirlist":dirlist, "color":"blue","alpha":0.5}
        single_tracks_to_iters(label_to_data=label_to_data,title=title,savefile=savefile,
                                    input_filename="q_values.txt",
                                    ylabel = r"$Q_2$",
                                    ylog=False,
                                    ylim=[0.55,1.02],
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
                                ylim = [0,4])
            hist.ax.vlines(x = 5, ymin = 0, ymax = 2.6, linestyle= "dashed",color = "black", label = r"$s_0^{(0)}$",linewidth =15)
            hist.ax.legend()
            hist.save_fig(savefile)

def s0_histogram_periodic_inc_dec():
    input_filename = "final_s0.txt"
    for l in [4,5,6]:
        n = 1
        label_to_data = {}
        savefile = "Panel_1/s0_histogram_periodic_l_{}_n_{}_inc_dec.png".format(l,n)
        label_to_data[r"$\sigma_T=\bar{\sigma}_{(0)}+2 \sigma^{(SD)}$"]={"dir":"data/kv_10_l_{}_n_{}_increase/".format(l,n), "color":"red", "bins":50,"linewidth":5}
        label_to_data[r"$\sigma_T=\bar{\sigma}_{(0)}-2 \sigma^{(SD)}$"]={"dir":"data/kv_10_l_{}_n_{}_decrease/".format(l,n), "color":"blue", "bins":20,"linewidth":5}

        title = r"$n_{total} = $"+"{}".format(l**3)
        hist = histogram(   label_to_data=label_to_data,
                            input_filename=input_filename,
                            title=title,
                            xlim = [4.5,5.5],
                            ylim = [0,20],
                            xticks=[4.5,5,5.5],
                            yticks = [0,10,20])
        hist.ax.vlines(x = 5, ymin = 0, ymax = 10, linestyle= "dashed",color = "black", label = r"$s_0^{(0)}$",linewidth =15)
        hist.ax.legend()
        hist.save_fig(savefile)

def stress_histogram_periodic():
    input_filename = "final_stresses.txt"
    for l in [6]:
        for n in [2,4]:
            label_to_data = {}
            dir = "data/kv_10_l_{}_n_{}/".format(l,n)
            savefile = "Panel_1/stress_histogram_periodic_l_{}_n_{}.png".format(l,n)
            label_to_data[r"$\sigma_{(trained)}$"] = {"dir": "data/kv_10_l_{}_n_{}/".format(l,n), "color":"red", "bins":20}

            title = r"$n_{total} = $"+"{}".format(l**3)
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
        savefile = "Panel_1/stress_histogram_periodic_l_{}_n_{}_inc_dec.png".format(l,n)
        # label_to_data = {r"$\sigma_T = \bar{\sigma}_{(0)}>$":{"dir":dir, "color":"red", "bins":20}}
        label_to_data[r"$\sigma_T=\bar{\sigma}_{(0)}+2 \sigma^{(SD)}$"]={"dir":"data/kv_10_l_{}_n_{}_increase/".format(l,n), "color":"red", "bins":50,"linewidth":5}
        label_to_data[r"$\sigma_T=\bar{\sigma}_{(0)}-2 \sigma^{(SD)}$"]={"dir":"data/kv_10_l_{}_n_{}_decrease/".format(l,n), "color":"blue", "bins":20,"linewidth":5}

        title = r"$n_{total} = $"+"{}".format(l**3)
        hist = histogram(   label_to_data=label_to_data,
                            input_filename=input_filename,
                            title = title,
                            xlabel = r"$\sigma$",
                            xlim = [0,0.8],
                            xticks=[0.2*i for i in range(5)],
                            ylim = [0,7])
        init_stresses = np.loadtxt("/Users/shabeebameen/Projects/tvm-fire/init/kv_10_l_{}/stresses.txt".format(l))
        hist.ax.vlines(x = np.mean(init_stresses)+2*np.std(init_stresses), ymin = 0, ymax = 5, linestyle= "dotted",color = "red",linewidth =25,alpha =1)
        hist.ax.vlines(x = np.mean(init_stresses)-2*np.std(init_stresses), ymin = 0, ymax = 5, linestyle= "dotted",color = "blue",linewidth =25)

        hist.ax.legend()
        hist.save_fig(savefile)

def main():
    os.makedirs("Panel_1",exist_ok=True)
    # write_stresses_periodic()
    write_final_s0_periodic()
    error_to_iters_single_pattern()
    error_to_iters_single_pattern_inc_dec()
    overlap_to_iters_single_pattern_inc_dec()
    overlap_to_iters_single_pattern()
    s0_histogram_periodic()
    s0_histogram_periodic_inc_dec()
    stress_histogram_periodic()
    # stress_histogram_periodic_inc_dec()
    
if __name__ == "__main__":
    # for i in range(100):
    #     os.system("cp init/kv_10_l_6/000/conf data/kv_10_l_6_n_1_increase/{:03d}/".format(i))
    #     os.system("cp init/kv_10_l_6/000/conf data/kv_10_l_6_n_1_decrease/{:03d}/".format(i))
    main()