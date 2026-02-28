from plots import *
import os


def error_to_iters_single_pattern():
    for l in [4,5,6]:
        label_to_data = {}

        for n in [1,2,4]:
            savefile = "Panel_2/error_to_iters_periodic_l_{}.png".format(l)
            title = r"$n_{total} = $"+"{}".format(l**3)
            dirlist = find_complete_runs(["data/kv_10_l_{}_n_{}/{:03d}/".format(l,n,i) for i in range(100)])
            label_to_data[r"$n_T=$"+"{}".format(n)]={"dirlist":dirlist, "color":n_cells_color_map[n], "alpha":0.6}
          
        single_tracks_to_iters(label_to_data=label_to_data,title=title,savefile=savefile,
                               xlim=[1,5000],yticks=[5*i for i in range(1,5000)])

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

def final_iteration_to_final_overlap_periodic():
    for l in [4,5,6]:  
        savefile = "Panel_2/final_iterations_to_final_overlap_l_{}.png".format(l)
        label_to_dir = {}
        # for n in [1,2,3,4,5,6]:
        for n in [1,2,4]:

            dirlist = find_complete_runs(["data/kv_10_l_{}_n_{}/{:03d}/".format(l,n,i) for i in range(100)])
            x_array = [len(np.loadtxt(dir+"costs.txt"))for dir in dirlist]
            y_array = [np.loadtxt(dir+"q_values.txt")[-1] for dir in dirlist]
            label_to_dir[r"$n_T=$"+"{}".format(n)]= {"x_array":x_array,"y_array":y_array,"marker":"o","color":n_cells_color_map[n], "s":3000, "alpha":0.6}
        plotter = scatter_plot(label_to_dir,
                               xlim = [10,5000],
                               ylim=[0.55,1.05],
                                title = r"$n_{total}=$"+"{}".format(l**3))
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
        label_to_dir= {"_none": {"x_array":x_array,"y_array":y_array,"marker":"o","color":"black","s":4000}}
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
        x_array = [1,2,3,4,5,6]
        y_array = [np.std(np.loadtxt("data/kv_10_l_{}_n_{}/final_s0.txt".format(l,n))) for n in x_array]
        label_to_dir[r"$n_{total} = $"+"{}".format(l**3)] = {"x_array":x_array,"y_array":y_array,"marker":l_to_marker[l],"color":"black","alpha": l_to_alpha[l], "color":l_to_color[l], "s":4000}
    plotter = scatter_plot(label_to_dir,
                    xlim =[0.5,6.5],
                    ylim=[0,0.5],
                    # title = r"$n_{total}=$"+"{}".format(l**3),
                    xticks= x_array,
                    yticks= [0,0.2,0.4],
                    xlog=False,
                    xlabel=r"$n_T$",
                    ylabel=r"$SD(s_0)$")
    plotter.ax.legend()
    plotter.save_fig(savefile)

def main():
    os.makedirs("Panel_2",exist_ok=True)
    error_to_iters_single_pattern()
    final_iteration_to_final_overlap_periodic()
    SD_s0_scatter_periodic()
    s0_histogram_periodic_combined()
    SD_s0_scatter_periodic_all_sizes()

if __name__ == "__main__":
    main()