from plots import *
import os
import copy
title = None

def error_to_iters_single_pattern():
    for l in [4,5,6]:
        label_to_data = {}
        for n in [1,2,4]:
            savefile = "Panels/Panel_2/error_to_iters_periodic_l_{}.png".format(l)
            # title = r"$n_{total} = $"+"{}".format(l**3)
            dirlist = find_complete_runs(["data/kv_10_l_{}_n_{}/{:03d}/".format(l,n,i) for i in range(100)])
            label_to_data[r"$N_T=$"+"{}".format(n)]={"dirlist":dirlist, "color":n_cells_color_map[n], "alpha":0.6, "s":4000}
          
        single_tracks_to_iters(label_to_data=label_to_data,title=title,savefile=savefile,
                               xlim=[1,5000],yticks=[5*i for i in range(1,5000)])

def s0_histogram_periodic_combined():
    input_filename = "final_s0.txt"
    for l in [4,5,6]:
        label_to_data = {}
        savefile = "Panels/Panel_2/s0_histogram_periodic_l_{}.png".format(l)

        for n in [4]:
            dir = "data/kv_10_l_{}_n_{}/".format(l,n)
            label_to_data[r"$N_T=$"+"{}".format(n)] = {"dir":dir, "color":n_cells_color_map[n], "bins":30}
        # title = r"$n_{total} = $"+"{}".format(l**3)
        hist = histogram(   label_to_data=label_to_data,
                            input_filename=input_filename,
                            title=title,
                            xlim = [3.5,6.5],
                            ylim = [0,3.5])
        hist.ax.vlines(x = 5, ymin = 0, ymax = 3.2, linestyle= "dashed",color = "black", label = r"$s_0^{(init)}$",linewidth =15)

        hist.ax.legend()
        hist.save_fig(savefile)

def s0_histogram_periodic_for_inset():
    input_filename = "final_s0.txt"
    l = 4
    n = 4        
    label_to_data = {}

    dir = "data/kv_10_l_{}_n_{}/".format(l,n)
    label_to_data[r"$N_T=$"+"{}".format(n)] = {"dir":dir, "color":n_cells_color_map[n], "bins":30}
        
    hist = histogram(   label_to_data=label_to_data,
                        input_filename=input_filename,
                        xlim = [3.5,6.5],
                        ylim = [0,2],
                        xticks= [4,5,6],
                        yticks= [0,1])
    # hist.ax.vlines(x = 5, ymin = 0, ymax = 3.2, linestyle= "dashed",color = "black", label = r"$s_0^{(init)}$",linewidth =15)
    # hist.ax.legend()
    return hist.ax


def final_iteration_to_final_overlap_periodic():
    for l in [4,5,6]:  
        savefile = "Panels/Panel_2/final_iterations_to_final_overlap_l_{}.png".format(l)
        label_to_dir = {}
        # for n in [1,2,3,4,5,6]:
        for n in [1,2,4]:

            dirlist = find_complete_runs(["data/kv_10_l_{}_n_{}/{:03d}/".format(l,n,i) for i in range(100)])
            x_array = [len(np.loadtxt(dir+"costs.txt"))for dir in dirlist]
            y_array = [np.loadtxt(dir+"q_values.txt")[-1] for dir in dirlist]
            label_to_dir[r"$N_T=$"+"{}".format(n)]= {"x_array":x_array,"y_array":y_array,"marker":"o","color":n_cells_color_map[n], "s":4000, "alpha":0.6}
        plotter = scatter_plot(label_to_dir,
                               xlim = [10,5000],
                               ylim=[0.55,1.05],
                               title=title)
        plotter.ax.legend()
        plotter.save_fig(savefile)

def SD_s0_scatter_periodic():
    for l in [4,5,6]:  
        savefile = "Panels/Panel_2/SD_s0_scatter_l_{}.png".format(l)
        label_to_dir = {}
        # for n in [1,2,3,4,5,6]:
        #     dirlist = find_complete_runs(["data/kv_10_l_{}_n_{}/{:03d}/".format(l,n,i) for i in range(100)])
        x_array = [1,2,3,4,5,6]
        y_array = [np.std(np.loadtxt("data/kv_10_l_{}_n_{}/final_s0.txt".format(l,n))) for n in x_array]
        label_to_dir= {"_none": {"x_array":x_array,"y_array":y_array,"marker":"o","color":"black","s":4000}}
        plotter = scatter_plot(label_to_dir,
                     xlim =[0.5,6.5],
                     ylim=[0,0.5],
                     title=title,
                    #  title = r"$n_{total}=$"+"{}".format(l**3),
                     xticks= x_array,
                     xlog=False,
                     xlabel=r"$N_T$",
                     ylabel=r"$SD(s_0)$")
        plotter.save_fig(savefile)

def SD_s0_scatter_periodic_for_inset():
    l = 4
    label_to_dir = {}
    x_array = [1,2,3,4,5,6]
    y_array = [np.std(np.loadtxt("data/kv_10_l_{}_n_{}/final_s0.txt".format(l,n))) for n in x_array]
    label_to_dir= {"_none": {"x_array":x_array,"y_array":y_array,"marker":"o","color":"black","s":4000, "alpha":0.6}}
    plotter = scatter_plot(label_to_dir,
                    xlim =[0.5,6.5],
                    ylim=[0,0.5],
                    title = title,
                    # title = r"$n_{total}=$"+"{}".format(l**3),
                    xticks= x_array,
                    xlog=False,
                    xlabel=r"$N_T$",
                    ylabel=r"$SD(s_0)$",
                
        )
    return plotter
def SD_s0_scatter_periodic_all_sizes():
    savefile = "Panels/Panel_2/SD_s0_scatter.png"
    label_to_dir = {}
    for l in [4,5,6]:  
        x_array = [1,2,3,4,5,6]
        y_array = [np.std(np.loadtxt("data/kv_10_l_{}_n_{}/final_s0.txt".format(l,n))) for n in x_array]
        label_to_dir[r"$N = $"+"{}".format(l**3)] = {"x_array":x_array,"y_array":y_array,"marker":l_to_marker[l],"color":"black","alpha": 0.8, "color":l_to_color[l], "s":4000}
    plotter = scatter_plot(label_to_dir,
                    xlim =[0.5,6.5],
                    ylim=[0,0.5],
                    # title = r"$n_{total}=$"+"{}".format(l**3),
                    xticks= x_array,
                    yticks= [0,0.2,0.4],
                    xlog=False,
                    xlabel=r"$N_T$",
                    ylabel=r"$SD(s_0)$")
    plotter.ax.legend()
    plotter.save_fig(savefile)

# def create_inset_SD_s0_plot():
#     plotter_main = SD_s0_scatter_periodic_for_inset()
#     fig,ax = s0_histogram_periodic_for_inset()
#     fig.savefig("test_inset.png")
#     # axins = inset_axes(plotter_main.ax, width="45%", height="45%", 
#     #                    bbox_to_anchor=(0.45,0.15,1,1),
#     #                bbox_transform=plotter_main.ax.transAxes,
#     #                loc="lower left")
#     # # axins.set_position([0.55, 0., 0.5, 0.5])

#     # # copy artists from ax2 to axins
#     # for line in plotter_inset.ax.get_lines():
#     #     axins.plot(
#     #         line.get_xdata(),
#     #         line.get_ydata(),
#     #         linestyle=line.get_linestyle(),
#     #         linewidth=line.get_linewidth(),
#     #         color=line.get_color(),
#     #         marker=line.get_marker(),
#     #     )
#     # for patch in plotter_inset.ax.patches:
#     #     new_patch = copy.copy(patch)
#     #     new_patch.set_transform(axins.transData)  # attach to inset axes coordinates
#     #     axins.add_patch(new_patch)
#     #     # copy limits
#     #     axins.set_xlim(plotter_inset.ax.get_xlim())
#     #     axins.set_ylim(plotter_inset.ax.get_ylim())
#     # # optionally copy labels
#     # axins.set_xticks(plotter_inset.ax.get_xticks())
#     # axins.set_yticks(plotter_inset.ax.get_yticks())
#     # # optionally copy labels
#     # axins.set_xlabel(plotter_inset.ax.get_xlabel())
#     # axins.set_ylabel(plotter_inset.ax.get_ylabel())
#     # plotter_main.save_fig("test.png")

    
#     axins = inset_axes(
#         plotter_main.ax,
#         width="45%",
#         height="45%",
#         bbox_to_anchor=(0.45,0.15,1,1),
#         bbox_transform=plotter_main.ax.transAxes,
#         loc="lower left"
#     )

#     # copy lines (your vertical dashed line)
#     for line in ax.lines:
#         new_line = copy.copy(line)
#         new_line.axes = None
#         new_line.figure = None
#         axins.add_line(new_line)

#     # copy histogram bars
#     for patch in ax.patches:
#         new_patch = copy.copy(patch)
#         new_patch.axes = None
#         new_patch.figure = None
#         axins.add_patch(new_patch)

#     # limits and ticks
#     # axins.set_xlim(plotter_inset.ax.get_xlim())
#     # axins.set_ylim(plotter_inset.ax.get_ylim())
#     # axins.set_xticks(plotter_inset.ax.get_xticks())
#     # axins.set_yticks(plotter_inset.ax.get_yticks())
#     # axins.set_xlabel(plotter_inset.ax.get_xlabel())
#     # axins.set_ylabel(plotter_inset.ax.get_ylabel())
#     plotter_main.save_fig("test.png")





def main():
    os.makedirs("Panels/Panel_2",exist_ok=True)
    SD_s0_scatter_periodic()
    create_inset_SD_s0_plot(SD_s0_scatter_periodic_for_inset(), s0_histogram_periodic_for_inset())
    error_to_iters_single_pattern()
    final_iteration_to_final_overlap_periodic()
    s0_histogram_periodic_combined()
    SD_s0_scatter_periodic_all_sizes()

if __name__ == "__main__":
    main()