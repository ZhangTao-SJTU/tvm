from plots import *
import os
colors = ['#543005','#dfc27d','#80cdc1','#003c30']
colors = ['#3E215D','#284E78','#5D8233','#ECD662']
colors = ['#573583','#2c6d0d','#c36001','#f3a850'

]
colors = ['#573583','#400602','#c36001','#f3a850']
colors = ['#7f266b','#2c6d0d','#904220','#f3a850']
colors = ['purple','green',"orange"]

#7f266b


n_cells = [10,20,40][::-1]
n_spheroid_color_map = dict(zip(n_cells,colors))
alphas = [1,0.7,0.2]
bins = [35,15,20]

def error_to_iters_spheroid():
    for l in [5,6]:
        savefile = "Panels/Panel_5/error_to_iters_spheroid_l_{}.png".format(l)
        label_to_data = {}
        title = r"$n_{total} = $"+"{}".format(l**3)

        for i,n in enumerate([10,20,40]):
            dirlist = find_complete_runs(["data/kv_10_l_{}_n_sp_{:03d}/{:03d}/".format(l,n,i) for i in range(100)])
            label_to_data[r"$n_{sp}=$"+"{}".format(n)]={"dirlist":dirlist, "color":n_spheroid_color_map[n], "alpha":alphas[i]}
        single_tracks_to_iters(label_to_data=label_to_data,title=title,savefile=savefile, )

def s0_histogram_spheroid_combined():
    input_filename = "final_s0.txt"
    for l in [5,6]:
        label_to_data = {}
        savefile = "Panels/Panel_5/s0_histogram_spheroid_l_{}.png".format(l)
        for i,n in enumerate([10,40]):
            dir = "data/kv_10_l_{}_n_sp_{:03d}/".format(l,n)
            label_to_data[r"$n_{sp}=$"+"{}".format(n)] = {"dir":dir, "color":n_spheroid_color_map[n], "bins":20,"alpha":0.8}
        title = r"$n_{total} = $"+"{}".format(l**3)
        hist = histogram(   label_to_data=label_to_data,
                            input_filename=input_filename,
                            title=title,
                            xlim = [3,7],
                            ylim = [0,3.5])
        hist.ax.vlines(x = 5, ymin = 0, ymax = 3.2, linestyle= "dashed",color = "black", label = r"$s_0^{(init)}$",linewidth =15)
        hist.ax.legend()
        hist.save_fig(savefile)

def s0_histogram_spheroid_for_inset():
    input_filename = "final_s0.txt"
    l = 6
    n = 20
    label_to_data = {}
    dir = "data/kv_10_l_{}_n_sp_{:03d}/".format(l,n)
    label_to_data[r"$n_{sp}=$"+"{}".format(n)] = {"dir":dir, "color":n_spheroid_color_map[n], "bins":20,"alpha":0.8}
    hist = histogram(   label_to_data=label_to_data,
                        input_filename=input_filename,
                        xlim = [3.5,6.5],
                        ylim = [0,2],
                        xticks= [4,5,6],
                        yticks= [0,1])
    # hist.ax.vlines(x = 5, ymin = 0, ymax = 3.2, linestyle= "dashed",color = "black", label = r"$s_0^{(init)}$",linewidth =15)
    return hist.ax
def SD_s0_scatter_spheroid():
    for l in [5,6]:  
        savefile = "Panels/Panel_5/SD_s0_scatter_l_{}.png".format(l)
        label_to_dir = {}
        # for n in [1,2,3,4,5,6]:
        #     dirlist = find_complete_runs(["data/kv_10_l_{}_n_{}/{:03d}/".format(l,n,i) for i in range(100)])
        x_array = [5,10,20,40]
        y_array = [np.std(np.loadtxt("data/kv_10_l_{}_n_sp_{:03d}/final_s0.txt".format(l,n))) for n in x_array]
        label_to_dir= {"_none": {"x_array":x_array,"y_array":y_array,"marker":"o","color":"black","s":3000,"alpha":0.5}}
        plotter = scatter_plot(label_to_dir,
                     xlim =[0,45],
                     ylim=[0,1.3],
                     title = r"$n_{total}=$"+"{}".format(l**3),
                     xticks= x_array,
                     xlog=False,
                     xlabel=r"$n_{sp}$",
                     ylabel=r"$SD(s_0)$")
        plotter.save_fig(savefile)

def SD_s0_scatter_spheroid_for_inset():
    l = 6  
    label_to_dir = {}
    # for n in [1,2,3,4,5,6]:
    #     dirlist = find_complete_runs(["data/kv_10_l_{}_n_{}/{:03d}/".format(l,n,i) for i in range(100)])
    x_array = [5,10,20,40]
    y_array = [np.std(np.loadtxt("data/kv_10_l_{}_n_sp_{:03d}/final_s0.txt".format(l,n))) for n in x_array]
    label_to_dir= {"_none": {"x_array":x_array,"y_array":y_array,"marker":"o","color":"black","s":3000,"alpha":0.5}}
    plotter = scatter_plot(label_to_dir,
                    xlim =[0,45],
                    ylim=[0,1.3],
                     title = r"$n_{total}=$"+"{}".format(l**3),
                     xticks= x_array,
                     xlog=False,
                     xlabel=r"$n_{sp}$",
                     ylabel=r"$SD(s_0)$")
    return plotter
def final_iteration_to_final_overlap_spheroid():
    for l in [5,6]:  
        savefile = "Panels/Panel_5/final_iterations_to_final_overlap_l_{}.png".format(l)
        label_to_dir = {}
        # for n in [1,2,3,4,5,6]:
        for n in [10,20,40]:

            dirlist = find_complete_runs(["data/kv_10_l_{}_n_sp_{:03d}/{:03d}/".format(l,n,i) for i in range(100)])
            x_array = [len(np.loadtxt(dir+"costs.txt"))for dir in dirlist]
            y_array = [np.loadtxt(dir+"q_values.txt")[-1] for dir in dirlist]
            label_to_dir[r"$n_{sp}=$"+"{}".format(n)]= {"x_array":x_array,"y_array":y_array,"marker":"o","color":n_spheroid_color_map[n], "s":2000, "alpha":0.5}
        plotter = scatter_plot(label_to_dir,
                               ylim=[0.6,1.05],
                               xlim = [50,20000],
                                title = r"$n_{total}=$"+"{}".format(l**3),)
        plotter.ax.legend(ncols = 2)
        plotter.save_fig(savefile)

def area_change_to_stress_change():
    for l in [5,6]:  
        savefile = "Panels/Panel_5/area_change_to_stress_change_l_{}.png".format(l)
        label_to_dir = {}
        # for n in [1,2,3,4,5,6]:
        for n in [10,20,40]:

            dirlist = find_complete_runs(["data/kv_10_l_{}_n_sp_{:03d}/{:03d}/".format(l,n,i) for i in range(100)])
            
            x_array = [len(np.loadtxt(dir+"costs.txt"))for dir in dirlist]
            y_array = [np.loadtxt(dir+"q_values.txt")[-1] for dir in dirlist]
            for dir in dirlist:
                if not os.path.isfile(dir+"areas.txt"):
                    continue
                areas = np.loadtxt(dir+"areas.txt")
                stresses = pd.read_csv(dir+"0000000.stresses.csv")
                y_array.append(areas[1]/areas[0]-1)
                x_array.append(stresses["Target"].to_list()[0]/stresses["Current"].to_list()[0] - 1)
   
                                   
            label_to_dir[r"$n_{sp}=$"+"{}".format(n)]= {"x_array":x_array,"y_array":y_array,"marker":"o","color":n_spheroid_color_map[n], "s":3000, "alpha":0.5}
        plotter = scatter_plot(label_to_dir,
                               ylim=[-0.12,0.1],
                               xlim = [-1,1],
                                title = r"$n_{total}=$"+"{}".format(l**3),
                                xlabel=r"$(\sigma_T - \sigma_{0})/\sigma_{0}$",
                                ylabel=r"$(S_{T} - S_{0})/S_{0}$",
                                xlog=False,
                                xticks=[0.5*i for i in range(-2,3)],
                                yticks=[0.1*i for i in range(-2,3)])

        plotter.ax.legend(ncol=2)
        plotter.save_fig(savefile)
def inset_s0_plot():
    for l in [5,6]:
        main_plot_file = "Panel_5/SD_s0_scatter_l_{}.png".format(l)
        inset_plot_file = "Panel_5/s0_histogram_spheroid_l_{}.png".format(l)
        filename ="Panel_5/SD_s0_with_inset_l_{}.png".format(l)
        width = 0.55
        inset_pos =[0.38, 0.32,width,width]
        create_inset_plot(main_plot_file,inset_plot_file,inset_pos,filename)
def main():
    os.makedirs("Panels/Panel_5",exist_ok=True)
    # error_to_iters_spheroid()
    # s0_histogram_spheroid_combined()
    # SD_s0_scatter_spheroid()
    # final_iteration_to_final_overlap_spheroid()
    # area_change_to_stress_change()
    create_inset_SD_s0_plot(SD_s0_scatter_spheroid_for_inset(),s0_histogram_spheroid_for_inset(),bbox_to_anchor=(0.4, 0.4, 1.2, 1.2),filename="Panels/Panel_5/SD_s0_with_inset.png")
if __name__ == "__main__":
    main()