from toolbox import manuscriptPlots
from toolbox.periodic import PeriodicTissue
from toolbox.minimization import FIREminimization
from toolbox import stress
import pandas as pd
import numpy as np
import os

# def s0_histogram(dir, final_iter):
#     plotter = manuscriptPlots.plot()
#     plotter.set_ylim(0,40)
#     plotter.set_xlim(4.5,5.4)
#     plotter.set_xticks([0.25*i for i in range(150)])
#     plotter.set_yticks([5*i for i in range(1,100)])
#     plotter.set_xlabel(r"$s_0$")
#     plotter.set_yScaled()
#     plotter.initialize_figure()
#     df = pd.read_csv("{}cellParameters.{}.input".format(dir,final_iter), sep=" ",header=None)
#     plotter.histogram_from_dataframe(df[2],color = "blue", bins = 20, label = "_final")
#     plotter.save_fig("{}/s0_histogram.png".format(dir))

def s0_histogram(dir, final_iter,filename = "s0_histogram.png"):
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(0,40)
    plotter.set_xlim(4.5,5.4)
    plotter.set_xticks([0.25*i for i in range(150)])
    plotter.set_yticks([5*i for i in range(1,100)])
    plotter.set_xlabel(r"$s_0$")
    plotter.set_yScaled()
    plotter.initialize_figure()
    df = pd.read_csv("{}{:04d}.cellParameters.input".format(dir,final_iter), sep=" ",header=None)
    plotter.histogram_from_dataframe(df[2],color = "blue", bins = 20, label = "_final")
    plotter.save_fig("{}{}".format(dir,filename))

# def combined_stress_histogram(dir, final_iter):
#     plotter = manuscriptPlots.plot()
#     plotter.set_ylim(0,7.5)
#     plotter.set_xlim(0,0.75)
#     plotter.set_xticks([0.25*i for i in range(150)])
#     plotter.set_yticks([5*i for i in range(1,100)])
#     plotter.set_xlabel(r"$\sigma_{shear}$")
#     plotter.set_yScaled()
#     plotter.initialize_figure()

#     #init
#     file = "init_config.txt"
#     sample = PeriodicTissue.from_config(dir, file)
#     min = FIREminimization.periodic_tissue(sample)
#     min.load_cell_parameters("cellParameters.init.input")
#     stresses = []
#     for cellID,cell in sample.cells_.items():
#         cell.max_shear_stress_ = stress.calculate_max_shear_stress(sample,cellID)
#         stresses.append(cell.max_shear_stress_)
#     df = pd.DataFrame({"stress": stresses})
#     plotter.histogram_from_dataframe(df["stress"],color = "red", bins = 20, label = "initial")

#     #final
#     file = "{}.bulk.txt".format(final_iter)
#     sample = PeriodicTissue.from_config(dir, file)
#     min = FIREminimization.periodic_tissue(sample)
#     min.load_cell_parameters("cellParameters.{}.input".format(final_iter))
#     stresses = []
#     for cellID,cell in sample.cells_.items():
#         cell.max_shear_stress_ = stress.calculate_max_shear_stress(sample,cellID)
#         stresses.append(cell.max_shear_stress_)
#     df = pd.DataFrame({"stress": stresses})
#     plotter.histogram_from_dataframe(df["stress"],color = "blue", bins = 20, label = "final")
#     plotter.save_fig("{}stress_histogram.png".format(dir))

def combined_stress_histogram(dir, final_iter):
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(0,7.5)
    plotter.set_xlim(0,0.75)
    plotter.set_xticks([0.25*i for i in range(150)])
    plotter.set_yticks([5*i for i in range(1,100)])
    plotter.set_xlabel(r"$\sigma_{shear}$")
    plotter.set_yScaled()
    plotter.initialize_figure()

    #init
    file = "init_config.txt"
    sample = PeriodicTissue.from_config(dir, file)
    min = FIREminimization.periodic_tissue(sample)
    # min.load_cell_parameters("0.cellParameters.input")
    stresses = []
    for cellID,cell in sample.cells_.items():
        cell.max_shear_stress_ = stress.calculate_max_shear_stress(sample,cellID)
        stresses.append(cell.max_shear_stress_)
    df = pd.DataFrame({"stress": stresses})
    plotter.histogram_from_dataframe(df["stress"],color = "red", bins = 20, label = "initial")

    #final
    file = "{:04d}.bulk.txt".format(final_iter)
    sample = PeriodicTissue.from_config(dir, file)
    min = FIREminimization.periodic_tissue(sample)
    min.load_cell_parameters("{:04d}.cellParameters.input".format(final_iter))
    stresses = []
    for cellID,cell in sample.cells_.items():
        cell.max_shear_stress_ = stress.calculate_max_shear_stress(sample,cellID)
        stresses.append(cell.max_shear_stress_)
    df = pd.DataFrame({"stress": stresses})
    plotter.histogram_from_dataframe(df["stress"],color = "blue", bins = 20, label = "final")
    plotter.save_fig("{}stress_histogram.png".format(dir))

def cost_plot(dir):
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(0,1)
    plotter.set_xlim(0,np.ceil(len(np.loadtxt("{}costs.txt".format(dir)))/10)*10)
    plotter.set_xticks([50*i for i in range(1000)])
    plotter.set_yticks([.5*i for i in range(1,100)])
    plotter.set_xlabel("Epochs")
    plotter.set_ylabel("Normalized Cost")
    plotter.set_yScaled()
    plotter.initialize_figure()
    array = np.loadtxt("{}costs.txt".format(dir))
    array = list(array)
    init_cost = np.loadtxt("{}initial_cost.txt".format(dir))
    array.insert(0, init_cost)
    array/=max(array)
    iters = np.arange(len(array))
    plotter.plot_xy(iters, array, color = "blue", label = "_final")
    plotter.save_fig("{}cost.png".format(dir))

def cost_plot(dir):
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(0,1)
    plotter.set_xlim(0,np.ceil(len(np.loadtxt("{}costs.txt".format(dir)))/10)*10)
    plotter.set_xticks([50*i for i in range(1000)])
    plotter.set_yticks([.5*i for i in range(1,100)])
    plotter.set_xlabel("Epochs")
    plotter.set_ylabel("Normalized Cost")
    plotter.set_yScaled()
    plotter.initialize_figure()
    array = np.loadtxt("{}costs.txt".format(dir))
    array = list(array)
    init_cost = np.loadtxt("{}initial_cost.txt".format(dir))
    array.insert(0, init_cost)
    array/=max(array)
    iters = np.arange(len(array))
    plotter.plot_xy(iters, array, color = "blue", label = "_final")
    plotter.save_fig("{}cost.png".format(dir))



def distance_plot(dir):
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(0,1)
    plotter.set_xlim(0,np.ceil(len(np.loadtxt("{}distances.txt".format(dir)))/10)*10)
    plotter.set_xticks([50*i for i in range(1000)])
    plotter.set_yticks([.5*i for i in range(1,100)])
    plotter.set_xlabel("Epochs")
    plotter.set_ylabel("Parameter Space Distance (normalized)")
    plotter.set_yScaled()
    plotter.initialize_figure()
    array = np.loadtxt("{}distances.txt".format(dir))
    array = list(array)
    init_distance = np.loadtxt("{}initial_distances.txt".format(dir))
    array.insert(0, init_distance)
    array/=max(array)
    iters = np.arange(len(array))
    plotter.plot_xy(iters, array, color = "blue", label = "_final")
    plotter.save_fig("{}pattern_distances.png".format(dir))

def calculate_parameter_space_distance(cellParametersA, cellParametersB):
    df = pd.read_csv(cellParametersA, sep=" ",header=None)
    cellID_to_s0 = {int(i): [float(s0)] for i, s0 in zip(df[0].to_numpy(), df[2].to_numpy())}
    df = pd.read_csv(cellParametersB, sep=" ",header=None)
    for i, row in df.iterrows():
        cellID = row[0]
        s0 = row[2]
        if cellID in cellID_to_s0:
            cellID_to_s0[cellID].append(s0)
        else:
            raise ValueError("CellID {} not found in first pattern".format(cellID))
    distance = 0
    for cellID, s0s in cellID_to_s0.items():
        # print("CellID: {}, s0s: {}".format(cellID, s0s))
        distance += (s0s[0] - s0s[1])**2
    distance = distance**0.5
    # print(distance)
    return distance
def self_distance_plot(dir):
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(0,1.1)
    plotter.set_xlim(0,np.ceil(len(np.loadtxt("{}costs.txt".format(dir)))/10)*10)
    plotter.set_xticks([50*i for i in range(150)])
    plotter.set_yticks([.5*i for i in range(1,100)])
    plotter.set_xlabel("Epochs")
    plotter.set_ylabel("Parameter Space Distance (normalized)")
    plotter.set_yScaled()
    plotter.initialize_figure()
    distances = []
    iters = []
    initial_file = "{:04d}.cellParameters.input".format(0)
    for i in range(len(np.loadtxt("{}costs.txt".format(dir)))):
        final_file = "{:04d}.cellParameters.input".format(i)
        distance = calculate_parameter_space_distance("{}/{}".format(dir,initial_file), "{}/{}".format(dir,final_file))
        distances.append(distance)
        iters.append(i)
    distances = np.array(distances)
    distances /= max(distances)
    iters = np.array(iters)
    plotter.plot_xy(iters, distances, color = "blue", label = "_final")
    plotter.save_fig("{}self_distance.png".format(dir))
def main():
    # dir = "2_cells_bidisperse_0.05/"
    # dir = "two_cells_test/"
    # dir = "patterns_final/patternA/"
    # dir = "patterns_final/patternB/"
    source_dir = "patterns_final/"
    if not os.path.isdir(source_dir+"graphs/"):
        os.makedirs(source_dir+"graphs/")
    distance_plot(source_dir)
    for dir in [source_dir+"patternA/", source_dir+"patternB/"]:
        print("Processing directory:", dir)
        cost = np.loadtxt("{}costs.txt".format(dir))
        final_iter = len(cost)-1
        self_distance_plot(dir)
        print("Final iteration:", final_iter)
        s0_histogram(dir, final_iter)
        combined_stress_histogram(dir, final_iter)
        cost_plot(dir)


    # cost = np.loadtxt("{}costs.txt".format(dir))
    # # distances = np.loadtxt("{}distances.txt".format(dir))
    # # final_iter = len(distances)-1
    # final_iter = len(cost)-1
    # self_distance_plot(dir)
    # print("Final iteration:", final_iter)
    # s0_histogram(dir, final_iter)
    # combined_stress_histogram(dir, final_iter)
    # cost_plot(dir)
    # distance_plot(dir)
def single_cell_combined_s0_histogram(inc_dir, dec_dir, filename = "single_cell_combined_s0_histogram.jpg"):
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(0,14)
    plotter.set_xlim(3.5,5.5)
    plotter.set_xticks([0.5*i for i in range(150)])
    plotter.set_yticks([10*i for i in range(1,100)])
    plotter.set_xlabel(r"$s_0$")
    plotter.set_yScaled()
    plotter.initialize_figure()
    df = pd.read_csv("{}histogram_data.csv".format(inc_dir))["Final_s0"]
    plotter.histogram_from_dataframe(df,color = "red", bins = 100, alpha = 0.7, label = r"$\sigma_{target} = \sigma_\mu+2\sigma_s$")
    df = pd.read_csv("{}histogram_data.csv".format(dec_dir))["Final_s0"]
    plotter.histogram_from_dataframe(df,color = "blue", bins = 20, alpha = 0.7, label = r"$\sigma_{target} = \sigma_\mu-2\sigma_s$")
    plotter.ax.vlines(5, ymin = 0, ymax = 10,linestyles= "--", color = "black",alpha = 1,label = r"$s_0^{(initial)}$")
    plotter.ax.legend(loc='upper left')
    plotter.save_fig("graphs/{}".format(filename))


def single_cell_combined_stress_histogram(inc_dir, dec_dir, filename = "single_cell_combined_stress_histogram.jpg"):
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(0,20)
    plotter.set_xlim(0,1)
    plotter.set_xticks([0.5*i for i in range(150)])
    plotter.set_yticks([20*i for i in range(1,100)])
    plotter.set_xlabel(r"$\sigma_{hidden}$")
    plotter.set_yScaled()
    plotter.initialize_figure()
    df = pd.read_csv("{}histogram_data.csv".format(inc_dir))["Final_Stress"]
    plotter.histogram_from_dataframe(df,color = "red", bins = 20, alpha = 0.7, label = r"$\sigma_{target} = \sigma_\mu+2\sigma_s$")
    df = pd.read_csv("{}histogram_data.csv".format(dec_dir))["Final_Stress"]
    plotter.histogram_from_dataframe(df,color = "blue", bins = 20, alpha = 0.7, label = r"$\sigma_{target} = \sigma_\mu-2\sigma_s$")
    plotter.ax.legend(loc='upper right')
    plotter.save_fig("graphs/{}".format(filename))

def single_cell_initial_stress_histogram(array,filename = "initial_stress_histogram.jpg"):
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(0,8)
    plotter.set_xlim(0,0.5)
    plotter.set_xticks([0.25*i for i in range(150)])
    plotter.set_yticks([5*i for i in range(1,100)])
    plotter.set_xlabel(r"$\sigma_{hidden}$")
    plotter.set_yScaled()
    plotter.initialize_figure()
    plotter.histogram_from_array(array,color = "black", bins = 20, alpha = 0.4, label = r"$\sigma_{hidden}^{(initial)}$")
    mean = np.mean(array)
    std = np.std(array)
    plotter.ax.vlines(mean+2*std, ymin = 0, ymax = 4,linestyles= "--", color = "red",alpha = 1,label = r"$\sigma_{target} = \sigma_\mu+2\sigma_s$")
    plotter.ax.vlines(mean-2*std, ymin = 0, ymax = 4,linestyles= "--", color = "blue",alpha = 1,label = r"$\sigma_{target} = \sigma_\mu-2\sigma_s$")
    plotter.ax.legend(loc='upper right')
    plotter.save_fig("graphs/{}".format(filename))


if __name__ == "__main__":
    main()