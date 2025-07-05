from toolbox import manuscriptPlots
from toolbox.periodic import PeriodicTissue
from toolbox.minimization import FIREminimization
from toolbox import stress
import pandas as pd
import numpy as np

def s0_histogram(dir, final_iter):
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(0,15)
    plotter.set_xlim(4,6)
    plotter.set_xticks([0.5*i for i in range(150)])
    plotter.set_yticks([5*i for i in range(1,100)])
    plotter.set_xlabel(r"$s_0$")
    plotter.set_yScaled()
    plotter.initialize_figure()
    df = pd.read_csv("{}cellParameters.{}.input".format(dir,final_iter), sep=" ",header=None)
    plotter.histogram_from_dataframe(df[2],color = "blue", bins = 20, label = "_final")
    plotter.save_fig("{}/s0_histogram.png".format(dir))

def combined_stress_histogram(dir, final_iter):
    plotter = manuscriptPlots.plot()
    plotter.set_ylim(0,10)
    plotter.set_xlim(0,1)
    plotter.set_xticks([0.5*i for i in range(150)])
    plotter.set_yticks([5*i for i in range(1,100)])
    plotter.set_xlabel(r"$\sigma_{shear}$")
    plotter.set_yScaled()
    plotter.initialize_figure()

    #init
    file = "init_config.txt"
    sample = PeriodicTissue.from_config(dir, file)
    min = FIREminimization.periodic_tissue(sample)
    min.load_cell_parameters("cellParameters.init.input")
    stresses = []
    for cellID,cell in sample.cells_.items():
        cell.max_shear_stress_ = stress.calculate_max_shear_stress(sample,cellID)
        stresses.append(cell.max_shear_stress_)
    df = pd.DataFrame({"stress": stresses})
    plotter.histogram_from_dataframe(df["stress"],color = "red", bins = 20, label = "initial")

    #final
    file = "{}.bulk.txt".format(final_iter)
    sample = PeriodicTissue.from_config(dir, file)
    min = FIREminimization.periodic_tissue(sample)
    min.load_cell_parameters("cellParameters.{}.input".format(final_iter))
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
    plotter.set_xlim(0,40)
    plotter.set_xticks([5*i for i in range(150)])
    plotter.set_yticks([.5*i for i in range(1,100)])
    plotter.set_xlabel("Epochs")
    plotter.set_ylabel("Normalized Cost")
    plotter.set_yScaled()
    plotter.initialize_figure()
    array = np.loadtxt("{}costs.txt".format(dir))
    array/=array[0]
    iters = np.arange(len(array))
    plotter.plot_xy(iters, array, color = "blue", label = "_final")
    plotter.save_fig("{}/cost.png".format(dir))

def main():
    dir = "7_bidisperse_5_4.9/"
    final_iter = 9
    s0_histogram(dir, final_iter)
    combined_stress_histogram(dir, final_iter)
    cost_plot(dir)


if __name__ == "__main__":
    main()