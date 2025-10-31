#plotting class
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

plt.style.use("toolbox/manuscript.mplstyle")
# plt.rcParams['text.usetex'] = True
import numpy as np
import pandas as pd

class plot:
    def __init__(self):
        self.title = None
        self.neutralColor = "#7d878a"
        self.leftMargin = 0.2
        self.bottomMargin = 0.13
        self.width = 0.9-self.leftMargin
        self.height = 0.87-self.bottomMargin
        self.xlabel = 'X-axis'
        self.ylabel = 'Y-axis'
        self.x0 = 0
        self.x1 = 1
        self.y0 = 0
        self.y1 = 1
        self.xticks = [0.5 * i for i in range(1,5)]
        self.yticks = [1+0.2*i for i in range(4)]
        self.xScaled= False
        self.yScaled = False
        self.fig = None
        self.ax = None
    def set_xLog(self):
        self.ax.xscale("log")
    def set_yLog(self):
        self.ax.yscale("log")
    def set_xScaled(self,yeah = True):
        self.xScaled = yeah
    def set_yScaled(self,yeah = True):
        self.yScaled = yeah
    def set_xlabel(self,xlabel):
        self.xlabel = xlabel
    def set_ylabel(self,ylabel):
        self.ylabel = ylabel
    def set_xlim(self,x0,x1):
        self.x0 = x0
        self.x1 = x1
    def set_ylim(self,y0,y1):
        self.y0 = y0
        self.y1 = y1
    def set_xticks(self,xticks):
        self.xticks = xticks
    def set_yticks(self,yticks):
        self.yticks = yticks
    def set_title(self,title):
        self.title = title
    def initialize_figure(self):
        self.fig, self.ax = plt.subplots()
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((0,0))
        self.ax.set_position([self.leftMargin, self.bottomMargin, self.width, self.height]) 
        if self.title is not None:
            self.ax.set_title(self.title)
        self.ax.set_xticks(self.xticks)
        self.ax.set_yticks(self.yticks)
        if self.xScaled:
            self.ax.xaxis.set_major_formatter(formatter)
        if self.yScaled:
            self.ax.yaxis.set_major_formatter(formatter)
        self.ax.set_xlim(self.x0, self.x1)
        self.ax.set_ylim(self.y0, self.y1)
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)

    def plot_xy(self,x_array,y_array, color = "#7d878a", label = "plot",alpha = 0.8):
        self.ax.plot(x_array,y_array, color = color, label = label,alpha = alpha)
    def plot_errorfill(self,x_array,y_array, err_array,color = "#7d878a", label = "plot",alpha = 0.3):
        self.ax.fill_between(x_array, y_array-err_array, y_array+err_array,color = color, label = label,alpha = alpha)
    def plot_max_min_fill(self,x_array,y_min_array, y_max_array,color = "#7d878a", label = "_plot",alpha = 0.3):
        self.ax.fill_between(x_array, y_min_array, y_max_array,color = color, label = label,alpha = alpha)
    def histogram_from_dataframe(self, data, from_array = False,  bins = 50, color = "#7d878a", alpha = 0.8,label = "plot"):
        if from_array:
            data = pd.DataFrame(data, columns = ["data"])
        data.plot(
            ax = self.ax,
            kind = "hist",
            density = True,
            linewidth  = 2,
            edgecolor = color,
            bins = bins,
            color = color,
            label = "_hidden",
            xlabel=self.xlabel,
            alpha = alpha)
        data.plot(ax = self.ax, kind = "kde",color = color, label = label, alpha = alpha, linewidth = 5)

    def histogram_from_array(self, array, bins = 50, color = "#7d878a", label = "_plot",alpha = 0.8):
        # df = pd.DataFrame(array, columns = ["data"])
        df = pd.DataFrame({"data":array})
        self.histogram_from_dataframe(df["data"], bins = bins, color = color, label = label,alpha = alpha)

    def save_fig(self,filename = "test.png"):
        # self.ax.legend()
        self.fig.savefig(fname=filename, dpi=300, transparent=False)