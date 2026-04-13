#plotting class
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
plt.style.use("toolbox/manuscript.mplstyle")
# plt.rcParams['text.usetex'] = True
from scipy import stats
import numpy as np
import pandas as pd

class plot:
    def __init__(self):
        self.title = None
        self.neutralColor = "#7d878a"
        # self.leftMargin = 0.2
        # self.bottomMargin = 0.13
        # self.width = 0.9-self.leftMargin
        # self.height = 0.87-self.bottomMargin
        self.xlabel = None
        self.ylabel = None
        self.xlim = None
        self.ylim = None
        self.xticks = None
        self.yticks = None
        self.xScaled= False
        self.yScaled = False
        self.xLog = False
        self.yLog = False
        self.fig = None
        self.ax = None
        self.transparent = True
    def set_xLog(self, isLog:bool = True):
        self.xLog = isLog
    def set_yLog(self, isLog:bool = True):
        self.yLog = isLog
    def set_xScaled(self,isScaled:bool = True):
        self.xScaled = isScaled
    def set_yScaled(self,isScaled:bool = True):
        self.yScaled = isScaled
    def set_xlabel(self,xlabel):
        self.xlabel = xlabel
    def set_ylabel(self,ylabel):
        self.ylabel = ylabel

    def set_xlim(self,xlim):
        self.xlim = xlim
    def set_ylim(self,ylim):
        self.ylim = ylim
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
        # self.ax.set_position([self.leftMargin, self.bottomMargin, self.width, self.height]) 
        if self.title is not None:
            self.ax.set_title(self.title)
        if self.xticks is not None:
            self.ax.set_xticks(self.xticks)
        if self.yticks is not None:
            self.ax.set_yticks(self.yticks)
        if self.xScaled:
            self.ax.xaxis.set_major_formatter(formatter)
        if self.yScaled:
            self.ax.yaxis.set_major_formatter(formatter)
        self.ax.set_xlim(self.xlim[0],self.xlim[1])
        self.ax.set_ylim(self.ylim[0],self.ylim[1])
        if self.xLog:
            self.ax.set_xscale("log")
        if self.yLog:
            self.ax.set_yscale("log")
        self.fig.supxlabel(self.xlabel)
        self.fig.supylabel(self.ylabel)

    def plot_xy(self,x_array,y_array, **kwargs):
        self.ax.plot(x_array,y_array, **kwargs)

        
    def plot_scatter(self, x_array,y_array,**kwargs):
        self.ax.scatter(x_array,y_array,**kwargs)
    def plot_errorbar(self,x_array,y_array, err_array, **kwargs):
        self.ax.errorbar(x_array, y_array, err_array,**kwargs)

    def plot_errorfill(self,x_array,y_array, err_array,color = "#7d878a", label = "plot",alpha = 0.3):
        self.ax.fill_between(x_array, y_array-err_array, y_array+err_array,color = color, label = label,alpha = alpha,)
    def plot_max_min_fill(self,x_array,y_min_array, y_max_array,color = "#7d878a", label = "_plot",alpha = 0.3):
        self.ax.fill_between(x_array, y_min_array, y_max_array,color = color, label = label,alpha = alpha)
    def histogram_from_dataframe(self, data, from_array = False, fit_type = None, **kwargs):
        if from_array:
            data = pd.DataFrame(data, columns = ["data"])
        data.plot(
            ax = self.ax,
            kind = "hist",
            density = True,
            # xlabel=self.xlabel,
            **kwargs)
        if fit_type is None:
            return
        if fit_type == "kde":
            data.plot(ax = self.ax, kind = "kde", label = "_hidden", color = kwargs["color"])
        elif fit_type == "gamma":
            a, loc, scale = stats.gamma.fit(data, floc=0)
            x = np.linspace(self.xlim[0], self.xlim[1], 500)
            pdf = stats.gamma.pdf(x, a, loc=loc, scale=scale)
            # gamma_label = r'$\alpha={:.2f}, \theta={:.2e}$'.format(a, scale)
            gamma_label = "_None"
            self.ax.plot(x, pdf, **kwargs, label = gamma_label, linewidth = 5)
    def histogram_from_array(self, array, fit_type = None, label = "_plot",**kwargs):
        # df = pd.DataFrame(array, columns = ["data"])
        df = pd.DataFrame({"data":array})
        self.histogram_from_dataframe(df["data"], fit_type = fit_type,  label = label, **kwargs)

    def save_fig(self,filename = "test.png"):
        # self.ax.legend()
        self.fig.savefig(fname=filename, transparent=self.transparent)
        plt.close(self.fig)

class plot_shared_x_axis(plot):
    def __init__(self):
        super().__init__()
        self.ax_top = None
        self.ax_bottom = None
        self.ylabel_top = None
        self.ylabel_bottom = None
        self.yticks_top = None
        self.yticks_bottom = None
        self.ylim_top = None
        self.ylim_bottom = None
        self.yLog_top = False
        self.yLog_bottom = False
        ## right y axis
        self.ax_bottom_right = None
        self.ylabel_bottom_right = None
        self.yticks_bottom_right = None
        self.ylim_bottom_right = None
        self.yLog_bottom_right = False
    def set_ylabel_top(self,ylabel):
        self.ylabel_top = ylabel
    def set_ylabel_bottom(self,ylabel):
        self.ylabel_bottom = ylabel
    def set_ylim_top(self,ylim):
        self.ylim_top = ylim
    def set_yticks_top(self,ticks):
        self.yticks_top = ticks
    def set_yticks_bottom(self,ticks):
        self.yticks_bottom = ticks
    def set_ylim_bottom(self,ylim):
        self.ylim_bottom = ylim
    def set_yLog_top(self, isLog:bool = True):
        self.yLog_top = isLog
    def set_yLog_bottom(self, isLog:bool = True):
        self.yLog_bottom = isLog
    def set_ylabel_bottom_right(self, ylabel):
        self.ylabel_bottom_right = ylabel

    def set_ylim_bottom_right(self, ylim):
        self.ylim_bottom_right = ylim

    def set_yticks_bottom_right(self, ticks):
        self.yticks_bottom_right = ticks

    def set_yLog_bottom_right(self, isLog: bool = True):
        self.yLog_bottom_right = isLog

    def initialize_sharedx_figure(self, height_ratios=[1,1]):
        self.fig, (self.ax_top, self.ax_bottom) = plt.subplots(
            2, 1,
            sharex=True,
            height_ratios=height_ratios
        )
        # self.fig.subplots_adjust(
        #     left=self.leftMargin,
        #     bottom=self.bottomMargin,
        #     right=self.leftMargin + self.width,
        #     top=self.bottomMargin + self.height,
        #     hspace=0
        # )
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((0,0))
        # Remove touching spines for clean look
        self.ax_top.spines["bottom"].set_visible(True)
        self.ax_bottom.spines["top"].set_visible(False)
        self.ax_top.tick_params(labelbottom=False)
        self.ax_top.set_xlim(self.xlim[0],self.xlim[1])
        self.ax_bottom.set_xlim(self.xlim[0],self.xlim[1])

        self.ax_top.set_ylim(self.ylim_top[0],self.ylim_top[1])
        self.ax_bottom.set_ylim(self.ylim_bottom[0],self.ylim_bottom[1])
        self.fig.supxlabel(self.xlabel)
        self.ax_top.set_ylabel(self.ylabel_top)
        self.ax_bottom.set_ylabel(self.ylabel_bottom)
        if self.xticks is not None:
            self.ax_top.set_xticks(self.xticks)
        if self.yticks_top is not None:
            self.ax_top.set_yticks(self.yticks_top)
        if self.yticks_bottom is not None:
            self.ax_bottom.set_yticks(self.yticks_bottom)
        if self.xLog:
            self.ax_top.set_xscale("log")
            self.ax_bottom.set_xscale("log")
        if self.yLog_top:
            self.ax_top.set_yscale("log")
        if self.yLog_bottom:
            self.ax_bottom.set_yscale("log")
        if self.title is not None:
            self.ax_top.set_title(self.title)
        self.ax_bottom_right = self.ax_bottom.twinx()
        
        if self.ylim_bottom_right is not None:
            self.ax_bottom_right.set_ylim(self.ylim_bottom_right[0], self.ylim_bottom_right[1])

        if self.yticks_bottom_right is not None:
            self.ax_bottom_right.set_yticks(self.yticks_bottom_right)

        if self.ylabel_bottom_right is not None:
            self.ax_bottom_right.set_ylabel(self.ylabel_bottom_right)

        if self.yLog_bottom_right:
            self.ax_bottom_right.set_yscale("log")
# after setting ylabels
        offset = -0.18
        self.ax_top.yaxis.set_label_coords(offset, 0.5)
        self.ax_bottom.yaxis.set_label_coords(offset, 0.5)

    def plot_xy_top(self,x_array,y_array, **kwargs):
        self.ax_top.plot(x_array,y_array, **kwargs)
    def plot_xy_bottom(self,x_array,y_array, **kwargs):
        self.ax_bottom.plot(x_array,y_array, **kwargs)
    def plot_xy_bottom_right(self, x_array, y_array, **kwargs):
        self.ax_bottom_right.plot(x_array, y_array, **kwargs)
    def plot_scatter_top(self,x_array,y_array, **kwargs):
        self.ax_top.scatter(x_array,y_array, **kwargs)
    def plot_scatter_bottom(self,x_array,y_array, **kwargs):
        self.ax_bottom.scatter(x_array,y_array, **kwargs)
    def plot_scatter_bottom_right(self, x_array, y_array, **kwargs):
        self.ax_bottom_right.scatter(x_array, y_array, **kwargs)
    def legend_top(self, **kwargs):
        handles_top, labels_top = self.ax_top.get_legend_handles_labels()
        handles_bottom, labels_bottom = self.ax_bottom.get_legend_handles_labels()
        handles_bottom_right, labels_bottom_right = self.ax_bottom_right.get_legend_handles_labels()

        handles = handles_top + handles_bottom + handles_bottom_right
        labels = labels_top + labels_bottom + labels_bottom_right

        unique = dict(zip(labels, handles))

        self.ax_top.legend(unique.values(), unique.keys(), **kwargs)