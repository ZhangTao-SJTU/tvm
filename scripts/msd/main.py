#!/usr/bin/env python3
"""
plot MSD vs. time
"""

import os
import numpy as np

def main(runids, n):
    """ main function """
    tauList = [400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000]
    MSD = {}
    for runid in runids:
        MSD[runid] = run(runid, n, tauList)
    plotMSDTime(tauList, MSD)

def run(runid_start, n, tauList):
    # look for number of cells
    runDir = "{:06d}".format(runid_start)
    with open(os.path.join(runDir, "cellCenter.txt"), "r") as file:
        for line in file:
            if len(line) <= 1:
                break
            lineSplit = line.split()
            if lineSplit[0] == "time":
                N = 0
            else:
                N += 1

    MSD_n = np.zeros((n, len(tauList)), dtype=np.float32)
    for runid in range(runid_start, runid_start + n):
        centers_dict = loadRun(runid, N)
        timestamps = list(centers_dict.keys())
        print("loaded cell center coordinates of run {:06d}, time interval: {:d}~{:d}".format(runid, timestamps[0], timestamps[-1]))
        MSD_n[runid - runid_start] = computeMSD(centers_dict, tauList)
    MSD = np.mean(MSD_n, axis=0)
    return MSD

def loadRun(runid, N):
    centers_dict = {}
    runDir = "{:06d}".format(runid)
    with open(os.path.join(runDir, "cellCenter.txt"), "r") as file:
        for line in file:
            if len(line) <= 1:
                continue
            lineSplit = line.split()
            if lineSplit[0] == "time":
                timestamp = int(float(lineSplit[1]))
                centers = np.zeros((N, 3), dtype=np.float32)
                centers_dict[timestamp] = centers
            else:
                cellID = int(lineSplit[0])
                x = float(lineSplit[1])
                y = float(lineSplit[2])
                z = float(lineSplit[3])
                centers[cellID, :] = (x, y, z)
    return centers_dict

def computeMSD(centers_dict, tauList):
    Lx, Ly, Lz = (8., 8., 8.)
    timestamps = list(centers_dict.keys())
    if tauList[-1] > timestamps[-1]:
        print("tauList is over range {:d}~{:d}".format(timestamps[0], timestamps[-1]))
    N = centers_dict[timestamps[0]].shape[0]
    dt = timestamps[1] - timestamps[0]
    MSD_N = np.zeros((len(tauList), N), dtype=np.float32)
    for i in range(len(tauList)):
        tau = tauList[i]
        t0List = []
        for time in range(0, timestamps[-1] - tau + 1, dt):
            t0List.append(time)
        dr2 = np.zeros((len(t0List), N), dtype=np.float32)
        for j in range(len(t0List)):
            t0 = t0List[j]
            t1 = t0 + tau
            centers0 = centers_dict[t0]
            centers1 = centers_dict[t1]
            for k in range(N):
                r0 = centers0[k, :]
                r1 = centers1[k, :]
                dx = (r1[0]-r0[0] + Lx/2.0)%Lx - Lx/2.0
                dy = (r1[1]-r0[1] + Ly/2.0)%Ly - Ly/2.0
                dz = (r1[2]-r0[2] + Lz/2.0)%Lz - Lz/2.0
                dr2[j, k] = dx*dx + dy*dy + dz*dz
        MSD_N[i, :] = np.mean(dr2, axis=0)
    MSD = np.mean(MSD_N, axis=1)
    return MSD

def plotMSDTime(tauList, MSD):
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    plt.rcdefaults()

    colors = {}
    colors[51] = "r"
    colors[56] = "g"
    colors[61] = "b"
    labels = {}
    labels[51] = "s0=5.00"
    labels[56] = "s0=5.40"
    labels[61] = "s0=5.80"

    fig = plt.figure(figsize=(6.4, 4.8))
    gridspec.GridSpec(12, 9)
    ax = plt.subplot2grid((12, 9), (1, 1), colspan=9, rowspan=12)

    for runid in MSD:
        ax.loglog(tauList, MSD[runid], linewidth=4, color=colors[runid], ls="-", label=labels[runid])

    # ax.axis([0, 1400, 0, 0.6])
    # ax.set_xlim([0, 40000])

    # plt.xticks([0, 200, 400, 600, 800, 1000, 1200, 1400])

    ax.tick_params(labelsize=12, width=4, length=5)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(2.0)

    legend = ax.legend(loc='best', frameon=False)
    for label in legend.get_texts():
        label.set_fontsize(12)
    for label in legend.get_lines():
        label.set_linewidth(3.5)  # the legend line width

    fig.savefig('../figures/{:06d}_MSD_time.png'.format(runid))
    plt.close(fig)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--runids", action="store", type=int, nargs='+', help="the id of the runs")
    parser.add_argument("-n", "--n", action="store", type=int, default=1, help="number of runs with the same parameters")
    args = parser.parse_args()

    if args.runids is None:
        print("no runids")
        exit(1)

    main(args.runids, args.n)
