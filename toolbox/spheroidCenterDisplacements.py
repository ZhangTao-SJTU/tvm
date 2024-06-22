from scipy import stats
import numpy as np
import pandas as pd
import os
from toolbox import functions

def main():
    if not os.path.isdir("sounok/spheroidCenterDisplacements/"):
        os.mkdir("sounok/spheroidCenterDisplacements/")
    timevals=[500*i for i in range(51)]
    for s0 in ["52","54","56","57","58"]:
        for gamma in ["025","100"]:
            spheroidCenterDisplacements = {i:[] for i in timevals}
            spheroidCenterDisplacements[0] = 0
            for run in range(30):
                
                dir = "ECM64_5/s0_{}_gamma_{}_run_{}/".format(s0,gamma,run)
                if not os.path.isdir(dir): continue
                thisSpheroidCentersDict={i:[] for i in timevals}
                for time in timevals:
                    if not os.path.isfile(dir+"{:07d}.cellInfo.txt".format(time)):
                        functions.make_time_cellInfo(dir,time)
                    thisSpheroidCenterNow = []
                    with open(dir+"{:07d}.cellInfo.txt".format(time),"r") as f:
                        lines = f.readlines()
                        for i, line in enumerate(lines):
                            if i==0: continue
                            thisSpheroidCenterNow.append([float(j) for j in lines[i].split()[1:4]])
                    thisSpheroidCenterNow = np.mean(thisSpheroidCenterNow,axis=0)
                    thisSpheroidCentersDict[time]=thisSpheroidCenterNow
                for i,time in enumerate(timevals):
                    if i==0: continue
                    spheroidCenterDisplacements[time].append(np.linalg.norm(
                        np.subtract(thisSpheroidCentersDict[time],
                                    thisSpheroidCentersDict[0])))
            t=[0]
            mean=[0]
            sem=[0]
            for time in timevals[1:]:
                t.append(time)
                mean.append(np.mean(spheroidCenterDisplacements[time]))
                sem.append(stats.sem(spheroidCenterDisplacements[time]))
            df = pd.DataFrame({"time":t,"mean":mean,"sem":sem})
            df.to_csv("sounok/spheroidCenterDisplacements/{}_{}.csv".format(s0,gamma),index=False)

    return

if __name__ == '__main__':
    main()