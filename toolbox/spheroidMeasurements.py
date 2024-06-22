import os
from toolbox import tissueSample
from toolbox import functions
from toolbox import stressTensor
import numpy as np
from scipy import stats
class SpheroidMeasurements:
    def __init__(self, dirList:list, testTimeVal:int = 25000):
        self._testTimeVal:int = testTimeVal
        self._dirList:list = None
        self._timevals:list[int] = None
        self._outputDir:str = None
        self._dirToSpheroids:dict[str:dict[int:tissueSample.Sample]] = None
        self._EvalDirList(dirList)
        return
    
    def _EvalDirList(self,dirList:list) -> None:
        print("Initial list of directories: ", dirList)
        self._dirList = []
        for testDir in dirList:
            if not os.path.isdir(testDir):
                continue
            if os.path.isfile(testDir + "{:07d}.sample.vtk".format(self._testTimeVal)):
                with open (testDir + "cellCenter.txt") as file:
                    lines = file.readlines()
                    test_value = lines[-2].split()[1]
                    if (test_value == "nan" or test_value == "-nan"):
                        continue
                    else: self._dirList.append(testDir)
        print("Validated list of directories: ", self._dirList)
        return
    
    def EvalSpheroids(self) -> None:
        self._dirToSpheroids = {dir:{} for dir in self._dirList}
        for dir in self._dirList:
            for time in self._timevals:
                self._dirToSpheroids[dir][time] = tissueSample.Sample(
                    configDir = dir,
                    simulationTime = time)
        return
    
    def EvalStressTensors(self) -> None:
        for dir in self._dirList:
            for time in self._timevals:
                sample = self._dirToSpheroids[dir][time]
                for cellID,cell in sample.cells_.items():
                    if cell.type_:
                        stressTensor.calculate_stress_tensor(sample,cellID)
        return
    
    def setTimevals(self,timevals:list[int]) -> None:
        self._timevals = timevals
        return
    
    def setOutputDir(self,outputDir:str) -> None:
        self._outputDir = outputDir
        if not os.path.isdir(outputDir):
            os.mkdir(outputDir)
        return
    
    def writeCellAttributes(self,timeArray) -> None:
        print("Writing (in separate files) cell attributes for times: ", timeArray)
        for time in timeArray:
            print("Processing for time value: ", time)
            filename = self._outputDir + "CellAttributes_{}.csv".format(time)
            print("results to be written to: ", filename)
            with open(filename,"w") as file:
                file.write("cellVolume,cellShape\n")
                for dir in self._dirList:
                    if not os.path.isfile(dir + "{:07d}.cellInfo.txt".format(time)):
                        functions.writeTimeCellInfo(dir,time)
                    with open(dir + "{:07d}.cellInfo.txt".format(time),"r") as cellInfoFile:
                        lines = cellInfoFile.readlines()
                        for i, line in enumerate(lines):
                            if i == 0 or len(line.split()) == 0:
                                continue
                            cellVolume = float(line.split()[-2])
                            cellShape = float(line.split()[-1])
                            file.write("{},{}\n".format(cellVolume,cellShape))
        return
    
    def writeCellStressTensor(self) -> None:
        for time in self._timevals:
            filename = self._outputDir + "CellStressTensor_{}.csv".format(time)
            with open(filename,"w") as file:
                file.write("isSurface,NormalStress\n")
                for dir in self._dirList:
                    sample = self._dirToSpheroids[dir][time]
                    for cellID,cell in sample.cells_.items():
                        if cell.type_:
                            normal = np.subtract(cell.center_,sample.sample_center_)
                            normal = normal / np.linalg.norm(normal)
                            normalStress = np.dot(cell.stress_tensor_, normal)
                            file.write("{},{}\n".format(int(cell.is_surface_),np.linalg.norm(normalStress)))
        return
    
    def calculateAverageSpheroidShapes(self,timeArray) -> None:
        timeToSpheroidShapes = {time:[] for time in timeArray}
        for dir in self._dirList:
            print("Processing for directory: ", dir)
            with open(dir + "spheroidShape.txt") as file:
                lines = file.readlines()
                for line in lines:
                    time = int(float(line.split()[0]))
                    if time in timeArray:
                        timeToSpheroidShapes[time].append(float(line.split()[-1]))

        filename = self._outputDir + "AverageSpheroidShapes.csv"
        print("results to be written to: ", filename)
        with open(filename,"w") as file:
            file.write("time,averageSpheroidShape,semSpheroidShape\n")
            for time,spArray in timeToSpheroidShapes.items():
                file.write("{},{},{}\n".format(time,np.mean(spArray),stats.sem(spArray)))
        return

    def calculateAverageCellDisplacements(self,timeArray) -> None:
        timeToDisplacements = {time:[] for time in timeArray}
        filename = self._outputDir + "AverageCellDisplacements.csv"
        print("results to be written to: ", filename)
        for i,time in enumerate(timeArray):
            if i == 0:
                continue
            currentTime = timeArray[i]
            previousTime = timeArray[i-1]
            for dir in self._dirList:
                cellIDToCurrentCenter = {}
                if not os.path.isfile(dir + "{:07d}.cellInfo.txt".format(currentTime)):
                    functions.writeTimeCellInfo(dir,currentTime)
                with open(dir + "{:07d}.cellInfo.txt".format(currentTime)) as file:
                    lines = file.readlines()
                    for linenum, line in enumerate(lines):
                        if linenum == 0 or len(line.split()) == 0:
                            continue
                        cellID = int(line.split()[0])
                        cellCenter = np.array([float(line.split()[1]),
                                               float(line.split()[2]),
                                               float(line.split()[3])])
                        cellIDToCurrentCenter[cellID] = cellCenter
                if not os.path.isfile(dir + "{:07d}.cellInfo.txt".format(previousTime)):
                    functions.writeTimeCellInfo(dir,previousTime)
                with open(dir + "{:07d}.cellInfo.txt".format(previousTime)) as file:
                    lines = file.readlines()
                    for linenum, line in enumerate(lines):
                        if linenum == 0 or len(line.split()) == 0:
                            continue
                        cellID = int(line.split()[0])
                        cellCenter = np.array([float(line.split()[1]),
                                               float(line.split()[2]),
                                               float(line.split()[3])])
                        displacement = np.linalg.norm(
                            np.subtract(cellIDToCurrentCenter[cellID],cellCenter))
                        timeToDisplacements[currentTime].append(displacement)

        with open(filename,"w") as file:
            file.write("time,averageDisplacement,semDisplacement\n")
            for time,dispArray in timeToDisplacements.items():
                if len(dispArray) == 0 or len(dispArray) == 1:
                    continue
                file.write("{},{},{}\n".format(time,np.mean(dispArray),stats.sem(dispArray)))            
        return