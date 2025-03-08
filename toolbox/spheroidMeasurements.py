import os
import numpy as np
from scipy import stats
from toolbox.spheroid import Spheroid
from toolbox import functions
from toolbox import stressTensor
from toolbox import cellAspectRatio
from toolbox import momentOfInertia
from toolbox import overlap

class SpheroidMeasurements:
    def __init__(self, dirList:list, testTimeVal:int = 25000):
        self._testTimeVal:int = testTimeVal
        self._dirList:list = []
        self._timevals:list[int] = []
        self._outputDir:str = None
        self._dirToSpheroids:dict[str:dict[int:Spheroid]] = None
        self._EvalDirList(dirList)
        self._spheroidsEvaluated = False
        return
    def setTimevals(self, timevals:list[int]) -> None:
        self._timevals = timevals
        return
    def appendTimevals(self,timevals:list[int]) -> None:
        for time in timevals:
            self._timevals.append(time)
        self._timevals = list(np.unique(self._timevals))
        return
    def setOutputDir(self,outputDir:str) -> None:
        self._outputDir = outputDir
        if not os.path.isdir(outputDir):
            os.mkdir(outputDir)
        return
    def _EvalDirList(self,dirList:list) -> None:
        print("Initial list of directories: ", dirList)
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
        if not len(self._dirList):
            raise ValueError("SpheroidMeasurements:_dirList is empty")
        else:
            print("Validated list of directories: ", self._dirList)
    def EvalSpheroids(self) -> None:
        self._dirToSpheroids = {dir:{} for dir in self._dirList}
        for dir in self._dirList:
            for time in self._timevals:
                self._dirToSpheroids[dir][time] = Spheroid.from_config(config_dir = dir, simulation_time = time)
        self._spheroidsEvaluated = True
        print("Spheroids evaluated.")
    def EvalCellNeighbors(self) -> None:
        for dir in self._dirList:
            for time in self._timevals:
                sample = self._dirToSpheroids[dir][time]
                if len(sample.cell_neighbors_):
                    continue
                sample.evaluate_cell_neighbors(sample)
        print("Cell neighbors evaluated.")
    def writeCellAttributes(self) -> None:
        if not len(self._timevals):
            print("Use SpheroidMeasurements.setTimevals(timevals:list[int])")
            raise ValueError("self._timearray is empty.")
        else:
            print("Writing (in separate files) cell volumes and shape indices for times: ", self._timevals)
        for time in self._timevals:
            print("Processing for time value: ", time)
            filename = self._outputDir + "CellAttributes_{}.csv".format(time)
            print("results to be written to: ", filename)
            with open(filename,"w") as file:
                file.write("{},{}\n".format("cellVolume","cellShape"))
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
        # Write stress tensor in voigt form 
        if not self._timevals:
            print("Use SpheroidMeasurements.setTimevals(timevals:list[int])")
            raise ValueError("self._timearray is empty.")
        elif not self._spheroidsEvaluated:
            print("Use SpheroidMeasurements.EvalSpheroids() to evaluate spheroids")
            raise ValueError("self._spheroidsEvaluated == False")
        else:
            print("Writing (separate files) cell.is_surface_, |sigma.rhat| for times: ",self._timevals)

        for time in self._timevals:
            filename = self._outputDir + "CellStressTensor_{}.csv".format(time)
            with open(filename,"w") as file:
                file.write("isSurface,xx,yy,zz\n")
                for dir in self._dirList:
                    sample = self._dirToSpheroids[dir][time]
                    for cellID,cell in sample.cells_.items():
                        if cell.type_:
                            cell.stress_tensor_ = stressTensor.calculate_stress_tensor(sample,cellID)
                            normal = np.subtract(cell.center_,sample.sample_center_)
                            normal = normal / np.linalg.norm(normal)
                            normalStress = np.dot(cell.stress_tensor_, normal)
                            file.write("{},{}\n".format(int(cell.is_surface_),np.linalg.norm(normalStress)))
        return
    
    def writeCellAspectRatios(self) -> None:
        if not self._timevals:
            print("Use SpheroidMeasurements.setTimevals(timevals:list[int])")
            raise ValueError("self._timearray is empty.")
        if not self._spheroidsEvaluated:
            print("Use SpheroidMeasurements.EvalSpheroids() to evaluate spheroids")
            raise ValueError("self._spheroidsEvaluated == False")
        for time in self._timevals:
            filename = self._outputDir + "AspectRatios_{}.csv".format(time)
            with open(filename,"w") as file:
                file.write("aspectRatio\n")
                for dir in self._dirList:
                    sample = self._dirToSpheroids[dir][time]
                    for cellID,cell in sample.cells_.items():
                        if cell.type_:
                            file.write("{}\n".format(
                                cellAspectRatio.calculate_aspect_ratio(sample,cellID)))
        return
    
    def writeCellAspectRatios_shape_tensor(self) -> None:
        if not self._timevals:
            print("Use SpheroidMeasurements.setTimevals(timevals:list[int])")
            raise ValueError("self._timearray is empty.")
        if not self._spheroidsEvaluated:
            print("Use SpheroidMeasurements.EvalSpheroids() to evaluate spheroids")
            raise ValueError("self._spheroidsEvaluated == False")
        for time in self._timevals:
            filename = self._outputDir + "AllShapeAspectRatios_{}.csv".format(time)
            with open(filename,"w") as file:
                file.write("aspectRatio\n")
                for dir in self._dirList:
                    sample = self._dirToSpheroids[dir][time]
                    for cellID,cell in sample.cells_.items():
                        if cell.type_:
                        # if cell.type_ and cell.is_in_chain_:
                            file.write("{}\n".format(
                                cellAspectRatio.calculate_aspect_ratio_from_shape_tensor(
                                    sample,cellID)))
        return
    
    def writeStressShapeProjections(self) -> None:
        if not self._timevals:
            print("Use SpheroidMeasurements.setTimevals(timevals:list[int])")
            raise ValueError("self._timearray is empty.")
        if not self._spheroidsEvaluated:
            print("Use SpheroidMeasurements.EvalSpheroids() to evaluate spheroids")
            raise ValueError("self._spheroidsEvaluated == False")
        for time in self._timevals:
            filename = self._outputDir + "StressShapeProjections_{}.csv".format(time)
            with open(filename,"w") as file:
                file.write("is_surface,max_shear,shapeProjection\n")
                for dir in self._dirList:
                    sample = self._dirToSpheroids[dir][time]
                    for cellID,cell in sample.cells_.items():
                        if cell.type_:
                            shape = cellAspectRatio.calculate_shape_tensor(sample,cellID)
                            stress = stressTensor.calculate_stress_tensor(sample,cellID)
                            _, shape_egvecs = np.linalg.eigh(shape)
                            stress_egvals, stress_egvecs = np.linalg.eigh(stress)
                            file.write("{},{},{}\n".format(
                                cell.is_surface_,
                                0.5 * (stress_egvals[-1] - stress_egvals[0]),
                                np.dot(shape_egvecs[-1],stress_egvecs[-1])))
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
            file.write("{},{},{}\n".format(
                "time",
                "averageSpheroidShape",
                "semSpheroidShape"))
            for time,spArray in timeToSpheroidShapes.items():
                file.write("{},{},{}\n".format(
                    time,
                    np.mean(spArray),
                    stats.sem(spArray)))
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
                file.write("{},{},{}\n".format(
                    time,
                    np.mean(dispArray),
                    stats.sem(dispArray)))            
        return


    def calculateAverageOverlap(self,is_chain_overlap = False) -> None:
        if not self._timevals:
            print("Use SpheroidMeasurements.setTimevals(timevals:list[int])")
            raise ValueError("self._timearray is empty.")
        if not self._spheroidsEvaluated:
            print("Use SpheroidMeasurements.EvalSpheroids() to evaluate spheroids")
            raise ValueError("self._spheroidsEvaluated == False")
        time_to_Qn = {time:[] for time in self._timevals}
        for i,time in enumerate(self._timevals):
            currentTime = time
            previousTime = self._timevals[i-1]
            for dir in self._dirList:
                if i == 0:
                    time_to_Qn[time].append(1)
                    continue
                if is_chain_overlap:
                    d1 = overlap.edit_cell_neighbors(
                        self._dirToSpheroids[dir][previousTime])
                else:    
                    d1 = self._dirToSpheroids[dir][previousTime].cell_neighbors_
                d2 = self._dirToSpheroids[dir][currentTime].cell_neighbors_
                time_to_Qn[currentTime].append(overlap.calculate_Q(d1,d2))
        if is_chain_overlap:
            filename = self._outputDir + "AverageChainOverlap.csv"
        else:
            filename = self._outputDir + "AverageOverlap.csv"
        with open(filename,"w") as file:
            file.write("time,mean,sem\n")
            for time,QnArray in time_to_Qn.items():
                # if len(QnArray) == 0 or len(QnArray) == 1:
                #     continue
                file.write("{},{},{}\n".format(
                    time,
                    np.mean(QnArray),
                    stats.sem(QnArray)))
        return

    def mark_chain_cells(self):
        cutoff = 0.8
        for dir in self._dirList:
            for time in self._timevals:
                sample = self._dirToSpheroids[dir][time]
                for cellID,cell in sample.cells_.items():
                    if not cell.type_:
                        continue
                    stress = stressTensor.calculate_stress_tensor(sample,cellID)
                    shape = cellAspectRatio.calculate_shape_tensor(sample,cellID)
                    _, shape_egvecs = np.linalg.eigh(shape)
                    _, stress_egvecs = np.linalg.eigh(stress)
                    if abs(np.dot(shape_egvecs[-1],stress_egvecs[-1])) > cutoff:
                        cell.is_in_chain_ = True
        return
    
    def writeMaxShearStress(self):
        time_to_max_shear_stress = {time:[] for time in self._timevals}
        for dir in self._dirList:
            for time in self._timevals:
                sample = self._dirToSpheroids[dir][time]
                for cellID,cell in sample.cells_.items():
                    if not cell.type_:
                        continue
                    if cell.is_surface_:
                        continue
                    if not cell.max_shear_stress_:
                        stress = stressTensor.calculate_stress_tensor(sample,cellID)
                        egvals = np.linalg.eigvalsh(stress)
                        cell.max_shear_stress_ = 0.5 * (egvals[-1] - egvals[0])
                    time_to_max_shear_stress[time].append(cell.max_shear_stress_)
        for time in self._timevals:
            filename = self._outputDir + "MaxShearStress_{}.csv".format(time)
            with open(filename,"w") as file:
                file.write("maxShearStress\n")
                for stress in time_to_max_shear_stress[time]:
                    file.write("{}\n".format(stress))
    
    def writeDistanceToShear(self):
        time = self._timevals[-1]
        filename = self._outputDir + "DistanceToShear.csv"
        rBins = np.linspace(0,4,17)
        distanceToStress = {r:[] for r in rBins}
        for dir in self._dirList:
            sample = self._dirToSpheroids[dir][time]
            for cellID,cell in sample.cells_.items():
                if not cell.type_:
                    continue
                if cell.is_surface_:
                    continue
                radial_distance = np.linalg.norm(np.subtract(cell.center_,sample.spheroid_center_))
                bin = min(rBins, key = lambda x: abs(x - radial_distance))
                if not cell.max_shear_stress_:
                    stress = stressTensor.calculate_stress_tensor(sample,cellID)
                    egvals = np.linalg.eigvalsh(stress)
                    cell.max_shear_stress_ = 0.5 * (egvals[-1] - egvals[0])
                distanceToStress[bin].append(cell.max_shear_stress_)
        with open(filename,"w") as file:
            file.write("radius,avgMaxShearStress,sem,samplesize\n")
            for r,stresses in distanceToStress.items():
                if len(stresses) < 2:
                    continue
                file.write("{},{},{},{}\n".format(r,np.mean(stresses),stats.sem(stresses),len(stresses)))
                
    def write_high_stress_cell_overlaps(self):
        shear_stresses = []
        for dir in self._dirList:
            for time in self._timevals:
                sample = self._dirToSpheroids[dir][time]
                for cellID,cell in sample.cells_.items():
                    if not cell.type_:
                        continue
                    if cell.is_surface_:
                        continue
                    if not cell.max_shear_stress_:
                        stress = stressTensor.calculate_stress_tensor(sample,cellID)
                        egvals = np.linalg.eigvalsh(stress)
                        cell.max_shear_stress_ = 0.5 * (egvals[-1] - egvals[0])
                    shear_stresses.append(cell.max_shear_stress_)
        cutoff = np.percentile(shear_stresses,90)

        time_to_Qn = {time:[] for time in self._timevals}
        time_to_high_stress_cell_Qn = {time:[] for time in self._timevals}
        for i,time in enumerate(self._timevals):
            currentTime = time
            previousTime = self._timevals[i-1]
            for dir in self._dirList:
                if i == 0:
                    time_to_Qn[time].append(1)
                    time_to_high_stress_cell_Qn[time].append(1)
                    continue
                if not len(self._dirToSpheroids[dir][previousTime].cell_neighbors_):
                    self._dirToSpheroids[dir][previousTime].evaluate_cell_neighbors()
                if not len(self._dirToSpheroids[dir][currentTime].cell_neighbors_):
                    self._dirToSpheroids[dir][currentTime].evaluate_cell_neighbors()
                neighbors_previous = self._dirToSpheroids[dir][previousTime].cell_neighbors_
                neighbors_current = self._dirToSpheroids[dir][currentTime].cell_neighbors_
                time_to_Qn[currentTime].append(overlap.calculate_Q(neighbors_previous,neighbors_current))
                previous_high_stress_cell_neighbors = {}
                for cellID,cell in self._dirToSpheroids[dir][previousTime].cells_.items():
                    if not cell.type_:
                        continue
                    if cell.is_surface_:
                        continue
                    if cell.max_shear_stress_ > cutoff:
                        previous_high_stress_cell_neighbors[cellID] = self._dirToSpheroids[dir][previousTime].cell_neighbors_[cellID]
                time_to_high_stress_cell_Qn[currentTime].append(
                    overlap.calculate_Q(previous_high_stress_cell_neighbors,neighbors_current))
        filename = self._outputDir + "overlap.csv"
        with open(filename,"w") as file:
            file.write("time,mean_total,sem_total\n")
            for time,QnArray in time_to_Qn.items():
                if len(QnArray) == 0 or len(QnArray) == 1:
                    continue
                file.write("{},{},{}\n".format(
                    time,
                    np.mean(QnArray),
                    stats.sem(QnArray)))
        
        filename = self._outputDir + "overlap_high_stress.csv"
        with open(filename,"w") as file:
            file.write("time,mean_total,sem_total\n")
            for time,QnArray in time_to_high_stress_cell_Qn.items():
                if len(QnArray) == 0 or len(QnArray) == 1:
                    continue
                file.write("{},{},{}\n".format(
                    time,
                    np.mean(QnArray),
                    stats.sem(QnArray)))