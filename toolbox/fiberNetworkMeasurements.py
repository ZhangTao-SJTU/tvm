import numpy as np
import pandas as pd
from scipy import stats
import os
from toolbox import fiberNetwork

def evalRbinsToNodeIDs(fiberNetwork,rBins,origin = None):
    if origin is None :
        origin = fiberNetwork.getOrigin()
    rBinsToNodeIDs = {i:[] for i in rBins}
    for nodeID in fiberNetwork.connectedNodes_:
        radialDistance = np.linalg.norm(
            np.subtract(
                fiberNetwork.nodeIDToCoordinates_[nodeID],
                origin))
        if (radialDistance>max(rBins) or 
            radialDistance<min(rBins)):
            continue
        bin = min(rBins, key = lambda x: abs(x - radialDistance))
        rBinsToNodeIDs[bin].append(nodeID)
    return rBinsToNodeIDs

def evalRbinsToEdgeIDs(fiberNetwork,rBins):
    rBinsToEdgeIDs = {i:[] for i in rBins}
    for edgeID,edge in fiberNetwork.edgeIDToEdge_.items():
        if (edge.radialDistance_>max(rBins) or 
            edge.radialDistance_<min(rBins)):
            continue
        bin = min(rBins, key = lambda x: abs(x - edge.radialDistance_))
        rBinsToEdgeIDs[bin].append(edgeID)
    return rBinsToEdgeIDs

class FiberNetworkMeasurements:
    def __init__(self, dirList:list, testTimeVal:int = 25000, linkerSpringNetwork:bool = False):
        self._linkerSpringNetwork:bool = linkerSpringNetwork
        self._testTimeVal:int = testTimeVal
        self._dirList:list = None
        self._timevals:list[int] = None
        self._outputDir:str = None
        self._dirToFiberNetworks:dict[str:dict[int:fiberNetwork.fiberNetwork]] = None
        self._EvalDirList(dirList)
        return
    
    def _EvalDirList(self,dirList:list) -> None:
        print("Initial list of directories: ", dirList)
        self._dirList = []
        for testDir in dirList:
            if not os.path.isdir(testDir):
                print("Directory does not exist: ", testDir)
                continue
            if self._linkerSpringNetwork:
                testFileExists = os.path.isfile(testDir + "{:07d}.link.vtk".format(self._testTimeVal))
            else:
                testFileExists = os.path.isfile(testDir + "{:07d}.ECM.vtk".format(self._testTimeVal))
            if testFileExists:
                with open (testDir + "cellCenter.txt") as file:
                    lines = file.readlines()
                    test_value = lines[-10].split()[1]
                    if (test_value == "nan" or test_value == "-nan"):
                        continue
                    else:
                        self._dirList.append(testDir)
        print("Validated list of directories: ", self._dirList)
        return
    def EvalFiberNetworks(self) -> None:
        self._dirToFiberNetworks = {dir:{} for dir in self._dirList}
        for dir in self._dirList:
            for time in self._timevals:
                self._dirToFiberNetworks[dir][time] = fiberNetwork.fiberNetwork(
                    configDir = dir,
                    time = time,
                    linkerSpringNetwork = self._linkerSpringNetwork)
                # if not self._linkerSpringNetwork:
                    # self._dirToFiberNetworks[dir][time].calculateEdgeAttributes()
        return
    def setTimevals(self,timevals:list[int]) -> None:
        self._timevals = timevals
        return
    def setOutputDir(self,outputDir:str) -> None:
        self._outputDir = outputDir
        if not os.path.isdir(outputDir):
            os.mkdir(outputDir)
        return
    
    def writeFiberEdgeAttributes(self) -> None:
        if self._linkerSpringNetwork:
            print("This function is not applicable for linker spring networks")
            print("Please use writeLinkerSpringAttributes() instead")
            return
        print("Writing (in separate files) fiber edge attributes for times: ", self._timevals)
        for time in self._timevals:
            print("Processing for time value: ", time)
            filename = self._outputDir + "FiberEdgeAttributes_{}.csv".format(time)
            print("results to be written to: ", filename)
            with open(filename,"w") as file:
                file.write("length,strain,tension\n")
                for dir in self._dirList:
                    network = self._dirToFiberNetworks[dir][time]
                    for _,edge in network.edgeIDToEdge_.items():
                        file.write("{},{},{}\n".format(edge.length_,edge.strain_,edge.tension_))
        return
    
    def calculateAverageLinkerSpringTensions(self) -> None:
        if not self._linkerSpringNetwork:
            print("This function is not applicable for fiber networks")
            print("Please use writeFiberEdgeAttributes() instead")
            return
        filename = self._outputDir + "AverageLinkerSpringTensions.csv"
        print("Writing average linker spring tensions to: ", filename)
        with open(filename,"w") as file:
            file.write("time,avgTension,semTension\n")
            for time in self._timevals:
                tensions = []
                for dir in self._dirList:
                    network = self._dirToFiberNetworks[dir][time]
                    ## IMPORTANT!~
                    ## we will consolidate these calculations later
                    ## strain is a more reasonable output in general.
                    for _,edge in network.edgeIDToEdge_.items():
                        if edge.tension_ is None:
                            length = np.linalg.norm(
                                np.subtract(
                                    network.nodeIDToCoordinates_[edge.nodes_[0]],
                                    network.nodeIDToCoordinates_[edge.nodes_[1]]))
                            l0 = (
                                network._linkerL0init
                                - network._linkerShrinkRate
                                * (network._time - network._linkerShrinkTimeInit)
                                /10)
                                #because the shrink is applied every 10 time steps
                            l0 = max(l0,network._linkerL0final)
                            tension = network._linkerK*(length-l0)
                            tensions.append(abs(tension))
                        else:
                            tensions.append(abs(edge.tension_))
                file.write("{},{},{}\n".format(time,np.mean(tensions),stats.sem(tensions)))
        return

    def calculateAverageDisplacementOnShells(
            self,
            output_filename:str = "AverageDisplacementOnShells.csv",
            initTimeval:int = 10000,
            finalTimeval:int = 25000,
            r_bins:list[float] = np.linspace(5,29,18)):
        print("Calculating average displacement on shells")
        filename = self._outputDir + output_filename
        print("Results to be written to: ", filename)
        shellRadiusToDisplacements = {i:[] for i in r_bins}
        for dir in self._dirList:
            initNetwork = self._dirToFiberNetworks[dir][initTimeval]
            finalNetwork = self._dirToFiberNetworks[dir][finalTimeval]  
            initial_rBinsToNodeIDs = evalRbinsToNodeIDs(
                fiberNetwork = initNetwork,
                rBins = r_bins)
            for rBin,nodeIDList in initial_rBinsToNodeIDs.items():
                if not len(nodeIDList):
                    continue
                for nodeID in nodeIDList:
                    if not nodeID in finalNetwork.connectedNodes_:
                        continue
                    displacement = np.subtract(
                                finalNetwork.nodeIDToCoordinates_[nodeID],
                                initNetwork.nodeIDToCoordinates_[nodeID])
                    # REMOVE TRANSLATIONAL MOTION
                    initOrigin = initNetwork.getOrigin()
                    finalOrigin = finalNetwork.getOrigin()
                    displacement = np.subtract(
                        displacement,
                        np.subtract(finalOrigin,initOrigin))
                    shellRadiusToDisplacements[rBin].append(np.linalg.norm(displacement))
            print("Processed ", dir)

        radius = []
        displacements = []
        displacements_samplesize = []
        displacements_sem = []

        for rBin, displacementList in shellRadiusToDisplacements.items():
            if len(displacementList)>1:
                radius.append(rBin)
                displacements.append(np.mean(displacementList))
                displacements_samplesize.append(len(displacementList))
                displacements_sem.append(stats.sem(displacementList))

        df = pd.DataFrame({
            "radius":radius,
            "displacements":displacements,
            "sem":displacements_sem,
            "samplesize":displacements_samplesize})
        df.to_csv(filename,index = False)
        return

    def calculateAverageRelativeDensityOnShells(
            self,
            initTimeval:int = 10000,
            finalTimeval:int = 25000,
            r_bins:list[float] = np.linspace(5,29,6)):
        print("Calculating average density change on shells")
        filename = self._outputDir + "AverageRelativeDensityOnShells.csv"
        print("Results to be written to: ", filename)
        shellRadiusToRelativeDensities = {i:[] for i in r_bins}
        for dir in self._dirList: 
            initNetwork = self._dirToFiberNetworks[dir][initTimeval]
            finalNetwork = self._dirToFiberNetworks[dir][finalTimeval]
            initOrigin = initNetwork.getOrigin() 
            initial_rBinsToNodeIDs = evalRbinsToNodeIDs(
                fiberNetwork = initNetwork,
                rBins = r_bins)
            final_rBinsToNodeIDs = evalRbinsToNodeIDs(
                fiberNetwork = finalNetwork,
                rBins = r_bins,
                origin = initOrigin)
            for rBin,nodeIDList in initial_rBinsToNodeIDs.items():
                if len(nodeIDList):
                    shellRadiusToRelativeDensities[rBin].append(
                        len(final_rBinsToNodeIDs[rBin])/len(nodeIDList))
            print("Processed ", dir)

        radius = []
        densities = []
        densitiesSem = []

        for rBin, densityList in shellRadiusToRelativeDensities.items():
            if len(densityList)>1:
                radius.append(rBin)
                densities.append(np.mean(densityList))
                densitiesSem.append(stats.sem(densityList))

        df = pd.DataFrame({
            "radius":radius,
            "densities":densities,
            "densitiesSem":densitiesSem})
        df.to_csv(filename,index = False)
        return
    
    def calculateAverageOrientationOnShells(
            self,
            r_bins:list[float] = np.linspace(5,29,8),
            time:int = 25000,
            strainCutoff:float = None):
        print("\n\nCalculating average orientation on shells...")
        if strainCutoff is None:
            print("No strain cutoff specified")       
            filename = self._outputDir + "AverageOrientationOnShells.csv"
        else:
            print("Strain cutoff specified: ", strainCutoff)
            filename = self._outputDir + "AverageOrientationOnShells_strainCutoff_{}.csv".format(strainCutoff)
        print("Results to be written to: ", filename)
        shellRadiusToOmegaList = {i:[] for i in r_bins}
        for dir in self._dirList:
            network = self._dirToFiberNetworks[dir][time] 
            rBinsToEdgeIDs = evalRbinsToEdgeIDs(
                fiberNetwork = network,
                rBins = r_bins)
            for rBin,edgeIDList in rBinsToEdgeIDs.items():
                if len(edgeIDList):
                    totalLength = 0
                    omega = np.zeros((3,3))
                    for edgeID in edgeIDList:
                        edge = network.edgeIDToEdge_[edgeID]
                        if strainCutoff is None or edge.strain_>strainCutoff:
                            omega = np.add(
                                omega,edge.length_
                                *np.outer(edge.etaSpherical_,edge.etaSpherical_))
                            totalLength += edge.length_
                    if totalLength > 0:
                        omega = np.multiply(1/totalLength,omega)
                        shellRadiusToOmegaList[rBin].append(omega)
                    else: continue     
            print("Processed ", dir)

        data = {"radius":[]}
        for i in range(3):
            for j in range(i,3):
                data["omega_{}{}".format(i,j)] = []
                data["omegaSem_{}{}".format(i,j)] = []
        
        for rBin, omegaList in shellRadiusToOmegaList.items():
            if len(omegaList)>1:
                omegaAverage = np.mean(omegaList,axis = 0)
                omegaSem = stats.sem(omegaList,axis = 0)
                data["radius"].append(rBin)
                for i in range(3):
                    for j in range(i,3):
                        data["omega_{}{}".format(i,j)].append(omegaAverage[i][j])
                        data["omegaSem_{}{}".format(i,j)].append(omegaSem[i][j])

        df = pd.DataFrame(data)
        df.to_csv(filename,index = False)
        return
    


                    