import numpy as np
import os
from toolbox import functions
from toolbox import topology

def coordinatesAndEdgesDict(file):
    with open(file,"r") as f:
        lines = f.readlines()
        POINTSFlag = False #flag for coordinates
        LINESFlag = False #flag for edges
        STRAINSFlag = False #flag for strains
        # list of coordinates [[x1,y1,z1],[x2,y2,z2],...
        coordinates = []
        # list of edges 
        # [[e1_coord1_ID, e1_coord2_ID],[e2_coord1_ID, e2_coord2_ID], ...]
        edges = []  
        # list of strains. This is read from the .vtk file
        strains = []
        for line in lines:
            if len(line.split()):
                if line.split()[0] == "POINTS":
                    POINTSFlag = True
                    continue
                if line.split()[0] == "LINES":
                    POINTSFlag = False
                    LINESFlag = True
                    continue
                if line.split()[0] == "CELL_DATA":
                    LINESFlag = False
                    continue
                if line.split()[0] == "LOOKUP_TABLE":
                    STRAINSFlag = True
                    continue
                if POINTSFlag:
                    coordinates.append(
                        [float(line.split()[0]),
                         float(line.split()[1]),
                         float(line.split()[2])])
                if LINESFlag:
                    edges.append(
                        [int(line.split()[1]),
                         int(line.split()[2])])
                if STRAINSFlag:
                    strains.append(float(line.split()[0]))

    coordinates_dict = {i: np.array(coordinates[i]) for i in range(len(coordinates))}
    #print("Time to make coordinates dict: ",time.time()-start)
    if len(edges) == len(strains):
        edges_dict = {i: topology.fiberEdge(edges[i],strains[i]) for i in range(len(edges))}
    else:
        edges_dict = {i: topology.fiberEdge(edges[i],None) for i in range(len(edges))}
    #print("Time to make coordinates and edges dict: ",time.time()-start)
    return coordinates_dict, edges_dict

class fiberNetwork:
    def __init__(
            self,
            configDir = "build/",
            time = 25000,
            linkerSpringNetwork:bool = False):
        
        if not configDir.endswith("/"):
            raise ValueError("configDir must end with a '/'")
        
        if not os.path.isdir(configDir):
            raise ValueError("configDir must be a valid directory")
        
        self._configDir = configDir
        self._linkerSpringNetwork = linkerSpringNetwork
        if self._linkerSpringNetwork: 
            self._file = configDir + "{:07d}".format(time) + ".link.vtk"
        else: 
            self._file = configDir + "{:07d}".format(time) + ".ECM.vtk"        
        self._time = time
        self._origin = functions.calculateOrigin(self._configDir,self._time)
        
        self._kStretch = None
        self._l0 = None
        self._kBend = None
        if self._linkerSpringNetwork:
            self._nLinks = None
            self._linkerK = None
            self._linkerL0init = None
            self._linkerL0final = None
            self._linkerShrinkRate = None
            self._linkerShrinkTimeInit = None
        self._loadConfFile()
        self._edgesInitialized = False
        self.nodeIDToCoordinates_,self.edgeIDToEdge_ = coordinatesAndEdgesDict(self._file)
        self.connectedNodes_ = np.unique([edge.nodes_ for _,edge in self.edgeIDToEdge_.items()])
        return

    def getOrigin(self):
        return self._origin
    def calculateEdgeAttributes(self):
        if self._edgesInitialized: return
        for edgeID, edge in self.edgeIDToEdge_.items():
            edge.nodes_ = sorted(
                edge.nodes_,
                key = lambda id: np.linalg.norm(
                    np.subtract(
                        self.nodeIDToCoordinates_[id],
                        self._origin)))
            
            # Adjust the origin for the node that is closest to the origin.
            node0_adjustedCoordinates = np.subtract(
                self.nodeIDToCoordinates_[edge.nodes_[0]],
                self._origin)
            # Calculate the radial distance of the edge.     
            edge.radialDistance_ = np.linalg.norm(node0_adjustedCoordinates)
            # Calculate the length, strain and unit vector for the edge.
            coordinates1 = self.nodeIDToCoordinates_[edge.nodes_[1]]
            coordinates0 = self.nodeIDToCoordinates_[edge.nodes_[0]]

            edge.length_ = self._l0 * (1 + edge.strain_)
            edge.tension_ = self._kStretch * (edge.length_ - self._l0)

            edge.etaCartesian_ = np.subtract(coordinates1,coordinates0)
            edge.etaCartesian_ = np.multiply(edge.etaCartesian_,1/np.linalg.norm(edge.etaCartesian_))

            cosTheta = (node0_adjustedCoordinates[2]
                        / edge.radialDistance_)
            sinTheta = (np.sqrt(node0_adjustedCoordinates[0]**2
                                + node0_adjustedCoordinates[1]**2)
                                / edge.radialDistance_)
            cosPhi = (node0_adjustedCoordinates[0]
                      / np.sqrt(node0_adjustedCoordinates[0]**2
                                + node0_adjustedCoordinates[1]**2))
            sinPhi = (node0_adjustedCoordinates[1]
                      / np.sqrt(node0_adjustedCoordinates[0]**2
                                + node0_adjustedCoordinates[1]**2))
                 
            edge.rHat_ = np.array(
                [sinTheta * cosPhi,
                 sinTheta * sinPhi,
                 cosTheta])
            edge.thetaHat_ = np.array(
                [cosTheta * cosPhi,
                 cosTheta * sinPhi,
                 - sinTheta])
            edge.phiHat_ = np.array(
                [- sinPhi,
                 cosPhi,
                 0])
            # Hence, calculate the unit vector for the edge in spherical coordinates. 
            edge.etaSpherical_ = np.array(
                [np.dot(edge.etaCartesian_, edge.rHat_),
                 np.dot(edge.etaCartesian_, edge.thetaHat_),
                 np.dot(edge.etaCartesian_,edge.phiHat_)])
            
        self._edgesInitialized = True
        return
    
    def _loadConfFile(self):
        with open(self._configDir + "conf","r") as file:
            lines = file.readlines()
            if self._linkerSpringNetwork:
                for line in lines:
                    if (line.split()[0]) == "link": 
                        self._nLinks = (int(line.split()[1]))
                        self._linkerK = (float(line.split()[2]))
                        self._linkerL0init = (float(line.split()[3]))
                        self._linkerL0final = (float(line.split()[4]))
                        self._linkerShrinkRate = (float(line.split()[5]))
                        self._linkerShrinkTimeInit = (float(line.split()[6]))
            else:
                for line in lines:                
                    if (line.split()[0]) == "fiber": 
                        self._kStretch = (float(line.split()[1]))
                        self._l0 = (float(line.split()[2]))
                        self._kBend = (float(line.split()[3]))
        return

