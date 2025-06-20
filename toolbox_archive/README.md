# Fiber Network Measurements:

Given a set of runs with the same parameters (i.e. the same conf file), we consider broadly two types of average measurements.

(a) The state of the fiber networks at some particular instance in time. For example:

    + Tension or stress histogram
    + Fiber network alignment
    + 

(b) Average differences between fiber networks 

1. Average Displacement on shells:
    
    __Step 0:__ The overall goal is to populate the following dictionary with all the runs: 
    
    shellRadiusToDisplacements:dict[double:list]
    
    That is, we break up the region around each spheroids into spherical shells (wrt the spheroid centers), and consider the displacement of nodes on each shell after a fixed time interval.

    (rBins = np.linspace(5,29,18) is a good shell resolution for the current systems)

    Algorithm:
    
    For each run:
    
    __Step 1:__ Evaluate the origin at time t_init=10000, and bin the node IDs of the fiber network by the shells they fall into (wrt this origin) at t_init.
    That is, use the evalRbinsToNodeIDs() to populate a dictionary initial_rBinsToNodeIDs:dict[radius:list] for this run.

    __Step 2:__ For each rBin in initial_rBinsToNodeIDs that has a non empty list, see if the nodeID exists on the connected network at t_final. If it does, append the displacement of this node (i.e. its shift in position between t_init and t_final) to shellRadiusToDisplacements[rBins]

    __Final Step:__ obtain averages sem and create a csv file...
      
1. CalculateAverageOrientationOnShells:

    __Step 0:__ Again, the overall goal is to populate the following dictionary with all the runs:

    shellRadiusToOmega