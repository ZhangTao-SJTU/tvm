import multiprocessing
from scipy import stats
import numpy as np
import pandas as pd
import os
from toolbox import tissue
from toolbox import stress

# This function takes a spheroid, 
# and returns a dict: {cellID: list of neighbor cellIDs}

# The neighbors of a cell are defined as cells that share a polygon with it.

def find_cell_neighbors(spheroid:tissue.Sample):
    cellID_to_neighbors:dict[int,list] = {}
    for cellID,cell in spheroid.cells_.items():
        if bool(cell.type_):
            cellID_to_neighbors[cellID] = []
    for i, cellID in enumerate(cellID_to_neighbors):
        for polygonID in spheroid.cells_[cellID].polygons_:
            for j in range(i+1, len(cellID_to_neighbors)):
                if polygonID in spheroid.cells_[list(cellID_to_neighbors.keys())[j]].polygons_:
                    cellID_to_neighbors[cellID].append(list(cellID_to_neighbors.keys())[j])
                    cellID_to_neighbors[list(cellID_to_neighbors.keys())[j]].append(cellID)

    return cellID_to_neighbors

# This function takes two dicts of the form {cellID: list of neighbor cellIDs}
# These are intended to be the states of a spheroid at two different times, 
# with d2 representing the later time.

# The function returns the average overlap between the two states.

def calculate_Q(d1:dict[int,list], d2:dict[int,list]):
    # Q = 1/N * sum_cells (w_cell)
    # w_cell = 1 if cell has changed less than 2 neighbors, 0 otherwise
    # N = number of cells

    # d1, d2 are dicts of the form:
    # {cellID: list of neighbor cellIDs}

    N = len(d1)
    Q = 0
    for cellID in d1:
        set1 = set(d1[cellID])
        set2 = set(d2[cellID])
        if len(set1.union(set2)) - len(set1.intersection(set2))<2:
            Q += 1/N
    return Q

# Given a list of equally prepared runs, located in the directiories in dir_list,
# this function calculates the average overlap between the states of the spheroids
# at time t-5000 and time t.

# The function returns a DataFrame with columns:
# time: the time t
# mean: the average overlap at time t for the runs in dir_list
# sem: the standard error of the mean of the overlap at time t

def calculate_average_overlap(dir_list:list[str],timevals=[5000*i for i in range(6)]):
    Qn={i:[] for i in timevals}
    Qn[0]=1
    
    for dir in dir_list:
        #equip spheroids
        spheroids={i:None for i in timevals}
        for time in timevals:
            spheroids[time] = tissue.Sample(config_dir=dir,
                            simulation_time=time)
        print("spheroids equipped for ",dir)
        for i,time in enumerate(timevals):
            if i==0: continue
            d1 = find_cell_neighbors(spheroids[timevals[i-1]])
            d2 = find_cell_neighbors(spheroids[time])
            Qn[time].append(calculate_Q(d1,d2))

    t=[0]
    mean=[1]
    sem=[0]

    for time in timevals[1:]:
        t.append(time)
        mean.append(np.mean(Qn[time]))
        sem.append(stats.sem(Qn[time]))
    return pd.DataFrame({"time":t,"mean":mean,"sem":sem})

def main():    
    quantity="overlap"
    output_dir = "sounok/{}_data/".format(quantity)
    if not os.path.exists("sounok/"): os.mkdir("sounok/")
    if not os.path.exists(output_dir): os.mkdir(output_dir)
    s0_vals=["52","54","56","57","58"]

    num_runs = 30

    for gamma in ["025","100"]:

        def process_iteration(s0):
            dir_list = []
            for i in range(num_runs):
                test_dir = "ECM64_5/s0_{}_gamma_{}_run_{}/".format(s0,gamma,i)
                if os.path.isdir(test_dir):
                    dir_list.append(test_dir)
            return calculate_average_overlap(dir_list)

        pool = multiprocessing.Pool()

        results = pool.map(process_iteration, s0_vals)

        for s, df in zip(s0_vals, results):
            df.to_csv(output_dir+"{}_{}.csv".format(s,gamma),index=False)
    return

if __name__ == "__main__":
    main()