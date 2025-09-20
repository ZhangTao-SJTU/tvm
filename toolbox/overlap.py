from scipy import stats
import numpy as np
import pandas as pd
from toolbox.spheroid import Spheroid
import copy

# This function takes a spheroid, 
# and returns a dict: {cellID: list of neighbor cellIDs}

# The neighbors of a cell are defined as cells that share a polygon with it.

'''
# def find_cell_neighbors(spheroid:Spheroid):
#     cellID_to_neighbors:dict[int,list] = {}
#     for cellID,cell in spheroid.cells_.items():
#         if not cell.type_:
#             continue
#         cellID_to_neighbors[cellID] = []
#     for i, cellID in enumerate(cellID_to_neighbors):
#         for polygonID in spheroid.cells_[cellID].polygons_:
#             for j in range(i+1, len(cellID_to_neighbors)):
#                 test_cellID = list(cellID_to_neighbors.keys())[j]
#                 if not polygonID in spheroid.cells_[test_cellID].polygons_:
#                     continue
#                 cellID_to_neighbors[cellID].append(test_cellID)
#                 cellID_to_neighbors[test_cellID].append(cellID)
#                 # if polygonID in spheroid.cells_[list(cellID_to_neighbors.keys())[j]].polygons_:
#                 #     cellID_to_neighbors[cellID].append(list(cellID_to_neighbors.keys())[j])
#                 #     cellID_to_neighbors[list(cellID_to_neighbors.keys())[j]].append(cellID)
#     spheroid.cell_neighbors_ = cellID_to_neighbors
#     return cellID_to_neighbors

# def edit_cell_neighbors(sample:Spheroid):
#     edited_cellID_to_neighbors = copy.deepcopy(sample.cell_neighbors_)
#     for cellID, neighbors in sample.cell_neighbors_.items():
#         if sample.cells_[cellID].is_in_chain_:
#             continue
#         else: del edited_cellID_to_neighbors[cellID]
#     return edited_cellID_to_neighbors
'''
# This function takes two dicts of the form {cellID: list of neighbor cellIDs}
# These are intended to be the states of a spheroid at two different times, 
# with d2 representing the later time.

# The function returns the average overlap between the two states.

def calculate_Q(d1:dict[int,list], d2:dict[int,list]) -> float:
    # Q = 1/N * sum_cells (w_cell)
    # w_cell = 1 if cell has changed less than 2 neighbors, 0 otherwise
    # N = number of cells
    # d1, d2 are dicts of the form:
    # {cellID: list of neighbor cellIDs}
    if not len(d1) or not len(d2):
        raise ValueError("Empty dict passed to calculate_Q")
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

# def calculate_average_overlap(dir_list:list[str]):
#     timevals=[5000*i for i in range(6)]
#     Qn={i:[] for i in timevals}
#     Qn[0]=1
    
#     for dir in dir_list:
#         #equip spheroids
#         spheroids={i:None for i in timevals}
#         for time in timevals:
#             spheroids[time] = tissueSample.Sample(config_dir=dir,
#                             simulation_time=time)
#         print("spheroids equipped for ",dir)
#         for i,time in enumerate(timevals):
#             if i==0: continue
#             d1 = find_cell_neighbors(spheroids[timevals[i-1]])
#             d2 = find_cell_neighbors(spheroids[time])
#             Qn[time].append(calculate_Q(d1,d2))

#     t=[0]
#     mean=[1]
#     sem=[0]

#     for time in timevals[1:]:
#         t.append(time)
#         mean.append(np.mean(Qn[time]))
#         sem.append(stats.sem(Qn[time]))
#     return pd.DataFrame({"time":t,"mean":mean,"sem":sem})

def calculate_average_overlap(dir_to_time_to_sample:dict[str,dict[int,Spheroid]]):
    timevals = list(dir_to_time_to_sample.values())[0].keys()
    Qn={i:[] for i in timevals}
    Qn[0] = 1
    for _,time_to_sample in dir_to_time_to_sample.items():
        for i,time in enumerate(timevals):
            if time==0:
                continue
            if not len(time_to_sample[timevals[i-1]].cell_neighbors_):
                time_to_sample[timevals[i-1]].evaluate_cell_neighbors()
            if not len(time_to_sample[time].cell_neighbors_):
                time_to_sample[time].evaluate_cell_neighbors
            d1 = time_to_sample[timevals[i-1]].cell_neighbors_
            d2 = time_to_sample[time].cell_neighbors_
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
    dir_list = ["samples/5.8_1"]
    timevals = [0,10000,25000]
    # timevals = [5000*i for i in range(6)]
    print("Timevals: ",timevals)
    dir_to_time_to_sample = {}
    for dir in dir_list:
        print("Loading samples for dir: ",dir)
        time_to_sample = {time:Spheroid.from_config(dir,time) for time in timevals}
        dir_to_time_to_sample[dir] = time_to_sample
    print("Samples loaded")
    df = calculate_average_overlap(dir_to_time_to_sample)
    df.to_csv("output.csv",index=False)
    return

# if __name__ == "__main__":
#     main()