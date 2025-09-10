from toolbox.cSection import makeSampleCrossSection
from toolbox.patterns import Patterns
from toolbox.periodic import PeriodicTissue
from toolbox.minimization import FIREminimization
from toolbox import stress
import matplotlib.pyplot as plt
from toolbox import stress
import os
import numpy as np
import pandas as pd
import glob
dir = "sp_5_cell_increase/run_13/"
df = pd.read_csv("{}initial_stress.csv".format(dir))
print(df["cellID"])
# ids = [int(os.path.basename(f).split(".")[0]) for f in glob.glob("{}*.bulk.txt".format(dir))]
# ids = sorted(ids)
# final_iter = 33
# for iter in range(final_iter+1):
file = "0000.bulk.txt"
sample = PeriodicTissue.from_config(dir,file)
sample.write_cell_collection_vtk(df["cellID"].to_numpy(),"target_cells_isolated.vtk",use_scalar=False)

# min = FIREminimization.periodic_tissue(sample)

# min.load_cell_parameters("cellParameters.{}.input".format(iter))
# for i, row in df.iterrows():
#     cellID = row["CellID"]
#     target_stress = row["Target"]
#     cell = sample.cells_[cellID]
#     # shear = stress.calculate_max_shear_stress(sample, cellID)
#     # cell.max_shear_stress_ = shear
#     # vtk_scalar = abs(cell.max_shear_stress_ - target_stress)/ target_stress
#     for polygonID in cell.polygons_:
#         polygon = sample.polygons_[polygonID]
#         # polygon.vtk_scalar_ = vtk_scalar
# for cellID in df["CellID"].to_numpy():
    # sample.write_cell_collection_vtk(df["CellID"].to_numpy(),"{}.target_cells.vtk".format(iter),use_scalar=True)
    # sample.write_cell_collection_vtk([cellID],"{}single.{}.vtk".format(dir,cellID),use_scalar=False)

