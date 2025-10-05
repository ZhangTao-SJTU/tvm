## create cross section of a periodic tissue
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

dir = "7_bidisperse_5_4.9_0.5_increase/"
# dir = "7_mono_0.5_decrease/"
df = pd.read_csv("{}stresses.csv".format(dir))
print(df["CellID"])
training_cell = df["CellID"].to_numpy()[0]
ids = [int(os.path.basename(f).split(".")[0]) for f in glob.glob("{}*.bulk.txt".format(dir))]
final_iter = max(ids)
file = "{}.bulk.txt".format(final_iter)
sample = PeriodicTissue.from_config(dir,file)
min = FIREminimization.periodic_tissue(sample)
min.load_cell_parameters("cellParameters.{}.input".format(final_iter))

for cellID, cell in sample.cells_.items():
    cell.vtk_scalar_ = cell.s0_
    cell.vtk_scalar_ = stress.calculate_max_shear_stress(sample, cellID)
normal = np.array([1,0,0])
center = sample.cells_[training_cell].center_
print(sample.cells_[training_cell].s0_)
makeSampleCrossSection(sample=min._config, center = center, normal = normal, filename = "cross_section.vtk")
single_cell = sample.extract_cell(training_cell)
makeSampleCrossSection(sample=single_cell, center=center, normal = normal, filename="single_cell_cross_section.vtk")

