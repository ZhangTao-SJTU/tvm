from toolbox.patterns import Patterns
from toolbox.periodic import PeriodicTissue
import matplotlib.pyplot as plt
import numpy as np
import os
if os.path.isdir("7_test"):
    os.system("rm -r 7_test")
os.system("cp -r init/7_test/ 7_test/")
dir = "7_test/"
file = "minimized.txt"

tissue = PeriodicTissue.from_config(dir,file)
# tr.minimize_config()
# tr._config.load_periodic_tissue_cell_properties()

# tr.set_target_cells_spheroid(spheroid_radius=0.65)
# tr.calculate_max_shear_stresses()
# for cellID in tr._target_cells:
#     cell = tr._config.cells_[cellID]
#     print(cellID, cell.max_shear_stress_)
# tr.set_target_stress(0.079)
# tr.run()


# for cellID,cell in tr._config.cells_.items():
#     if cell.crossBoundary_:
#         continue
#     shear.append(cell.max_shear_stress_)
#     for polygonID in cell.polygons_:
#         polygon = tr._config.polygons_[polygonID]
#         polygon.vtk_scalar_ = cell.max_shear_stress_
# tr._config.write_periodic_vtk("test.vtk",use_scalar=True)
# plt.hist(shear, bins = 10)
# print(tr._target_cells)
