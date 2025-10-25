# from toolbox.periodic import PeriodicTissue
# from toolbox.minimization import FIREminimization
# from toolbox.training import Training
# from toolbox.pattern_trainer import train_random_cells,resume_run
# from toolbox.periodic import PeriodicTissue
# from toolbox import stress
# import numpy as np

# dir = "tests/6_5.0_4_cells/"
# # sample = PeriodicTissue.from_config(dir,"sample.topo")

# # min_instance = Training.periodic_tissue(sample)
# # min_instance.set_cpp_executable_dir("/home/shabeeb/Projects/tvm-fire/build/")
# # min_instance.minimize_config()
# # sample =  PeriodicTissue.from_config(dir,"minimized.txt")
# # stresses = [stress.calculate_max_shear_stress(sample,cellID) for cellID in sample.cells_]
# # train_random_cells(dir, target_stress = np.mean(stresses), n_cells = 4,cpp_executable_dir = "/home/shabeeb/Projects/tvm-fire/build/")
# resume_run(dir,max_iters = 1000)

from scripts.data_processing import download_all
# download_all("4_cells_mean_l_4/")
# download_all("4_cells_mean_l_5/")
# download_all("4_cells_mean_l_6/")
download_all("6_cells_mean_l_6/")