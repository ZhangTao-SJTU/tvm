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

# from scripts.data_processing import download_all
# # download_all("4_cells_mean_l_4/")
# # download_all("4_cells_mean_l_5/")
# # download_all("4_cells_mean_l_6/")
# download_all("6_cells_mean_l_6/")

import os

for l in [4,5,6]:
    for a,b in [(1,1),(1,2),(2,2)]:
        local_folder = "data/{}_{}_l_{}/".format(a,b,l)
        os.makedirs(local_folder,exist_ok=True)
        for i in range(100):
            subfolder = local_folder + "{:03d}/".format(i)
            os.makedirs(subfolder,exist_ok=True)
            cluster_folder = "mameen@smatter-login.syr.edu:/home/mameen/{}_{}_l_{}/{:03d}/".format(a,b,l,i)
            for file in ["info.csv", "init_config.txt","initial_cost.txt","initial_stress_A.csv","initial_stress_B.csv","q_values.txt","costs.txt"]:
                if os.path.isfile("/Users/shabeebameen/Projects/tvm-fire/{}".format(subfolder+file)):
                    print("{} already exists".format(subfolder+file))
                    continue
                os.system("scp {} /Users/shabeebameen/Projects/tvm-fire/{}".format(cluster_folder+file,subfolder))
                print("{} downloaded".format(subfolder+file))
