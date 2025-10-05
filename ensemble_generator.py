# create voronoi tesellation, then minimize configuration

import os
from scripts.minimize_config import minimize_config_in_dir
dir = "/home/shabeeb/Projects/tvm-fire/init_homogeneous_4/"
os.makedirs(dir, exist_ok= True)
for run in range(100):
    run_dir = dir+"{:03d}/".format(run)
    os.makedirs("{}".format(run_dir), exist_ok= True)
    if os.path.isfile("{}minimized.txt".format(run_dir)):
        continue
    print(run_dir)
    os.system("cp conf_l_4  {}conf".format(run_dir))
    minimization_found = False
    while not minimization_found:
        try:
            os.system("cd {} && python /home/shabeeb/Projects/tvm-fire/scripts/tvm/main.py conf".format(run_dir))
            minimize_config_in_dir(run_dir,"/home/shabeeb/Projects/tvm-fire/build/")
            minimization_found = True
        except:
            print("Minimization failed, retrying with a new tesellation")
            continue