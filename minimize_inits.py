from toolbox.periodic import PeriodicTissue
from toolbox.training import Training
from toolbox.stress import calculate_max_shear_stress
import numpy as np

def main():
    dir = ["init/kv_10_l_{}/".format(l) for l in [4,5,6]]
    
    for d in dir:
        for run_dir in [d+"{:03d}/".format(i) for i in range(100)]:
            sample = PeriodicTissue.from_config(run_dir, "minimized.txt")
            training_instance = Training.from_sample(sample)
            training_instance.set_cpp_executable_dir("/home/shabeeb/Projects/tvm-fire/build/")
            training_instance.edit_conf(kv = 10, final_time = 100000, log = 1000)
            training_instance.minimize_config()

def fix_conf():
    dir = ["init/kv_10_l_{}/".format(l) for l in [4,5,6]]
    for d in dir:
        for run_dir in [d+"{:03d}/".format(i) for i in range(100)]:
            sample = PeriodicTissue.from_config(run_dir, "minimized.txt")
            training_instance = Training.from_sample(sample)
            training_instance.edit_conf(kv = 10, final_time = 10000, log = 100)
def calculate_stresses():
    dir = ["init/kv_10_l_{}/".format(l) for l in [4,5,6]]
    for d in dir:
        stresses = []
        for run_dir in [d+"{:03d}/".format(i) for i in range(100)]:
            sample = PeriodicTissue.from_config(run_dir, "minimized.txt")
            for cellID in sample.cells_:
                stress = calculate_max_shear_stress(sample, cellID)
                stresses.append(stress)
        np.savetxt(d+"stresses.txt", stresses)
if __name__ == "__main__":
    calculate_stresses()