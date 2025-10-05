from toolbox.periodic import PeriodicTissue
from toolbox.minimization import FIREminimization
import sys
def minimize_config_in_dir(dir, input_file = "sample.topo", cpp_executable_dir = "/home/shabeeb/Projects/tvm-fire/build/"):
    tissue = PeriodicTissue.from_config(dir,input_file)
    minimizer = FIREminimization.periodic_tissue(tissue)
    minimizer.set_cpp_executable_dir(cpp_executable_dir)
    minimizer.minimize_config()

def main():
    dir = sys.argv[1]
    minimize_config_in_dir(dir)

if __name__ == "__main__":
    main()