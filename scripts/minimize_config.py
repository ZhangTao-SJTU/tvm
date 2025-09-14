from toolbox.periodic import PeriodicTissue
from toolbox.minimization import FIREminimization
import sys
import os

def main():
    dir = sys.argv[1]
    # dir = os.getcwd()+"/"
    tissue = PeriodicTissue.from_config(dir,"sample.topo")
    minimizer = FIREminimization.periodic_tissue(tissue)
    minimizer.set_cpp_executable_dir("/home/mameen/tvm/build/")
    minimizer.minimize_config()

if __name__ == "__main__":
    main()