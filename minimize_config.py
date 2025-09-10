# from toolbox.periodic import PeriodicTissue
# from toolbox.patterns import Patterns
# import os
# if os.path.isdir("test/"):
#     os.system("rm -r test/")
# os.system("cp -r init/7_3/ test/")
# run_dir = "test/"
# file = "minimized.txt"
# cpp_executable_dir = "tvm/build/"
# tissue = PeriodicTissue.from_config(run_dir,file)
# training_instance = Patterns.periodic_tissue(tissue)
# training_instance.set_cpp_executable_dir(cpp_executable_dir)
# training_instance.minimize_config()
from toolbox.periodic import PeriodicTissue
from toolbox.minimization import FIREminimization
import sys
def main():
    dir = sys.argv[1]
    tissue = PeriodicTissue.from_config(dir,"sample.topo")
    minimizer = FIREminimization.periodic_tissue(tissue)
    minimizer.set_cpp_executable_dir("/Users/shabeebameen/Projects/tvm-fire/build/")

    # minimizer.set_cpp_executable_dir("/home/mameen/tvm/build/")
    minimizer.minimize_config()
if __name__ == "__main__":
    main()