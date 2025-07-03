#take kv to 10[786, 408, 508, 870, 900]
from toolbox.cSection import makeSampleCrossSection
from toolbox.patterns import Patterns
from toolbox.periodic import PeriodicTissue
import matplotlib.pyplot as plt
import os
import numpy as np

def make_bidisperse_sample():
    test_sample = "7_1/"
    if os.path.isdir(test_sample):
        os.system("rm -r {}".format(test_sample))
    os.system("cp -r init/{} {}".format(test_sample,test_sample))
    dir = test_sample
    file = "minimized.txt"

    tissue = PeriodicTissue.from_config(dir,file)
    training_instance = Patterns.periodic_tissue(tissue)
    for cellID,cell in training_instance._config.cells_.items():
        if cellID%4:
            cell.s0_ = 5
        else:
            cell.s0_ = 4.9
    training_instance.write_cell_parameters()
    training_instance.minimize_config()

if __name__ == "__main__":
    make_bidisperse_sample()