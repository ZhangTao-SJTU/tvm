from toolbox.cSection import makeSampleCrossSection
from toolbox.patterns import Patterns
from toolbox.periodic import PeriodicTissue
import matplotlib.pyplot as plt
import os
import numpy as np

test_sample = "8_0_5_1_sample/"
if os.path.isdir(test_sample):
    os.system("rm -r {}".format(test_sample))
os.system("cp -r init/{} {}".format(test_sample,test_sample))
dir = test_sample
file = "minimized.txt"
tissue = PeriodicTissue.from_config(dir,file)
training_instance = Patterns.periodic_tissue(tissue)
# training_instance.edit_conf(kv = 10)
training_instance.minimize_config()
training_instance.set_target_cell_to_stress({
    11:0.12,
    157:0.28,
    494:0.20,
    166:0.13,
    378:0.20})
training_instance.set_clamping_FIRE_only(False)
training_instance.set_clamping_correction_factor(20)
training_instance.set_clamping_max_iters(2)
training_instance.set_tolerance(1e-4)
training_instance.run()
