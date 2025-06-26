# Doesnt work.
from toolbox import pulseDrive
import os
if os.path.isdir("7_test"):
    os.system("rm -r 7_test")
os.system("cp -r init/7_test/ 7_test/")
dir = "7_test/"
file = "minimized.txt"
tr = pulseDrive.from_config(dir,file)
tr.edit_conf(s0 = 5.2, kv = 100)
for _ in range(5):
    tr.pick_random_modified_cell()
tr.run_s0_pulsing(min_s0 = 5.2,max_s0 = 5.6, iterations = 10)

# tr = pulseDrive.from_config("test/","8.bulk.txt")
# tr.set_iter_counter(9)
# tr.single_iteration()