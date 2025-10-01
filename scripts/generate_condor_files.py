import os
import sys
import numpy as np

def create_sub_file(script,dir):
    lines = []
    lines.append("executable = /home/mameen/examples/singularity_wrapper.sh\n")
    lines.append("arguments  = /home/mameen/examples/ubuntu18_povray_paula.img {}run_job.sh\n".format(dir))
    lines.append("transfer_input_files = /home/mameen/scripts/{}, /home/mameen/toolbox, {}\n".format(script,dir))
    lines.append("should_transfer_files = YES\n")
    lines.append("output     = output.txt\n")
    lines.append("error      = error.txt\n")
    lines.append("log        = log.txt\n")
    lines.append("getenv     = True\n")
    lines.append("request_cpus = 1\n")
    lines.append("request_memory = 200 MB\n")
    # lines.append('Requirements = TARGET.vm_name == "its-u20-nfs-20210413" && regexp("CRUSH", TARGET.name)\n')
    lines.append("queue \n")
    # with open("job_{}.sub".format(run_num), "w") as f:
    with open("{}run_job.sub".format(dir), "w") as f:
        for line in lines:
            f.write(line)

def create_exec_file(script,dir):
    lines = []
    lines.append("#!/bin/bash\n")
    lines.append("source /home/mameen/.bashrc\n")
    # lines.append("pip install --user -e .\n")
    lines.append("export PYTHONPATH=$PWD:$PYTHONPATH\n")
    lines.append("python {} {}\n".format(script,dir))
    # with open("run_job_{}.sh".format(run_num), "w") as f:
    with open("{}run_job.sh".format(dir), "w") as f:
        for line in lines:
            f.write(line)
    # os.system("chmod +x run_job_{}.sh".format(run_num))
    os.system("chmod +x {}run_job.sh".format(dir))
    
def main():
    experiment = "/home/mameen/1_cell_increase_2_sigma/" 
    # os.system("cp -r /home/mameen/init_homogenous/ {}".format(experiment))
    stresses = np.loadtxt(experiment + "stresses.txt")
    # script = "/home/mameen/scripts/single_pattern.py"
    script = "single_pattern.py"

    for i in range(100):
        run_dir = experiment + "{:03d}/".format(i)
        create_exec_file(script,run_dir)
        create_sub_file(script,run_dir)
        # os.system("echo '0\n{}' > {}stress_limits".format(np.mean(stresses),run_dir))
        # os.system("echo '{}\n100' > {}stress_limits".format(np.mean(stresses),run_dir))
        # os.system("echo '{}' > {}target".format(np.mean(stresses) + 2*np.std(stresses),run_dir))
        os.system("echo '{}' > {}target".format(np.mean(stresses),run_dir))
        os.system("echo '1e-13' > {}tolerance".format(run_dir))
        os.system("echo '4' > {}n_cells".format(run_dir))
        # os.system("rm {}error.txt {}output.txt {}log.txt".format(run_dir,run_dir,run_dir))
        #submit job
        os.system("cd {} && condor_submit {}run_job.sub".format(run_dir,run_dir))

def resubmit():
    for experiment in ["/home/mameen/1_cell_increase_2_sigma/","/home/mameen/1_cell_decrease_2_sigma/","2_cells_mean/"]: 
        for i in range(100):
            run_dir = experiment + "{:03d}/".format(i)
            os.system("cd {} && condor_submit {}run_job.sub".format(run_dir,run_dir))


if __name__ == "__main__":
    resubmit()
# import os
# # for experiment in ["sp_2_cell_increase",
# #                    "sp_2_cell_decrease",
# #                    "sp_5_cell_decrease"]:
# # for experiment in ["sp_5_cell_increase",
# #                    "sp_5_cell_decrease"]:
# # for experiment in ["mp_2_cells"]:
# for experiment
#     for i in range(20):
#         dir = "/home/mameen/{}/run_{}/".format(experiment,i)
#         print(dir)
#         if not os.path.isfile(dir + "distances.txt"):
#             continue
#         with open(dir + "distances.txt", "r") as f:
#             lines = f.readlines()
#             print(len(lines), lines[-1])
#         for pattern in ["patternA/", "patternB/"]:
#             if not os.path.isfile(dir + pattern + "costs.txt"):
#                 continue
#             with open(dir + pattern + "costs.txt", "r") as f:
#                 lines = f.readlines()
#                 if len(lines) < 2:
#                     continue
#                 print(pattern, len(lines), lines[-1])