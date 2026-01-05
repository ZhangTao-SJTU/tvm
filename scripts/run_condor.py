import os
import sys
import numpy as np
from toolbox.pattern_trainer import remove_last_iteration

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
    lines.append("request_memory = 500 MB\n")
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

def resubmit(experiment):
    for i in range(100):
        run_dir = experiment + "{:03d}/".format(i)
        if os.path.isfile(run_dir + "error.txt"):
            if os.path.isfile(run_dir + "costs.txt"):
                print(run_dir, "error file and costs exists, deleting last iteration and editing conf")
                remove_last_iteration(run_dir)
            os.system("rm {}error.txt".format(run_dir))
        os.system("echo '1e-12' > {}tolerance".format(run_dir))
        os.system("cd {} && condor_submit {}run_job.sub".format(run_dir,run_dir))
    
def main():
    stresses = np.loadtxt("init/init_homogeneous_6/stresses.txt")
    target = np.mean(stresses)
    tolerance = 1e-6
    # for experiment in ["/home/mameen/3_cells_mean_l_6/"]:
    #     n_cells = 3
    #     script = "single_pattern.py"
    #     for i in range(100):
    #         run_dir = experiment + "{:03d}/".format(i)
    #         create_exec_file(script,run_dir)
    #         create_sub_file(script,run_dir)
    #         # os.system("echo '0\n{}' > {}stress_limits".format(np.mean(stresses),run_dir))
    #         # os.system("echo '{}\n100' > {}stress_limits".format(np.mean(stresses),run_dir))
    #         # os.system("echo '{}' > {}target".format(np.mean(stresses) + 2*np.std(stresses),run_dir))
    #         # os.system("echo '{}' > {}target".format(np.mean(stresses),run_dir))
    #         os.system("echo '{}' > {}target".format(target,run_dir))
    #         os.system("echo '{}' > {}tolerance".format(tolerance,run_dir))
    #         os.system("echo '{}' > {}n_cells".format(n_cells,run_dir))
    #         # os.system("rm {}error.txt {}output.txt {}log.txt".format(run_dir,run_dir,run_dir))
    #         #submit job
    #         os.system("cd {} && condor_submit {}run_job.sub".format(run_dir,run_dir))
    for experiment in ["/home/mameen/6_cells_mean_l_6/"]:
        n_cells = 6
        script = "single_pattern.py"
        for i in range(100):
            run_dir = experiment + "{:03d}/".format(i)
            create_exec_file(script,run_dir)
            create_sub_file(script,run_dir)
            # os.system("echo '0\n{}' > {}stress_limits".format(np.mean(stresses),run_dir))
            # os.system("echo '{}\n100' > {}stress_limits".format(np.mean(stresses),run_dir))
            # os.system("echo '{}' > {}target".format(np.mean(stresses) + 2*np.std(stresses),run_dir))
            # os.system("echo '{}' > {}target".format(np.mean(stresses),run_dir))
            os.system("echo '{}' > {}target".format(target,run_dir))
            os.system("echo '{}' > {}tolerance".format(tolerance,run_dir))
            os.system("echo '{}' > {}n_cells".format(n_cells,run_dir))
            # os.system("rm {}error.txt {}output.txt {}log.txt".format(run_dir,run_dir,run_dir))
            #submit job
            os.system("cd {} && condor_submit {}run_job.sub".format(run_dir,run_dir))
    # resubmit("/home/mameen/2_cells_mean_l_4/")
    # resubmit("/home/mameen/2_cells_mean_l_5/")
    # resubmit("/home/mameen/2_cells_mean_l_6/")
    # resubmit("/home/mameen/4_cells_mean_l_4/")
    # resubmit("/home/mameen/4_cells_mean_l_5/")
    # resubmit("/home/mameen/4_cells_mean_l_6/")

def multiple_patterns():
    script = "multiple_patterns.py"
    stresses = np.loadtxt("init/init_homogeneous_6/stresses.txt")
    target = np.mean(stresses)
    tolerance = 1e-6
    max_iters = 500

    def submit_jobs(script,run_dir,target,tolerance,n_cells_A,n_cells_B,max_iters):
        create_exec_file(script,run_dir)
        create_sub_file(script,run_dir)
        os.system("echo '{}' > {}target".format(target,run_dir))
        os.system("echo '{}' > {}tolerance".format(tolerance,run_dir))
        os.system("echo '{}' > {}n_cells_A".format(n_cells_A,run_dir))
        os.system("echo '{}' > {}n_cells_B".format(n_cells_B,run_dir))
        os.system("echo '{}' > {}max_iters".format(max_iters,run_dir))
        #submit job
        os.system("cd {} && condor_submit {}run_job.sub".format(run_dir,run_dir))

    for experiment in ["/home/mameen/1_1_l_4/","/home/mameen/1_1_l_5/","/home/mameen/1_1_l_6/"]:
        n_cells_A = 1
        n_cells_B = 1
        for i in range(100):
            run_dir = experiment + "{:03d}/".format(i)
            submit_jobs(script,run_dir,target,tolerance,n_cells_A,n_cells_B,max_iters)
    for experiment in ["/home/mameen/1_2_l_4/","/home/mameen/1_2_l_5/","/home/mameen/1_2_l_6/"]:
        n_cells_A = 1
        n_cells_B = 2
        for i in range(100):
            run_dir = experiment + "{:03d}/".format(i)
            submit_jobs(script,run_dir,target,tolerance,n_cells_A,n_cells_B,max_iters)
    for experiment in ["/home/mameen/2_2_l_4/","/home/mameen/2_2_l_5/","/home/mameen/2_2_l_6/"]:
        n_cells_A = 2
        n_cells_B = 2
        for i in range(100):
            run_dir = experiment + "{:03d}/".format(i)
            submit_jobs(script,run_dir,target,tolerance,n_cells_A,n_cells_B,max_iters)

def single_pattern_by_cell():
    def submit_jobs(script,run_dir,target,tolerance,n_cells,max_iters):
        create_exec_file(script,run_dir)
        create_sub_file(script,run_dir)
        os.system("echo '{}' > {}target".format(target,run_dir))
        os.system("echo '{}' > {}tolerance".format(tolerance,run_dir))
        os.system("echo '{}' > {}n_cells".format(n_cells,run_dir))
        os.system("echo '{}' > {}max_iters".format(max_iters,run_dir))
        #submit job
        os.system("cd {} && condor_submit {}run_job.sub".format(run_dir,run_dir))
    
    script = "single_pattern_by_cell.py"
    l = 4
    stresses = np.loadtxt("init/init_homogeneous_{}/stresses.txt".format(l))
    target = np.mean(stresses)
    tolerance = 1e-6
    max_iters = 500
    for n_cells in [3,4,5,6,7]:
        experiment = "/home/mameen/new_cells_{}_l_{}/".format(n_cells,l)
        for i in range(10):
            run_dir = experiment + "{:03d}/".format(i)
            submit_jobs(script,run_dir,target,tolerance,n_cells,max_iters)
    
if __name__ == "__main__":
    single_pattern_by_cell()
    # resubmit()