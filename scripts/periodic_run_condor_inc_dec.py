import os
import sys
import numpy as np
from toolbox.periodic import PeriodicTissue
from toolbox.spheroid import Spheroid
from toolbox.patterns import Patterns
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

# def resubmit(experiment):
#     for i in range(100):
#         run_dir = experiment + "{:03d}/".format(i)
#         if os.path.isfile(run_dir + "error.txt"):
#             if os.path.isfile(run_dir + "costs.txt"):
#                 print(run_dir, "error file and costs exists, deleting last iteration and editing conf")
#                 remove_last_iteration(run_dir)
#             os.system("rm {}error.txt".format(run_dir))
#         os.system("echo '1e-12' > {}tolerance".format(run_dir))
#         os.system("cd {} && condor_submit {}run_job.sub".format(run_dir,run_dir))
    

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

def start_increase_runs():
    script = "single_pattern.py"
    tolerance = 1e-6
    system_sizes = [4,5,6]
    n_cells = [1]
    n_runs = 100
    max_iters = 50000
    clear_interval = 50
    learning_rate = 10
    for l in system_sizes:
        stresses = np.loadtxt("init/kv_10_l_{}/stresses.txt".format(l))
        target = np.mean(stresses) + 2* np.std(stresses)
        for n in n_cells:
            experiment = "/home/mameen/kv_10_l_{}_n_{}_increase/".format(l,n)
            for i in range(n_runs):
                run_dir = experiment + "{:03d}/".format(i)
                create_exec_file(script,run_dir)
                create_sub_file(script,run_dir)
                os.system("echo {} > {}target".format(target,run_dir))
                os.system("echo {} > {}tolerance".format(tolerance,run_dir))
                os.system("echo {} > {}max_iters".format(max_iters,run_dir))
                os.system("echo {} > {}clear_interval".format(clear_interval,run_dir))
                os.system("echo {} > {}learning_rate".format(learning_rate,run_dir))
                sample = PeriodicTissue.from_config(run_dir, "minimized.txt")
                Patterns.find_random_target_cells(sample, n_cells = n)
                #submit job
                print("Submitting job in dir:", run_dir)
                os.system("cd {} && condor_submit {}run_job.sub".format(run_dir,run_dir))

def start_decrease_runs():
    script = "single_pattern.py"
    tolerance = 1e-6
    system_sizes = [4,5,6]
    n_cells = [1]
    n_runs = 100
    max_iters = 50000
    clear_interval = 50
    learning_rate = 10
    for l in system_sizes:
        stresses = np.loadtxt("init/kv_10_l_{}/stresses.txt".format(l))
        target = np.mean(stresses) - 2* np.std(stresses)
        for n in n_cells:
            experiment = "/home/mameen/kv_10_l_{}_n_{}_decrease/".format(l,n)
            for i in range(n_runs):
                run_dir = experiment + "{:03d}/".format(i)
                create_exec_file(script,run_dir)
                create_sub_file(script,run_dir)
                os.system("echo {} > {}target".format(target,run_dir))
                os.system("echo {} > {}tolerance".format(tolerance,run_dir))
                os.system("echo {} > {}max_iters".format(max_iters,run_dir))
                os.system("echo {} > {}clear_interval".format(clear_interval,run_dir))
                os.system("echo {} > {}learning_rate".format(learning_rate,run_dir))
                sample = PeriodicTissue.from_config(run_dir, "minimized.txt")
                Patterns.find_random_target_cells(sample, n_cells = n)
                #submit job
                print("Submitting job in dir:", run_dir)
                os.system("cd {} && condor_submit {}run_job.sub".format(run_dir,run_dir))

def remove_bad_runs():
    system_sizes = [4,5,6]
    n_cells = [1,2,3,4,5,6]
    n_runs = 100
    to_remove = []
    for l in system_sizes:
        for n in n_cells:
            experiment = "/home/mameen/kv_10_l_{}_n_{}/".format(l,n)
            for i in range(n_runs):
                run_dir = experiment + "{:03d}/".format(i)
                if not os.path.isfile(run_dir + "error.txt"):
                    continue
                
                with open(run_dir + "error.txt", 'r') as f:
                    lines = f.readlines()
                    if lines[-1].split()[0] == 'WARNING:':
                        continue
                    print(run_dir)
                    if lines[-1].split()[0] == 'ValueError:':
                        print("known error")
                        to_remove.append(run_dir)
                    elif lines[-1].split()[0] == 'numpy.linalg.LinAlgError:':
                        print("numpy error")
                        to_remove.append(run_dir)
                    elif lines[-1].split()[0] == 'OSError:':
                        print("OSerror")
                        to_remove.append(run_dir)
                    elif lines[-1].split()[0] == 'IndexError:':
                        print("Index Error")
                        to_remove.append(run_dir)
                    elif lines[-1].split()[0] == 'pandas.errors.EmptyDataError:':
                        print("Pandas Error")
                        to_remove.append(run_dir)
                    else:
                        print("Other error:")
                        print(lines[-1].split())
    for dir in to_remove:
        os.system("rm -rf {}".format(dir))
    print("Removed bad runs, if any.")
def check_costs():
    system_sizes = [4,5,6]
    n_cells = [1,2,3,4,5,6]
    n_runs = 100
    for l in system_sizes:
        for n in n_cells:
            experiment = "/home/mameen/kv_10_l_{}_n_{}/".format(l,n)
            for i in range(n_runs):
                run_dir = experiment + "{:03d}/".format(i)
                if not os.path.isdir(run_dir):
                    continue
                if os.path.isfile(run_dir + "costs.txt"):
                    costs = np.loadtxt(run_dir + "costs.txt")
                    print(run_dir, costs[-1])
                else:
                    print(run_dir, " No costs file yet!")
def resubmit():
    script = "single_pattern.py"
    tolerance = 1e-6
    system_sizes = [4,5,6]
    n_cells = [1,2,3,4,5,6]
    n_runs = 100
    max_iters = 50000
    clear_interval = 50
    learning_rate = 10
    for l in system_sizes:
        for n in n_cells:
            experiment = "/home/mameen/kv_10_l_{}_n_{}/".format(l,n)
            for i in range(n_runs):
                run_dir = experiment + "{:03d}/".format(i)
                if os.path.isdir(run_dir):
                    continue
                stresses = np.loadtxt("init/kv_10_l_{}/stresses.txt".format(l))
                target = np.mean(stresses)
                os.system("cp -r /home/mameen/init/kv_10_l_{}/{:03d} {}".format(l,i, run_dir))
                create_exec_file(script,run_dir)
                create_sub_file(script,run_dir)
                os.system("echo {} > {}target".format(target,run_dir))
                os.system("echo {} > {}tolerance".format(tolerance,run_dir))
                os.system("echo {} > {}max_iters".format(max_iters,run_dir))
                os.system("echo {} > {}clear_interval".format(clear_interval,run_dir))
                os.system("echo {} > {}learning_rate".format(learning_rate,run_dir))
                sample = PeriodicTissue.from_config(run_dir, "minimized.txt")
                Patterns.find_random_target_cells(sample, n_cells = n)
                #submit job
                print("Submitting job in dir:", run_dir)
                os.system("cd {} && condor_submit {}run_job.sub".format(run_dir,run_dir))
if __name__ == "__main__":
    # remove_bad_runs()
    # check_costs()
    start_decrease_runs()
    start_increase_runs()