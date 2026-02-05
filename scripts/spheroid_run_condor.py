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
    lines.append("queue \n")
    with open("{}run_job.sub".format(dir), "w") as f:
        for line in lines:
            f.write(line)

def create_exec_file(script,dir):
    lines = []
    lines.append("#!/bin/bash\n")
    lines.append("source /home/mameen/.bashrc\n")
    lines.append("export PYTHONPATH=$PWD:$PYTHONPATH\n")
    lines.append("python {} {}\n".format(script,dir))
    with open("{}run_job.sh".format(dir), "w") as f:
        for line in lines:
            f.write(line)
    os.system("chmod +x {}run_job.sh".format(dir))

def main():
    script = "single_pattern.py"
    tolerance = 1e-6
    system_sizes = [5,6]
    n_spheroids = [5,10,20,40]
    n_targets = 1
    n_runs = 100
    max_iters = 50000
    clear_interval = 50
    learning_rate = 10
    for l in system_sizes:
        stresses = np.loadtxt("init/kv_10_l_{}/stresses.txt".format(l))
        target = np.mean(stresses)
        for n_s in n_spheroids:
            experiment = "/home/mameen/kv_10_l_{}_n_sp_{:03d}/".format(l,n_s)
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
                Patterns.find_target_cells_in_spheroid(sample, n_spheroid = n_s, n_cells = 1)
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
    main()