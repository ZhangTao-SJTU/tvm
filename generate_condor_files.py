import os
import sys

def create_sub_file(script,dir):
    lines = []
    lines.append("executable = /home/mameen/examples/singularity_wrapper.sh\n")
    lines.append("arguments  = /home/mameen/examples/ubuntu18_povray_paula.img {}run_job.sh\n".format(dir))
    lines.append("transfer_input_files = {}, /home/mameen/setup.py, /home/mameen/toolbox \n".format(script))
    lines.append("should_transfer_files = YES\n")
    lines.append("output     = output.txt\n")
    lines.append("error      = error.txt\n")
    lines.append("log        = log.txt\n")
    lines.append("getenv     = True\n")
    lines.append("request_cpus = 1\n")
    lines.append("request_memory = 1 GB\n")
    lines.append('Requirements = TARGET.vm_name == "its-u20-nfs-20210413" && regexp("CRUSH", TARGET.name)\n')


    lines.append("queue \n")
    # with open("job_{}.sub".format(run_num), "w") as f:
    with open("{}run_job.sub".format(dir), "w") as f:
        for line in lines:
            f.write(line)


def create_exec_file(script,init_dir, run_dir):
    lines = []
    lines.append("#!/bin/bash\n")
    lines.append("source /home/mameen/.bashrc\n")
    lines.append("pip install --user -e .\n")
    # lines.append("mkdir -p {}\n".format(run_dir))
    lines.append("python {} {} {}\n".format(script, init_dir, run_dir))
    # with open("run_job_{}.sh".format(run_num), "w") as f:
    with open("{}run_job.sh".format(run_dir), "w") as f:
        for line in lines:
            f.write(line)
    # os.system("chmod +x run_job_{}.sh".format(run_num))
    os.system("chmod +x {}run_job.sh".format(run_dir))
    
def main():
    for experiment in ["sp_2_cell_increase",
                       "sp_2_cell_decrease",
                       "sp_5_cell_decrease"]:

        head_dir ="/home/mameen/{}/".format(experiment)
        script = "/home/mameen/{}.py".format(experiment)
        os.makedirs(head_dir, exist_ok=True)
        for i in range(20):
            init_dir = "/home/mameen/init_homogeneous/7_{}/".format(i)
            run_dir = head_dir + "run_{}/".format(i)
            os.makedirs(run_dir, exist_ok=True)
            create_exec_file(script,init_dir, run_dir)
            create_sub_file(script,run_dir)
            os.system("cd {} && condor_submit {}run_job.sub".format(run_dir,run_dir))

if __name__ == "__main__":
    main()
import os
# for experiment in ["sp_2_cell_increase",
#                    "sp_2_cell_decrease",
#                    "sp_5_cell_decrease"]:
# for experiment in ["sp_5_cell_increase",
#                    "sp_5_cell_decrease"]:
for experiment in ["mp_2_cells"]:
    for i in range(20):
        dir = "/home/mameen/{}/run_{}/".format(experiment,i)
        print(dir)
        if not os.path.isfile(dir + "distances.txt"):
            continue
        with open(dir + "distances.txt", "r") as f:
            lines = f.readlines()
            print(len(lines), lines[-1])
        for pattern in ["patternA/", "patternB/"]:
            if not os.path.isfile(dir + pattern + "costs.txt"):
                continue
            with open(dir + pattern + "costs.txt", "r") as f:
                lines = f.readlines()
                if len(lines) < 2:
                    continue
                print(pattern, len(lines), lines[-1])