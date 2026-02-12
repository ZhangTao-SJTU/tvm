import os
import sys

def download_single_pattern_cellwise():
    for n in [6]:
        local_folder = "data/single_pattern_cellwise_n_{}/".format(n)
        os.makedirs(local_folder,exist_ok=True)
        for i in range(5,10):
            subfolder = local_folder + "{:03d}/".format(i)
            os.makedirs(subfolder,exist_ok=True)
            cluster_folder = "mameen@smatter-login.syr.edu:/home/mameen/single_pattern_cellwise_n_{}_l_4/{:03d}/".format(n,i)
            for file in ["info.csv", "init_config.txt","initial_cost.txt","q_values.txt","costs.txt"]:
                if os.path.isfile("/Users/shabeebameen/Projects/tvm-fire/{}".format(subfolder+file)):
                    print("{} already exists".format(subfolder+file))
                    continue
                os.system("scp {} /Users/shabeebameen/Projects/tvm-fire/{}".format(cluster_folder+file,subfolder))

def download_single_patterns(l,n):
    local_folder = "data/kv_10_l_{}_n_{}/".format(l,n)
    os.makedirs(local_folder,exist_ok=True)
    for i in range(100):
        subfolder = local_folder + "{:03d}/".format(i)
        os.makedirs(subfolder,exist_ok=True)
        cluster_folder = "mameen@smatter-login.syr.edu:/home/mameen/kv_10_l_{}_n_{}/{:03d}/".format(l,n,i)
        for file in ["q_values.txt","costs.txt","minimized.txt","cellParameters.input"]:
            # if os.path.isfile("/Users/shabeebameen/Projects/tvm-fire/{}".format(subfolder+file)):
            #     print("{} already exists".format(subfolder+file))
            #     continue
            # if not os.path.isfile(cluster_folder + "minimized.txt"):
            #     continue
            os.system("scp {} /Users/shabeebameen/Projects/tvm-fire/{}".format(cluster_folder+file,subfolder))

def download_single_patterns_spheroid(l,n):
    l = int(l)
    n = int(n)
    local_folder = "data/kv_10_l_{}_n_sp_{:03d}/".format(l,n)
    os.makedirs(local_folder,exist_ok=True)
    for i in range(100):
        subfolder = local_folder + "{:03d}/".format(i)
        os.makedirs(subfolder,exist_ok=True)
        cluster_folder = "mameen@smatter-login.syr.edu:/home/mameen/kv_10_l_{}_n_sp_{:03d}/{:03d}/".format(l,n,i)
        os.system("scp {}files/0000000.stresses.csv /Users/shabeebameen/Projects/tvm-fire/{}".format(cluster_folder,subfolder))

        # for file in ["q_values.txt","costs.txt","minimized.txt","cellParameters.input"]:
        # for file in ["frozen_cells.txt",
        #              "spheroid_cells.txt",
        #              "target_cells.txt"
        #              ]:
            # if os.path.isfile("/Users/shabeebameen/Projects/tvm-fire/{}".format(subfolder+file)):
            #     print("{} already exists".format(subfolder+file))
            #     continue
            # if not os.path.isfile(cluster_folder + "minimized.txt"):
            #     continue
            # os.system("scp {} /Users/shabeebameen/Projects/tvm-fire/{}".format(cluster_folder+file,subfolder))

if __name__ == "__main__":
    # l = sys.argv[1]
    # n = sys.argv[2]
    # download_single_patterns_spheroid(l,n)
    for l in [5,6]:
        for n in [5,10,20,40]:
            download_single_patterns_spheroid(l,n)