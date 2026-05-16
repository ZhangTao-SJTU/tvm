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
        # for file in ["q_values.txt","costs.txt","minimized.txt","cellParameters.input", "files/0000000.stresses.csv"]:
        for file in ["SD_s0.csv"]:
            os.system("scp {} /Users/shabeebameen/Projects/tvm-fire/{}".format(cluster_folder+file,subfolder))
def download_single_patterns_old(n):
    local_folder = f"data/{n}_cells_mean_l_6/"
    os.makedirs(local_folder,exist_ok=True)
    for i in range(100):
        subfolder = local_folder + "{:03d}/".format(i)
        os.makedirs(subfolder,exist_ok=True)
        cluster_folder = f"mameen@smatter-login.syr.edu:/home/mameen/{n}_cells_mean_l_6/{i:03d}/"
        # for file in ["q_values.txt","costs.txt","errors.txt","minimized.txt","cellParameters.input", "0000.stresses.csv","initial_stress.csv"]:
        for file in ["errors.txt"]:
            os.system(f"scp {cluster_folder+file} /Users/shabeebameen/Projects/tvm-fire/{subfolder}")

def download_single_patterns_old_inc_dec(inc_or_dec):
    local_folder = f"data/1_cell_{inc_or_dec}_2_sigma/"
    os.makedirs(local_folder,exist_ok=True)
    for i in range(100):
        subfolder = local_folder + "{:03d}/".format(i)
        os.makedirs(subfolder,exist_ok=True)
        cluster_folder = f"mameen@smatter-login.syr.edu:/home/mameen/1_cell_{inc_or_dec}_2_sigma/{i:03d}/"
        for file in ["q_values.txt","costs.txt","errors.txt","minimized.txt","cellParameters.input", "0000.stresses.csv","initial_stress.csv"]:
            os.system(f"scp {cluster_folder+file} /Users/shabeebameen/Projects/tvm-fire/{subfolder}")

def download_single_patterns_inc_dec(l,n=1, inc_or_dec="increase"):
    local_folder = "data/kv_10_l_{}_n_{}_{}/".format(l,n,inc_or_dec)
    os.makedirs(local_folder,exist_ok=True)
    for i in range(100):
        subfolder = local_folder + "{:03d}/".format(i)
        os.makedirs(subfolder,exist_ok=True)
        cluster_folder = "mameen@smatter-login.syr.edu:/home/mameen/kv_10_l_{}_n_{}_{}/{:03d}/".format(l,n,inc_or_dec,i)
        for file in ["q_values.txt","costs.txt","minimized.txt","cellParameters.input", "files/0000000.stresses.csv"]:
            # if os.path.isfile("/Users/shabeebameen/Projects/tvm-fire/{}".format(subfolder+file)):
            #     print("{} already exists".format(subfolder+file))
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
        # for file in ["q_values.txt","costs.txt","minimized.txt","cellParameters.input", "files/0000000.stresses.csv"]:
        for file in ["SD_s0.csv"]:

            # if os.path.isfile("/Users/shabeebameen/Projects/tvm-fire/{}".format(subfolder+file)):
            #     print("{} already exists".format(subfolder+file))
            #     continue
            os.system("scp {} /Users/shabeebameen/Projects/tvm-fire/{}".format(cluster_folder+file,subfolder))

def download_multiple_patterns(subpattern):
    l=4
    folder = "kv_10_l_{}_p".format(l)
    for p in subpattern:
        folder+="_{}".format(p)
    folder+="/"
    local_folder = "data/"+folder
    os.makedirs(local_folder,exist_ok=True)
    for i in range(100):
        subfolder = local_folder + "{:03d}/".format(i)
        os.makedirs(subfolder,exist_ok=True)
        cluster_folder = "mameen@smatter-login.syr.edu:/home/mameen/"+folder+"{:03d}/".format(i)
        # os.system("scp {}files/0000000.stresses.csv /Users/shabeebameen/Projects/tvm-fire/{}".format(cluster_folder,subfolder))

        for file in ["q_values.txt",
                     "costs.txt",
                     "minimized.txt",
                     "cellParameters.input",
                     "files/0000000.stresses.csv",
                     "info.csv"]:
            os.system("scp {} /Users/shabeebameen/Projects/tvm-fire/{}".format(cluster_folder+file,subfolder))


def download_multiple_patterns_new(subpattern):
    l=4
    folder = "new_kv_10_l_{}_p".format(l)
    for p in subpattern:
        folder+="_{}".format(p)
    folder+="/"
    local_folder = "data/"+folder
    os.makedirs(local_folder,exist_ok=True)
    for i in range(100):
        subfolder = local_folder + "{:03d}/".format(i)
        os.makedirs(subfolder,exist_ok=True)
        cluster_folder = "mameen@smatter-login.syr.edu:/home/mameen/"+folder+"{:03d}/".format(i)
        # os.system("scp {}files/0000000.stresses.csv /Users/shabeebameen/Projects/tvm-fire/{}".format(cluster_folder,subfolder))

        for file in ["q_values.txt",
                     "costs.txt",
                     "minimized.txt",
                     "cellParameters.input",
                     "files/0000000.stresses.csv",
                     "info.csv"]:
            os.system("scp {} /Users/shabeebameen/Projects/tvm-fire/{}".format(cluster_folder+file,subfolder))

if __name__ == "__main__":
    # subpattern = sys.argv[1:]
    # print(subpattern)
    # download_multiple_patterns_new(subpattern)

    # l = sys.argv[1]
    # inc_or_dec = sys.argv[2]
    # download_single_patterns_inc_dec(l,inc_or_dec=inc_or_dec)

    l = sys.argv[1]
    n = sys.argv[2]
    download_single_patterns(l,n)
    # download_single_patterns_spheroid(l,n)

    #download old data

    # n = sys.argv[1]
    # download_single_patterns_old(n)