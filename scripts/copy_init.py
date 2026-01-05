import os

def main():
    l = 4
    for n_cells in [3,4,5,6,7]:
        experiment = "/home/mameen/new_cells_{}_l_{}/".format(n_cells,l)
        os.makedirs(experiment, exist_ok=True)
        for i in range(10):
            os.system("cp -r /home/mameen/init/init_homogeneous_{}/{:03d} {}".format(l,i,experiment))
if __name__ == "__main__":
    main()