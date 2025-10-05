import os

def intialize_test_dir(input_dir,output_dir):
    os.system("cp -r {} {}".format(input_dir,output_dir))
    os.system("echo '1e-12' > {}tolerance".format(output_dir))
    os.system("echo '4' > {}n_cells".format(output_dir))
    os.system("echo '0.23295009545919323' > {}target".format(output_dir))
    os.system("echo '10' > {}learning_rate".format(output_dir))

if __name__ == "__main__":
    input_dir = "/home/shabeeb/Projects/tvm-fire/init_homogeneous_4/000/"
    output_dir = "/home/shabeeb/Projects/tvm-fire/test_4_cells/"
    intialize_test_dir(input_dir,output_dir)
