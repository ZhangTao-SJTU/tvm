from toolbox.patterns import Patterns
from toolbox.periodic import PeriodicTissue
from toolbox.stress import calculate_max_shear_stress
from toolbox.pattern_trainer import find_random_target_cells, train_target_cells, resume_run
import os
import numpy as np
import pandas as pd
import glob
import sys

class multiple_patterns:
    def __init__(self):
        self._run_dir = None
        self._target_cell_to_stress_A = None
        self._target_cell_to_stress_B = None
        self._target_stress = None
        self._tolerance = 1e-8
        self._n_cells_A = 2
        self._n_cells_B = 2
        self._max_iters = 10
        self._convergence_check_interval = 10
        self._learning_rate = 10
        self._net_error = None
        self._distance = None
        self._epoch = 0
        self._current_pattern = None

    @classmethod
    def from_dir(cls, dir):
        inst = cls()
        inst._run_dir = dir
        # set cpp_executable_dir based on existing directories
        if os.path.isdir("/Users/shabeebameen/Projects/tvm-fire/build/"):
            cls._cpp_executable_dir = "/Users/shabeebameen/Projects/tvm-fire/build/"
        elif os.path.isdir("/home/shabeeb/Projects/tvm-fire/build/"):
            cls._cpp_executable_dir = "/home/shabeeb/Projects/tvm-fire/build/"
        elif os.path.isdir("/home/mameen/tvm/build/"):
            cls._cpp_executable_dir = "/home/mameen/tvm/build/"

        if os.path.isfile("{}target".format(dir)):
            with open("{}target".format(dir), "r") as f:
                inst._target_stress = float(f.read().strip())
        if os.path.isfile("{}tolerance".format(dir)):
            with open("{}tolerance".format(dir), "r") as f:
                inst._tolerance = float(f.read().strip())
        if os.path.isfile("{}learning_rate".format(dir)):
            with open("{}learning_rate".format(dir), "r") as f:
                inst._learning_rate = float(f.read().strip())
        if os.path.isfile("{}max_iters".format(dir)):
            with open("{}max_iters".format(dir), "r") as f:
                inst._max_iters = int(f.read().strip())
        if os.path.isfile("{}n_cells_A".format(dir)):
            with open("{}n_cells_A".format(dir), "r") as f:
                inst._n_cells_A = int(f.read().strip())
        if os.path.isfile("{}n_cells_B".format(dir)):
            with open("{}n_cells_B".format(dir), "r") as f:
                inst._n_cells_B = int(f.read().strip())
        # create directory for saving end-of-epoch files
        os.makedirs(inst._run_dir + "files/", exist_ok=True)
        return inst
  
    # set self._target_cell_to_stress_A and self._target_cell_to_stress_B
    # requires initialization of self._target_stress (can be read from dir/target if using alternate constructor.
    # This will be the uniform target stress for all target cells in both patterns.
    # Save the corresponding vtks and initial_stresses.csv files in the run_dir
    def set_new_uniform_target_stress_patterns(self):
        target_cells_A = find_random_target_cells(self._run_dir, n_cells = self._n_cells_A,output_vtk_file = "target_cells_A.vtk")
        self._target_cell_to_stress_A ={i:self._target_stress for i in target_cells_A}
        target_cells_B = find_random_target_cells(self._run_dir, n_cells = self._n_cells_B, exclude_cells = target_cells_A, output_vtk_file = "target_cells_B.vtk")
        self._target_cell_to_stress_B ={i:self._target_stress for i in target_cells_B}
        # Save initial_stress_A/B.csv files
        for pattern, target_dict in {"A":self._target_cell_to_stress_A, "B":self._target_cell_to_stress_B}.items():
            tissue = PeriodicTissue.from_config(self._run_dir, "minimized.txt")
            initial_stress = {cellID: calculate_max_shear_stress(tissue,cellID)for cellID in target_dict}
            df = pd.DataFrame(list(initial_stress.items()), columns=['CellID', 'Current'])
            # additionally save the uniform target stress column
            df['Target'] = self._target_stress       
            df.to_csv("{}initial_stress_{}.csv".format(self._run_dir,pattern), index=False)
    
    def load_target_stress_patterns(self):
        for pattern, target_dict in {"A":self._target_cell_to_stress_A,"B":self._target_cell_to_stress_B}.items():
            df = pd.read_csv("{}initial_stress_{}.csv".format(self._run_dir,pattern))
            target_dict = dict(zip(df["CellID"].to_numpy().astype(int), df["Target"].to_numpy().astype(float)))
        
    def single_iteration(self):
        # start a new run if no costs.txt file exists
        # Otherwise, use resume_run
        # Either way, first train pattern A for self._max_iters
        self._current_pattern = "A"
        print("Training pattern {}".format(self._current_pattern))
        if not os.path.isfile("{}costs.txt".format(self._run_dir)):
            print("...Starting a new run.")
            train_target_cells(self._run_dir, self._target_cell_to_stress_A, learning_rate=self._learning_rate, max_iters=self._max_iters, cpp_executable_dir=self._cpp_executable_dir, tolerance=self._tolerance)
        else:
            resume_run(self._run_dir, target_cell_to_stress = self._target_cell_to_stress_A, learning_rate=self._learning_rate, max_iters=self._max_iters, cpp_executable_dir=self._cpp_executable_dir, tolerance=self._tolerance)
        self.write_info()
        # Now train pattern B for self._max_iters
        self._current_pattern = "B"
        print("Training pattern {}".format(self._current_pattern))
        resume_run(self._run_dir, target_cell_to_stress = self._target_cell_to_stress_B, learning_rate=self._learning_rate, max_iters=self._max_iters, cpp_executable_dir=self._cpp_executable_dir, tolerance=self._tolerance)
        self.write_info()
        self.clear_dir()
        self._epoch += 1

    def run(self, epochs = 100000):
        print(self._target_cell_to_stress_A)
        print(self._target_cell_to_stress_B)
        print("Tolerance:", self._tolerance)
        print("Max iterations:", self._max_iters)
        print("cpp_executable_dir:", self._cpp_executable_dir)
        for i in range(epochs):
            self.single_iteration()
            if self._net_error < self._tolerance:
                print("Converged with net error:", self._net_error)
                break
            # check if distance is not changing every convergence_check_interval epochs
            if i>0 and not i%self._convergence_check_interval:
                distances = pd.read_csv("{}info.csv".format(self._run_dir))["Distance"].to_numpy()
                if np.allclose(distances[-self._convergence_check_interval:], distances[-1]):
                    print("Parameter space distance did not change for the last {} epochs.".format(self._convergence_check_interval))
                    break
    def evaluate_net_error(self):
        tissue = PeriodicTissue.from_config(self._run_dir, "minimized.txt")
        trainer = Patterns.periodic_tissue(tissue)
        target_stress = {}
        for p in ["A","B"]:
            df = pd.read_csv("{}initial_stress_{}.csv".format(self._run_dir,p))
            for _, row in df.iterrows():
                cellID = int(row["CellID"])
                target_stress[cellID] = float(row["Target"])
        trainer.set_target_cell_to_stress(target_stress)
        trainer.load_cell_parameters()
        self._net_error = trainer.evaluate_cost()

    def evaluate_parameter_space_distance(self):
        fileA = "{}cellParameters.init.input".format(self._run_dir)
        fileB = "{}cellParameters.input".format(self._run_dir)
        df_A = pd.read_csv(fileA, sep=" ",header=None)
        df_B = pd.read_csv(fileB, sep=" ",header=None)
        if not (df_A[0].to_numpy() == df_B[0].to_numpy()).all():
            raise ValueError("Files have different number of cells {} vs {}".format(df_A.shape[0], df_B.shape[0]))
        self._distance = np.sqrt(np.sum((df_A[2].to_numpy()-df_B[2].to_numpy())**2))
        
    # Moves files to files/ directory and appends info to info.csv file

    def write_info(self):
        errors = np.loadtxt("{}costs.txt".format(self._run_dir))
        iter = len(errors)-1
        for file in ["cellParameters.input","bulk.txt","stresses.csv"]:
            filename = "{}{:07d}.{}".format(self._run_dir,iter,file)
            os.system("cp {} {}".format(filename, "{}files/".format(self._run_dir)))
        self.evaluate_net_error()
        self.evaluate_parameter_space_distance()
        if not os.path.isfile("{}info.csv".format(self._run_dir)):
            with open("{}info.csv".format(self._run_dir), "w") as f:
                f.write("Epoch,Iter,Pattern,Error,Distance\n")
        with open("{}info.csv".format(self._run_dir), "a") as f:
            f.write("{},{},{},{},{}\n".format(self._epoch, iter, self._current_pattern, self._net_error, self._distance))
    def clear_dir(self):
        for file in ["cellParameters.input","bulk.txt","stresses.csv"]:
            list_of_files = sorted(glob.glob("{}*.{}".format(self._run_dir,file)))
            # keep only the latest file
            for f in list_of_files[:-1]:
                os.remove(f)
def main():
    if len(sys.argv) < 2:
        print("Usage: python multiple_patterns.py <run_dir>")
        sys.exit(1)
    run_dir = sys.argv[1]
    trainer = multiple_patterns.from_dir(run_dir)
    if not os.path.isfile("{}costs.txt".format(run_dir)):
        trainer.set_new_uniform_target_stress_patterns()
    else:
        trainer.load_target_stress_patterns()
    trainer.run()

if __name__ == "__main__":
    main()