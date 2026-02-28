from toolbox.patterns import Patterns
from toolbox.periodic import PeriodicTissue
from toolbox.multiplePatterns import MultiplePatterns
import os
import numpy as np
import pandas as pd
import sys

def find_cpp_executable_dir():
    candidates = [
        "/Users/shabeebameen/Projects/tvm-fire/build/",
        "/home/shabeeb/Projects/tvm-fire/build/",
        "/home/mameen/tvm/build/",
    ]
    for path in candidates:
        if os.path.isdir(path):
            return path
    raise FileNotFoundError("No tvm-fire build directory found")

def read_int_list(filepath):
    with open(filepath, "r") as f:
        return [int(line.strip()) for line in f.readlines()]

def load_config(run_dir):
    config = {}
    if os.path.isfile("{}target".format(run_dir)):
        config["target_stress"] = float(np.loadtxt("{}target".format(run_dir)))
    if os.path.isfile("{}target_cells.txt".format(run_dir)):
        config["target_cells"] = read_int_list("{}target_cells.txt".format(run_dir))
    if os.path.isfile("{}frozen_cells.txt".format(run_dir)):
        config["frozen_cells"] = read_int_list("{}frozen_cells.txt".format(run_dir))
    if os.path.isfile("{}tolerance".format(run_dir)):
        config["tolerance"] = float(np.loadtxt("{}tolerance".format(run_dir)))
    if os.path.isfile("{}learning_rate".format(run_dir)):
        config["learning_rate"] = float(np.loadtxt("{}learning_rate".format(run_dir)))
    if os.path.isfile("{}max_iters".format(run_dir)):
        config["max_iters"] = int(np.loadtxt("{}max_iters".format(run_dir)))
    if os.path.isfile("{}clear_interval".format(run_dir)):
        config["clear_interval"] = int(np.loadtxt("{}clear_interval".format(run_dir)))
    if os.path.isfile("{}subpattern_n_cells".format(run_dir)):
        config["subpattern_n_cells"] = read_int_list("{}subpattern_n_cells".format(run_dir))
    return config

def main():
    if len(sys.argv) != 2:
        raise ValueError("Usage: python run_multiple_patterns.py <run_dir>")
    run_dir = sys.argv[1]
    cfg = load_config(run_dir)
    # Build the trainer via the inherited classmethod chain
    tissue = PeriodicTissue.from_config(run_dir, "minimized.txt")
    trainer = Patterns.from_sample(tissue)
    # Apply settings
    trainer.set_cpp_executable_dir(find_cpp_executable_dir())
    trainer.set_target_cell_to_stress({cellID:cfg["target_stress"] for cellID in cfg["target_cells"]})
    trainer.set_tolerance(cfg.get("tolerance", 1e-6))
    trainer.set_learning_rate(cfg.get("learning_rate", 10))
    trainer.set_clear_interval(cfg.get("clear_interval", 100))
    trainer.set_frozen_cells(cfg.get("frozen_cells", []))
    trainer.initialize()
    
    multiple_pattern_trainer = MultiplePatterns.from_pattern(trainer,cfg["subpattern_n_cells"])
    multiple_pattern_trainer.run()

if __name__ == "__main__":
    main()