# Single fixed cell pulse drive
from toolbox import training
tr = training.Training.from_config(config_dir = "10_0/", input_filename = "minimized.txt")
tr.set_expansion_factor(0.25)
tr.pick_fixed_cell()
tr.pulse_drive(n_iterations = 10)