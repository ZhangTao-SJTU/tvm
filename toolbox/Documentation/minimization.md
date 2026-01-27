Class Hierarchy (Parent -> Child)
FIREminimization -> Training -> Patterns 

Class structure:

1. minimization.FIREminimization:

Attributes and their initializations: 
    self._config = None
    self._dir = None
    self._cpp_executable_dir = None

@classmethods:
from_sample(cls, tissue):
    sample = cls()
    sample._config = tissue
    sample._dir = tissue.config_dir_
    return sample

Routines:
minimize_config(FIRE_only = False):
    i. writes self._config to "sample.topo"
        NOTE: This WILL erase any existing sample.topo in self._dir
    ii. Rename self._dir/minimized.txt to self._dir/minimized_old.txt. This is done because the cpp executable would append to any existing minimized.txt file.

    iii. Run the tvm executable. This will run using the new sample.topo and any existing cellParameters.txt in self._dir.If code doesnt crash, it produces a minimized.txt file.
    
    iv. set file to minimized.txt and load the minimized configuration to self._config

    v. load_cell_parameters() loads from cellParameters.input to self.config

    Result: by side effect, a successful completion of this function will have (i) the minimized.txt (ii) loaded minimized configuration in self._config 
     
write_configuration(filename = "sample.topo"):
    write configuration in self._dir+filename
    current support for periodic and spheroid

load_cell_parameters(filename = "cellParameters.input"):
    load associated cell attributes (v0_,s0_ is_fixed) from the input file