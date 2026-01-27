Class Hierarchy (Parent -> Child)
FIREminimization -> Training -> Patterns 

Class structure:

1. training.Training(FIREminimization):

Training is a child class of FIREminimization. This makes sense, because we are viewing training a special type of minimization involving parameter changes (along with inertial relaxation)

In practice, this class will primarily be used to package some attributes and routines that are useful for any training protocol. 

Of course, by being a child class of FIREminimization, it further has access to the minimization routines.

Attributes and initialization::
    self._iter_counter = 0
    self._tolerance = 1e-8
    self._cost = None
    self._cost_values = None
    self._q_values = None
    self._initial_config = None

@classmethod
def from_sample(cls,tissue):
    inst = super().from_sample(tissue)
    return inst  
    
Setters:
set_iter_counter
set_tolerance
set_initial_config:
    i. Store initial config(PeriodicTissue or Spheroid)
    ii. Evaluate cell neighbors (for calculating Q_2...)

Functions:
+ write_cell_parameters():
    Writes the file cellParameters.input, with current state of cells from self._config.cells_ (if spheroid, skip the empty cells)
    (That is, parameters are to be written to a new file cellParameters.input from the info in cell objects from self._config.cells_)
+ edit_conf(**kwargs):
    First, this function reads the current conf file in self._dir.
    Then it rewrites the conf file with new parameters specified in kwargs
    Usage:
+ edit_spheroid(**kwargs)