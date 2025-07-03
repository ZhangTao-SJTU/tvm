1. minimization.FIREminimization:

Attributes:
    self._config = None
        For now, supports objects of type periodicTissue.PeriodicTissue
        Results of the minimization is reloaded into this object.
    self._dir = None
    self._modified_cells = []
        Array of cellIDs of "modified cells".
        Modified cells can have a different
            cell.v0_
            cell.s0_
            cell.is_fixed_
        than the rest of the tissue.

@classmethods:
periodic_tissue(PeriodicTissue)

Routines:
minimize_config(FIRE_only = False):

Minimization concludes with loading the new attributes into the tissue, and also loading cellParameters.input for any modified cells during the minimization.

write_configuration(filename = "sample.topo"):

load_cell_parameters():
    updates self._modified_cells and associated cell attributes (v0_,s0_ is_fixed) directly from the file cellParameters.input.
    (If cellParameters doesnt exist: set self._modified_cells to None and return.)

1. training.Training(FIREminimization):

Training is a child class of FIREminimization. This makes sense, because we are viewing training a special type of minimization involving parameter changes (along with inertial relaxation)

In practice, this class will primarily be used to package some attributes and routines that are useful for any training protocol. 

Of course, by being a child class of FIREminimization, it further has access to the minimization routines.

Attributes:
    self._iter_counter = 0
    self._learning_rate = 5e-1
    self._lambda = 1.1
    self._tolerance = 1e-8
    self._cost = None
    self._cost_values = None
@classmethods:
periodic_tissue(tissue:periodic.PeriodicTissue)

Setters:
set_iter_counter
set_lambda
set_learning_rate
set_tolerance

Functions:

+ write_cell_parameters():

    Writes the file cellParameters.input, with current state of cells from self._modified_cells
    (That is, parameters are to be written to a new file cellParameters.input from the info in cell objects from self._config.cells_)

    This is a part of training because modified cells can be the hidden cells (their parameters change in the training process)


+ pick_random_modified_cell(v0 = 1, s0 = 5.2, is_fixed = False):

    pick a cell; cellID appended to self._modified_cells
    Criteria:
    1. No cross boundary cells
    2. No modified cells should share polygons
    3. Should not be too close to the boundary (for visualization).

+ edit_conf(**kwargs):

    First, this function reads the current conf file in self._dir.

    Then it rewrites the conf file with new parameters specified in kwargs

    Usage:

1. patterns.Patterns(Training):

classmethod:

periodic_tissue(tissue:periodic.PeriodicTissue)


Attributes:
    self._target_cells = None
    Target cells are the input/output cells (as opposed to hidden cells). 
    
    At parts of the iteration, these may or may not also be "modified cells" (i.e. in the list FIREminimization._modified cells, meaning they are in the cellParameters.input)
    self._target_stress = None

Setters:
set_target_stress(stress:double):

set_target_cells_spheroid(cellIDs:list):
    A tool for selecting target cells on a spherical shell
    For convenience a target_cells.vtk is also outputted.
    
    ##
    Edits: make target cells not share any polygons.
    ##

set_central_target_cells(n_cells = 1):
    sets single or multiple target cells near the center (for better visibility)
Functions:

calculate_max_shear_stresses():
    For each cellID in self._target_cells, this functions quips self._config.cells_[cellID].max_shear_stress_ with the max shear stress.

evaluate_cost():
clamp_target_cells(tol = 1e-3, max_iters = 10):
    This function clamps the self._target_cells to a desired value of max shear stress. This is achieved by successively changing the s0 for these cells, minimizing the configuration and checking the resulting target stress. These iterations are terminated after success to tolerance or when max iterations are reached.

    Ideas implemented:
    1. For each target cell, the target stress at a given clamping iteration should be an overshooting/undershooting of the final target_stress.

    However, we should limit the overshooting within a reasonable amount.
    
        if (self._target_stress-cell.max_shear_stress)>0:
        target stress should be greater than self.target_stress, and the amount should be reflective of large this difference is.
        else: target stress should be smaller

        implementation: 

single_iteration():
run():
    





