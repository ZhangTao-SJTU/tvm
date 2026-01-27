Class Hierarchy (Parent -> Child)
FIREminimization -> Training -> Patterns 

Class structure:

1. patterns.Patterns(Training):

Attributes and initialization:
    self._target_cell_to_stress = None
    self._frozen_cells = []
    self._clamping_max_iters = 10
    self._clamping_s0_lower_limit = 4
    self._clamping_s0_upper_limit = 6
    self._clamping_correction_factor = 0.1
    self._clamping_FIRE_only = True
    self._clamping_tolerance = 1e-15
    self._learning_rate = 10
    self._clear_interval = 100

@classmethod:
from_sample(tissue)

Setters:
set_clamping_FIRE_only(bool)
set_clamping_max_iters():
set_clamping_correction_factor():
set_target_cell_to_stress({cellID:stress}):
    (i) Sets self._target_cell_to_stress
    (ii) Produce vtk of target cells
set_frozen_cells()
set_learning_rate()
set_clear_interval()

Functions:

calculate_max_shear_stresses():
    store target cell stresses in cell.max_shear_stress

evaluate_cost():
    return average of absolute percentage difference of current cell stresses wrt target cell stresses

initialize():
    (i) Write cell parameters and check minimization
    (ii) store minimized.txt as 0000000.bulk.txt. store initial cell parameters as 0000000.cellParameters.input. 
    (iii) set_initial_config() (function in training class)
    (iv) self._cost_values = [initial_cost]
        self._q_values = [1]
        save to files
    (v) save initial stresses in 0000000.stresses.csv
    (vi) self._iter_counter += 1

solve_cell_s0_for_target_stress(cellID,target_stress):
    ONLY SIDE EFFECT: target cell.s0 gets target stress (in usage this will be slightly overdriven from actual target stress for the cell i.e. from self._target_cell_to_stress)


clamp_target_cells():
    This function clamps the self._target_cells to a desired value of max shear stress. This is achieved by successively changing the s0 for these cells, minimizing the configuration and checking the resulting target stress. These iterations are terminated after success to tolerance or when max iterations are reached.

    -> it is more effective to set max iterations low, because the clamping effects can be completely erased by reconnection.

    Ideas implemented:
    1. For each target cell, the target stress at a given clamping iteration should be an overshooting/undershooting of the final target_stress.

    However, we should limit the overshooting within a reasonable amount.
    
        if (self._target_stress-cell.max_shear_stress)>0:
        target stress should be greater than self.target_stress, and the amount should be reflective of large this difference is.
        else: target stress should be smaller
    
    Algorithm:
    Run clamping iterations upto convergence or max iters:
        i. Evaluate which cells need to be clamped (if none, clamping has converged!)
        ii. Record pre-clamping s0
        iii. On each clampable target cell:
            Solve s0 for (over)clamped target stress
        iv. Write cell parameters, now with clamped s0 for target cells. If there is no change in s0 from preclamped values, clamping is no longer useful -> exit clamping
        v. Minimize configuration (FIRE only).

single_iteration():
    i. Record free state areas of hidden cells
    ii. Clamp. At the end of clamping one should have new target cell s0s in cellParameters.input
    iii. Update hidden cell s0s
    iv. Unclamp: reset target cell s0, write cell parameters, minimize config (with overdamping!)
    v. Write configuration (bulk.txt), copy and store unclamped s0 (cellParameters.input)

run_to_max_iters():
    evaluate cost
    for max iters:
        if cost<tolerance: break
        i. single_iteration()
        ii. evaluate cost
        iii. write cost, qval, stresses
        iv. if self._iter_counter % self._clear_interval == 0:
                self.clear_directory()
        v. self._iter_counter ++

clear_directory()
    Move {:07d iter_counter}.bulk.txt, ... to self._dir/files
    Delete the rest of the files


set_random_target_cells(self, n_cells = 1, target_stress = 1, **kwargs)
    Kwargs:
        "stress_limits"
        "exclude_cells"
