from toolbox.training import Training

class pulseDrive(Training):
    def __init__(self):
        super().__init__()

    @classmethod
    def from_config(cls, config_dir, input_filename):
        inst = super().from_config(config_dir,input_filename)
        return inst
    
    # A single iteration, consisting of:
    # 1. writing current state of (driven) cell parameters
    # 2. Minimizing configuration (subject to current state of driven cells)
    # 3. Writing resulting configuration, bulk vtk and input cells vtk
    #   (labelled by iteration number)
    # 4. Incrementing the iteration counter (this is mainly helpful for resuming runs)
    def single_iteration(self):
        self.write_cell_parameters()
        self.minimize_config(FIRE_only = True)
        self.write_configuration(filename = "{}.bulk.txt".format(self._iter_counter))
        self._config.write_periodic_vtk(filename = "{}.bulk.vtk".format(self._iter_counter))
        self._config.write_cell_collection_vtk(
            cells_array = self._modified_cells,
            filename = "{}.input.vtk".format(self._iter_counter))
        self._iter_counter += 1
    
    # Pulse the s0 of driven cells between min_s0 and max_s0
    def run_s0_pulsing(self, min_s0 = 5, max_s0 = 5.6, iterations = 10):
        for iter in range(iterations):
            for cellID in self._modified_cells:
                cell = self._config.cells_[cellID]
                if iter%2:
                    cell.s0_ = max_s0
                else:
                    cell.s0_ = min_s0
            self.single_iteration()
    
