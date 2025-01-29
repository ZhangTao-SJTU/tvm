Broad Goal:

- To make v0 and s0 individual cell properties
- To make these input quantities at initialization.

--- Initial survey ---

Forces and energies are updated by routines in Volume.cpp and Interface.cpp.

1. In tvm.cpp, the new run->volume_ and run->interface_ objects are initialized

1. Further, in tvm.cpp, LoadConf() loads s0_ into run->interface_

    In other words, double s0_ in the scope of Interface.cpp is actually a global property of the tissue.

    (Minor update: remove the default initialization s0_ = 5.4 from Interface.cpp)

1. v0_ = 1 is hard coded into Energy/Volume.cpp, in the following functions:

        Volume::updatePressure()
        Volume::updateEnergy()

    (in fact,there is no v0_ variable, which we will introduce on a cell level)

    Volume::updatePressure() updates Cell::pressure_
    for every cell in the run.

    Volume::updateEnergy() only calculate value for Volume::Energy_

1. A summary of the algorithm in Volume::updateForces():
    - vertex->volumeForce_ = [0,0,0] for all vertices
    - Volume::updateVolume() 
    - Volume::updatePressure()
    - For each cell, and for each polygon of the cell, run Volume::updatePolygonForces(cell,polygon). This function distributes the cell pressure to vertex->volumeForce_ of the vertex.

    (minor updates: using list comprehension where helpful in Volume.cpp)


    --- Hence, introducing cell-level v0 in Volume::updatePressure() is enough. ---

1. A summary of the algorithm in Interface::updateForces():

    - vertex->interfaceForce_ = [0,0,0] for all vertices
    - for all polygons; polygon->updateArea();
    - Interface::updateTension()

        This function updates the attribute double polygon->tension_ for each polygon.

        First, this is set to 0 for each polygon. Then, for each cell, a tension of 2*(s-s0) is added to the tension of each polygon

    - for each polygon in the run, Interface::updatePolygonForces()
    Similar to Volume::updatePolygonForces(), this redistributes net polygon tension to vertex->interfaceForce.

    (minor updates: using list comprehension where helpful in Interface.cpp)

    --- in short, introducing cell-level s0 in Interface::updateTension() is enough. ---

/|\|/|/|\|/|/|\|/|/|\|/|/|\|/|/|\|/|/|\|/

--- Changes Log ---

1. Introduce double v0_ and s0_ as Cell attributes (in Cell.h and Cell.cpp)
    (why double? in the code v and s are calculated as doubles. just avoiding unexpected behavior when sutractracting float from double in the force/energy calculations.)
1. Set v0_ = 1 as a default initialization.
1. In tvm.cpp, make loadconf() initialize the s0 value for each cell 
    
    (instead of initializing run_->interface_->s0_)
1. In Interface.h and Interface.cpp:
    - Remove attribute Interface::s0_
    - Edit Interface::updateTension()
    - Edit Interface::updateEnergy()
1. In Volume.h and Volume.cpp
    - Edit Volume::updatePressure()
    - Edit Volume::updateEnergy()


    So far, we have only refactored the code to accept v0_ and s0_ as individual cell properties. The values for each cell are initialized as follows:

    - v0_ = 1 is initialized at the constructor for every new cell.
    - s0_ is set in tvm.cpp with LoadConf(), which sets 
    the the s0 value for every cell to the one specified in conf 

    In order to input these cell level values for some (or even all) cells we can take a similar approach as in the fixed cells.

1. Designing the input file

    We will have a file called cellParameters.input. Contents to be formatted as follows

        <cellID>    <v0>    <s0>    <int(is_fixed_)>
    
    For example

        7   0.9 5.4 1

    meaning, set cell with cellID 7 to have v0 = 0.9; s0_ = 5.4; is_fixed = True
    - In tvm.cpp this file should be read AFTER loadConf(). Otherwise loadConf will overwrite any new s0 values
    - s0 is a target surface area (and NOT a target shape index!) when changing v0 values from 1, please make reasonable adjustments for s0.
    
    Based on these considerations, we introduce a function loadCellProperties() in tvm.cpp

    - Notice that this function basically obsoletes InitializeFixed()
    - Default behavior: if cellParameters.input is not found or cannot be opened, skip fixing any cell parameters and minimize homogeneous configuration.

/|\|/|/|\|/|/|\|/|/|\|/|/|\|/|/|\|/|/|\|/

This concludes our introduction of individually addressable cell properties
    - v0
    - s0
    - whether the cell is "fixed" (i.e vertices of the cell do not move, though they may be changed by reconnection)
to the minimization, via a file "cellProperties.input" in the run folder.

So far, the format of this file is:

        <cellID>    <v0>    <s0>    <int(is_fixed_)>
    
but for enhanced readability we can introduce headers in the future or make a cascaded format
(so that all quantities do not need to be specified)


    
Next steps: command line arguments to avoid overdamped stage (because this introduces noise)
