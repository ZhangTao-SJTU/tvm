To Do List:

Broad Goal:

FIRE minimize the configuration.


Cookbook of changes (in the order in which I addressed them):

1. FIRE minimization setup: 

    For now lets FIRE minimize in terms of

        vertex -> volumeforce_
        vertex -> interfaceForce_


    These are calculated with volume->updateforces(), interface->updateforces() respectively. If replacement functions are written, we can just replace the respective function calls.

    Let's do this by introducing a new function in run:

        run -> FIREminimize()

    This function does NOT use t_init,t_f and dt from conf (i.e run->dt_). Instead, 
    
        FIRE_itermax
        FIRE_dt
        FIRE_dtmax

    are parameters in the scope of this function.

    Here are the steps in FIREminimize():

    +   Log, dump, and reconnect by iteration multiples instead of simulation_time_.
    However, still maintain the run->simulation_time_ (which increases by dt_ every iteration) in order to use the dump functions. 

    + Write FIRE versions of the velocities and positions functions:

            FIREupdateVerticesPosition(double& FIRE_dt)
            FIREupdateVerticesVelocity(double& FIRE_dt)


    Reasons:

    i. No overdamped motion

    ii. FIRE_dt is a variable in the scope of Run* run
    
1. FIRE minimization implementation:

    We need a dot product function that computes f.f, f.v and v.v. it is okay to calculate them at the same time, and just after velocities are updated.

    To facilitate this, we introduce in the scope of the Run object the following attributes:
    
        double FIRE_ff;
        double FIRE_fv; //power
        double FIRE_vv;

    and the function updates the value of these:

        int FIREupdateForceVelocityProjections();

1. Implement a COM polygon center. Rewrite the following function in Polygon.cpp:
        
        Polygon::updateCenter()

    Note: this routine is only called in Run::updateGeoinfo()

        Run::updateGeoinfo() 
    is called in every iteration of both these functions
        
        Run::overdampedMotion() 
        Run::FIREminimized()
    
    Reason: the resulting exact forces are MUCH easier to code.

1. Exact volume forces: rewrite (in Energy/Volume.cpp) the function:

        Volume::updateForces()
    
    Algorithm:

    Upon calling this function, vertex->volumeForce_ is first set to {0,0,0} for each vertex

    The routine the iterates over the cells. Each iteration adds the contribution from that
    cell to the vertex->volumeForce_.

    The routine requires the anticlockwise orientation of every "list" below:

        for (auto cell : run->cells_){
            for (auto polygon: cell->polygons_){
                list = polygon->vertices_;
            }
        }
    
    We note that the routine Cell::updatePolygonDirections() updates <long int,bool>cell->polygonDirections_().
    
    We further note that this polygon directions routine is called (by way of Volume::updatePolygonDirections) during:

    (a) Initialization in tvm.cpp
    (b) As part of processing reconnection events

    and that it includes a built in volume update for the cells.
    
    Implementation:

    We will use Cell::updateVolume() for inspiration since it uses similar concepts.
