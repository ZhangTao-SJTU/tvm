To Do List:

Broad Goal:

FIRE minimize the configuration.


Cookbook of changes (in the order in which I addressed them):

1. FIRE minimization setup: 

    For now lets FIRE minimize in terms of

        vertex -> volumeforce_
        vertex -> interfaceForce_


    The first two are calculated with volume->updateforces(), interface->updateforces() respectively. If replacement functions are written, we can just replace the respective function calls.

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

    To facilitate this, we introduce to the run class the following attributes:
    
        double FIRE_ff;
        double FIRE_fv; //power
        double FIRE_vv;

    and the function updates the value of these:

        int FIREupdateForceVelocityProjections();


