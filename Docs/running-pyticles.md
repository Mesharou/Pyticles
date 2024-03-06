# How to run Pyticles

- [Set up your CROCO simulation](#set-up-your-croco-simulation)

- [Set up Pyticles experiment](#setting-up-pyticles-experiment)

- [Run Pyticles](#run-pyticles)

## Set up your CROCO simulation

The first step before running a Pyticles experiment is to specify CROCO files path and metadata using the Class `files` from `Modules/R_files.py`. Pyticles will then be able to load CROCO variables from Netcdf files with the appropriate file name and time-step.

- Create a simulation name. For instance `my_croco_simul`.

  **NOTE** Ensure that this simulation name is not already used in R_files.py

- Edit `Modules/R_files.py` in the Class `files` *i.e* after `if polgyr in simul:`
  add the following block:

    ```Python
    elif 'my_croco_simul' in simul:
      self.realyear = True
      self.realyear_origin = datetime(1979, 1, 1)
      self.model = 'croco'
      self.digits = 5
      
      folder = '/path/to/croco_root_folder/'
      self.grd = folder + 'croco_grd.nc'
      self.his = folder + 'croco_his.'
      self.tfile = 24
      self.tstart = 0
      self.tend = 10000
    ```

  Here all variables are designed to retrieve CROCO file names
  
  - **self.realyear**: True for interannual files, False for climatology
  
  - **self.realyear_origin**: corresponding to first CROCO scrum_time

  - **self.model**:
    - `ucla`: UCLA ROMS model 
    - `croco` (format **croco_his_0200.nc**)
    - `croco_xios` (xios format **croco_his_2012-10-24-2012-10-28.nc**)
    - `croco_date` (format **roms_his_Y20M1.nc**)
  
    If your files have a different format, you can either change their names fo fit (using symbolic link for instance) or modify R_files.py by adding a `model`.

  - **self.digits:** The number of digits used to name CROCO files with `self.model = croco`.

    example for time step 200 of croco outputs:
    - with `self.digits = 4`,  **croco_his_0200.nc**
    - with `self.digits = 5`,  **croco_his_00200.nc**

  - **folder**: Path to croco_files

  - **self.his**: prefix name for history files

  - **self.tfile**: number of output time steps in a CROCO file

  - **self.tstart**: First CROCO time step to take for Pyticles Simulation. (Usually **0**). We will see later that there are other variables chose CROCO output

  - **self.tend**: Last Pyticles time step. just take a large integer. As for tstart, initial and final steps of Pyticles are determined in `Input/input_file.py`

  **TODO**
  - self.real_year_tstart
  - add this doc to class files

Once your simulation is set in R_files.py it is safe to test it first in a Jupyter notebook.

To initialize a simulation we need:

- A simulation name: `my_croco_simul`
- (Optionally) coordinates to zoom: [10500,11300,4750,5550,[1,300,1]]
- (Optionally) a time to load:  0 (here)

```Python
    # --> project related modules
    sys.path.append('pyticles/Modules/')
    from R_files import load

    # --> load my_simul
    parameters = "[10500,11300,4750,5550,[1,300,1]] 0"
    simul = load(simul=parameters, output=False)
```

See `path/to/notebooks.ipynb` for practical case

## Setting up Pyticles experiment

All Pyticles options are defined in `Inputs/input_file.py`.
There are different sections, due to some constraint some variables are defined in sections they do not belong to. Note that some sections are hardcoded and shall not be edited by users.

For more detail information see directly in `Inputs/input_file.py`.

- CROCO inputs: metadata about CROCO files
- Particles dynamics: numerical schemes for Pyticles, output frequency and CROCO fields
- Pyticles Outputs: path for netcdf output file, variables to save, figures
- Particles seeding: Specify seeding frequency and location

Here we present only some important parameters.

- **dfile:** Pyticles output time step in units of CROCO outputs time step.
  - dfile = 1, Pyticles output frequency = CROCO output frequency
  - dfile = 1/2, Pyticles output frequency = 2 x CROCO output frequency
  - dfile = -1, backward experiment

- Integers, **start_file** and **end_file**: CROCO output time step for initial and
   last Pyticles time step.

- (optional) boolean, **restart**: If True, restart a Pyticles simulation (typically to re-run a simulation that crashed)

- (optional) integer > 0, **restart_time**: number of CROCO time steps since initialization
  used to restart. Note that if simulation backward the sign of dfile is automatically applied.

- (optional) **itime_trap** and **trap_file**: use to initialize a Pyticles simulation from a former Pyticles simulation. Very specific

- string, **my_simul**: name of CROCO simulation specified in R_files.py to retrieve CROCO files

In terms of particles dynamics there are many possibilities.
Particles can be advected in 3D, both passively or with a sedimentation velocity.
2D advection is also an option either at surface or at constant depth.

Numerically, the default option is a bi-linear in space linear in time
interpolation using RK4 time stepping. A small time step is used by Pyticles to
respect CFL condition. Either statically by setting up maximum horizontal and 
vertical maximum velocities

```Python
    umax = 2
    vmax = 2
    wmax = 2*1e-3
```

Or dynamically using `ìnline_cfl = True` (much slower).

- Time stepping options are

  ```
  timestep = 'RK4' # Choices are
            # FE (forward-Euler)
            # RK2, RK4 (Runge-Kutta 2nd and 4th order)
            # AB2, AB3, AB4 (Adams-Bashforth 2,3,4th order)
            # ABM4 (Adams-Bashforth 4th order + Adams-Moulton corrector).
  ```

- Spatial interpolation: Default is linear, Available :
  - `#define CUBIC_INTERPOLATION`
  - `#define CRSPL_INTERPOLATION`
  - `#define WENO_INTERPOLATION`

  Beware these higher order schemes have not been rigorously tested
  To define them, in Modules/interp_3d_for_pyticles.F.

  **To Activate them Activate ccp keys : NEW_VERSION and chosen numerical scheme
  Compile cpp keys use make command**

  **TODO**
  Clarifier les clés CCP 
  - advdepth (can do anything )
  - advsurf (read only 2D field straight from netcdf file)
  routines dans update_largemem_xyz
  renommer les clés CPP
  
## Run Pyticles

A feature of Pyticles is the ability to run OMP parallel.
To run the code you only need to 

- Activate your Python environment Pyticles

- Ensure that pyticles Fortran routines have been built with the same env. 
[see install for more detail](install.md)

- Use one of the following commands

### Debug mode

  ```Bash
  python -i Pyticles.py  1 
  ```

### Production, using 4 processors

```Bash
python Pyticles.py 4 &> my_run.out
```

### On a cluster using PBS

```
qsub run_Pyticles.pbs
```

### TIPS

- To keep track of your experiments it can be nice to copy `input_file.py`
  to a file `my_experiment.py` edit the latter and force copy before running
  Pyticles

```Bash
cp -f Inputs/my_experiment.py Inputs/input_file.py
```
