# How to run Pyticles

## Set up your CROCO simulation

The first step before running a Pyticles experiment is to specifiy CROCO files path and metadata using the class load from `Modules/R_files.py`. Pyticles will then be able to load CROCO variables from netcdf files with the appropriate file name and timestep. 

- Create a simualtion name. For instance `my_croco_simul`.

  **NOTE** Ensure that this simulation name does not exist

- Edit `Modules/R_files.py` in the class files section (ie **after line 992**)

  add the following block

    ```Python 
    elif 'my_croco_simul' in simul:
            self.realyear = True
            self.realyear_origin = datetime(1979,1,1)
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
  
  - **self.realyear**: if True takes into account real calendar (see CROCO CPP keys)
  
  - **self.realyear_origin**: corresponding to first CROCO srum_time

  - **self.model**:
    - `croco` (format **croco_his_0200.nc**)
    - `croco_xios` (xios format **croco_his_2012-10-24-2012-10-28.nc**)
    
    For interannual AGRIF files like **croco_his_Y2005M03.nc** Maybe check `croco_agrif_jc` or `croco_lionel` model options (to be tested).

    If your files have a different format, you can either change their names fo fit (using symbolic link for instance) or modify R_files.py by adding a `model`.

  - **self.digits:** The number of digits used to name CROCO files with `self.model = croco`.

    example for time step 200 of croco outputs:
    - with `self.digits = 4`,  **croco_his_0200.nc**
    - with `self.digits = 5`,  **croco_his_00200.nc**

  - **folder**: Path to croco_files.

  - **self.his**: prefix name for history files.

  - **self.tfile**: number of output time steps in a CROCO file.

  - **self.tstart**: First CROCO time step to take for Pyticles Simuation. (Usually **0**). We will see later that there are other variables chose CROCO output 

  - **self.tend**: Last Pyticles time step. just take a large integer. As for tstart, initial and final steps of Pyticles are determined in `Input/input_file.py`

  **TODO**
  - self.real_year_tsart
  - add this doc to class files

Once your simulation is set in R_files.py it is safe to test it first in a Jupyter notebook.

To initialize a simulation we need:
- a simulation name: my_croco_simul
- (optionnally) coordinates to zoom: [10500,11300,4750,5550,[1,300,1]]
- (optionnaly) a time to load:  0 (here)

``` Python
    # --> project related modules
    sys.path.append('pyticles/Modules/')
    from R_files import load

    # --> load my_simul
    parameters = " [10500,11300,4750,5550,[1,300,1]] 0"
    simul = load(simul=parameters, output=False)

```