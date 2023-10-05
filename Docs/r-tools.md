# R_tools

R_tools is an optional toolbox to load CROCO variables from Netcdf files and
compute others that are missing. It is quite fast as it is written in
Fortran but it comes with the cost of building a new environment dedicated to
postprocessing and to build the toolbox.

`R_vars.py`` is a higher level module using R_tools that makes life easier when it
comes to access CROCO variables.

## Variables accessible using R_vars

    'temp': ['Temperature', r'$T\,(^{\circ}C)$', [0,0,1]],\
    'salt': ['Salinity', 'PSU', [0,0,1]],\
    'u': ['u', 'm/s', [1,0,1]],\
    'v': ['v', 'm/s', [0,1,1]],\
    'ubar': ['ubar', 'm/s', [1,0,-1]],\
    'vbar': ['vbar', 'm/s', [0,1,-1]],\
    'zeta': ['SSH', r'$\eta\,(m)$' ,[0,0,-1]],\
    'hbls': ['Thickness of KPP surface boundary layer', 'm', [0,0,-1]],\
    'hbbls': ['Thickness of KPP bottom boundary layer', 'm', [0,0,-1]],\
    'hbls_rho': ['Surface mixed-layer (based on rho = rhos+ 0.03)', 'm', [0,0,-1]],\
    'hbls_t': ['Surface mixed-layer (based on t = ts - 0.2)', 'm', [0,0,-1]],\
    'AKt': ['Temperature vertical diffusion coef', 'm2/s', [0,0,0]],\
    'AKv': ['Momentum vertical diffusion coef', 'm2/s', [0,0,0]],\
    'omega': ['S-coordinate vertical velocity', 'm/s ?', [0,0,0]],\
    \
    'psi': ['psi', 'Streamfunction for depth integrated flow' ,[1,1,-1]],\
    'psi_surf': ['psi_surf', 'Streamfunction of surface flow' ,[1,1,-1]],\
    \
    'rho': ['in-situ density', 'kg.m-3', [0,0,1]],\
    'rho1': ['in-situ density standard pressure', 'kg.m-3', [0,0,1]],\
    'rhop': ['potential density', 'kg.m-3', [0,0,1]],\
    'bvf': ['Brunt-Vaisala Frequency squared: N2', 's-2', [0,0,0]],\
    'buoy': ['buoyancy', 'm/s-2', [0,0,1]],\
    \
    'w': ['Vertical velocity', r'$w\,(m\,s^{-1})$', [0,0,1]],\
    \
    'absvrt': ['Absolute Vorticity', 's-1' ,[1,1,1]],\
    'vrt': ['Relative Vorticity', r'$\frac{\zeta}{f}$' ,[1,1,1]],\
    'pv': ['Potential Vorticity', 'PVU' ,[1,1,0]],\
    'pvr': ['Potential Vorticity on rho levels', 'PVU' ,[1,1,1]]

## Files description

- `R_vars.py`: top-level module calling `R_tools` and `R_tools_fort`
- `R_tools.py`: upper level module API for Fortran routines
- `R_tools_fort.F`: main Fortan code that we be pre-compiled in C to include Fortran routines
- `R_tools_fort_routines/`: folder with all Fortran routines
- `R_tools_fort.f`: pre-compiled Fortran script (intermediate file after make receipe)
- `R_tools_fort_xxxxx.so`: compiled C module with f2py (imported in `R_tools.py`)

## Python environment for R_tools

You can create your own Python environment but it is not guaranteed that build
will work. Some `yaml` files that have been tested are present in `Install/pyticles-pp-xxx.yaml`. 

```YAML
name: pyticles-pp-3.11
channels:
  - conda-forge
dependencies:
  - python=3.11
  - numpy
  - scipy
  - matplotlib
  - netcdf4
  - basemap
  - cartopy
  - xarray
  - dask
  - ipython
  - jupyter
```

To create postprocessing tools environment using
 [conda](https://docs.conda.io/en/latest/) or
[micromamba](https://mamba.readthedocs.io/en/latest/micromamba-installation.html)

- Conda

    ```Bash
    conda env create -f pyticles-pp-xxx.yaml
    ```

- Micromamba

    ```Bash
    micromamba env create -f pyticles-pp-xxx.yaml
    ```

## Build R_tools

- Conda

    ```Bash
    source activate /path/to/postproc/env
    cd Modules
    make R_tools
    ```

- Micromamba

    ```Bash
    micromamba activate /path/to/postproc/env
    cd Modules
    make R_tools
    ```

To Check whether the build actually worked you shall see a `R_tools_fort_xxx.so` file in `Modules/`

[Troubleshooting](install.md#troubleshooting)

## Datarmor

- A Pyticles postprocessing environment is pre-built at
 `/home2/datawork/jcollin/conda/envs/pyticles-pp-3.11`.
  Additionnal packages are:

  - cartopy
  - xarray
  - dask
  - ipython
  - jupyter

- Build R_tools

  ```Bash
  cd Install
  qsub datarmor_build_tools.pbs
  ```

  A new file with python version name `Module/R_tools_fort.cpython-$python-version-x86_64-linux-gnu.so` should appear.
  If not check `datarmor_build.pbs.oxxxx` for errors.

## Tutorial

See [notebooks](../Postprocessing/tutorial/croco-variable-manipulation.ipynb)