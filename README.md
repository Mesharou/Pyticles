# Pyticles

## What is it?


Pyticles is a Python/Fortran hybrid parallelized code for Lagrangian particles 3D advection from high-resolution ROMS/CROCO model data.



![image info](./Figures/2d_advection.png)


## Publications using the Pyticles code

- Gula, J., M.J. Molemaker & J.C. McWilliams, **2014** : Submesoscale cold filaments in the Gulf Stream, J. Phys. Oceanogr., 44, 2617-2643.
- Vic, C., G. Roullet, X. Carton, X. Capet, M.J. Molemaker & J. Gula, **2015** : Eddy-topography interactions and the fate of the Persian Gulf Outflow, J. Geophys. Res. Ocean., 120, doi :10.1002/2015JC011033.
- Gula, J., M.J. Molemaker & J.C. McWilliams, **2016** : Submesoscale dynamics of a Gulf Stream frontal eddy in the South Atlantic Bight, J. Phys. Oceanogr., 46, 305-325.
- Klymak, J.M., R.K. Shearman, J. Gula, C.M. Lee, E.A. D’Asaro, L. Thomas, R. Harcourt, A. Shcherbina, M.A. Sundermeyer, M.J. Molemaker & J.C. McWilliams, **2016** : Submesoscale streamers exchange water on the North Wall of the Gulf Stream, Geophys. Res. Lett., 43, 1226–1233.
- Ramanantsoa, J.D., M. Krug, P. Penven, M. Rouault & J. Gula, **2018** : Coastal upwelling south of Madagascar : temporal and spatial variability, Journal of Marine System, 178, 29-37.
- Ramanantsoa, J.D., P. Penven, M. Krug, J. Gula & M. Rouault, **2018** : Uncovering a new current : the South-west MAdagascar Coastal Current (SMACC), Geophys. Res. Lett., 45, 1930-1938, doi :10.1002/2017GL075900.
- Ragoasha, N., S. Herbette, G. Cambon, J. Veitch, C. Reason, C. Roy, 2019: Lagrangian pathways in the southern Benguela upwelling system, Journal of Marine Systems, 195, 50-66, https://doi.org/10.1016/j.jmarsys.2019.03.008.
- Gula, J., T. Blacic, & R.E. Todd, **2019** : Submesoscale coherent vortices in the Gulf Stream, Geophys. Res. Lett., 46, 2704-2714. https ://doi.org/10.1029/2019GL081919
- Smilenova, A., J. Gula, M. Le Corre, L. Houpert & Y. Reecht, **2020** : A persistent deep anticyclonic vortex in the Rockall Trough sustained by anticyclonic vortices shed from the slope current and wintertime convection, J. Geophys. Res. Ocean, https ://doi.org/10.1029/2019JC015905.
- Vic, C., S. Hascoet, J. Gula, T. Huck, & C. Maes, **2022** : Oceanic mesoscale cyclones cluster surface Lagrangian material, Geophys. Res. Lett, 49, e2021GL097488. https ://doi.org/10.1029/2021GL097488.
- Wang, L., J. Gula, J. Collin & L. Memery, **2022** : Effects of mesoscale dynamics on the distribution and transport of gravitationally sinking particles to the deep ocean, submitted.
- Qu, L., L. Thomas, A. Wienkers, R. Hetland, D. Kobashi, J. Taylor, F. Hsu, J. MacKinnon, R. Shearman, J. Nash, **2022** : Rapid Vertical Exchange at Fronts in the Northern Gulf of Mexico, https://doi.org/10.21203/rs.3.rs-1620471/v1
- Tagliabue, A., A. Lough, C. Vic, V. Roussenov, J. Gula, M. Lohan, J. Resing, & R. Williams, **2022** : Mechanisms driving the dispersal of hydrothermal iron from the northern Mid Atlantic Ridge, Geophys. Res. Lett., https://doi.org/10.1002/essoar.10512044.1

## DOI

The latest release has been archived on zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4973786.svg)](https://doi.org/10.5281/zenodo.4973786)



## License

[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg)](http://creativecommons.org/licenses/by-nc-sa/4.0/)

## Installation guide

### Requirements:
- netcdf4
- gfortran
- conda

### clone:
    
    git clone https://github.com/Mesharou/Pyticles.git

### Install
#### Classic:

    cd Install/
    conda env create -f conda_install.yml

#### DATARMOR specific:
tcsh users: 
    
    module load anaconda-py3.6/4.7.12

create new env
    
    conda create --name pyticles

    conda init tcsh
    conda activate pyticles

install modules

    conda install numpy
    conda install scipy
    conda install matplotlib
    conda install netCDF4

bash users:
    
    module load conda/latest
    module load NETCDF/4.3.3.1-mpt-intel2016

    cd Install/
    conda env create -f conda_install.yml

### Build modules
    cd Modules/
    make all

### Getting started
All Pyticles parameters and options are defined in Inputs/input\_file.py

You must define a simulation name and modify Modules/R\_files.py accordingly
to set the croco file's path 

To run Pyticles sequentially (interactively)
    
    python -i Pyticles.py 1 

In parallel with $np OpenMP procs
    
    python Pyticles.py $np > out


