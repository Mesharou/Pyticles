# Pyticles installation guide

Tested on ubuntu 20.04 with gcc >= 9

## Create Python environment for Pyticles

Pyticles uses Fortan routines compiled using f2py. Note that compilation is very sensitive to Numpy version. It is strongly recommended to use those provided in `Install/pyticles-xxx.yml` folder. However, depending on your architecture you may need to change either compilation options in `Modules/Makefile` and `Modules/Make_tools` or Numpy version.

Also user are encouraged to separate the environment they use for postprocessing from Pyticles tool. Adding new libraries to your Python environment may modify Numpy version and break Pyticles.

Before creating your own Pyticles environment you may check whether they are available on your local server like Styx and Datarmor.

### Micromamba

[Micromamba](https://mamba.readthedocs.io/en/latest/micromamba-installation.html) is faster and more efficient than conda when it comes to build Python environment. Once downloaded, edit `Install/mamba_create.sh`

- `mamba_exe`: path to micromamba binary
- `mamba_root_prefix` : path where to install Pyticles environment
- `pyticles_env` : Python version for Pyticles

Then run it

```Bash
cd Install
chmod +x mamba_create.sh
./mamba_create.sh
```

### Conda

For conda users the solution is similar. Once conda is installed you can create Pyticles from yaml file of your choice.

```Bash
cd Install
conda env create -f pyticles-xxx.yml
```

## Build Pyticles

When your Pyticles environment is ready Fortran routines can be built under
Pyticles environment. Note that R_tools for postrpocessing shall not be built
under Pyticles env as Pyticles core is very

- Activate Pyticles with **micromamba**

    ```Bash
    micromamba activate /path/to/pyticles/env
    ```

- Activate Pyticles with **conda**

    ```Bash
    source activate /path/to/pyticles/env/
    ```

- Build

    ```Bash
    cd Modules
    make pyticles
    ```

    After successful build you should have two python modules `pyticles_3d_sig_sa.cpython-xxx-arch.so`

### Troubleshooting

If compilation fails there may be some things to check

- Ensure that gcc version is >= 9

    ```Bash
    gcc --version
    ```

- Check Makefile `FCFLAGS` point to right path using the following shell command:

    ```Bash
    echo | gcc -E -Wp,-v -
    ```

    And edit Makefile accordingly

    ```Bash
    FCFLAGS = -I/path/to/gcc -c --f90flags="-extend_source -O1 -g -check all -CA -CB -CS"
    ```

- Verify that **Numpy version correspond to your code version**
  - Python-3.11 Numpy-1.26 for Pyticles-60bc00a
  - Python-3.8 for Pyticles-5fce563f141d907f336effebe9ccdc09caba760d

- If still stuck [create an issue](https://github.com/Mesharou/Pyticles/issues)

### Datarmor specific

- Pyticles 3-11 is already built at `/home2/datawork/jcollin/conda/envs/pyticles-3.11`

    Pyticles environment may be activated using the following commands in a computing node:

    ```Bash
    module purge 
    modue load conda/latest
    source activate `/home2/datawork/jcollin/conda/envs/pyticles-3.11`
    ```

- Build Pyticles

    ```Bash
    cd Install
    qsub datarmor_build_pyticles.pbs
    ```

- Create a user specific Python environment anyway using FTP queue

    ```Bash
    qsub -q ftp -I -l mem=10g -l walltime=00:30:00 -S /bin/bash
    ```

    Then you can download and install micromamba or use conda module with `module load conda/latest`
    and install interactively packages.

- Note that there is a toolbox for postprocessing 
  See [R_tools](./r-tools.md#datarmor) for detailed information
