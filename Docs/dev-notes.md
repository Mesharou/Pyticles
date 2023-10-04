# Development notes

## Tutorial feature

### Discussion

@jeremy:

- Installer Pyticles
  - build les tools (pyticles (part.F + tools R_vars etc...)
  - format markdown dans un dossier docs
  - question d'ajouter des outils de postprocessing au projet comme Rvars ? R_tools_gigatl ?

- Mettre en place sa configuration
  - CROCO: R_files.py
  - particules: input_file.py
  - format jupyter notebook

- Analyser et visualiser les sorties de Pyticles
  - méthode simple avec R_tools
  - méthode avancée avec xarray
  - format jupyter notebook
  - quelles diagnostiques type ?

@jonathan:

- Ce serait bien d'avoir des une doc un peu plus propre
pour l'installation
- quelques exemples de notebooks pour le post-traitement et la visualisation.
- Après effectivement je ne sais pas s'il vaut mieux développer les outils R_tools.py ou xarray.

### Code structure and methodology

- Create `branch` for `new-feature` - merge party to master - tag and release notes ?
- How to deal with Jupyter notebooks (clear outputs)
- Documentation in `Docs/xxx.md`
- Dev notes in `Docs/dev-notes.md` used to maintain release notes
- Tutorials in `?`: `Tutorials/*ipynb` or `Postprocessing/*ipynb`

### Install tools

- Do we keep separate environment for postprocessing and Pyticles or not ?
  - ipython
  - xarray
  - dask
  - cartopy
  - jupyter
  - "+" pyticles librairies => provide yaml + pre-installed env (datarmor/styx/lops ?)

- Do we merge Makefile and Make_tools
- Do we need to add CFLAGS to R_tools as well ?
- Include path is not consistent in Makefiles but compilation works anyway. Is this a problem to consider ?
- compile with gcc-9 ok ?

### postprocessing tools

- R_tools files were excluded from project. Check with Jonathan if we keep them or not.

  - Modules/Make_tools
  - Modules/R_smooth.py
  - Modules/R_tools.py
  - Modules/R_tools_fort.F
  - Modules/R_tools/
  - Modules/R_vars.py

- **Why use R_vars and not R_vars_gula ?**

- Issue with R_files.py gigatl1_1h_surf (2 days shift (rings a bell))


### Postprocessing and visualization tutorials

What to add exactly ?

- section of CROCO variables + particles position (ZY, XY)
- particles trajectory ? more specific ?
- interpolation of CROCO variables onto particles position (already there no ?)
- 