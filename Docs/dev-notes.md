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

### Install tools

- Include path is not consistent in Makefiles but compilation works anyway. Is this a problem to consider ?

- Do we keep separate environment for postprocessing and Pyticles or not ?

- Do we merge Makefile and Make_tools

- Do we need to add CFLAGS to R_tools as well ?

### postprocessing tools

- R_tools files were excluded from project. Check with Jonathan if we keep them or not.

  - Modules/Make_tools
  - Modules/R_smooth.py
  - Modules/R_tools.py
  - Modules/R_tools_fort.F
  - Modules/R_tools/
  - Modules/R_vars.py


- Issue with R_files.py gigatl1_1h_surf (2 days shift (rings a bell))
