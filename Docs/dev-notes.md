# Development notes

## Tutorial feature

### FIXME
- Issue with R_files.py gigatl1_1h_surf (2 days shift (rings a bell))
il y a une varialble shift dans R_files.py (cf swot R_files si besoin)

    files
    __init__()
    self.init = 0

    if 'tides' not in simul:
        self.shift = 2 * 24

- clarifier les clés CPP de adv3d etc...
  à noter iso = isopycne par rapport à rho0
  détailler des cas simples

- rename keys `agrif_jc` and `croco_lionel`
  - agrif_jc: agrif ? interannuel ?
  - croco_lionel ? maybe this could come with an example in R_files

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
- How to deal with Jupyter notebooks: clear outputs
  
  To clear jupyter notebook outputs, create a file `.git/hooks/pre-commit`
  containing the following code:
  ```Bash
  #!/bin/bash
  for f in $(git diff --name-only --cached); do
      if [[ $f == *.ipynb ]]; then
          jupyter nbconvert --clear-output --inplace $f
          git add $f
      fi
  done

  if git diff --name-only --cached --exit-code
  then
      echo "No changes detected after removing notebook output"
      exit 1
  fi
  ```
  pre-commit file must be executable `chmod +x .git/pre-commit`

  Also in order to function, **jupyter shall be installed in development
  Python env**.


- Documentation in `Docs/xxx.md`
- Dev notes in `Docs/dev-notes.md` used to maintain release notes
- Tutorials in `?`: `Tutorials/*ipynb`

### Install tools

- Created separate environment for postprocessing and Pyticles.
Postprocessing environment has the following additional packages:
  - ipython
  - xarray
  - dask
  - cartopy
  - jupyter

  Environment files are in `Ìnstall/*.yaml`. Postprocessing and Pyticles are
  already available on Datarmor. 

### postprocessing tools

- R_tools files were excluded from project. Now they are back.

  - Modules/Make_tools
  - Modules/R_smooth.py
  - Modules/R_tools.py
  - Modules/R_tools_fort.F
  - Modules/R_tools/
  - Modules/R_vars.py

  Note this not the development branch that is used (no there are some 
  important missing variables: Okubo-weiss etc...)

### Postprocessing and visualization tutorials

What to add exactly ?

- variable 2D + position des particules
- démo interpolation var croco interp part => plot ou netcdf
- section verticale
- example xarray (pas xgcm)

## Documentation


