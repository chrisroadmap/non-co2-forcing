# non-co2-forcing

Calculate the split of CO2 and non-CO2 forcing from Earth System Models.

## installation

### requirements
- `anaconda` for `python3`
- `python>=3.7`
- for the Cummins calibration, `R>=4.4.0` and `cmake>=3.2`

### set up environments for Python and R

First, create the `conda` environment for `python`. You will require an installation of `anaconda`. We recommend [miniconda](https://docs.anaconda.com/miniconda/) or [miniforge](https://github.com/conda-forge/miniforge). After installing `conda`, run the following commands

 
```
conda env create -f environment.yml
conda activate fair-calibrate
```

As of around June 2024, the `R` environments no longer install properly from `conda`, so we'll need to do this manually.

First, you'll require a version of `R` at least as new as 4.1.1. Grab this from https://www.r-project.org/. Second, you'll need Donald Cummins' `EBM` package binary, which can be downloaded from [here](https://github.com/donaldcummins/EBM/archive/refs/tags/v1.1.0.tar.gz). Download this file into the `r_scripts` directory. Then

```
cd r_scripts
R
> install.packages(c("expm", "FKF", "nloptr", "numDeriv"))
> install.packages("./EBM_1.1.0.tar.gz", repos = NULL)