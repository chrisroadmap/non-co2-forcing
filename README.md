# non-co2-forcing

Calculate the split of CO2 and non-CO2 forcing from Earth System Models.

## installation

### requirements
- `anaconda` for `python3`
- `python>=3.7`
- for the Cummins calibration, `R>=4.1.1` and `cmake>=3.2`

### set up environments for Python and R

First, create the `conda` environment for `python`. You will require an installation of `anaconda`. We recommend [miniconda](https://docs.anaconda.com/miniconda/) or [miniforge](https://github.com/conda-forge/miniforge). After installing `conda`, run the following commands

 
```
conda env create -f non-co2-forcing
conda activate fair-calibrate
```

As of around June 2024, the `R` environments no longer install properly from `conda`, so we'll need to do this manually.

First, you'll require a version of `R` at least as new as 4.1.1 and `cmake` at least as new as 3.2. Depending on your Linux distribution this might do the trick

```
sudo apt install r-base-core
sudo apt install cmake
```

If not, grab `R` from https://www.r-project.org/. `cmake` is probably also freely distributed somewhere.

Second, you'll need Donald Cummins' `EBM` package binary, which can be downloaded from [here](https://github.com/donaldcummins/EBM/archive/refs/tags/v1.1.0.tar.gz). Download this file into base directory. Then

```
R
> install.packages(c("expm", "FKF", "nloptr", "numDeriv"))
> install.packages("./EBM-1.1.0.tar.gz", repos = NULL)
```

## Reproduction

`cd` to the `scripts` directory and run everything in order. Scripts number 1 and 2 are executable, after this you should open them in a `jupyter` notebook.
