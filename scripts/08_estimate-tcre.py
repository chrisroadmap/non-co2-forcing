# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
import copy

from tqdm.auto import tqdm
import numpy as np
import pickle
import pandas as pd
import matplotlib.pyplot as pl

# %%
with open('../output/results.pickle', 'rb') as handle:
    data = pickle.load(handle)

# %%
# has a strange offset issue in 1pctCO2 and historical simulation also starts a bit warm - better to not trust
del data['KIOST-ESM']

# %%
non_co2_fraction = {}
non_co2_absolute = {}
all_absolute = {}
for model in data.keys():
    all_absolute[model] = {}
    non_co2_absolute[model] = {}
    non_co2_fraction[model] = {}
    n_runs = len(data[model]['historical'])
    tempsum_non_co2_absolute = 0
    tempsum_all_absolute = 0
    for irun, run in enumerate(data[model]['historical']):
        print(model, run)
        all_absolute[model][run] = data[model]['historical'][run]['forcing_fitted']
        non_co2_absolute[model][run] = data[model]['historical'][run]['forcing_nonco2']
        if not np.sum(np.isnan(non_co2_absolute[model][run])) and len(non_co2_absolute[model][run])>=165:
            tempsum_all_absolute = tempsum_all_absolute + data[model]['historical'][run]['forcing_fitted'][:165]
            tempsum_non_co2_absolute = tempsum_non_co2_absolute + data[model]['historical'][run]['forcing_nonco2'][:165]
        else:
            n_runs = n_runs - 1
    non_co2_fraction[model]['mean'] = tempsum_non_co2_absolute / tempsum_all_absolute / n_runs
    non_co2_absolute[model]['mean'] = tempsum_non_co2_absolute / n_runs
    all_absolute[model]['mean'] = tempsum_all_absolute / n_runs

# %%
for model in data.keys():
    pl.plot(np.arange(1850, 2015), non_co2_absolute[model]['mean']);
pl.xlim(1850, 1870);

# %%
plotted = copy.deepcopy(non_co2_fraction)

# %%
for model in data.keys():
    plotted[model]['mean'][plotted[model]['mean'] > 0.5] = np.nan
    plotted[model]['mean'][plotted[model]['mean'] < -0.5] = np.nan

# %%
for model in data.keys():
    pl.plot(np.arange(1850, 2015), plotted[model]['mean']);
pl.ylim(-0.5, 0.5)

# %%
pl.plot(np.arange(1850, 2015), non_co2_absolute['CanESM5']['mean']);

# %%
pl.plot(np.arange(1850, 2015), non_co2_fraction['CanESM5']['mean']);
pl.ylim(-1, 1)

# %%
# Damon Matthews estimate

# TCRE = \delta T / E * (1-f_nc)
