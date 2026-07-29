# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.3
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
def moving_average(a, window=11):
    ret = np.cumsum(a, dtype=float)
    ret[window:] = ret[window:] - ret[:-window]
    return ret[window - 1:] / window


# %%
for model in data.keys():
    pl.plot(np.arange(1855, 2010), moving_average(non_co2_absolute[model]['mean'])/moving_average(all_absolute[model]['mean']))
pl.ylim(-1, 1);

# %%
for model in data.keys():
    plotted = moving_average(non_co2_absolute[model]['mean'])/moving_average(all_absolute[model]['mean'])
    plotted[plotted<-5] = np.nan
    plotted[plotted>5] = np.nan
    pl.plot(np.arange(1855, 2010), plotted);

# %%
for model in data.keys():
    pl.plot(np.arange(1850, 2015), non_co2_absolute[model]['mean']);
pl.xlim(1850, 1900);

# %%
for model in data.keys():
    pl.plot(np.arange(1850, 2015), non_co2_absolute[model]['mean']);

# %%
for model in data.keys():
    pl.plot(np.arange(1860, 2005), moving_average(non_co2_absolute[model]['mean'], 21));

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
# since CMIP includes natural forcing too, don't do anthro only

# %%
df_gcb = pd.read_csv('../data/gcb_2025_v1.0.csv', index_col=0)

# %%
df_gcb

# %%
df_gcb.cumsum().sum(axis=1)

# %%
df_gcb_rebased = df_gcb - df_gcb.loc[1850:1900].mean()
df_gcb_rebased

# %%
df_gcb_rebased.cumsum().loc[2024]

# %%
df_gcb.cumsum().loc[2024] - df_gcb.cumsum().loc[1900]

# %%
df_gcb.cumsum().loc[2024] - df_gcb.cumsum().loc[1876]

# %%
df_hadcrut = pd.read_csv('../data/HadCRUT.5.1.0.0.analysis.ensemble_series.global.annual.csv', index_col=0)
df_hadcrut_rebased = pd.DataFrame()
for real in range(1, 201):
    df_hadcrut_rebased[f'Realization {real}'] = df_hadcrut[f'Realization {real}'] - df_hadcrut.loc[1850:1899, f'Realization {real}'].mean()

# %%
df_hadcrut_rebased

# %%
df_hadcrut_rebased.loc[1900].mean()

# %%
delta_T = df_hadcrut_rebased.loc[2005:2014].mean().mean()

# %%
delta_T

# %%
cumul_E = (df_gcb.cumsum().loc[2009] - df_gcb.cumsum().loc[1874]).sum()

# %%
df_gcb.cumsum().loc[2000:]

# %%
cumul_E

# %%
non_co2_fraction_20052014 = {}
non_co2_absolute_20052014 = {}
for model in data.keys():
    non_co2_fraction_20052014[model] = {}
    non_co2_absolute_20052014[model] = {}
    n_runs = len(data[model]['historical'])
    tempsum_fraction = 0
    tempsum_absolute = 0
    for irun, run in enumerate(data[model]['historical']):
        non_co2_fraction_20052014[model][run] = (
            data[model]['historical'][run]['forcing_nonco2'][155:165].mean() / 
            data[model]['historical'][run]['forcing_fitted'][155:165].mean()
        )
        non_co2_absolute_20052014[model][run] = data[model]['historical'][run]['forcing_nonco2'][155:165].mean()
        if not np.isnan(non_co2_fraction_20052014[model][run]):
            tempsum_fraction = tempsum_fraction + non_co2_fraction_20052014[model][run]
            tempsum_absolute = tempsum_absolute + data[model]['historical'][run]['forcing_nonco2'][155:165].mean()
        else:
            n_runs = n_runs - 1
    non_co2_fraction_20052014[model]['mean'] = tempsum_fraction / n_runs
    non_co2_absolute_20052014[model]['mean'] = tempsum_absolute / n_runs

# %%
non_co2_fraction_mean_20052014 = {}
non_co2_absolute_mean_20052014 = {}
for model in data.keys():
    non_co2_fraction_mean_20052014[model] = non_co2_fraction_20052014[model]['mean']
    non_co2_absolute_mean_20052014[model] = non_co2_absolute_20052014[model]['mean']

non_co2_fraction_mean_20052014_df = pd.DataFrame(non_co2_fraction_mean_20052014, index = ['mean']).T
non_co2_absolute_mean_20052014_df = pd.DataFrame(non_co2_absolute_mean_20052014, index = ['mean']).T

# %%
pl.hist(non_co2_fraction_mean_20052014_df)

# %%
non_co2_fraction_mean_20052014_df.mean()

# %%
delta_T / cumul_E * (1-non_co2_fraction_mean_20052014_df.mean())  # K GtC-1
1000 * delta_T / cumul_E * (1-non_co2_fraction_mean_20052014_df.mean())  # K (1000 GtC)-1
#1000 * delta_T / (44.009 / 12.011 * cumul_E) * (1-non_co2_fraction_mean_df.mean())  # K (1000 GtCO2)-1
# since CMIP includes natural forcing too, don't do anthro only

# %%
