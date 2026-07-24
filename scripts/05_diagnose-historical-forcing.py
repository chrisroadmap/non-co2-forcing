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
import numpy as np
import pickle
import pandas as pd
import matplotlib.pyplot as pl

# %%
pl.style.use('../defaults.mplstyle')
pl.rcParams['xtick.minor.visible'] = False
pl.rcParams['xtick.major.top'] = False
pl.rcParams['xtick.top'] = False
pl.rcParams['xtick.major.size'] = 0
#pl.rcParams['xtick.bottom'] = False

# %%
with open('../output/results.pickle', 'rb') as handle:
    data = pickle.load(handle)

# %%
len(data.keys())

# %%
# has a strange offset issue in 1pctCO2 and historical simulation also starts a bit warm - better to not trust
del data['KIOST-ESM']

# %%
len(data.keys())

# %%
ebm_df = pd.read_csv('../output/calibrations/4xCO2_cummins_ebm2_cmip6.csv')

# %%
non_co2_fraction = {}
non_co2_absolute = {}
for model in data.keys():
    non_co2_fraction[model] = {}
    non_co2_absolute[model] = {}
    n_runs = len(data[model]['historical'])
    tempsum_fraction = 0
    tempsum_absolute = 0
    for irun, run in enumerate(data[model]['historical']):
        non_co2_fraction[model][run] = (
            data[model]['historical'][run]['forcing_nonco2'][155:165].mean() / 
            data[model]['historical'][run]['forcing_fitted'][155:165].mean()
        )
        non_co2_absolute[model][run] = data[model]['historical'][run]['forcing_nonco2'][155:165].mean()
        if not np.isnan(non_co2_fraction[model][run]):
            tempsum_fraction = tempsum_fraction + non_co2_fraction[model][run]
            tempsum_absolute = tempsum_absolute + data[model]['historical'][run]['forcing_nonco2'][155:165].mean()
        else:
            n_runs = n_runs - 1
    non_co2_fraction[model]['mean'] = tempsum_fraction / n_runs
    non_co2_absolute[model]['mean'] = tempsum_absolute / n_runs

# %%
non_co2_fraction_mean = {}
non_co2_absolute_mean = {}
for model in data.keys():
    non_co2_fraction_mean[model] = non_co2_fraction[model]['mean']
    non_co2_absolute_mean[model] = non_co2_absolute[model]['mean']

non_co2_fraction_mean_df = pd.DataFrame(non_co2_fraction_mean, index = ['mean']).T
non_co2_absolute_mean_df = pd.DataFrame(non_co2_absolute_mean, index = ['mean']).T

# %%
non_co2_fraction_mean_df.sort_values('mean')

# %%
fig, ax = pl.subplots(figsize=(18/2.54, 10/2.54))
imodel = 0
for fillbounds in np.arange(-0.5, len(data.keys()), 2):
    ax.fill_between([fillbounds, fillbounds+1], -0.75, 0.4, color='0.95')
for model, row in non_co2_fraction_mean_df.sort_values('mean').iterrows():
    # print(row.model, non_co2_fraction[row.model]['mean'])
    for run in data[model]['historical']:
        ax.scatter(imodel, non_co2_fraction[model][run], marker='x', color='0.6', s=9)
    ax.scatter(imodel, row['mean'], marker='_', color='k', s=49)
    imodel=imodel+1
ax.axhline(0, ls=':', color='k')
ax.set_xticks(np.arange(len(non_co2_fraction_mean_df)));
ax.set_xticklabels(non_co2_fraction_mean_df.sort_values('mean').index, rotation=90);
ax.set_xlim(-0.8, len(data.keys())-0.2)
ax.set_ylim(-0.75, 0.4)
ax.set_title('non-CO2 forcing fraction in CMIP6 historical simulations, 2005-14 relative to 1850')
fig.tight_layout()
pl.savefig('../plots/non-co2-fraction-historical.png')

# %%
fig, ax = pl.subplots(figsize=(18/2.54, 10/2.54))
imodel = 0
for fillbounds in np.arange(-0.5, len(data.keys()), 2):
    ax.fill_between([fillbounds, fillbounds+1], -0.8, 1, color='0.95')
for model, row in non_co2_absolute_mean_df.sort_values('mean').iterrows():
    # print(row.model, non_co2_fraction[row.model]['mean'])
    for run in data[model]['historical']:
        ax.scatter(imodel, non_co2_absolute[model][run], marker='x', color='0.6', s=9)
    ax.scatter(imodel, row['mean'], marker='_', color='k', s=49)
    imodel=imodel+1
ax.axhline(0, ls=':', color='k')
ax.set_xticks(np.arange(len(non_co2_absolute_mean_df)));
ax.set_xticklabels(non_co2_absolute_mean_df.sort_values('mean').index, rotation=90);
ax.set_xlim(-0.8, len(data.keys())-0.2)
ax.set_ylim(-0.8, 1.0)
ax.set_title('non-CO2 forcing in CMIP6 historical simulations, 2005-14 relative to 1850')
ax.set_ylabel('W m$^{-2}$')
fig.tight_layout()
pl.savefig('../plots/non-co2-absolute-historical.png')

# %%
