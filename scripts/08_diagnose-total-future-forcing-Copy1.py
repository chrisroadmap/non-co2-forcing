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
import numpy as np
import pickle
import pandas as pd
import matplotlib.pyplot as pl
from scipy.stats import linregress

# %%
pl.style.use('../defaults.mplstyle')

# %%
with open('../output/results_scenarios.pickle', 'rb') as handle:
    data = pickle.load(handle)

# %%
len(data.keys())

# %%
ebm_df = pd.read_csv('../output/calibrations/4xCO2_cummins_ebm2_cmip6.csv')

# %%
scenarios_future = ['ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585']
scenarios = ['historical'] + scenarios_future

# %%
start_index = {
    'historical': 155,
    'ssp119': 240,
    'ssp126': 240,
    'ssp245': 240,
    'ssp370': 240,
    'ssp585': 240,
}

# %%
end_index = {
    'historical': 165,
    'ssp119': 250,
    'ssp126': 250,
    'ssp245': 250,
    'ssp370': 250,
    'ssp585': 250,
}

# %%
non_co2_fraction = {}
non_co2_absolute = {}
for model in data.keys():
    non_co2_fraction[model] = {}
    non_co2_absolute[model] = {}
    for scenario in scenarios:
        # print(model, scenario)
        if scenario not in data[model]:
            continue
        # exclude models not running to end of century
        if model in ['MPI-ESM-1-2-HAM', 'BCC-ESM1'] and scenario != 'historical':
            continue
        n_runs = len(data[model][scenario])
        # if n_runs == 0:
        #     print('zero!')
        #     continue
        non_co2_fraction[model][scenario] = {}
        non_co2_absolute[model][scenario] = {}
        tempsum_fraction = 0
        tempsum_absolute = 0
        for irun, run in enumerate(data[model][scenario]):
            non_co2_fraction[model][scenario][run] = (
                data[model][scenario][run]['forcing_nonco2'][start_index[scenario]:end_index[scenario]].mean() / 
                data[model][scenario][run]['forcing_fitted'][start_index[scenario]:end_index[scenario]].mean()
            )
            non_co2_absolute[model][scenario][run] = data[model][scenario][run]['forcing_nonco2'][start_index[scenario]:end_index[scenario]].mean()
            if not np.isnan(non_co2_fraction[model][scenario][run]):
                tempsum_fraction = tempsum_fraction + non_co2_fraction[model][scenario][run]
                tempsum_absolute = tempsum_absolute + data[model][scenario][run]['forcing_nonco2'][start_index[scenario]:end_index[scenario]].mean()
            else:
                n_runs = n_runs - 1
        non_co2_fraction[model][scenario]['mean'] = tempsum_fraction / n_runs
        non_co2_absolute[model][scenario]['mean'] = tempsum_absolute / n_runs

# %%
non_co2_fraction_mean = {}
non_co2_absolute_mean = {}
for model in non_co2_fraction.keys():
    non_co2_fraction_mean[model] = {}
    non_co2_absolute_mean[model] = {}
    for scenario in scenarios:
        if scenario not in non_co2_fraction[model]:
            continue
        non_co2_fraction_mean[model][scenario] = non_co2_fraction[model][scenario]['mean']
        non_co2_absolute_mean[model][scenario] = non_co2_absolute[model][scenario]['mean']

non_co2_fraction_mean_df = pd.DataFrame(non_co2_fraction_mean).T
non_co2_absolute_mean_df = pd.DataFrame(non_co2_absolute_mean).T

# %%
non_co2_fraction_mean_df

# %%
non_co2_fraction_mean_df.mean()

# %%
non_co2_fraction_mean_df.sort_values('historical')

# %%
colors = {
    'historical': '#000000',
    'ssp119'    : '#00a9cf',
    'ssp126'    : '#003466',
    'ssp245'    : '#f69320',
    'ssp370'    : '#df0000',
    'ssp585'    : '#980002',
}

# %%
fullnames = {
    'historical': 'Historical',
    'ssp119'    : 'SSP1-1.9',
    'ssp126'    : 'SSP1-2.6',
    'ssp245'    : 'SSP2-4.5',
    'ssp370'    : 'SSP3-7.0',
    'ssp585'    : 'SSP5-8.5',
}

# %%
xoffset = {
    'historical': 0,
    'ssp119'    : -0.3,
    'ssp126'    : -0.15,
    'ssp245'    : 0,
    'ssp370'    : +0.15,
    'ssp585'    : +0.30,
}

# %%
markersize = {
    'historical': 49,
    'ssp119'    : 25,
    'ssp126'    : 25,
    'ssp245'    : 25,
    'ssp370'    : 25,
    'ssp585'    : 25,
}

# %%
pl.rcParams['xtick.minor.visible'] = False
pl.rcParams['xtick.major.top'] = False
pl.rcParams['xtick.top'] = False
pl.rcParams['xtick.major.size'] = 0

fig, ax = pl.subplots(figsize=(18/2.54, 10/2.54))
imodel = 0
for fillbounds in np.arange(-0.5, len(data.keys()), 2):
    ax.fill_between([fillbounds, fillbounds+1], -0.75, 0.45, color='0.95')
for model, row in non_co2_fraction_mean_df.sort_values('historical').iterrows():
    # print(row.model, non_co2_fraction[row.model]['mean'])     
    for iscen, scenario in enumerate(scenarios):
        if model=='EC-Earth3-Veg':
            label=fullnames[scenario]
        else:
            label=None
        if scenario not in data[model]:
            continue
        # temporary flag for either missing model or wrong output data in dict: will be fixed
        if model in ['MPI-ESM-1-2-HAM', 'BCC-ESM1'] and scenario != "historical":
            continue
        for run in data[model][scenario]:
            ax.scatter(imodel+xoffset[scenario], non_co2_fraction[model][scenario][run], marker='x', color=colors[scenario], s=9, alpha=0.75/len(data[model][scenario]))
        ax.scatter(imodel+xoffset[scenario], row[scenario], marker='_', color=colors[scenario], s=markersize[scenario], zorder=7, label=label)
    imodel=imodel+1
ax.axhline(0, ls=':', color='k')
ax.set_xticks(np.arange(len(non_co2_fraction_mean_df)));
ax.set_xticklabels(non_co2_fraction_mean_df.sort_values('historical').index, rotation=90);
ax.set_xlim(-0.8, len(non_co2_fraction_mean_df)-0.2)
ax.set_ylim(-0.75, 0.45)
ax.set_title('non-CO$_2$ radiative forcing fraction in CMIP6 historical (2005-14) and future (2090-99) simulations')
ax.legend(loc='lower right')
fig.tight_layout()
pl.savefig('../plots/non-co2-fraction-historical-future.png')

# %%
pl.rcParams['xtick.minor.visible'] = True
pl.rcParams['xtick.major.top'] = True
pl.rcParams['xtick.top'] = True
pl.rcParams['xtick.major.size'] = 3.5

fig, ax = pl.subplots(figsize=(9/2.54, 9/2.54))
for scenario in ['ssp370']:
    ax.scatter(non_co2_fraction_mean_df['historical'], non_co2_fraction_mean_df[scenario], facecolor='None', edgecolors=colors[scenario]) 
    mask = ~np.isnan(non_co2_fraction_mean_df['historical']) & ~np.isnan(non_co2_fraction_mean_df[scenario])
    lr = linregress(non_co2_fraction_mean_df['historical'][mask], non_co2_fraction_mean_df[scenario][mask])
    xmin = non_co2_fraction_mean_df['historical'].min()
    xmax = non_co2_fraction_mean_df['historical'].max()
    pl.plot(np.linspace(xmin, xmax), np.linspace(xmin, xmax) * lr.slope + lr.intercept, color=colors[scenario]) 
    print(scenario, lr)
