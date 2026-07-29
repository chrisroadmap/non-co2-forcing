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

# %% [markdown]
# # Calculate total, CO2 and non-CO2 forcing in each model
#
# We take the EBM calibrations from 4xCO2, solve them iteratatively with N and T, to estimate F.

# %%
import os
import glob
import pickle
from pathlib import PurePath

import matplotlib.pyplot as pl
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from fair.energy_balance_model import EnergyBalanceModel
from scipy.optimize import minimize

# %%
ebm_df = pd.read_csv('../output/calibrations/4xCO2_cummins_ebm2_cmip6.csv')


# %%
# function we want to minimise
def cost_function(forcing, tas_truth, rndt_truth, C1, C2, kappa1, kappa2, epsilon):
    ebm = EnergyBalanceModel(
        ocean_heat_capacity=[C1, C2],
        ocean_heat_transfer=[kappa1, kappa2],
        deep_ocean_efficacy=epsilon,
        gamma_autocorrelation=1000,
    )
    ebm.add_forcing(forcing, timestep=1)
    ebm.run()
    tas_model = ebm.temperature[:, 0]
    rndt_model = ebm.toa_imbalance

    return np.sum((tas_model - tas_truth) ** 2 + (rndt_model - rndt_truth) ** 2)


# %%
overrides = {}
for idx, row in ebm_df.iterrows():
    overrides[row.model] = 'r*' + row.run[-6:]

# %%
overrides

# %%
overrides['GISS-E2-1-G'] = 'r*i1p1f1'
overrides['GISS-E2-1-H'] = 'r*i1p1f2'
overrides['CanESM5'] = 'r*i1p2f1'


# %%
data = {}

for idx, row in tqdm(ebm_df.iterrows()):
    # pattern_match = 'r*' + row.run[-6:]
    available_files = glob.glob(f"../data/cmip6-hbf/cmip_data/{row.model}/historical/*_{overrides[row.model]}_*")
    if len(available_files) == 0:
        continue
    data[row.model] = {}
    data[row.model]['historical'] = {}
    for file in available_files:
        run = PurePath(file).parts[6].split("_")[2]
        
        # remove short and fast warming runs from EC-Earth3
        if row.model=='EC-Earth3':
            if int(run.split('i')[0][1:]) >= 100:
                continue
        data[row.model]['historical'][run] = {}
        hist_df = pd.read_csv(file, index_col=0)
        hist_df['rndt'] = hist_df['rsdt'] - hist_df['rsut'] - hist_df['rlut']
        
        # AWI has NaN tas in 1850 - relatively safe to call it zero
        if row.model=='AWI-CM-1-1-MR':
            hist_df.loc[1850, 'tas'] = 0

        # TODO: remove short run 1 and run 4 from E3SM-1-0

        # TODO: remove r10 from EC-Earth3-Veg - missing 1880

        # EBM parameters
        C1 = ebm_df.loc[ebm_df['model']==row.model, 'C1'].values[0]
        C2 = ebm_df.loc[ebm_df['model']==row.model, 'C2'].values[0]
        kappa1 = ebm_df.loc[ebm_df['model']==row.model, 'kappa1'].values[0]
        kappa2 = ebm_df.loc[ebm_df['model']==row.model, 'kappa2'].values[0]
        epsilon = ebm_df.loc[ebm_df['model']==row.model, 'epsilon'].values[0]
    
        # initialise and get convergence
        forcing_initial_guess = hist_df.rndt.values + kappa1 * hist_df.tas.values - (epsilon - 1)*kappa2*(hist_df.tas.values - 0)
        res = minimize(cost_function, forcing_initial_guess, args=(hist_df.tas.values, hist_df.rndt.values, C1, C2, kappa1, kappa2, epsilon))

        # run EBM one more time in forward mode with forcing
        ebm = EnergyBalanceModel(
            ocean_heat_capacity=[C1, C2],
            ocean_heat_transfer=[kappa1, kappa2],
            deep_ocean_efficacy=epsilon,
            gamma_autocorrelation=1000
        )
        ebm.add_forcing(res.x, timestep=1)
        ebm.run()
        
        data[row.model]['historical'][run]['temperature_fitted'] = ebm.temperature[:, 0]
        data[row.model]['historical'][run]['temperature_truth'] = hist_df.tas.values
        data[row.model]['historical'][run]['toa_fitted'] = ebm.toa_imbalance
        data[row.model]['historical'][run]['toa_truth'] = hist_df.rndt.values
        data[row.model]['historical'][run]['forcing_fitted'] = res.x

# %%
# # manual trimming - EC-Earth3 simulations with r>100 are too short and warm very quickly
# ecearth3_runs = list(data['EC-Earth3']['historical'].keys())
# for run in ecearth3_runs:
#     if int(run.split('i')[0][1:]) >= 100:
#         del data['EC-Earth3']['historical'][run]

# %%
os.makedirs('../plots/diagnostics/temperature', exist_ok=True)
for idx, row in ebm_df.iterrows():
    if row.model not in data:
        continue
    n_runs = len(data[row.model]['historical'])
    grid_size = int(np.ceil(np.sqrt(n_runs)))
    fig, ax = pl.subplots(grid_size, grid_size, figsize=(16, 12), squeeze=False)
    for irun, run in enumerate(data[row.model]['historical']):
        ax[irun//grid_size, irun%grid_size].plot(np.arange(1850, 1850+len(data[row.model]['historical'][run]['temperature_fitted'])), data[row.model]['historical'][run]['temperature_fitted'], color='k')
        ax[irun//grid_size, irun%grid_size].plot(np.arange(1850, 1850+len(data[row.model]['historical'][run]['temperature_fitted'])), data[row.model]['historical'][run]['temperature_truth'], color='b')
        ax[irun//grid_size, irun%grid_size].set_title(f'{row.model}, {run}')
    fig.tight_layout()
    pl.savefig(f'../plots/diagnostics/temperature/{row.model}.png')

# %%
os.makedirs('../plots/diagnostics/toa', exist_ok=True)
for idx, row in ebm_df.iterrows():
    if row.model not in data:
        continue
    n_runs = len(data[row.model]['historical'])
    grid_size = int(np.ceil(np.sqrt(n_runs)))
    fig, ax = pl.subplots(grid_size, grid_size, figsize=(16, 12), squeeze=False)
    for irun, run in enumerate(data[row.model]['historical']):
        ax[irun//grid_size, irun%grid_size].plot(np.arange(1850, 1850+len(data[row.model]['historical'][run]['toa_fitted'])), data[row.model]['historical'][run]['toa_fitted'], color='k')
        ax[irun//grid_size, irun%grid_size].plot(np.arange(1850, 1850+len(data[row.model]['historical'][run]['toa_truth'])), data[row.model]['historical'][run]['toa_truth'], color='b')
        ax[irun//grid_size, irun%grid_size].set_title(f'{row.model}, {run}')
    fig.tight_layout()
    pl.savefig(f'../plots/diagnostics/toa/{row.model}.png')

# %%
os.makedirs('../plots/diagnostics/forcing', exist_ok=True)
for idx, row in ebm_df.iterrows():
    if row.model not in data:
        continue
    n_runs = len(data[row.model]['historical'])
    grid_size = int(np.ceil(np.sqrt(n_runs)))
    fig, ax = pl.subplots(grid_size, grid_size, figsize=(16, 12), squeeze=False)
    for irun, run in enumerate(data[row.model]['historical']):
        ax[irun//grid_size, irun%grid_size].plot(np.arange(1850, 1850+len(data[row.model]['historical'][run]['forcing_fitted'])), data[row.model]['historical'][run]['forcing_fitted'], color='k')
        ax[irun//grid_size, irun%grid_size].set_title(f'{row.model}, {run}')
    fig.tight_layout()
    pl.savefig(f'../plots/diagnostics/forcing/{row.model}.png')

# %%
# estimate CO2 and non-CO2 forcing
df_co2_conc_ssp = pd.read_csv('../data/ssp_co2_concentration.csv', index_col=0)
alpha = ebm_df['F_4xCO2'] / np.log(4)

x0=284.316999854786

def myhre(x, alpha):
    return alpha * np.log(x/x0)

co2_conc = df_co2_conc_ssp.loc[:, 'ssp245'].values

for idx, row in ebm_df.iterrows():
    if row.model not in data:
        continue
    for irun, run in enumerate(data[row.model]['historical']):
        sim_length = len(data[row.model]['historical'][run]['forcing_fitted'])
        data[row.model]['historical'][run]['forcing_co2'] = myhre(co2_conc[:sim_length], ebm_df.loc[ebm_df.model==row.model, ['F_4xCO2']].values[0,0] / np.log(4))
        data[row.model]['historical'][run]['forcing_nonco2'] = data[row.model]['historical'][run]['forcing_fitted'] - data[row.model]['historical'][run]['forcing_co2']

# %%
os.makedirs('../plots/diagnostics/forcing_nonco2', exist_ok=True)
for idx, row in ebm_df.iterrows():
    if row.model not in data:
        continue
    n_runs = len(data[row.model]['historical'])
    grid_size = int(np.ceil(np.sqrt(n_runs)))
    fig, ax = pl.subplots(grid_size, grid_size, figsize=(16, 12), squeeze=False)
    for irun, run in enumerate(data[row.model]['historical']):
        ax[irun//grid_size, irun%grid_size].plot(np.arange(1850, 1850+len(data[row.model]['historical'][run]['forcing_nonco2'])), data[row.model]['historical'][run]['forcing_nonco2'], color='k')
        ax[irun//grid_size, irun%grid_size].set_title(f'{row.model}, {run}')
    fig.tight_layout()
    pl.savefig(f'../plots/diagnostics/forcing_nonco2/{row.model}.png')

# %%
with open('../output/results.pickle', 'wb') as handle:
    pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)

# %%
