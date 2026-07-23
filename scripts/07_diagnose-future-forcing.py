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
import os
from pathlib import PurePath
import glob

from tqdm.auto import tqdm
import numpy as np
import pickle
import pandas as pd
import matplotlib.pyplot as pl

from fair.energy_balance_model import EnergyBalanceModel
from scipy.optimize import minimize

# %%
with open('../output/results.pickle', 'rb') as handle:
    data = pickle.load(handle)

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
overrides = {model: '*' for model in data}
overrides['GISS-E2-1-H'] = 'r*i1p1f2'
overrides['CanESM5'] = 'r*i1p2f1'

# %%
# has a strange offset issue in 1pctCO2 and historical simulation also starts a bit warm - better to not trust
del data['KIOST-ESM']

# %%
for model in tqdm(data):
    for scenario in ['ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585']:
        available_files = glob.glob(f"../data/cmip6-hbf/cmip_data/{model}/{scenario}/*{overrides[model]}*.csv")
        if len(available_files) > 0:
            data[model][scenario] = {}
        for file in available_files:
            run = PurePath(file).parts[6].split("_")[2]
        
            # remove short and fast warming runs from EC-Earth3
            if model=='EC-Earth3':
                if int(run.split('i')[0][1:]) >= 100:
                    continue

            # remove dodgy r10
            if model=='EC-Earth3-Veg':
                if int(run.split('i')[0][1:]) == 10:
                    continue
            
            if run in data[model]['historical']:
                data[model][scenario][run] = {}
                scen_df = pd.read_csv(file, index_col=0)
                scen_df['rndt'] = scen_df['rsdt'] - scen_df['rsut'] - scen_df['rlut']
                data[model][scenario][run]['toa_truth'] = np.zeros(165 + len(scen_df['rndt']))
                data[model][scenario][run]['toa_truth'][:165] = data[model]['historical'][run]['toa_truth'][:165]
                data[model][scenario][run]['toa_truth'][165:] = scen_df['rndt']
                data[model][scenario][run]['temperature_truth'] = np.zeros(165 + len(scen_df['tas']))
                data[model][scenario][run]['temperature_truth'][:165] = data[model]['historical'][run]['temperature_truth'][:165]
                data[model][scenario][run]['temperature_truth'][165:] = scen_df['tas']
            

                # EBM parameters
                C1 = ebm_df.loc[ebm_df['model']==model, 'C1'].values[0]
                C2 = ebm_df.loc[ebm_df['model']==model, 'C2'].values[0]
                kappa1 = ebm_df.loc[ebm_df['model']==model, 'kappa1'].values[0]
                kappa2 = ebm_df.loc[ebm_df['model']==model, 'kappa2'].values[0]
                epsilon = ebm_df.loc[ebm_df['model']==model, 'epsilon'].values[0]
    
                # initialise and get convergence
                forcing_initial_guess = (
                    data[model][scenario][run]['toa_truth'] + kappa1 * 
                    data[model][scenario][run]['temperature_truth'] - (epsilon - 1)*kappa2*(data[model][scenario][run]['temperature_truth'] - 0)
                )
                res = minimize(
                    cost_function, 
                    forcing_initial_guess, 
                    args=(
                        data[model][scenario][run]['temperature_truth'], 
                        data[model][scenario][run]['toa_truth'], C1, C2, kappa1, kappa2, epsilon
                    )
                )

                # run EBM one more time in forward mode with forcing
                ebm = EnergyBalanceModel(
                    ocean_heat_capacity=[C1, C2],
                    ocean_heat_transfer=[kappa1, kappa2],
                    deep_ocean_efficacy=epsilon,
                    gamma_autocorrelation=1000
                )
                ebm.add_forcing(res.x, timestep=1)
                ebm.run()
        
                data[model][scenario][run]['temperature_fitted'] = ebm.temperature[:, 0]
                data[model][scenario][run]['toa_fitted'] = ebm.toa_imbalance
                data[model][scenario][run]['forcing_fitted'] = res.x

# %%
# estimate CO2 and non-CO2 forcing
df_co2_conc_ssp = pd.read_csv('../data/ssp_co2_concentration.csv', index_col=0)
alpha = ebm_df['F_4xCO2'] / np.log(4)

x0=284.316999854786

def myhre(x, alpha):
    return alpha * np.log(x/x0)

for scenario in ['ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585']:
    co2_conc = df_co2_conc_ssp.loc[:, scenario].values

    for model in data:
        if scenario not in data[model]:
            continue
        n_runs = len(data[model][scenario])
        if n_runs == 0:
            continue
        for irun, run in enumerate(data[model][scenario]):
            sim_length = len(data[model][scenario][run]['forcing_fitted'])
            data[model][scenario][run]['forcing_co2'] = myhre(co2_conc[:sim_length], ebm_df.loc[ebm_df.model==model, ['F_4xCO2']].values[0,0] / np.log(4))
            data[model][scenario][run]['forcing_nonco2'] = data[model][scenario][run]['forcing_fitted'] - data[model][scenario][run]['forcing_co2']

# %%
with open('../output/results_scenarios.pickle', 'wb') as handle:
    pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)

# %%
for scenario in ['ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585']:
    os.makedirs(f'../plots/diagnostics/temperature/{scenario}', exist_ok=True)
    for model in data:
        if scenario not in data[model]:
            continue
        n_runs = len(data[model][scenario])
        if n_runs == 0:
            continue
        grid_size = int(np.ceil(np.sqrt(n_runs)))
        fig, ax = pl.subplots(grid_size, grid_size, figsize=(16, 12), squeeze=False)
        for irun, run in enumerate(data[model][scenario]):
            ax[irun//grid_size, irun%grid_size].plot(np.arange(1850, 1850+len(data[model][scenario][run]['temperature_fitted'])), data[model][scenario][run]['temperature_fitted'], color='k')
            ax[irun//grid_size, irun%grid_size].plot(np.arange(1850, 1850+len(data[model][scenario][run]['temperature_truth'])), data[model][scenario][run]['temperature_truth'], color='b')
            ax[irun//grid_size, irun%grid_size].set_title(f'{model}, {scenario}, {run}')
        fig.tight_layout()
        pl.savefig(f'../plots/diagnostics/temperature/{scenario}/{model}.png')

# %%
for scenario in ['ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585']:
    os.makedirs(f'../plots/diagnostics/forcing_nonco2/{scenario}', exist_ok=True)
    for model in data:
        if scenario not in data[model]:
            continue
        n_runs = len(data[model][scenario])
        if n_runs == 0:
            continue
        grid_size = int(np.ceil(np.sqrt(n_runs)))
        fig, ax = pl.subplots(grid_size, grid_size, figsize=(16, 12), squeeze=False)
        for irun, run in enumerate(data[model][scenario]):
            ax[irun//grid_size, irun%grid_size].plot(np.arange(1850, 1850+len(data[model][scenario][run]['forcing_nonco2'])), data[model][scenario][run]['forcing_nonco2'], color='k')
            ax[irun//grid_size, irun%grid_size].set_title(f'{model}, {scenario}, {run}')
        fig.tight_layout()
        pl.savefig(f'../plots/diagnostics/forcing_nonco2/{scenario}/{model}.png')
