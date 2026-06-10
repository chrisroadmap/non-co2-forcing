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
# # Calculate total forcing in each model
#
# We take the EBM calibrations and invert them.

# %%
import os
import glob
from pathlib import PurePath

import matplotlib.pyplot as pl
import numpy as np
import pandas as pd

from fair.energy_balance_model import EnergyBalanceModel
from scipy.optimize import minimize

# %%
ebm_df = pd.read_csv('../output/calibrations/4xCO2_cummins_ebm2_cmip6.csv')


# %%
def cost_function(forcing, tas_truth, rndt_truth, C1, C2, kappa1, kappa2, epsilon):
    ebm = EnergyBalanceModel(
        ocean_heat_capacity=[C1, C2],
        ocean_heat_transfer=[kappa1, kappa2],
        deep_ocean_efficacy=epsilon
    )
    ebm.add_forcing(forcing, timestep=1)
    ebm.run()
    tas_model = ebm.temperature[:, 0]
    rndt_model = ebm.toa_imbalance

    return np.sum((tas_model - tas_truth) ** 2 + (rndt_model - rndt_truth) ** 2)


def run_model(test_df, model):
    test_df['rndt'] = test_df['rsdt'] - test_df['rsut'] - test_df['rlut']
    C1 = ebm_df.loc[ebm_df['model']==model, 'C1'].values[0]
    C2 = ebm_df.loc[ebm_df['model']==model, 'C2'].values[0]
    kappa1 = ebm_df.loc[ebm_df['model']==model, 'kappa1'].values[0]
    kappa2 = ebm_df.loc[ebm_df['model']==model, 'kappa2'].values[0]
    epsilon = ebm_df.loc[ebm_df['model']==model, 'epsilon'].values[0]

    forcing_initial_guess = test_df.rndt.values + kappa1 * test_df.tas.values - (epsilon - 1)*kappa2*(test_df.tas.values - 0)
    res = minimize(cost_function, forcing_initial_guess, args=(test_df.tas.values, test_df.rndt.values, C1, C2, kappa1, kappa2, epsilon))

    return res, ebm


# %%
for idx, row in ebm_df.iterrows():
    available_files = glob.glob(f"../data/cmip6-hbf/cmip_data/{row.model}/historical/*")
    print(available_files)

# %%
for idx, row in ebm_df.iterrows():
    try:
        hist_df = pd.read_csv(f'../data/cmip6-hbf/cmip_data/{row.model}/1pctCO2/{row.model}_1pctCO2_{row.run}_anomalies.csv', index_col=0)
        # res, ebm = run_model(test_df, row.model)
        test_df['rndt'] = test_df['rsdt'] - test_df['rsut'] - test_df['rlut']
        C1 = row.C1
        C2 = row.C2
        kappa1 = row.kappa1
        kappa2 = row.kappa2
        epsilon = row.epsilon
        
        ebm = EnergyBalanceModel(
            ocean_heat_capacity=[C1, C2],
            ocean_heat_transfer=[kappa1, kappa2],
            deep_ocean_efficacy=epsilon
        )

        ebm.add_forcing(myhre(1.01 ** np.arange(len(test_df.tas.values)) * x0, row.alpha), timestep=1)
        ebm.run()
    
        ax[idx//7, idx%7].plot(ebm.temperature[:, 0], color='k')
        ax[idx//7, idx%7].plot(test_df.tas.values, color='k', ls='--')
        ax[idx//7, idx%7].plot(ebm.toa_imbalance, color='r')
        ax[idx//7, idx%7].plot(test_df.rndt.values, color='r', ls='--')
        ax[idx//7, idx%7].set_title(row.model)
    except:
        print(row.model)

fig.tight_layout()

# %%
test_df = pd.read_csv('../data/cmip6-hbf/cmip_data/UKESM1-0-LL/historical/UKESM1-0-LL_historical_r1i1p1f2_anomalies.csv', index_col=0)

# %%
test_df2 = pd.read_csv('../data/cmip6-hbf/cmip_data/UKESM1-0-LL/ssp126/UKESM1-0-LL_ssp126_r1i1p1f2_anomalies.csv', index_col=0)

# %%
test_df

# %%
test_df2

# %%
test_df = pd.concat((test_df, test_df2))

# %%
test_df['rndt'] = test_df['rsdt'] - test_df['rsut'] - test_df['rlut']

# %%
test_df

# %%
C1 = ebm_df.loc[ebm_df['model']=='UKESM1-0-LL', 'C1'].values[0]
C2 = ebm_df.loc[ebm_df['model']=='UKESM1-0-LL', 'C2'].values[0]
kappa1 = ebm_df.loc[ebm_df['model']=='UKESM1-0-LL', 'kappa1'].values[0]
kappa2 = ebm_df.loc[ebm_df['model']=='UKESM1-0-LL', 'kappa2'].values[0]
epsilon = ebm_df.loc[ebm_df['model']=='UKESM1-0-LL', 'epsilon'].values[0]


# %%
def cost_function(forcing, tas_truth, rndt_truth, C1, C2, kappa1, kappa2, epsilon):
    ebm = EnergyBalanceModel(
        ocean_heat_capacity=[C1, C2],
        ocean_heat_transfer=[kappa1, kappa2],
        deep_ocean_efficacy=epsilon
    )
    ebm.add_forcing(forcing, timestep=1)
    ebm.run()
    tas_model = ebm.temperature[:, 0]
    rndt_model = ebm.toa_imbalance

    return np.sum((tas_model - tas_truth) ** 2 + (rndt_model - rndt_truth) ** 2)


# %%
forcing_initial_guess = test_df.rndt.values + kappa1 * test_df.tas.values - (epsilon - 1)*kappa2*(test_df.tas.values - 0)

# %%
res = minimize(cost_function, forcing_initial_guess, args=(test_df.tas.values, test_df.rndt.values, C1, C2, kappa1, kappa2, epsilon))

# %%
pl.plot(res.x)
pl.plot(forcing_initial_guess)

# %%
ebm = EnergyBalanceModel(
    ocean_heat_capacity=[C1, C2],
    ocean_heat_transfer=[kappa1, kappa2],
    deep_ocean_efficacy=epsilon
)
ebm.add_forcing(res.x, timestep=1)
ebm.run()

# %%
pl.plot(ebm.temperature[:, 0])
pl.plot(test_df.tas.values)
pl.plot(ebm.toa_imbalance)
pl.plot(test_df.rndt.values)

# %%
pl.plot(ebm.temperature[:, 0] - test_df.tas.values)
pl.plot(ebm.toa_imbalance - test_df.rndt.values)

# %%
