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

# %% [markdown]
# # Check calibration and EBM against 1pctCO2 run

# %%
import os
import glob

import matplotlib.ticker
import matplotlib.pyplot as pl
import numpy as np
import pandas as pd

from fair.energy_balance_model import EnergyBalanceModel
from scipy.optimize import minimize

# %%
pl.style.use('../defaults.mplstyle')

# %%
x0=284.316999854786

def myhre(x, alpha):
    return alpha * np.log(x/x0)

co2_conc_1pct = 1.01 ** np.arange(151) * x0

# %%
ebm_df = pd.read_csv('../output/calibrations/4xCO2_cummins_ebm2_cmip6.csv')

# %%
alpha = ebm_df['F_4xCO2'] / np.log(4)

# %%
ebm_df['alpha'] = alpha

# %%
ebm_df

# %%
myhre(1.01 ** np.arange(151) * x0, ebm_df.loc[ebm_df['model']=='UKESM1-0-LL', 'alpha'].values[0])

# %%
# we don't have 1pctCO2 or scenario projections from KACE or FIO, so delete them now
ebm_df.drop(ebm_df["model"].isin(["KACE-1-0-G", "FIO-ESM-2-0"]).index)

# %%
ebm_df_sorted = ebm_df.sort_values(by=['model']).reset_index(drop=True)

# %%
fig, ax = pl.subplots(7, 7, figsize=(18/2.54, 22/2.54))

for idx, row in ebm_df_sorted.iterrows():
    try:
        test_df = pd.read_csv(f'../data/cmip6-hbf/cmip_data/{row.model}/1pctCO2/{row.model}_1pctCO2_{row.run}_anomalies.csv', index_col=0)
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

        ax[idx//7, idx%7].plot(test_df.tas.values, color='r', ls='--', label='ESM', lw=1)
        ax[idx//7, idx%7].plot(ebm.temperature[:, 0], color='k', label='EBM fit', lw=1)
        
        ax[idx//7, idx%7].set_title(row.model, fontsize=7)
    except:
        ax[idx//7, idx%7].set_title(row.model, fontsize=7)
        ax[idx//7, idx%7].axis('off')
        print(row.model)

    if idx==1:
        ax[idx//7, idx%7].legend(frameon=False)
    if idx%7==0:
        ax[idx//7, idx%7].set_ylabel('°C')
    ax[idx//7, idx%7].set_xlim(0, 150)
    ax[idx//7, idx%7].set_ylim(0, 8)
    ax[idx//7, idx%7].set_xticks([0, 75, 150])
    xTick_objects = ax[idx//7, idx%7].xaxis.get_major_ticks()
    xTick_objects[0].label1.set_horizontalalignment('left')
    xTick_objects[-1].label1.set_horizontalalignment('right')
    ax[idx//7, idx%7].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(5))
   

fig.tight_layout()
pl.savefig('../plots/1pctCO2-validation-tas.png')

# %% [markdown]
# Notes:
#
# exclude the following scenarios
#
# - FIO-ESM-2-0 (missing 1pctCO2)
# - KACE-1-0-G (missing 1pctCO2)
# - GISS-E2-1-G (looks like 1pctCO2 is held after 70 years)
# - KIOST-ESM (doesn't seem to be in balance in year zero)

# %%
fig, ax = pl.subplots(7, 7, figsize=(18/2.54, 22/2.54))

for idx, row in ebm_df_sorted.iterrows():
    try:
        test_df = pd.read_csv(f'../data/cmip6-hbf/cmip_data/{row.model}/1pctCO2/{row.model}_1pctCO2_{row.run}_anomalies.csv', index_col=0)
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
    
        ax[idx//7, idx%7].plot(test_df.rndt.values, color='g', ls='--', label='ESM', lw=1)
        ax[idx//7, idx%7].plot(ebm.toa_imbalance, color='k', label='EBM fit', lw=1)
        
        ax[idx//7, idx%7].set_title(row.model, fontsize=7)
    except:
        ax[idx//7, idx%7].set_title(row.model, fontsize=7)
        ax[idx//7, idx%7].axis('off')
        print(row.model)

    if idx==3:
        ax[idx//7, idx%7].legend(frameon=False)
    if idx%7==0:
        ax[idx//7, idx%7].set_ylabel('W m$^{-2}$')
    ax[idx//7, idx%7].set_xlim(0, 150)
    ax[idx//7, idx%7].set_ylim(-1, 4)
    ax[idx//7, idx%7].set_xticks([0, 75, 150])
    xTick_objects = ax[idx//7, idx%7].xaxis.get_major_ticks()
    xTick_objects[0].label1.set_horizontalalignment('left')
    xTick_objects[-1].label1.set_horizontalalignment('right')
    ax[idx//7, idx%7].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(5))
   

fig.tight_layout()
pl.savefig('../plots/1pctCO2-validation-rndt.png')

# %%
