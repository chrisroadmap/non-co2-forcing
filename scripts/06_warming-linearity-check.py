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
# # Are CO2 and non-CO2 components of warming linear?
#
# yes they are. Move to notebook 4. Also move the scenario stuff to notebook 4

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
model = 'CESM2'

# %%
data[model]['historical']['r1i1p1f1'].keys()

# %%
for model in data:
    for run in data[model]['historical']:
        C1 = ebm_df.loc[ebm_df['model']==model, 'C1'].values[0]
        C2 = ebm_df.loc[ebm_df['model']==model, 'C2'].values[0]
        kappa1 = ebm_df.loc[ebm_df['model']==model, 'kappa1'].values[0]
        kappa2 = ebm_df.loc[ebm_df['model']==model, 'kappa2'].values[0]
        epsilon = ebm_df.loc[ebm_df['model']==model, 'epsilon'].values[0]
        
        ebm = EnergyBalanceModel(
            ocean_heat_capacity=[C1, C2],
            ocean_heat_transfer=[kappa1, kappa2],
            deep_ocean_efficacy=epsilon,
            gamma_autocorrelation=1000
        )
        ebm.add_forcing(data[model]['historical'][run]['forcing_co2'], timestep=1)
        ebm.run()
        data[model]['historical'][run]['temperature_co2'] = ebm.temperature[:, 0]
    
        ebm = EnergyBalanceModel(
            ocean_heat_capacity=[C1, C2],
            ocean_heat_transfer=[kappa1, kappa2],
            deep_ocean_efficacy=epsilon,
            gamma_autocorrelation=1000
        )
        ebm.add_forcing(data[model]['historical'][run]['forcing_nonco2'], timestep=1)
        ebm.run()
        data[model]['historical'][run]['temperature_nonco2'] = ebm.temperature[:, 0]

# %%
for model in data:
    for run in data[model]['historical']:
        print(
            (
                data[model]['historical'][run]['temperature_fitted'] - 
                data[model]['historical'][run]['temperature_nonco2'] -  
                data[model]['historical'][run]['temperature_co2']
            ).sum()
        )
