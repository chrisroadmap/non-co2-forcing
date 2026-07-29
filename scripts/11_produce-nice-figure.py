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
import pickle

import matplotlib
import matplotlib.pyplot as pl
import numpy as np
import pandas as pd

# %%
pl.style.use('../defaults.mplstyle')

# %%
with open('../output/results_scenarios.pickle', 'rb') as handle:
    data = pickle.load(handle)

# %%
# delete problematic scenarios
del data['E3SM-1-0']['historical']['r1i1p1f1']
del data['E3SM-1-0']['historical']['r4i1p1f1']
del data['EC-Earth3-Veg']['historical']['r10i1p1f1']
del data['MPI-ESM-1-2-HAM']['ssp370']
del data['BCC-ESM1']['ssp370']

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
for scenario in ['ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585']:
    for run in data['CanESM5'][scenario]:
        pl.plot(data['CanESM5'][scenario][run]['temperature_fitted'], color=colors[scenario], ls='--', alpha=1)
        pl.plot(data['CanESM5'][scenario][run]['temperature_truth'], color=colors[scenario], ls='-', alpha=1)

# %%
df_can = pd.DataFrame(index=np.arange(1850,2101))
#for run in data['CanESM5'][scenario]:
df_can = pd.concat([pd.Series(data['CanESM5'][scenario][run]['forcing_nonco2'], name=run, index=np.arange(1850,2101)) for run in data['CanESM5'][scenario]], axis=1)

# %%
np.min(df_can, axis=1)

# %%
df_can = pd.concat([pd.Series(data['CanESM5']['historical'][run]['forcing_nonco2'], name=run, index=np.arange(1850,2015)) for run in data['CanESM5']['historical']], axis=1)
pl.fill_between(
    np.arange(1850, 2015),
    np.min(df_can, axis=1),
    np.max(df_can, axis=1),
    color=colors['historical'],
    alpha=0.15,
    lw=0,
)
pl.plot(
    np.arange(1850, 2015),
    np.mean(df_can, axis=1),
    color=colors['historical'],
)

for scenario in ['ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585']:
    df_can = pd.concat([pd.Series(data['CanESM5'][scenario][run]['forcing_nonco2'], name=run, index=np.arange(1850,2101)) for run in data['CanESM5'][scenario]], axis=1)
    pl.fill_between(
        np.arange(2014, 2101),
        np.min(df_can, axis=1)[164:],
        np.max(df_can, axis=1)[164:],
        color=colors[scenario],
        alpha=0.15,
        lw=0,
    )
    pl.plot(
        np.arange(2014, 2101),
        np.mean(df_can, axis=1)[164:],
        color=colors[scenario],
    )
pl.ylim(-2.5, 2.5)
pl.xlim(1850, 2100)

# %%
sorted(data)

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
fig, ax = pl.subplots(7, 7, figsize=(18/2.54, 22/2.54))

for idx, model in enumerate(sorted(data)):
    ax[idx//7, idx%7].set_title(model, fontsize=7)
    df_hist = pd.concat([pd.Series(data[model]['historical'][run]['forcing_nonco2'][:165], name=run, index=np.arange(1850,2015)) for run in data[model]['historical']], axis=1)
    if model=='CanESM5':
        label = fullnames['historical']
    else:
        label = None
    ax[idx//7, idx%7].fill_between(
        np.arange(1850, 2015),
        np.min(df_hist, axis=1),
        np.max(df_hist, axis=1),
        color=colors['historical'],
        alpha=0.15,
        lw=0,
    )
    ax[idx//7, idx%7].plot(
        np.arange(1850, 2015),
        np.mean(df_hist, axis=1),
        color=colors['historical'],
        lw=0.5,
        label=label
    )

    for scenario in ['ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585']:
        if scenario not in data[model]:
            continue
        if model=='CanESM5':
            label = fullnames[scenario]
        else:
            label = None
        df_future = pd.concat([pd.Series(data[model][scenario][run]['forcing_nonco2'], name=run, index=np.arange(1850,len(data[model][scenario][run]['forcing_nonco2'])+1850)) for run in data[model][scenario]], axis=1)
        ax[idx//7, idx%7].fill_between(
            np.arange(2014, 1850+len(df_future)),
            np.min(df_future, axis=1)[164:],
            np.max(df_future, axis=1)[164:],
            color=colors[scenario],
            alpha=0.15,
            lw=0,
        )
        ax[idx//7, idx%7].plot(
            np.arange(2014, 1850+len(df_future)),
            np.mean(df_future, axis=1)[164:],
            color=colors[scenario],
            lw=0.5,
            label=label
        )
    ax[idx//7, idx%7].set_ylim(-3, 3)
    ax[idx//7, idx%7].set_xlim(1850, 2100)


    if idx%7==0:
        ax[idx//7, idx%7].set_ylabel('W m$^{-2}$')
#     ax[idx//7, idx%7].set_xlim(0, 150)
#     ax[idx//7, idx%7].set_ylim(0, 8)
    ax[idx//7, idx%7].set_xticks([1850, 1975, 2100])
    xTick_objects = ax[idx//7, idx%7].xaxis.get_major_ticks()
    xTick_objects[0].label1.set_horizontalalignment('left')
    xTick_objects[-1].label1.set_horizontalalignment('right')
    ax[idx//7, idx%7].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(5))

ax[6, 3].axis('off')
ax[6, 4].axis('off')
ax[6, 5].axis('off')
ax[6, 6].axis('off')

l = {}
for scenario in ['historical', 'ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585']:
    l[scenario] = ax[6, 3].plot(0, 0, color=colors[scenario], lw=0.5, label=label)

ax[6, 3].legend([val[0] for val in l.values()], [fullnames[key] for key in l.keys()], loc='upper left', frameon=False)

fig.tight_layout()
pl.savefig('../plots/non-co2-timeseries.png')

# %%
