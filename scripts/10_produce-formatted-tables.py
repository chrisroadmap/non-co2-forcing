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
import pandas as pd

# %%
df = pd.read_csv('../output/calibrations/4xCO2_cummins_ebm2_cmip6.csv')

# %%
df_sorted = df.sort_values(by=["model"])
df_sorted


# %%
def significant(x, n=2):
    if x == 0:
        res = x
    else:
        res = round(x, -int(np.floor(np.log10(abs(x)))) + (n - 1))
    if res >= 10**n:
        res = int(res)
    return res


# %%
round(8983, 2)

# %%
print(r"\begin{tabular}{r|rrrrrrr}")
print(r"Model & $C_1$ & $C_2$ & $\kappa_1$ & $\kappa_2$ & $\epsilon$ & $F_{4\times\mathrm{CO}_2}$ & $\alpha$ \\")
print(r"\midrule")

for idx, row in df_sorted.iterrows():
    print(rf"{row.model} & {row.C1:.1f} & {row.C2:.0f} & {row.kappa1:.2f} & {row.kappa2:.2f} & {row.epsilon:.2f} & {row.F_4xCO2:.2f} & {row.F_4xCO2/np.log(4):.2f} \\") 

print(r"\end{tabular}")

# %%
