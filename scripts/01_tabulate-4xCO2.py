#!/usr/bin/env python
# coding: utf-8

"""Combine Hege-Beate's data into one file."""

import glob
import os
from pathlib import PurePath

import pandas as pd

print("Making nice 4xCO2 data...")

available_files = glob.glob(
    "../data/cmip6-hbf/cmip_data/*/abrupt-4xCO2/"
    "*_abrupt-4xCO2_*_anomalies.txt"
)

maxlen = 0

# where we have more than one ensemble member in a model, there is usually only one well-
# behaving ensemble member
multi_runs = {
    "GISS-E2-1-G": "r1i1p1f1",
    "GISS-E2-1-H": "r1i1p1f1",
    "MRI-ESM2-0": "r1i1p1f1",
    "EC-Earth3": "r3i1p1f1",
    "FIO-ESM-2-0": "r1i1p1f1",
    "CanESM5": "r1i1p2f1",
    "FGOALS-f3-L": "r1i1p1f1",
    "CNRM-ESM2-1": "r1i1p1f2",
}

models = []
runs = []
lines = []
for file in available_files:
    model = PurePath(file).parts[4]
    run = PurePath(file).parts[6].split("_")[2]
    if model in multi_runs:
        if run != multi_runs[model]:
            continue
    models.append(model)
    runs.append(run)
    df = pd.read_csv(file, index_col=0)
    vars = {}
    for var in ["tas", "rlut", "rsut", "rsdt"]:
        vars[var] = df[var].values.squeeze()
        if len(vars[var]) > maxlen:
            maxlen = len(vars[var])
        line = [
            "CMIP",
            model,
            run,
            "CMIP6",
            "unspecified",
            "World",
            "abrupt-4xCO2",
            "W m^-2",
            var,
        ]
        line.extend(vars[var][:150])
        lines.append(line)
    vars["rndt"] = vars["rsdt"] - vars["rsut"] - vars["rlut"]
    line = [
        "CMIP",
        model,
        run,
        "CMIP6",
        "unspecified",
        "World",
        "abrupt-4xCO2",
        "W m^-2",
        "rndt",
    ]
    line.extend(vars["rndt"][:150])
    lines.append(line)

df = pd.DataFrame(
    lines,
    columns=(
        [
            "activity_id",
            "climate_model",
            "member_id",
            "mip_era",
            "model",
            "region",
            "scenario",
            "unit",
            "variable",
        ]
        + ["X%d" % year for year in range(1850, 2000)]
    ),
)

to_remove = [f"X{year}" for year in range(1850, 2000)]
df.dropna(subset=to_remove, inplace=True)

os.makedirs(
    "../output/calibrations/",
    exist_ok=True,
)

df.to_csv(
    "../output/calibrations/4xCO2_cmip6.csv",
    index=False,
)
