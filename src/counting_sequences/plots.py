#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""plots.py: Contains plots to visualize counter results."""
import matplotlib
import matplotlib.pyplot as plt
import pypipegraph as ppg
import pandas as pd
import pypipegraph as ppg
from typing import List
from pathlib import Path

__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


def plot_count_corr(csv1: Path, csv2: Path, dependencies: List[ppg.Job] = []):
    outfile = Path("results") / f"{csv1.name} vs {csv2.name}.png"
    print(outfile)

    def __plot():
        df1 = pd.read_csv(csv1, sep="\t")
        df2 = pd.read_csv(csv2, sep="\t")
        df1["index"] = [
            f"{sample} {name}" for sample, name in zip(df1["Sample"], df1["Name"])
        ]
        df2["index"] = [
            f"{sample} {name}" for sample, name in zip(df2["Sample"], df2["Name"])
        ]
        df1 = df1.set_index("index")
        df2 = df2.set_index("index")
        f = plt.figure()
        plt.plot(df1["Read Count"], df2.loc[df1.index]["Read Count"], marker=".", ls="")
        plt.title(f"Correlation of counts {csv1.name} vs {csv2.name}")
        plt.xlabel(f"Count {csv1.name}")
        plt.xlabel(f"Count {csv2.name}")
        plt.tight_layout()
        f.savefig(outfile)

    return ppg.FileGeneratingJob(outfile, __plot).depends_on(dependencies)
