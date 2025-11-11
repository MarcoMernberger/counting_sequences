#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""plots.py: Contains plots to visualize counter results."""
import matplotlib
import matplotlib.pyplot as plt

# import pypipegraph as ppg
import pandas as pd
import numpy as np
from typing import List
from pathlib import Path
from pypipegraph import Job, MultiFileGeneratingJob, FileGeneratingJob
from mbf.align import Sample
from .util import get_reads_for_lanes_callable

# from mplots import MPPlotJob


__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


def plot_count_corr(
    csv1: Path,
    csv2: Path,
    dependencies: List[Job] = [],
    index_column: str = "Name",
    count_Column: str = "Read Count",
):
    outfile = Path("results") / f"{csv1.name} vs {csv2.name}.png"
    print(outfile)

    def __plot():
        df1 = pd.read_csv(csv1, sep="\t")
        df2 = pd.read_csv(csv2, sep="\t")
        df1["index"] = [
            f"{sample} {name}"
            for sample, name in zip(df1["Sample"], df1[index_column])
        ]
        df2["index"] = [
            f"{sample} {name}"
            for sample, name in zip(df2["Sample"], df2[index_column])
        ]
        df1 = df1.set_index("index")
        df2 = df2.set_index("index")
        print(df1.shape, df2.shape)
        print(len(df1[count_Column]), len(df2.loc[df1.index][count_Column]))
        print(set(df1.index).difference(set(df2.index)))
        print(set(df2.index).difference(set(df1.index)))
        print(len(set(df1.index).intersection(set(df2.index))))
        print(df1[count_Column])
        print(df2.loc[df1.index][count_Column])
        f = plt.figure()
        plt.plot(
            df1[count_Column],
            df2.loc[df1.index][count_Column],
            marker=".",
            ls="",
        )
        plt.title(f"Correlation of counts {csv1.name} vs {csv2.name}")
        plt.xlabel(f"Count {csv1.name}")
        plt.xlabel(f"Count {csv2.name}")
        plt.tight_layout()
        f.savefig(outfile)

    return FileGeneratingJob(outfile, __plot).depends_on(dependencies)


def plot_count_corr_m(
    csv1: Path,
    csv2: Path,
    dependencies: List[Job] = [],
    index_column: str = "Name",
    count_Column: str = "Read Count",
    sample_column: str = "Sample",
):
    outfile = Path("results") / f"{csv1.name} vs {csv2.name}.png"
    print(outfile)

    def __plot():
        df1 = pd.read_csv(csv1, sep="\t")
        df2 = pd.read_csv(csv2, sep="\t")
        df1["index"] = [
            f"{sample} {name}"
            for sample, name in zip(df1["Sample"], df1[index_column])
        ]
        df2["index"] = [
            f"{sample} {name}"
            for sample, name in zip(df2["Sample"], df2[index_column])
        ]
        df1 = df1.set_index("index")
        df2 = df2.set_index("index")
        print(df1.shape, df2.shape)
        raise ValueError()
        f = plt.figure()
        plt.plot(
            df1[count_Column],
            df2.loc[df1.index][count_Column],
            marker=".",
            ls="",
        )
        plt.title(f"Correlation of counts {csv1.name} vs {csv2.name}")
        plt.xlabel(f"Count {csv1.name}")
        plt.xlabel(f"Count {csv2.name}")
        plt.tight_layout()
        f.savefig(outfile)

    return MultiFileGeneratingJob(outfile, __plot).depends_on(dependencies)


def plot_reads_for_lanes(
    outfile: Path, lanes: List[Sample], dependencies: List[Job] = []
) -> Job:
    df_calc_func = get_reads_for_lanes_callable(lanes, dependencies)

    def __plot(df):
        df_plot = df.set_index("Sample")
        f = plt.figure(figsize=(10, 10))
        x = np.arange(df_plot.shape[0])
        plt.bar(x, df_plot["Reads"])
        plt.xticks(x, df_plot.index, rotation=75)
        plt.ylabel("Read count")
        plt.xlabel("Lane")
        plt.title("Read counts for lanes")
        plt.tight_layout()
        return f

    return MPPlotJob(
        outfile, calc_function=df_calc_func, plot_function=__plot
    ).depends_on(dependencies)
