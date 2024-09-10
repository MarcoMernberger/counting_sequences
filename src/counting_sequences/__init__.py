# -*- coding: utf-8 -*-

__version__ = '1.0.0'


from .counting import SequenceCounter, SequenceCounter2
from .adapter import CutadaptMatch, Paired_Filtered_Trimmed_From_Job
from .plots import plot_count_corr, plot_reads_for_lanes
from .util import (
    get_reads_for_lanes_callable,
    get_reads_for_lanes_df,
)
from .ngmerge import ngmerge
