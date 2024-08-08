# -*- coding: utf-8 -*-
from pkg_resources import get_distribution, DistributionNotFound

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = __name__
    __version__ = get_distribution(dist_name).version
except DistributionNotFound:
    __version__ = "unknown"
finally:
    del get_distribution, DistributionNotFound

from .counting import SequenceCounter
from .adapter import CutadaptMatch, Paired_Filtered_Trimmed_From_Job
from .plots import plot_count_corr, plot_reads_for_lanes
from .util import (
    generate_stitched_fastq,
    get_reads_for_lanes_callable,
    get_reads_for_lanes_df,
)
from ngmerge import ngmerge
