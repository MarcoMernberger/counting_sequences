#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""ngmerge.py: Contains a Wrapper for NGmerge."""

import pypipegraph2 as ppg
import subprocess
from pathlib import Path
from typing import Optional, List, Dict
from pypipegraph2 import Job


__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


def ngmerge(
    merged_output_file: Path,
    r1: Path,
    r2: Path,
    merged_failed_prefix: Optional[str] = None,
    dependencies: List[Job] = [],
    options: Dict[str, str] = {},
) -> Job:
    """
    Wrapper for NGmerge.

    Returns a pipegraph job that creates a new merged fastq file from r1 and r2
    input fastqs. In addition, non-overlapping reads or ambiguously merged
    reads are discarded and written to a separate file according to the NGmerge
    manual.


    Parameters
    ----------
    merged_output_file : Path
        Output file path for the new merged single end fastq file.
    r1 : Path
        Path to R1 file.
    r2 : Path
        Path to R2 file.
    merged_failed_prefix : str, optional
        Prefix for fastq files containing reads that failed to merge, by default None.
    dependencies : List[Job], optional
        List of dependencies, by default [].
    options : Dict[str, str], optional
        Additional options to pass to ngmerge, by default {}.

    Returns
    -------
    Job
        MultiFileGeneratingJob that creates merged/failed fastq and log files.
    """
    merged_output_file.parent.mkdir(parents=True, exist_ok=True)
    deps = dependencies
    deps.append(ppg.ParameterInvariant(f"PI_{merged_output_file}", list(options)))
    merged_prefix = merged_failed_prefix
    if merged_failed_prefix is None:
        merged_prefix = merged_output_file.stem + "_failed"

    def merge(output_file):
        cmd = [
            "NGmerge",
            "-1",
            str(r1),
            "-2",
            str(r2),
            "-s",
            "-o",
            str(merged_output_file),
            "-f",
            merged_prefix,
        ]
        for k, v in options.items():
            if v == "":
                cmd.append(k)
            else:
                cmd.extend([k, v])
        try:
            subprocess.check_call(cmd)
        except subprocess.CalledProcessError:
            print(" ".join(cmd))
            raise

    job = ppg.FileGeneratingJob(merged_output_file, merge).depends_on(deps)
    return job
