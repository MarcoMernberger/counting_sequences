#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""ngmerge.py: Contains a Wrapper for NGmerge."""

import pypipegraph2 as ppg
import subprocess
from pathlib import Path
from typing import List, Dict
from pypipegraph2 import Job


__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


def ngmerge(
    output_file: Path,
    r1: Path,
    r2: Path,
    dependencies: List[Job] = [],
    options: Dict[str, str] = {},
):
    """
    generate_stitched_fastq wrapper for ngmerge.

    Returns a pipegraph job that creates a new merged fastq file from r1 and r2
    input fastqs. In addition, non-overlapping reads or ambiguously merged
    reads are discarded and written to a separate file according to the NGmerge
    manual.

    Parameters
    ----------
    output_file : Path
        Output file path for the new fastq file.
    r1 : Path
        Path to R1 file.
    r2 : Path
        Path to R2 file.
    dependencies : List[Job], optional
        List of dependencies, by default [].
    options : Dict[str, str], optional
        Additional options to pass to ngmerge, by default {}.

    Returns
    -------
    Job
        FileGeneratingJob that creates the merged bam file.
    """
    output_file.parent.mkdir(parents=True, exist_ok=True)
    deps = dependencies
    deps.append(ppg.ParameterInvariant(f"PI_{output_file}", list(options)))

    def __dump(output_file):
        if not output_file.exists():
            cmd = [
                "NGmerge",
                "-1",
                str(r1),
                "-2",
                str(r2),
                "-s",
                "-o",
                str(output_file),
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

    job = ppg.FileGeneratingJob(output_file, __dump).depends_on(deps)
    return job
