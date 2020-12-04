#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Dummy conftest.py for counting_sequences.

    If you don't know what this is for, just leave it empty.
    Read more about conftest.py under:
    https://pytest.org/latest/plugins.html
"""
import sys
import subprocess
import pathlib
import pytest

from pypipegraph.testing.fixtures import (  # noqa:F401
    new_pipegraph,
    pytest_runtest_makereport,
)
from mbf_qualitycontrol.testing.fixtures import new_pipegraph_no_qc  # noqa:F401
from pypipegraph.testing import force_load
from pathlib import Path
# import pytest
root = pathlib.Path(__file__).parent.parent
print("root", root)
sys.path.append(str(root / "src"))
print("the path is", sys.path)
subprocess.check_call(["python3", "setup.py", "build_ext", "-i"], cwd=root)


class DummyLane:

    def __init__(self, name, filename):
        self.name = name
        self.filename = filename

    def get_aligner_input_filenames(self):
        return [self.filename, None]        

    def prepare_input(self):
        return []

        
@pytest.fixture
def raw_lane():
    return DummyLane("TestSample", Path(__file__).parent.parent / "data" / "test.fastq")
