#!/usr/bin/env python
# -*- coding: utf-8 -*-
from pathlib import Path
from counting_sequences import SequenceCounter, CutadaptMatch
import pypipegraph as ppg
import pandas as pd
import pytest
import collections
from counting_sequences import CutadaptMatch

__author__ = "MarcoMernberger"
__copyright__ = "MarcoMernberger"
__license__ = "mit"


def test_cut_adapt_filter():
    forward_read = b"TTGCTTTATCTGTTCACTTGTGCCCTGACTTTCAAthisiswhatshouldremainCTCTGTNTCCTTCCTCTTCCTCCTTTCCTCTTCTTCCAGCTCCCTNATGTTTTCTTCATCGGTACGGATTGCACAGTTCCTCCNNTCTTCCGCTTCCTCGCCACCGCCCTCAACCCGCAGCCACCGCACATCACTTATCCGTTGATTCCCGACCCTTACTTTGATTGTTCCTCCTATACCTACTATTACCCGCTGGGGCCGTAGAGTCTATTCTTC"
    reverse_read = b"CAGTGAGGAATCAGAGGNCTCCGGACCCTGGGCAACCAGCCCTGthisiswhatshouldremainTCGTCTCTCCAGCCCCAGCTGCTCACCATCTCTCTCTCAGCAGCGCTCATGGTGGTGGCCGCGCCTCACACCCTCCTTCTTTTTCTGTTACTGCTTGTAGTTTGCCATTGCGCGGTTGCTTGTTCCGTGCGTTTTTGTTTTATCCTCCCCCCTCCTTTTTTCTTCTGTTTTTTTTGTTTTTCCTCTTTTTTTTTTTTGTTTGTTTTT"
    forward_adapter = "CCTGACTTTCAA"
    reverse_adapter = "CAGGGCTGGTTG"
    forward_adapter_rc = "TTGAAAGTCAGG"
    reverse_adapter_rc = "CAACCAGCCCTG"
    min_overlap = 5
    cutadapt = CutadaptMatch(forward_adapter, reverse_adapter, min_overlap=min_overlap)
    filter_func = cutadapt.filter()

    # test example
    filtered = filter_func(
        forward_read, forward_read, "r1", reverse_read, reverse_read, "r2"
    )
    r1_trimmed, r2_trimmed = filtered[0].decode(), filtered[3].decode()
    assert r1_trimmed.startswith("thisiswhatshouldremain")
    assert r2_trimmed.startswith("thisiswhatshouldremain")

    # test full adapters flanking
    forward_read = b"TTGCTTTATCTGTTCACTTGTGCCCTGACTTTCAAthisiswhatshouldremainCAGGGCTGGTTGCTCTGTNTCCTTCCTCTTCCTCCTTTCCTCTTCTTCCAGCTCCCTNATGTTTTCTTCATCGGTACGGATTGCACAGTTCCTCCNNTCTTCCGCTTCCTCGCCACCGCCCTCAACCCGCAGCCACCGCACATCACTTATCCGTTGATTCCCGACCCTTACTTTGATTGTTCCTCCTATACCTACTATTACCCGCTGGGGCCGTAGAGTCTATTCTTC"
    reverse_read = b"CAGTGAGGAATCAGAGGNCTCCGGACCCTGGGCAACCAGCCCTGthisiswhatshouldremainTTGAAAGTCAGGTCGTCTCTCCAGCCCCAGCTGCTCACCATCTCTCTCTCAGCAGCGCTCATGGTGGTGGCCGCGCCTCACACCCTCCTTCTTTTTCTGTTACTGCTTGTAGTTTGCCATTGCGCGGTTGCTTGTTCCGTGCGTTTTTGTTTTATCCTCCCCCCTCCTTTTTTCTTCTGTTTTTTTTGTTTTTCCTCTTTTTTTTTTTTGTTTGTTTTT"
    filtered = filter_func(
        forward_read, forward_read, "r1", reverse_read, reverse_read, "r2"
    )
    print(filtered)
    r1_trimmed, r2_trimmed = filtered[0].decode(), filtered[3].decode()
    assert r1_trimmed == r2_trimmed == "thisiswhatshouldremain"
    # test partial adapters flanking forward
    forward_read = b"ACTTTCAAthisiswhatshouldremainCAGGGCTGGTTGCTCTGTNTCCTTCCTCTTCCTCCTTTCCTCTTCTTCCAGCTCCCTNATGTTTTCTTCATCGGTACGGATTGCACAGTTCCTCCNNTCTTCCGCTTCCTCGCCACCGCCCTCAACCCGCAGCCACCGCACATCACTTATCCGTTGATTCCCGACCCTTACTTTGATTGTTCCTCCTATACCTACTATTACCCGCTGGGGCCGTAGAGTCTATTCTTC"
    reverse_read = (
        b"CAGTGAGGAATCAGAGGNCTCCGGACCCTGGGCAACCAGCCCTGthisiswhatshouldremainTTGAAAGT"
    )

    filtered = filter_func(
        forward_read, forward_read, "r1", reverse_read, reverse_read, "r2"
    )
    print(filtered)
    r1_trimmed, r2_trimmed = filtered[0].decode(), filtered[3].decode()
    assert r1_trimmed == r2_trimmed == "thisiswhatshouldremain"

    # test partial adapters flanking both
    forward_read = b"ACTTTCAAthisiswhatshouldremainCAGGGC"
    reverse_read = b"AGCCCTGthisiswhatshouldremainTTGAAAGT"
    filtered = filter_func(
        forward_read, forward_read, "r1", reverse_read, reverse_read, "r2"
    )
    print(filtered)
    r1_trimmed, r2_trimmed = filtered[0].decode(), filtered[3].decode()
    assert r1_trimmed == r2_trimmed == "thisiswhatshouldremain"

    # test partial adapters too short
    forward_read = b"TCAAthisiswhatshouldremainCAGGGC"
    reverse_read = b"AGCCCTGthisiswhatshouldremainTTG"
    filtered = filter_func(
        forward_read, forward_read, "r1", reverse_read, reverse_read, "r2"
    )
    assert filtered is None

    # test end adapter missing
    forward_read = b"ACTTTCAAthisiswhatshouldremain"
    reverse_read = b"AGCCCTGthisiswhatshouldremain"
    filtered = filter_func(
        forward_read, forward_read, "r1", reverse_read, reverse_read, "r2"
    )
    print(filtered)
    r1_trimmed, r2_trimmed = filtered[0].decode(), filtered[3].decode()
    assert r1_trimmed == r2_trimmed == "thisiswhatshouldremain"

    # test start adapter missing in r1
    forward_read = b"thisiswhatshouldremainCAGGGCTGGTTGCTCTGTNTCCTTCCTCTTCCTCCTTTCCTCTTCTTCCAGCTCCCTNATGTTTTCTTCATCGGTACGGATTGCACAGTTCCTCCNNTCTTCCGCTTCCTCGCCACCGCCCTCAACCCGCAGCCACCGCACATCACTTATCCGTTGATTCCCGACCCTTACTTTGATTGTTCCTCCTATACCTACTATTACCCGCTGGGGCCGTAGAGTCTATTCTTC"
    reverse_read = b"GCAACCAGCCCTGthisiswhatshouldremainTTGAAAGT"
    filtered = filter_func(
        forward_read, forward_read, "r1", reverse_read, reverse_read, "r2"
    )
    assert filtered is None

    # test start adapter missing in r2
    forward_read = b"CCTGACTTTCAAthisiswhatshouldremainCAGGGCTGGTTGCTCTGTNTCCTTCCTCTTCCTCCTTTCCTCTTCTTCCAGCTCCCTNATGTTTTCTTCATCGGTACGGATTGCACAGTTCCTCCNNTCTTCCGCTTCCTCGCCACCGCCCTCAACCCGCAGCCACCGCACATCACTTATCCGTTGATTCCCGACCCTTACTTTGATTGTTCCTCCTATACCTACTATTACCCGCTGGGGCCGTAGAGTCTATTCTTC"
    reverse_read = b"thisiswhatshouldremainTTGAAAGT"
    filtered = filter_func(
        forward_read, forward_read, "r1", reverse_read, reverse_read, "r2"
    )
    assert filtered is None
    print("anoter")
    # test no empty strings
    forward_read = b"CCTGACTTTCAACAGGGCTGGTTGCTCTGTNTCCTTCCTCTTCCTCCTTTCCTCTTCTTCCAGCTCCCTNATGTTTTCTTCATCGGTACGGATTGCACAGTTCCTCCNNTCTTCCGCTTCCTCGCCACCGCCCTCAACCCGCAGCCACCGCACATCACTTATCCGTTGATTCCCGACCCTTACTTTGATTGTTCCTCCTATACCTACTATTACCCGCTGGGGCCGTAGAGTCTATTCTTC"
    reverse_read = b"CAGTGAGGAATCAGAGGNCTCCGGACCCTGGGCAACCAGCCCTGTTGAAAGTCAGG"
    filtered = filter_func(
        forward_read, forward_read, "r1", reverse_read, reverse_read, "r2"
    )
    print(filtered)
    assert filtered is None
