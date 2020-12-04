#!/usr/bin/env python
# -*- coding: utf-8 -*-
from pathlib import Path
from counting_sequences import SequenceCounter
import pypipegraph as ppg
import pandas as pd
import pytest
import collections

__author__ = "MarcoMernberger"
__copyright__ = "MarcoMernberger"
__license__ = "mit"


@pytest.fixture
def scouter(tmpdir):
    scouter = SequenceCounter(
        sequence_file_path=Path(__file__).parent.parent / "data" / "seq_in.csv",
        name=None,
        start_seq_to_trim="CCTCTT",
        trimmed_length=150,
        result_folder=Path(tmpdir) / "test"
    )
    return scouter


def test_init(scouter, tmpdir):
    assert scouter.name == "SC_CCTCTT_150"
    assert scouter.start_seq_to_trim == "CCTCTT"
    assert scouter.start_seq_len == 6
    assert scouter.trimmed_length == 150
    assert scouter.result_dir.exists()
    assert scouter.sequence_file_path.name == "seq_in.csv"


@pytest.mark.usefixtures("new_pipegraph_no_qc")
def test_write_trim_predefines(tmpdir, scouter):

    scouter.write_predefined_sequences()
    outputfile = scouter.result_dir / "predefined_sequences.tsv"
    ppg.run_pipegraph()
    df = pd.read_csv(outputfile, sep="\t")
    scouter.assert_predefined(df["Full Sequence"].values, df["Sequence"].values)
    assert outputfile.exists()
    df_new = pd.read_csv(outputfile, sep="\t")
    df_new.index = df_new["Name"]
    print(df.head())
    assert df_new.loc["1>A_test3"]["Duplicate"]
    assert df_new.loc["1>A_test4"]["Duplicate"]
    assert df_new.loc["1>A_test3"]["Deduplicated"]
    assert not df_new.loc["1>A_test4"]["Deduplicated"]
    assert df_new.loc["1>A_test3"]["Duplicate Entries"] == "1>A_test3;1>A_test4"


@pytest.mark.usefixtures("new_pipegraph_no_qc")
def test_write_trim_predefine_alert1(tmpdir, scouter):
    scouter = SequenceCounter(
        sequence_file_path=Path(__file__).parent.parent / "data" / "seq_fail1.csv",
        name=None,
        start_seq_to_trim="TTGCTTTACCTCCTTTT",
        trimmed_length=133,
        result_folder=Path(tmpdir) / "test"
    )
    df_sequence_df = pd.read_csv(scouter.sequence_file_path, sep="\t")
    df_sequence_df["Name"] = [
        f"{a}_{b}"
        for a, b in zip(
            df_sequence_df["Alteration"].astype(str),
            df_sequence_df["Effect"].astype(str),
        )
    ]
    index_function = scouter._find_index
    df_sequence_df = df_sequence_df.rename(columns={"Sequence": "Full Sequence"})
    sequences = []
    seen = collections.defaultdict(set)
    for _, row in df_sequence_df.iterrows():
        fullseq = row["Full Sequence"]
        seen[fullseq].add(row["Name"]) 
        index1 = index_function(fullseq)
        index2 = index1 + min(len(fullseq), scouter.trimmed_length)
        seq = fullseq[index1:index2]
        sequences.append(seq)
    df_sequence_df["Sequence"] = sequences
    duplicates = []
    dup_names = []
    for seq in df_sequence_df["Full Sequence"].values:
        dups = seen[seq]
        duplicates.append(len(dups) > 1)
        dup_names.append(",".join(list(dups)))
    df_sequence_df["Duplicate"] = duplicates
    df_sequence_df["Duplicate Entries"] = dup_names
    with pytest.raises(Exception):
        scouter.assert_predefined(df_sequence_df["Full Sequence"].values, df_sequence_df["Sequence"].values)


@pytest.mark.usefixtures("new_pipegraph_no_qc")
def test_write_trim_predefine_alert2(tmpdir, scouter):
    scouter = SequenceCounter(
        sequence_file_path=Path(__file__).parent.parent / "data" / "seq_fail2.csv",
        name=None,
        start_seq_to_trim="TTGCTTTACCTCCTTTT",
        trimmed_length=133,
        result_folder=Path(tmpdir) / "test"
    )
    df_sequence_df = pd.read_csv(scouter.sequence_file_path, sep="\t")
    df_sequence_df["Name"] = [
        f"{a}_{b}"
        for a, b in zip(
            df_sequence_df["Alteration"].astype(str),
            df_sequence_df["Effect"].astype(str),
        )
    ]
    index_function = scouter._find_index
    df_sequence_df = df_sequence_df.rename(columns={"Sequence": "Full Sequence"})
    sequences = []
    seen = collections.defaultdict(set)
    for _, row in df_sequence_df.iterrows():
        fullseq = row["Full Sequence"]
        seen[fullseq].add(row["Name"]) 
        index1 = index_function(fullseq)
        index2 = index1 + min(len(fullseq), scouter.trimmed_length)
        seq = fullseq[index1:index2]
        sequences.append(seq)
    df_sequence_df["Sequence"] = sequences
    duplicates = []
    dup_names = []
    for seq in df_sequence_df["Full Sequence"].values:
        dups = seen[seq]
        duplicates.append(len(dups) > 1)
        dup_names.append(",".join(list(dups)))
    df_sequence_df["Duplicate"] = duplicates
    df_sequence_df["Duplicate Entries"] = dup_names
    with pytest.raises(Exception):
        scouter.assert_predefined(df_sequence_df["Full Sequence"].values, df_sequence_df["Sequence"].values)


@pytest.mark.usefixtures("new_pipegraph_no_qc")
def test_write_fastq(tmpdir, raw_lane, scouter):
    scouter.write_fastq_count(raw_lane)
    ppg.run_pipegraph()
    output_file = scouter.result_dir / f"{raw_lane.name}_{scouter.name}_all_reads.tsv"
    assert output_file.exists()


@pytest.mark.usefixtures("new_pipegraph_no_qc")
def test_write_fastq_trimmed(tmpdir, raw_lane, scouter):
    scouter.write_fastq_count_trimmed(raw_lane)
    ppg.run_pipegraph()
    output_file = scouter.result_dir / f"{raw_lane.name}_{scouter.name}_all_reads.tsv"
    output_file2 = scouter.result_dir / f"{raw_lane.name}_{scouter.name}_all_reads_trimmed.tsv"
    assert output_file.exists()
    assert output_file2.exists()


@pytest.mark.usefixtures("new_pipegraph_no_qc")
def test_write_count_table(tmpdir, raw_lane, scouter):
    scouter.write_count_table(raw_lane, row_order=None)
    ppg.run_pipegraph()
    output_file = scouter.result_dir / f"{raw_lane.name}_{scouter.name}_sequence_count.tsv"
    output_file2 = scouter.result_dir / f"{raw_lane.name}_{scouter.name}_sequence_count_unmatched.tsv"
    assert output_file.exists()
    assert output_file2.exists()


@pytest.mark.usefixtures("new_pipegraph_no_qc")
def test_count_fastq(tmpdir, raw_lane, scouter):
    scouter.write_fastq_count(raw_lane)
    ppg.run_pipegraph()
    output_file = scouter.result_dir / f"{raw_lane.name}_{scouter.name}_all_reads.tsv"
    df = pd.read_csv(output_file, sep="\t")
    expected = {
        "CGTACAAGAGACAAGCAATCAGTGAGGAATCAGAGGCCTCCGGACCCTGGGCAACCAGCCCTGTCGTCTCTCCAGCCCCAGCTGCTCACCATCGCTATCTGAGCGCCACTCTTGTTGGGGCCAGCGCCTCCCACCTTCCCTCTTTTGCTTTACCTCCTTTTAGTTGGCCTTGCCCCGGCCCCGGTCCCTTGCCAAAATGTCTTGTTTAGCCCCGGGTGCTCCTGTCGGGTCTTGACTGATTCACACTTGATATTCTTGTCTTCTGGTTCTTGCTCTGATGAGCACACGTCTGCACTCTTGT":
        2,
        "CGTACAAGAGACAAGCAATCAGTGAGGAATCCCTCTTTTGCTTTACCTCCTTTTAGTTGGCCTTGCCCCGGCCCCGGTCCCTTGCCAAAATGTCTTGTTTAGCCCCGGGGTGCTCCTGTCGGGTCTTGACTGATTCACACTTGATATTCTTGTCTTCTGGTTCTTGCTCTGATGAGCACACGTCTGC":
        1,
        "CGTACAAGAGACAAGCAATCAGTGAGGAATCCCTCTTTTGCTTTACCTCCTTTTAGTTGGCCTTGCCCCGGCCCCGGTCCCTTGCCAAAATGTCTTGTTTAGCCCCGGGGTGCTCCTGTCGGGTCTTGACTGATTCACACTTGATATTCTTGTCTTCTGGTTCTTGCTCTGATGAGCACACGTCTGCACACACACA":
        1,
        "CGTACAAGAGACAAGCAATCAGTGAGGAATCCCTCTTTTGCTTTACCTCCTTTTAGTTGGCCTTGCCCCGGCCCCGGTCCCTTGCCAAAATGTCTTGTTTAGCCCCGGGGTGCTCCTGTCGGGTCTTGACTGATTCACACTTGATATTCTTGTCTTCTGGTTCTTGCTCTGATGAGCACACGTCTG":
        1,
        "AGGAATCCCTCTTTTGCTTTACCTCCTTTTAGCCTCTTTTGCCCCGGCCCCGGTCCCTTGCCAAAATGTCTTGTTTAGCCCCGGGTGCTCCTGTCGGGTCTTGACTGATTCACACTTGATATTCTTGTCTTCTGGTTCTTGCTCTGATGAGCACACGTCTGCAACCA":
        2,
        "AGGAATCGCTTTACCTCCTTTTAGTTGAAATTGCCCCGGCCCCGGTCCCTTGCCAAAATGTCTTGTTTAGCCCCGGGTGCTCCTGTCGGGTCTTGACTGATTCACACTTGATATTCTTGTCTTCTGGTTCTTGCTCTGATGAGCACACGTCTGCAACCA":
        1
    }
    for _, row in df.iterrows():
        assert expected[row["Sequence"]] == row["Count"]


@pytest.mark.usefixtures("new_pipegraph_no_qc")
def test_count_fastq_trimmed(tmpdir, raw_lane, scouter):
    scouter.write_fastq_count_trimmed(raw_lane)
    ppg.run_pipegraph()
    output_file = scouter.result_dir / f"{raw_lane.name}_{scouter.name}_all_reads_trimmed.tsv"
    df = pd.read_csv(output_file, sep="\t")
    expected = {
        "TTGCTTTACCTCCTTTTAGTTGGCCTTGCCCCGGCCCCGGTCCCTTGCCAAAATGTCTTGTTTAGCCCCGGGTGCTCCTGTCGGGTCTTGACTGATTCACACTTGATATTCTTGTCTTCTGGTTCTTGCTCTGATGAGCACACGTCTGCA":
        2,
        "TTGCTTTACCTCCTTTTAGTTGGCCTTGCCCCGGCCCCGGTCCCTTGCCAAAATGTCTTGTTTAGCCCCGGGGTGCTCCTGTCGGGTCTTGACTGATTCACACTTGATATTCTTGTCTTCTGGTTCTTGCTCTGATGAGCACACGTCTGC":
        2,
        "TTGCTTTACCTCCTTTTAGTTGGCCTTGCCCCGGCCCCGGTCCCTTGCCAAAATGTCTTGTTTAGCCCCGGGGTGCTCCTGTCGGGTCTTGACTGATTCACACTTGATATTCTTGTCTTCTGGTTCTTGCTCTGATGAGCACACGTCTG":
        1,
        "TTGCTTTACCTCCTTTTAGCCTCTTTTGCCCCGGCCCCGGTCCCTTGCCAAAATGTCTTGTTTAGCCCCGGGTGCTCCTGTCGGGTCTTGACTGATTCACACTTGATATTCTTGTCTTCTGGTTCTTGCTCTGATGAGCACACGTCTGCA":
        2,
        "AGGAATCGCTTTACCTCCTTTTAGTTGAAATTGCCCCGGCCCCGGTCCCTTGCCAAAATGTCTTGTTTAGCCCCGGGTGCTCCTGTCGGGTCTTGACTGATTCACACTTGATATTCTTGTCTTCTGGTTCTTGCTCTGATGAGCACACGT":
        1
    }
    for _, row in df.iterrows():
        assert expected[row["Sequence"]] == row["Count"]


@pytest.mark.usefixtures("new_pipegraph_no_qc")
def test_count_samples_fast(tmpdir, raw_lane, scouter):
    scouter.write_count_table(raw_lane, row_order=None)
    ppg.run_pipegraph()
    output_file = scouter.result_dir / f"{raw_lane.name}_{scouter.name}_sequence_count.tsv"
    df = pd.read_csv(output_file, sep="\t")
    df.index = df.Name
    df_predefined = pd.read_csv(scouter.sequence_file_path, sep="\t")
    df_predefined["Name"] = [
        f"{a}_{b}"
        for a, b in zip(
            df_predefined["Alteration"].astype(str),
            df_predefined["Effect"].astype(str),
        )
    ]
    expected = {
        "1>A_test1": 2,
        "1>A_test2": 2,
        "1>A_test3": 2,
        "1>A_test4": 2,
        "1>A_test5": 0,
    }
    for _, row in df_predefined.iterrows():
        assert df.loc[row["Name"]]["Read Count"] == expected[row["Name"]]
    