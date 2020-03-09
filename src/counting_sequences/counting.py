#!/usr/bin/env python
# -*- coding: utf-8 -*-
from pathlib import Path
from mbf_align import Sample, fastq2
from typing import List, Optional, Dict
from pypipegraph import FileGeneratingJob, Job
import collections
import pypipegraph as ppg
import pandas as pd
import numpy as np
import subprocess

__author__ = "MarcoMernberger"
__copyright__ = "MarcoMernberger"
__license__ = "mit"


def count_raw_input_reads(gz_filename1):
    if gz_filename1.endswith(".gz"):
        p1 = subprocess.Popen(["gunzip", "-c", gz_filename1], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["wc", "-l"], stdin=p1.stdout, stdout=subprocess.PIPE)
        x = int(p2.communicate()[0][:-1]) / 4
    else:
        p2 = subprocess.Popen(["wc", "-l", gz_filename1], stdout=subprocess.PIPE)
        t = p2.communicate()[0]
        try:
            x = int(t.split()[0])
            x = x / 4
        except ValueError:
            print("Error ...")
            print(t)
            raise
    return x


class SequenceCounter():
    """
    Trims reads from fastqs and compares them to a set of predefined sequences
    for identity.

    Reads in the fastqs are trimmed by searching for the first occurence of
    a given start sequence and taking the bases right of the first occurence
    up to specified length.
    If no such sequence is given, all reads start at index 0.
    If no length is specified, take the rest of the sequence.

    Returns
    -------
    [type]
        [description]

    Raises
    ------
    ValueError
        [description]
    ValueError
        [description]
    ValueError
        [description]
    ValueError
        [description]
    """
    def __init__(self, sequence_file_path: str, name: str = None, 
            start_seq_to_trim: Optional[str] = None, trimmed_length: int = None, 
            result_folder: str = "results/counts", dependencies = []
    ):
        """
        __init__ [summary]

        [extended_summary]

        Parameters
        ----------
        sequence_file_path : str
            Path to .tsv file containing the predefined sequences to count.
        name : str, optional
            A unique name, by default None.
        start_seq_to_trim : Optional[str], optional
            Start sequence for trimming, by default None.
        trimmed_length : int, optional
            Maximum length of trimmed reads, by default None.
        result_folder : str, optional
            Output folder to save the count results, by default "results/counts".
        """
        self.name = name if name is not None else f"SC_{start_seq_to_trim}_{trimmed_length}"
        self.start_seq_to_trim = start_seq_to_trim
        self.start_seq_len = 0 if start_seq_to_trim is None else len(start_seq_to_trim)
        self.trimmed_length = 5000 if trimmed_length is None else trimmed_length
        self.result_dir = Path(result_folder)
        self.result_dir.mkdir(parents=True, exist_ok=True)
        self.sequence_file_path = sequence_file_path
        self.dependencies = dependencies

    def get_dependencies(self):
        return self.dependencies + [ppg.FileTimeInvariant(self.sequence_file_path)]

    def write_predefined_sequences(self):
        outputfile = self.result_dir / "predefined_sequences.tsv"

        def __dump():
            df_sequence_df = pd.read_csv(self.sequence_file_path, sep="\t")
            df_sequence_df["Name"] = [
                f"{a}_{b}"
                for a, b in zip(
                    df_sequence_df["Alteration"].astype(str),
                    df_sequence_df["Effect"].astype(str),
                )
            ]
            index_function = self._get_index_function()
            df_sequence_df = df_sequence_df.rename(columns={"Sequence": "Full Sequence"})
            sequences = []
            seen = collections.defaultdict(set)
            for _, row in df_sequence_df.iterrows():
                fullseq = row["Full Sequence"]
                seen[fullseq].add(row["Name"])
                index1 = index_function(fullseq)
                index2 = index1 + min(len(fullseq), self.trimmed_length)
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
            df_sequence_df.to_csv(outputfile, sep="\t", index=False)
            self.assert_predefined(df_sequence_df["Full Sequence"].values, df_sequence_df["Sequence"].values)
        return ppg.FileGeneratingJob(outputfile, __dump).depends_on(self.dependencies)

    def write_fastq_count(self, raw_lane: Sample, dependencies: List[ppg.Job] = []) -> FileGeneratingJob:
        """
        Counts all occuring sequences in a fastq file.


        Very basic counter to count all occuring sequences in a fastq file.
        In addition, offers the option to trim the reads from the first occurence
        of a start (adapter) kmer up to a certain length.
        trim the reads first to a certain length.
        Found sequences and counts are dumped as a csv table.

        Returns
        -------
        FileGeneratingJob
            A job that generates the sequence count table.
        """
        output_file = self.result_dir / f"{raw_lane.name}_{self.name}_all_reads.tsv"

        def __write():
            df_counter = self.count_fastq(raw_lane)
            df_counter.to_csv(output_file, sep="\t", index=False)

        return ppg.FileGeneratingJob(output_file, __write).depends_on(dependencies).depends_on(self.dependencies)

    def count_fastq(self, raw_lane: Sample):
        """
        Counts all occuring sequences in a fastq file.


        Very basic counter to count all occuring sequences in a fastq file.
        In addition, the reads are trimmed in an equal manner to the predefined
        sequences to speed up counting by eliminating the need for substring
        matching.

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing the trimmed sequences and corresponding counts.
        """
        counter_to_df: Dict[str, List] = {"Sequence": [], "Count": []}
        fn1 = raw_lane.get_aligner_input_filenames()[0]
        counter: Dict[str, int] = collections.Counter()    
        for seq1, _, _ in fastq2.iterate_fastq(str(fn1), False):
            seq1 = seq1.decode()
            counter[seq1] += 1
        for seqs in counter:
            counter_to_df["Sequence"].append(seqs)
            counter_to_df["Count"].append(counter[seqs])
        df_counter = pd.DataFrame(counter_to_df)
        return df_counter

    def write_count_table(
        self,
        raw_lane: Sample,
        row_order: Optional[List[str]] = None,
        dependencies: List[Job] = [],
    ):
        output_file = self.result_dir / f"{raw_lane.name}_{self.name}_sequence_count.tsv"
        output_file2 = self.result_dir / f"{raw_lane.name}_{self.name}_sequence_count_unmatched.tsv"

        def __write():
            df_ret, df_unmatched = self.count_samples_fast(raw_lane, row_order)
            df_ret.to_csv(output_file, sep="\t", index=False)
            df_unmatched.to_csv(output_file2, sep="\t", index=False)

        return ppg.MultiFileGeneratingJob([output_file, output_file2], __write).depends_on(
            self.write_fastq_count_trimmed(raw_lane)).depends_on(
                self.write_fastq_count(raw_lane)).depends_on(self.dependencies).depends_on(
                    self.write_predefined_sequences())

    def write_fastq_count_trimmed(self, raw_lane: Sample):
        """
        Counts all occuring sequences in a fastq file.


        Very basic counter to count all occuring sequences in a fastq file.
        In addition, the reads are trimmed in an equal manner to the predefined
        sequences to speed up counting by eliminating the need for substring
        matching.

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing the trimmed sequences and corresponding counts.
        """
        output_file = self.result_dir / f"{raw_lane.name}_{self.name}_all_reads_trimmed.tsv"

        def __write():
            read_counter_file = self.result_dir / f"{raw_lane.name}_{self.name}_all_reads.tsv"
            df_read_counter = pd.read_csv(read_counter_file, sep="\t")
            counter: Dict[str, int] = collections.Counter()
            index_function = self._get_index_function()
            for _, row in df_read_counter.iterrows():
                seq1 = row["Sequence"]
                count = row["Count"]
                index1 = index_function(seq1)
                index2 = index1 + min(len(seq1), self.trimmed_length)
                seq = seq1[index1:index2]
                counter[seq] += count
            df = pd.DataFrame.from_dict(
                {"Sequence": list(counter.keys()), "Count": list(counter.values())})
            df.to_csv(output_file, sep="\t")

        return FileGeneratingJob(output_file, __write).depends_on(self.write_fastq_count(raw_lane)).depends_on(self.dependencies)

    def _find_index(self, sequence):
        index1 = sequence.find(self.start_seq_to_trim)
        if index1 == -1:
            return 0
        else:
            return index1 + len(self.start_seq_to_trim)

    def _get_index_function(self):
        if self.start_seq_to_trim is None:
            return lambda x: 0 
        else:
            return self._find_index

    def count_samples_fast(self, raw_lane, row_order=None):
        """
        Counts the occurence of all the predefined sequences given in
        self.sequence_file_path file and returns a DataFrame with additional
        information.

        The counting is done by simply reading a file with counted and trimmed
        reads from the input fastq file and trimming the given sequences
        in an equal manner to the predefined sequences.
        Then it just iterates over the predefined trimmed sequences and adds
        the already counted reads by dict lookup.

        Parameters
        ----------
        raw_lane : Sample
            The sample to count.
        row_order : Optional[List[str]], optional
            order to sort the rows by, by default None

        Returns
        -------
        pandas.DataFrame
            Dataframe containing the counted sequences and additional 
            information.
        """

        def get_pos(row_pos):
            poses = [int(alt) for alt in row_pos.split(",")]
            return np.max(poses)

        # get the counted reads from the fastq
        read_counter_file = self.result_dir / f"{raw_lane.name}_{self.name}_all_reads_trimmed.tsv"
        df_read_counter = pd.read_csv(read_counter_file, sep="\t")
        counter = collections.Counter(pd.Series(df_read_counter.Count.values, index=df_read_counter.Sequence).to_dict())
        to_df = {
            "Read Count": [],
            "Sequence": [],
            "Name": [],
            "Frequency": [],
            "Position": [],
            "Alteration": [],
            "Location": [],
            "Type": [],
            "AA Position": [],
            "Effect": [],
            "Mutant Codon": [],
            "WT Codon": [],
            "New AA": [],
            "Original AA": [],
        }
        fn1 = raw_lane.get_aligner_input_filenames()[0]
        total_counts = count_raw_input_reads(str(fn1))
        df_sequence_df = pd.read_csv(self.result_dir / "predefined_sequences.tsv", sep="\t")
        for _, row in df_sequence_df.iterrows():
            seq = row["Sequence"]
            pos = get_pos(row["Position"])
            count = counter[seq]
            seq_name = str(row["Name"])
            to_df["Name"].append(seq_name)
            to_df["Read Count"].append(count)
            to_df["Frequency"].append(count / total_counts)
            to_df["Sequence"].append(seq)
            to_df["Position"].append(pos)
            to_df["Alteration"].append(row["Alteration"]),
            to_df["Location"].append(row["Location"]),
            to_df["Type"].append(row["Type"]),
            to_df["AA Position"].append(row["AA Position"]),
            to_df["Effect"].append(row["Effect"]),
            to_df["Mutant Codon"].append(row["Mutant Codon"]),
            to_df["WT Codon"].append(row["WT Codon"]),
            to_df["New AA"].append(row["New AA"]),
            to_df["Original AA"].append(row["Original AA"])
        df_ret = pd.DataFrame(to_df)
        df_ret["Name"] = pd.Categorical(df_ret["Name"], categories=row_order)
        columns_in_order = [
            "Name",
            "Read Count",
            "Frequency",
            "Position",
            "Alteration",
            "Location",
            "Type",
            "Sequence",
            "AA Position",
            "Effect",
            "Mutant Codon",
            "WT Codon",
            "New AA",
            "Original AA",
        ]
        df_ret = df_ret[columns_in_order]
        df_ret = df_ret.sort_values("Name")
        df_ret.index = df_ret["Name"]
        # now the non_matched
        for seq in df_ret["Sequence"].unique():
            del counter[seq]
        to_df = {
            "Read Count": [],
            "Sequence": [],
            "Name": [],
            "Frequency": [],
        }
        for seq in counter:
            # pass over the remaining seqs in counter
            count = counter[seq]
            to_df["Name"].append("no_match")
            to_df["Read Count"].append(count)
            to_df["Frequency"].append(count / total_counts)
            to_df["Sequence"].append(seq)

        df_unmatched = pd.DataFrame(to_df)
        return df_ret, df_unmatched

    def assert_predefined(self, predefines, predefined_trimmed):
        """
        Checks wether the specified trimmed_length and start_sequence_to_trim
        does lead to ambigous trimmed reads.
        """
        predefined = set(predefines)
        trimmed = set(predefined_trimmed)
        assert len(predefined) == len(trimmed)


def grep_unmatched_reads(output_file, raw_lane, df_file, dependencies=[]):
    if isinstance(output_file, str):
        outfile = Path(output_file)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    def __dump():
        full_fasta_file = raw_lane.get_aligner_input_filenames()[0]
        df = pd.read_csv(df_file, sep="\t")
        check = {}
        for s in df["Sequence"]:
            check[s] = True
        with outfile.open("w") as output:
            for seq1, qual1, name1 in fastq2.iterate_fastq(str(full_fasta_file), False):
                seq1, qual1, name1 = seq1.decode(), qual1.decode(), name1.decode()
                if seq1 not in check:
                    output.write(f"@{name1}\n{seq1}\n+\n{qual1}\n")

    return ppg.FileGeneratingJob(outfile, __dump).depends_on(dependencies)
