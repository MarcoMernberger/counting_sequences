#!/usr/bin/env python
# -*- coding: utf-8 -*-
import collections
import pypipegraph as ppg
import pandas as pd
import numpy as np
from pathlib import Path
from mbf_align import Sample, fastq2
from typing import List, Optional, Dict, Union, Callable, Tuple
from pypipegraph import FileGeneratingJob, Job
from .util import count_raw_input_reads
from pypipegraph import Job
from mbf_align import Sample
from pandas import DataFrame

__author__ = "MarcoMernberger"
__copyright__ = "MarcoMernberger"
__license__ = "mit"


class SequenceCounter:
    """
    Trims reads from fastqs and compares them to a set of predefined sequences
    for identity.

    Reads in the fastqs are trimmed by searching for the first occurence of
    a given start (stop) sequence. In addition, resulting sequences can be
    trimmed to a predefined length.
    If no such sequence is given, all reads start at index 0.
    If no length is specified, take the rest of the sequence.
    The set of predefined sequences are trimmed is the same manner.
    Counting the sequences is then reduced to an identity comparison for speed.

    Parameters
    ----------
    sequence_file_path : str
        [description]
    name : str, optional
        [description], by default None
    seqs_to_trim_reads : Optional[Tuple[str]], optional
        [description], by default None
    seqs_to_trim_predefined : Optional[Tuple[str]], optional
        [description], by default None
    trimmed_length : int, optional
        [description], by default None
    result_folder : str, optional
        [description], by default "results/counts"
    dependencies : list, optional
        [description], by default []
    sequence_df_filter : Union[Callable, None], optional
        [description], by default None
    """

    def __init__(
        self,
        sequence_file_path: str,
        name: str = None,
        seqs_to_trim_reads: Optional[Tuple[str]] = None,
        seqs_to_trim_predefined: Optional[Tuple[str]] = None,
        trimmed_length: int = None,
        result_folder: str = "results/counts",
        dependencies=[],
        sequence_df_filter: Union[Callable, None] = None,
    ):
        """Contructor"""
        self.name = (
            name
            if name is not None
            else f"SC_{seqs_to_trim_reads}_{trimmed_length}"
        )
        self.seqs_to_trim_reads = seqs_to_trim_reads
        self.seqs_to_trim_predefined = seqs_to_trim_predefined
        self.trimmed_length = trimmed_length
        self.result_dir = Path(result_folder)
        self.result_dir.mkdir(parents=True, exist_ok=True)
        self.sequence_file_path = sequence_file_path
        self.dependencies = dependencies
        self.sequence_df_filter = sequence_df_filter

    @staticmethod
    def combine(inputfiles: Dict[str, Path], outfile: Path, dependencies: List[Job]) -> Job:
        """
        combine combines result files as generated by self.write_count_table().

        This collects a number of result count tables and concatenates them.

        Parameters
        ----------
        inputfiles : Dict[str, Path]
            dictionary of result files to combine.
        outfile : Path
            Path of the output table.
        dependencies : List[Job]
            Dependency jobs on which this depends.

        Returns
        -------
        Job
            FileGeneratingJob that creates the table.
        """
        outfile.parent.mkdir(parents=True, exist_ok=True)

        def dump():
            concats = []
            for sample_name in inputfiles:
                fileobj = inputfiles[sample_name]
                df = pd.read_csv(str(fileobj), sep="\t")
                df.insert(0, "Sample", [sample_name] * len(df))
                concats.append(df)
            df_ret = pd.concat(concats)
            df_ret.to_csv(outfile, sep="\t", index=False)

        return ppg.FileGeneratingJob(outfile, dump).depends_on(dependencies)

    def get_dependencies(self) -> List[Job]:
        """
        get_dependencies returns basic dependencies of the counter.

        Returns
        -------
        List[Job]
            List of dependencies.
        """
        return self.dependencies + [
            ppg.FileTimeInvariant(self.sequence_file_path),
            ppg.FunctionInvariant(
                f"{self.name}_write_predefined_sequences",
                self.write_predefined_sequences,
            ),
            ppg.ParameterInvariant(
                f"{self.name}_parameter",
                [self.seqs_to_trim_reads, self.seqs_to_trim_predefined, self.trimmed_length],
            ),
            ppg.FunctionInvariant(f"{self.name}_sequence_filter", self.sequence_df_filter),
            ppg.FunctionInvariant(f"{self.name}_get_trim_sequence_function", self.get_trim_sequence_function),
            ppg.FunctionInvariant(f"{self.name}_get_index_function", self._get_index_function),
            ppg.ParameterInvariant(f"{self.name}_params", [self.trimmed_length, self.seqs_to_trim_reads, self.seqs_to_trim_predefined])
        ]

    def get_trim_sequence_function(self, index_function: Callable) -> Callable:
        """
        get_trim_sequence_function returns a function that trims a given
        sequence according to the parameters.

        This assumes that the first index sequence must be present in a valid
        sequence. If not, an empty string is returned.

        Parameters
        ----------
        index_function : Callable
            Callable that returns a start/stop index for the sequence given.

        Returns
        -------
        Callable
            Callable that trims a given sequence and returns the trimmed string.
        """
        def _trim(fullseq):
            index1, index2 = index_function(fullseq)
            if index1 == -1:
                return ""
            return fullseq[index1:index2]

        def _trim_to_length(fullseq):
            trimmed_seq = _trim(fullseq)
            return trimmed_seq[:min(len(trimmed_seq), self.trimmed_length)]

        if self.trimmed_length is None:
            return _trim
        else:
            return _trim_to_length

    def write_predefined_sequences(self) -> Job:
        """
        write_predefined_sequences returns a job that creates a table
        of predefined sequences to search for.

        This job takes a table of input sequences that can be filtered
        via self.sequence_df_filter. Susequently, these sequences can be trimmed
        to a fixed length and / or by removing given start/end sequences.

        Returns
        -------
        Job
            Job that creates the table.
        """
        outputfile = self.result_dir / f"{self.name}_predefined_sequences.tsv"

        def __dump():
            df_sequence_df = pd.read_csv(self.sequence_file_path, sep="\t")
            if self.sequence_df_filter is not None:
                df_sequence_df = self.sequence_df_filter(df_sequence_df)
            if "Name" not in df_sequence_df.columns:
                try:
                    df_sequence_df["Name"] = [
                        f"{a}_{b}"
                        for a, b in zip(
                            df_sequence_df["Alteration"].astype(str),
                            df_sequence_df["Effect"].astype(str),
                        )
                    ]
                except KeyError:
                    print("You need a Name column to identify all the sequences.")
                    raise
            df_sequence_df = df_sequence_df.rename(
                columns={"Sequence": "Full Sequence"}
            )
            sequences = []
            deduplicated = []
            seen = collections.defaultdict(set)
            index_function = self._get_index_function(self.seqs_to_trim_predefined)
            trim_function = self.get_trim_sequence_function(index_function)
            df_sequence_df["Sequence"] = df_sequence_df["Full Sequence"].apply(
                trim_function
            )
            df_sequence_df["Duplicate"] = df_sequence_df.duplicated(
                subset="Sequence", keep="first"
            )
            df_sequence_df["Duplicate Full"] = df_sequence_df.duplicated(
                subset="Full Sequence", keep="first"
            )
            for _, row in df_sequence_df.iterrows():
                fullseq = row["Full Sequence"]
                deduplicated.append(fullseq not in seen)
                seen[fullseq].add(row["Name"])
                seq = trim_function(fullseq)
                sequences.append(seq)
            df_sequence_df["Duplicate Count"] = 1
            df_sequence_df["Duplicate Entries"] = ""
            for _, grouped_duplicates in df_sequence_df.groupby("Full Sequence"):
                if len(grouped_duplicates) > 1:
                    df_sequence_df.loc[
                        grouped_duplicates.index, "Duplicate Entries"
                    ] = ",".join(grouped_duplicates["Name"])
                    df_sequence_df.loc[
                        grouped_duplicates.index, "Duplicate Count"
                    ] = len(grouped_duplicates)
            df_sequence_df.to_csv(outputfile, sep="\t", index=False)
            try:
                pd.testing.assert_frame_equal(
                    df_sequence_df[["Duplicate"]], df_sequence_df[["Duplicate Full"]]
                )
            except AssertionError:
                duplicates = df_sequence_df[
                    df_sequence_df["Duplicate"] != df_sequence_df["Duplicate Full"]
                ]
                print(duplicates)
                print(duplicates[["Name", "Sequence", "Full Sequence"]])

        return ppg.FileGeneratingJob(outputfile, __dump).depends_on(self.get_dependencies())

    def write_fastq_count(
        self, raw_lane: Sample, dependencies: List[ppg.Job] = []
    ) -> FileGeneratingJob:
        """
        write_fastq_count counts all occuring sequences in a fastq file.

        Very basic counter to count all occuring sequences in a fastq file.
        In addition, offers the option to trim the reads from the first occurence
        of a start (adapter) kmer up to a certain length.
        trim the reads first to a certain length.
        Found sequences and counts are dumped as a csv table.

        Parameters
        ----------
        raw_lane : Sample
            Sample to count.
        dependencies : List[ppg.Job], optional
            List of dependency jobs, by default [].

        Returns
        -------
        FileGeneratingJob
            A job that generates the sequence count table.
        """
        output_file = self.result_dir / f"{raw_lane.name}_{self.name}_all_reads.tsv"

        def __write():
            df_counter = self.count_fastq(raw_lane)
            df_counter = df_counter.sort_values("Count", ascending=False)
            df_counter.to_csv(output_file, sep="\t", index=False)

        return (
            ppg.FileGeneratingJob(output_file, __write, empty_ok=True)
            .depends_on(dependencies)
            .depends_on(self.dependencies)
            .depends_on(raw_lane.prepare_input())
        )

    def count_fastq(self, raw_lane: Sample) -> DataFrame:
        """
        count_fastq counts all occuring sequences in a fastq file.

        Very basic counter to count all occuring sequences in a fastq file.
        In addition, the reads are trimmed in an equal manner to the predefined
        sequences to speed up counting by eliminating the need for substring
        matching.

        Parameters
        ----------
        raw_lane : Sample
            Sample to count.

        Returns
        -------
        DataFrame
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
    ) -> Job:
        """
        write_count_table writes a count table that contains the absolute and
        relative counts of predefined sequences contained in a sample.

        This is the main result we generate.

        Parameters
        ----------
        raw_lane : Sample
            Sample to count.
        row_order : Optional[List[str]], optional
            row order to sort by, by default None.
        dependencies : List[Job], optional
            List of dependency jobs, by default [].

        Returns
        -------
        Job
            Job that creates the output table.
        """
        output_file = (
            self.result_dir / f"{raw_lane.name}_{self.name}_sequence_count.tsv"
        )
        output_file2 = (
            self.result_dir
            / f"{raw_lane.name}_{self.name}_sequence_count_unmatched.tsv"
        )

        def __write():
            df_ret, df_unmatched = self.count_samples_fast(raw_lane, row_order)
            df_ret.to_csv(output_file, sep="\t", index=False)
            df_unmatched.to_csv(output_file2, sep="\t", index=False)

        return (
            ppg.MultiFileGeneratingJob([output_file, output_file2], __write)
            .depends_on(ppg.FunctionInvariant(f"{self.name}_count_fastq", self.count_fastq))
            .depends_on(self.write_fastq_count_trimmed(raw_lane))
            .depends_on(self.write_fastq_count(raw_lane))
            .depends_on(self.dependencies)
            .depends_on(self.write_predefined_sequences())
            .depends_on(raw_lane.prepare_input())
            .depends_on(self.get_dependencies())
            .depends_on(
                ppg.FunctionInvariant(
                    f"{self.name}_count_samples_fast", self.count_samples_fast
                )
            )
        )

    def write_fastq_count_trimmed(self, raw_lane: Sample) -> Job:
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
        output_file = (
            self.result_dir / f"{raw_lane.name}_{self.name}_all_reads_trimmed.tsv"
        )

        def __write():
            read_counter_file = (
                self.result_dir / f"{raw_lane.name}_{self.name}_all_reads.tsv"
            )
            df_read_counter = pd.read_csv(read_counter_file, sep="\t")
            counter: Dict[str, int] = collections.Counter()
            index_function = self._get_index_function(self.seqs_to_trim_reads)
            trim_function = self.get_trim_sequence_function(index_function)
            for _, row in df_read_counter.iterrows():
                seq1 = row["Sequence"]
                count = row["Count"]
                seq = trim_function(seq1)
                if len(seq) == 0:
                    continue
                counter[seq] += count
            df = pd.DataFrame.from_dict(
                {"Sequence": list(counter.keys()), "Count": list(counter.values())}
            )
            df.to_csv(output_file, sep="\t")

        return (
            FileGeneratingJob(output_file, __write)
            .depends_on(self.write_fastq_count(raw_lane))
            .depends_on(self.dependencies)
            .depends_on(raw_lane.prepare_input())
        )

    def _get_index_function(self, start_stop_seq_to_trim: Union[Tuple[str, str], None]) -> Callable:
        """
        _get_index_function returns an index function that returns an index
        to trim a sequnece by.

        This returns an index corrresponding to the first occurence of the start
        sequence supplied and the last occurence of the end sequence. Default
        ist the complete sequnece.

        Parameters
        ----------
        start_stop_seq_to_trim : Union[Tuple[str, str], None]
            Tuple of start/end sequences to search and trim by.

        Returns
        -------
        Callable
            Index function that returns the index to trim by.
        """
        start, stop = None, None
        if start_stop_seq_to_trim is not None:
            start, stop = start_stop_seq_to_trim

        def __find_index_from_start(seq):
            return seq.find(start) + len(start)

        def __find_index_from_end(seq):
            return seq.rfind(stop)

        def __return_length(seq):
            return len(seq)

        def __return_null(seq):
            return 0

        if start_stop_seq_to_trim is None:
            first_index_function = __return_null
            second_index_function = __return_length
        else:
            if start is None or start == "":
                first_index_function = __return_null
            else:
                first_index_function = __find_index_from_start
            if stop is None or stop == "":
                second_index_function = __return_length
            else:
                second_index_function = __find_index_from_end

        def _find_index(sequence):
            index1 = first_index_function(sequence)
            index2 = second_index_function(sequence)
            if index1 == -1:
                index1 = 0
            if index2 == -1:
                index2 = len(sequence)
            return index1, index2

        return _find_index

    def count_samples_fast(self, raw_lane: Sample, row_order=Union[None, List[str]]) -> Tuple[DataFrame, DataFrame]:
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
        read_counter_file = (
            self.result_dir / f"{raw_lane.name}_{self.name}_all_reads_trimmed.tsv"
        )
        df_read_counter = pd.read_csv(read_counter_file, sep="\t")
        counter = collections.Counter(
            pd.Series(
                df_read_counter.Count.values, index=df_read_counter.Sequence
            ).to_dict()
        )
        fn1 = raw_lane.get_aligner_input_filenames()[0]
        total_counts = count_raw_input_reads(str(fn1))
        df_sequence_df = pd.read_csv(
            self.result_dir / f"{self.name}_predefined_sequences.tsv", sep="\t"
        )
        counts = [counter[seq] for seq in df_sequence_df["Sequence"]]
        result = df_sequence_df.copy()
        result["Read Count"] = counts
        result["Frequency"] = result["Read Count"] / total_counts
        result["Name"] = pd.Categorical(result["Name"], categories=row_order)
        columns_in_order = [
            "Name",
            "Read Count",
            "Frequency",
        ]
        if "Position" in result.columns:
            result["Position"] = result["Position"].apply(get_pos)
            columns_in_order.append("Position")
        for col in df_sequence_df.columns:
            if col not in columns_in_order:
                columns_in_order.append(col)
        result = result[columns_in_order]
        result = result.sort_values("Name")
        result.index = result["Name"]
        # now the non_matched
        for seq in result["Sequence"].unique():
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
        return result, df_unmatched

    def generate_fastq_from_unmatched(self, raw_lane: Sample) -> FileGeneratingJob:
        """
        Generates a fastq for the unmatched files.

        The resulting fastq contains each unmatched read once. This can then be aligned and
        viewed in a browser to check what is exactly there.

        Parameters
        ----------
        raw_lane : Sample
            The sample for which we run the unmatched sequences.

        Returns
        -------
        FileGeneratingJob
            The job generating the fastq.
        """
        infile = (
            self.result_dir
            / f"{raw_lane.name}_{self.name}_sequence_count_unmatched.tsv"
        )
        outfile = (
            self.result_dir
            / f"{raw_lane.name}_{self.name}_sequence_count_unmatched.fastq"
        )

        def __dump_fastq():
            df = pd.read_csv(infile, sep="\t")
            df.index = df["Name"]
            with outfile.open("w") as outp:
                for name, row in df.iterrows():
                    outp.write(
                        f"@{name}_{row['Read Count']}\n{row['Sequence']}\n+\n{'F'*len(row['Sequence'])}\n"
                    )

        return (
            ppg.FileGeneratingJob(outfile, __dump_fastq)
            .depends_on(self.write_count_table(raw_lane))
            .depends_on(self.dependencies)
            .depends_on(self.get_dependencies())
        )


def grep_unmatched_reads(outfile: Path, raw_lane: Sample, df_file: DataFrame, dependencies=[]) -> Job:
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
