import pypipegraph as ppg
import mbf_genomes
import collections
import cutadapt
import cutadapt.align
import mbf_align
import pandas as pd
from pathlib import Path
from .util import read_fastq_iterator, get_fastq_iterator, reverse_complement
from typing import Callable, List, Union, Optional
from pypipegraph import Job
from mbf_align.fastq2 import iterate_fastq
from mbf_align import Sample


AdapterMatch = collections.namedtuple(
    "AdapterMatch", ["astart", "astop", "rstart", "rstop", "matches", "errors"]
)


class Paired_Filtered_Trimmed_From_Job(mbf_align.fastq2.Straight):
    """Filter reads with a callback func that takes seq1,qual1, name1,
    seq2, qual2, name2 and returns truncated reads/qualities
    """

    def __init__(
        self,
        sample: Sample,
        exon: str,
        wt: str,
        threshold: float = 0.75,
        sampling: int = 10000,
        dependencies: List = [],
        adapter_max_len: Optional[int] = None
    ):
        self.dependencies = dependencies
        mbf_align.fastq2.Straight.__init__(self)
        self.sample = sample
        self.exon = exon
        self.threshold = threshold
        self.sampling = sampling
        self.wt = wt
        self.name = f"Paired_Filtered_Trimmed_FromJob_{sample.name}"
        self.log_file = (
            self.sample.cache_dir / f"{self.sample.name}_adapter_sequences.log"
        )
        self.adapter_max_len = adapter_max_len

    def get_dependencies(self, output_filenames):
        return [
            self.generate_adapters_for_trimming(),
            ppg.ParameterInvariant(
                f"{self.name}_params",
                [self.exon, self.sampling, self.threshold, self.wt],
            ),
            self.set_adapter_matcher(),
        ]

    def generate_aligner_input_paired(
        self,
        output_filename1,
        output_filename2,
        list_of_fastq_pairs,
        reverse_reads,
        read_creator="fastq",
    ):
        if read_creator == "fastq":
            our_iter = iterate_fastq
        else:
            raise ValueError("Invalid read creator")  # pragma: no cover
        counter = 0
        seen = 0
        with open(output_filename1, "wb") as op1:
            with open(output_filename2, "wb") as op2:
                for fn1, fn2 in list_of_fastq_pairs:
                    for tup in zip(
                        our_iter(fn1, reverse_reads), our_iter(fn2, reverse_reads)
                    ):
                        seq1, qual1, name1 = tup[0]
                        seq2, qual2, name2 = tup[1]
                        seen += 1
                        filtered = self.filter_func(
                            seq1, qual1, name1, seq2, qual2, name2
                        )
                        if filtered is not None:
                            s1, q1, n1, s2, q2, n2 = filtered
                            if s1 is not None and s2 is not None:
                                op1.write(
                                    (b"@" + n1 + b"\n" + s1 + b"\n+\n" + q1 + b"\n")
                                )
                                op2.write(
                                    (b"@" + n2 + b"\n" + s2 + b"\n+\n" + q2 + b"\n")
                                )
                                counter += 1

    def generate_aligner_input(
        self, output_filename, list_of_fastqs, reverse_reads, read_creator="fastq"
    ):
        """This allows to see both mate pairs and select one of them"""

        if read_creator == "fastq":
            our_iter = iterate_fastq
        else:
            raise ValueError("Invalid read creator")  # pragma: no cover
        counter = 0
        seen = 0
        with open(output_filename, "wb") as op:
            for fn1, fn2 in list_of_fastqs:
                for tup in zip(
                    our_iter(fn1, reverse_reads), our_iter(fn2, reverse_reads)
                ):
                    seq1, qual1, name1 = tup[0]
                    seq2, qual2, name2 = tup[1]
                    seen += 1
                    filtered = self.filter_func(seq1, qual1, name1, seq2, qual2, name2)
                    if filtered is not None:
                        s1, q1, n1 = filtered
                        if s1 is not None:
                            op.write((b"@" + n1 + b"\n" + s1 + b"\n+\n" + q1 + b"\n"))
                            counter += 1

    def generate_adapters_for_trimming(
        self,
    ):
        # find start of wt in r1
        # find end.reverse of wt in r2
        # check lengths and raise if we do not encounter that in threshold*100 %
        r1, r2 = self.sample.get_aligner_input_filenames()

        def check():
            forwards = collections.Counter()
            reverses = collections.Counter()
            start = self.wt[:10].encode()
            end = self.wt[-10:].encode()
            start_reverse = reverse_complement(self.wt[:10]).encode()
            end_reverse = reverse_complement(self.wt[-10:]).encode()
            # check forward
            with self.log_file.open("w") as outp:
                total, found = 0, 0
                # for inp in [r1, r2]:
                for tup in zip(iterate_fastq(r1, False), iterate_fastq(r2, False)):
                    seq1, _, _ = tup[0]
                    seq2, _, _ = tup[1]

                    index1_end = seq1.find(start)
                    index2_end = seq2.find(end_reverse)
                    if index1_end != -1 and index2_end != -1:
                        index1_start = 0
                        index2_start = 0
                        if self.adapter_max_len is not None:
                            index1_start = max(0, index1_end - self.adapter_max_len)
                        if self.adapter_max_len is not None:
                            index2_start = max(0, index2_end - self.adapter_max_len)
                        forwards[seq1[index1_start:index1_end]] += 1  # this will take the whole beginning of the read. since we apparently cannot expect any kind of rational design, we have to assume that the relevant part will start at the very end of the read
                        reverses[seq2[index2_start:index2_end]] += 1
                        found += 1
                    else:
                        index1_end = seq2.find(start)
                        index2_end = seq1.find(end_reverse)
                        if index1_end != -1 and index2_end != -1:
                            index1_start = 0
                            index2_start = 0
                            if self.adapter_max_len is not None:
                                index1_start = max(0, index1_end - self.adapter_max_len)
                            if self.adapter_max_len is not None:
                                index2_start = max(0, index2_end - self.adapter_max_len)
                            forwards[seq2[index1_start:index1_end]] += 1
                            reverses[seq1[index2_start:index2_end]] += 1
                            found += 1
                        else:
                            # no flanking sequences found
                            pass
                    total += 1
                    if total >= self.sampling:
                        break
                freq = float(found) / total
                if freq < self.threshold:
                    raise ValueError(
                        f"The retrieved forward/reverse adapters occur in only {freq*100:.2f}% reads. You may discard to many reads. Check manually."
                    )
                forward_adapter = forwards.most_common(1)[0][0].decode()
                outp.write(
                    f"Most common start sequences:\n{forwards.most_common(10)}\n"
                )
                outp.write(f"Estimated forward adapter freqency = {freq}\n")
                reverse_adapter = reverse_complement(
                    reverses.most_common(1)[0][0].decode()
                )
                outp.write(f"Most common end sequences:\n{reverses.most_common(10)}\n")
                outp.write(f"Estimated reverse adapter freqency = {freq}\n")
                outp.write(f"{forward_adapter}\n{reverse_adapter}")

        return ppg.FileGeneratingJob(self.log_file, check).depends_on(self.dependencies)

    def set_adapter_matcher(self):
        def __load():
            with Path(self.log_file).open("r") as inp:
                readlines = inp.readlines()
            adapter_start = readlines[-2][:-1]
            adapter_end = readlines[-1]
            adapter_matcher = CutadaptMatch(
                adapter_start, adapter_end, cutadapt_error_rate=0
            )
            return adapter_matcher.filter()

        return ppg.AttributeLoadingJob(
            f"{self.name}_load", self, "filter_func", __load
        ).depends_on(self.generate_adapters_for_trimming())


class CutadaptMatch:
    def __init__(
        self,
        adapter_start: str,
        adapter_end: str,
        cutadapt_error_rate: int = 0,
        min_overlap: int = 3,
    ) -> None:
        self.adapter_sequence_begin = adapter_start
        self.adapter_sequence_end = adapter_end
        self.maximal_error_rate = cutadapt_error_rate
        self.min_overlap = min_overlap
        self.adapter_sequence_begin_reverse = mbf_genomes.common.reverse_complement(
            self.adapter_sequence_begin
        )
        self.adapter_sequence_end_reverse = mbf_genomes.common.reverse_complement(
            self.adapter_sequence_end
        )
        self.where_fwd = (
            cutadapt.align.START_WITHIN_SEQ2
            | cutadapt.align.START_WITHIN_SEQ1
            # | cutadapt.align.STOP_WITHIN_SEQ2
            | cutadapt.align.STOP_WITHIN_SEQ1
        )
        self.where_rev = (
            cutadapt.align.START_WITHIN_SEQ2
            | cutadapt.align.STOP_WITHIN_SEQ2
            | cutadapt.align.STOP_WITHIN_SEQ1
        )
        self.adapters = {}
        for aname, asequence, where in [
            ("adapter_sequence_begin", self.adapter_sequence_begin, self.where_fwd),
            (
                "adapter_sequence_begin_reverse",
                self.adapter_sequence_begin_reverse,
                self.where_fwd,
            ),
            ("adapter_sequence_end", self.adapter_sequence_end, self.where_rev),
            (
                "adapter_sequence_end_reverse",
                self.adapter_sequence_end_reverse,
                self.where_rev,
            ),
        ]:
            alen = len(asequence)
            if isinstance(asequence, str) and alen > 0:
                adapter = cutadapt.align.Aligner(
                    asequence if asequence else "",
                    self.maximal_error_rate / len(asequence),
                    # where,
                    wildcard_ref=True,
                    wildcard_query=False,
                    min_overlap=self.min_overlap,
                )
            else:
                adapter = None
            self.adapters[aname] = adapter

    def match(self, adapter, seq):
        alignment = adapter.locate(seq)
        if alignment is None:
            return None
        _match = AdapterMatch(*alignment)
        return _match

    def filter(self):
        def filter_func(seq1, qual1, name1, seq2, qual2, name2):
            seq1 = seq1.decode()
            seq2 = seq2.decode()
            match_begin_fwd1 = self.match(self.adapters["adapter_sequence_begin"], seq1)
            match_begin_fwd2 = self.match(self.adapters["adapter_sequence_begin"], seq2)
            if match_begin_fwd1 is None and match_begin_fwd2 is None:
                return None
            elif match_begin_fwd1 is None:
                # forward adapter is in read2
                seq1, qual1, seq2, qual2 = seq2, qual2, seq1, qual1
                match_begin_fwd = match_begin_fwd2
            elif match_begin_fwd2 is None:
                match_begin_fwd = match_begin_fwd1
            else:
                # both not none, take the best fit
                if match_begin_fwd1.errors <= match_begin_fwd2.errors:
                    match_begin_fwd = match_begin_fwd1
                else:
                    match_begin_fwd = match_begin_fwd2
                    seq1, qual1, seq2, qual2 = seq2, qual2, seq1, qual1
            i1 = match_begin_fwd.rstop
            # adapter_begin forward found
            match_end_fwd = self.match(
                self.adapters["adapter_sequence_end_reverse"],
                mbf_genomes.common.reverse_complement(seq1),
            )
            if match_end_fwd is not None:
                # adapter_end forward found
                i2 = len(seq1) - match_end_fwd.rstop
            else:
                i2 = len(seq1)
            # now the second read must have the reverse adapters
            match_end_rev = self.match(
                self.adapters["adapter_sequence_end_reverse"], seq2
            )
            if match_end_rev is not None:
                j1 = match_end_rev.rstop
                match_begin_rev = self.match(
                    self.adapters["adapter_sequence_begin"],
                    mbf_genomes.common.reverse_complement(seq2),
                )

                if match_begin_rev is None:
                    j2 = len(seq2)
                else:
                    j2 = len(seq2) - match_begin_rev.rstop  # match_begin_rev.rstart
            else:
                return None
            s1 = seq1[i1:i2]
            q1 = qual1[i1:i2]
            s2 = seq2[j1:j2]
            q2 = qual2[j1:j2]
            if s1 == "" or s2 == "":
                return None
            return (s1.encode(), q1, name1, s2.encode(), q2, name2)

        return filter_func

    def trim_and_sort(self):
        def trim_func(seq1, qual1, name1, seq2, qual2, name2):
            seq1 = seq1.decode()
            seq2 = seq2.decode()
            match_begin_fwd1 = None
            match_begin_fwd2 = None

            match_begin_fwd1
            match_begin_fwd1
            match_begin_fwd1 = self.match(self.adapters["adapter_sequence_begin"], seq1)
            match_begin_fwd2 = self.match(self.adapters["adapter_sequence_begin"], seq2)
            if match_begin_fwd1 is None and match_begin_fwd2 is None:
                # forward adapter nowhere to be found, discard
                print(self.adapters["adapter_sequence_begin"])
                print(seq1)
                print("<<match_begin_fwd1<<", match_begin_fwd1)
                print("here0")
                return None
            elif match_begin_fwd1 is None:
                # forward adapter is in read2
                seq1, qual1, seq2, qual2 = seq2, qual2, seq1, qual1
                match_begin_fwd = match_begin_fwd2
            elif match_begin_fwd2 is None:
                match_begin_fwd = match_begin_fwd1
            else:
                # both not none, take the best fit
                if match_begin_fwd1.errors <= match_begin_fwd2.errors:
                    match_begin_fwd = match_begin_fwd1
                else:
                    match_begin_fwd = match_begin_fwd2
                    seq1, qual1, seq2, qual2 = seq2, qual2, seq1, qual1
            i1 = match_begin_fwd.rstop
            # adapter_begin forward found
            match_end_fwd = self.match(
                self.adapters["adapter_sequence_end_reverse"],
                mbf_genomes.common.reverse_complement(seq1),
            )
            if match_end_fwd is not None:
                # adapter_end forward found
                print("<<match_end_fwd<<", match_end_fwd)
                i2 = len(seq1) - match_end_fwd.rstop
            else:
                i2 = len(seq1)
            # now the second read must have the reverse adapters
            match_end_rev = self.match(
                self.adapters["adapter_sequence_end_reverse"], seq2
            )
            # print("j1", self.adapter_sequence_end_reverse, match_end_rev)
            if match_end_rev is not None:
                j1 = match_end_rev.rstop
                match_begin_rev = self.match(
                    self.adapters["adapter_sequence_begin"],
                    mbf_genomes.common.reverse_complement(seq2),
                )

                if match_begin_rev is None:
                    j2 = len(seq2)
                else:
                    j2 = len(seq2) - match_begin_rev.rstop  # match_begin_rev.rstart
            else:
                # reverse read is not matching, discard
                print("here")
                return None
            s1 = seq1[i1:i2]
            q1 = qual1[i1:i2]
            s2 = seq2[j1:j2]
            q2 = qual2[j1:j2]
            print(seq1)
            print(seq2)
            print("----------")
            print(s1)
            print(s2)
            print(i1, i2, j1, j2)
            # raise ValueError()
            if s1 == "" or s2 == "":
                print("here2")
                return None
            print(i1, i2, j1, j2)
            print((s1.encode(), q1, name1, s2.encode(), q2, name2))
            raise ValueError()
            return (s1.encode(), q1, name1, s2.encode(), q2, name2)

        return filter_func

    def filter_single(self) -> Callable:
        """
        filter_single returns a filter function that filters paired end reads to a single read.

        filter_single returns a filter function that takes paired end reads,
        identifies the read that contains forward/reverse adapter sequences, trims thoser sequences
        and returns this single mate pair trimmed at the adapter occurences.

        Returns
        -------
        Callable
            Filter function for fastq processor.
        """

        def filter_func(seq1, qual1, name1, seq2, qual2, name2):
            seq1 = seq1.decode()
            seq2 = seq2.decode()
            match_begin_fwd1 = self.match(self.adapters["adapter_sequence_begin"], seq1)
            match_begin_fwd2 = self.match(self.adapters["adapter_sequence_begin"], seq2)
            if match_begin_fwd1 is None and match_begin_fwd2 is None:
                # forward adapter nowhere to be found, discard
                return None
            elif match_begin_fwd1 is None:
                # forward adapter is in read2, switch the reads
                seq1, qual1, name1, seq2, qual2, name2 = (
                    seq2,
                    qual2,
                    name2,
                    seq1,
                    qual1,
                    name1,
                )
                match_begin_fwd = match_begin_fwd2
                match_end_fwd = self.match(self.adapters["adapter_sequence_end"], seq2)
            elif match_begin_fwd2 is None:
                # forward adapter is in first read
                match_begin_fwd = match_begin_fwd1
                match_end_fwd = self.match(self.adapters["adapter_sequence_end"], seq1)
            else:
                # both not none, check the end
                match_end_fwd1 = self.match(self.adapters["adapter_sequence_end"], seq1)
                match_end_fwd2 = self.match(self.adapters["adapter_sequence_end"], seq2)
                if match_end_fwd1 is None and match_end_fwd2 is not None:
                    # take the second, switch the reads
                    seq1, qual1, name1, seq2, qual2, name2 = (
                        seq2,
                        qual2,
                        name2,
                        seq1,
                        qual1,
                        name1,
                    )
                    match_end_fwd = match_end_fwd2
                    match_begin_fwd = match_begin_fwd2
                else:
                    # take the first, as it contains the end adapter
                    match_end_fwd = match_end_fwd1
                    match_begin_fwd = match_begin_fwd1
            i1 = match_begin_fwd.rstop
            # adapter_begin forward found
            if match_end_fwd is not None:
                # adapter_end forward found
                i2 = match_end_fwd.rstart
            else:
                i2 = len(seq1)
            s1 = seq1[i1:i2]
            q1 = qual1[i1:i2]
            if s1 == "":
                return None
            return (s1.encode(), q1, name1)

        return filter_func

    def count_adapter_occurences(
        self,
        r1,
        r2,
        output_file,
        max=100000,
        dependencies=[],
    ):
        if isinstance(output_file, str):
            outfile = Path(output_file)
        outfile.parent.mkdir(parents=True, exist_ok=True)

        def __dump():
            iter1 = get_fastq_iterator(r1)
            iter2 = get_fastq_iterator(r2)
            counter = collections.Counter()
            examples = {}
            count = 0
            for tup in zip(iter1, iter2):
                seq1, name1, _ = tup[0]
                seq2, name2, _ = tup[1]
                # seq1 = seq1.decode()
                # seq2 = seq2.decode()
                match_begin_fwd1 = (
                    self.match(self.adapters["adapter_sequence_begin"], seq1)
                    is not None
                )
                match_begin_fwd2 = (
                    self.match(self.adapters["adapter_sequence_begin"], seq2)
                    is not None
                )
                match_end_fwd1 = (
                    self.match(self.adapters["adapter_sequence_end"], seq1) is not None
                )
                match_end_fwd2 = (
                    self.match(self.adapters["adapter_sequence_end"], seq2) is not None
                )
                match_begin_rev1 = (
                    self.match(self.adapters["adapter_sequence_begin_reverse"], seq1)
                    is not None
                )
                match_begin_rev2 = (
                    self.match(self.adapters["adapter_sequence_begin_reverse"], seq2)
                    is not None
                )
                match_end_rev1 = (
                    self.match(self.adapters["adapter_sequence_end_reverse"], seq1)
                    is not None
                )
                match_end_rev2 = (
                    self.match(self.adapters["adapter_sequence_end_reverse"], seq2)
                    is not None
                )
                key = (
                    match_begin_fwd1,
                    match_begin_fwd2,
                    match_end_fwd1,
                    match_end_fwd2,
                    match_begin_rev1,
                    match_begin_rev2,
                    match_end_rev1,
                    match_end_rev2,
                )
                counter[key] += 1
                if key not in examples:
                    examples[key] = (seq1, seq2, name1, name2)
                count += 1
                if count >= max:
                    break
            to_df = {
                "Begin_forward_1": [],
                "Begin_forward_2": [],
                "End_forward_1": [],
                "End_forward_2": [],
                "Begin_reverse_1": [],
                "Begin_reverse_2": [],
                "End_reverse_1": [],
                "End_reverse_2": [],
                "Count": [],
                "Example": [],
            }
            for key in counter:
                to_df["Begin_forward_1"].append(key[0])
                to_df["Begin_forward_2"].append(key[1])
                to_df["End_forward_1"].append(key[2])
                to_df["End_forward_2"].append(key[3])
                to_df["Begin_reverse_1"].append(key[4])
                to_df["Begin_reverse_2"].append(key[5])
                to_df["End_reverse_1"].append(key[6])
                to_df["End_reverse_2"].append(key[7])
                to_df["Count"].append(counter[key])
                to_df["Example"].append(examples[key])
            df = pd.DataFrame(to_df)
            df.to_csv(output_file, sep="\t", index=False)

        return ppg.FileGeneratingJob(outfile, __dump).depends_on(dependencies)

    @staticmethod
    def count_most_common_sequences(
        r1: Union[str, Path],
        r2: Optional[Union[str, Path]],
        output_file: Union[str, Path],
        max: int = 100000,
        dependencies: List[Job] = [],
        index: Optional[int] = None,
    ):
        if isinstance(output_file, str):
            outfile = Path(output_file)
        else:
            outfile = output_file
        outfile.parent.mkdir(parents=True, exist_ok=True)

        def __dump():
            iter1 = get_fastq_iterator(r1)
            iterlist = [iter1]
            if r2 is not None:
                iter2 = get_fastq_iterator(r2)
                iterlist.append(iter2)
            counter = collections.Counter()
            examples = {}
            count = 0
            for tup in zip(*iterlist):
                seqs = []
                names = []
                for read in tup:
                    s, n, _ = read
                    s = s[:index]
                    seqs.append(s)
                    names.append(n)
                key = tuple(seqs)
                counter[key] += 1
                if key not in examples:
                    examples[key] = tuple(names)
                count += 1
                if count >= max:
                    break
            if r2 is not None:
                to_df = {
                    "Seq1": [],
                    "Seq2": [],
                    "Count": [],
                    "Example": [],
                }
                for key in counter:
                    to_df["Seq1"].append(key[0])
                    to_df["Seq2"].append(key[1])
                    to_df["Count"].append(counter[key])
                    to_df["Example"].append(examples[key])
            else:
                to_df = {
                    "Seq": [],
                    "Count": [],
                    "Example": [],
                }
                for key in counter:
                    to_df["Seq"].append(key[0])
                    to_df["Count"].append(counter[key])
                    to_df["Example"].append(examples[key][0])

            df = pd.DataFrame(to_df)
            df = df.sort_values("Count", ascending=False)
            df.to_csv(output_file, sep="\t", index=False)

        pj = ppg.ParameterInvariant(
            f"most_common_{str(outfile)}",
            [index, max, str(output_file), str(r1), str(r2)],
        )
        return (
            ppg.FileGeneratingJob(outfile, __dump)
            .depends_on(dependencies)
            .depends_on(pj)
        )
