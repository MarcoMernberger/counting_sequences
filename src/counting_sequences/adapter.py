import pypipegraph as ppg
import mbf_genomes
import collections
import cutadapt
import cutadapt.align
import pandas as pd
import gzip
from pathlib import Path
from util import read_fastq_iterator


AdapterMatch = collections.namedtuple(
    "AdapterMatch", ["astart", "astop", "rstart", "rstop", "matches", "errors"]
)


class CutadaptMatch:
    def __init__(self, adapter_start: str, adapter_end: str, cutadapt_error_rate=0):
        self.adapter_sequence_begin = adapter_start
        self.adapter_sequence_end = adapter_end
        self.maximal_error_rate = cutadapt_error_rate
        self.adapter_sequence_begin_reverse = mbf_genomes.common.reverse_complement(
            self.adapter_sequence_begin
        )
        self.adapter_sequence_end_reverse = mbf_genomes.common.reverse_complement(
            self.adapter_sequence_end
        )
        self.where_fwd = (
            cutadapt.align.START_WITHIN_SEQ2 | cutadapt.align.STOP_WITHIN_SEQ2
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
            if isinstance(asequence, str):
                adapter = cutadapt.align.Aligner(
                    asequence if asequence else "",
                    self.maximal_error_rate / len(asequence),
                    where,
                    wildcard_ref=True,
                    wildcard_query=False,
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
                # forward adapter nowhere to be found, discard
                return None
            elif match_begin_fwd1 is None:
                # forward adapter is in read2
                seq1, qual1, seq2, qual2 = seq2, qual2, seq1, qual1
                match_begin_fwd = match_begin_fwd2
            elif match_begin_fwd2 is None:
                match_begin_fwd = match_begin_fwd1
            else:
                # both not none, take the best fit
                if match_begin_fwd1.errors < match_begin_fwd2.errors:
                    match_begin_fwd = match_begin_fwd1
                else:
                    match_begin_fwd = match_begin_fwd2
                    seq1, qual1, seq2, qual2 = seq2, qual2, seq1, qual1
            i1 = match_begin_fwd.rstop
            # adapter_begin forward found
            match_end_fwd = self.match(self.adapters["adapter_sequence_end"], seq1)
            if match_end_fwd is not None:
                # adapter_end forward found
                i2 = match_end_fwd.rstart
            else:
                i2 = len(seq1)
            # now the second read must have the reverse adapters
            match_end_rev = self.match(
                self.adapters["adapter_sequence_end_reverse"], seq2
            )
            if match_end_rev is not None:
                j1 = match_end_rev.rstop
                match_begin_rev = self.match(
                    self.adapters["adapter_sequence_begin_reverse"], seq2
                )
                if match_begin_rev is None:
                    j2 = len(seq2)
                else:
                    j2 = match_begin_rev.rstart
            else:
                # reverse read is not matching, discard
                return None
            s1 = seq1[i1:i2]
            q1 = qual1[i1:i2]
            s2 = seq2[j1:j2]
            q2 = qual2[j1:j2]
            if s1 == "" or s2 == "":
                return None
            return (s1.encode(), q1, name1, s2.encode(), q2, name2)

        return filter_func

    def filter_single(self):
        def filter_func(seq1, qual1, name1, seq2, qual2, name2):
            seq1 = seq1.decode()
            seq2 = seq2.decode()
            match_begin_fwd1 = self.match(self.adapters["adapter_sequence_begin"], seq1)
            match_begin_fwd2 = self.match(self.adapters["adapter_sequence_begin"], seq2)
            if match_begin_fwd1 is None and match_begin_fwd2 is None:
                # forward adapter nowhere to be found, discard
                return None
            elif match_begin_fwd1 is None:
                # forward adapter is in read2
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
                    # take the second
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
                    # take the first
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
        self, r1, r2, output_file, max=100000, dependencies=[],
    ):
        if isinstance(output_file, str):
            outfile = Path(output_file)
        outfile.parent.mkdir(parents=True, exist_ok=True)
        our_iter = read_fastq_iterator

        def __dump():
            with gzip.open(r1, "r") as op1:
                with gzip.open(r2, "r") as op2:
                    counter = collections.Counter()
                    examples = {}
                    count = 0
                    for tup in zip(our_iter(op1), our_iter(op2)):
                        print(tup)
                        seq1, name1, _ = tup[0]
                        seq2, name2, _ = tup[1]
                        seq1 = seq1.decode()
                        seq2 = seq2.decode()
                        match_begin_fwd1 = (
                            self.match(self.adapters["adapter_sequence_begin"], seq1)
                            is not None
                        )
                        match_begin_fwd2 = (
                            self.match(self.adapters["adapter_sequence_begin"], seq2)
                            is not None
                        )
                        match_end_fwd1 = (
                            self.match(self.adapters["adapter_sequence_end"], seq1)
                            is not None
                        )
                        match_end_fwd2 = (
                            self.match(self.adapters["adapter_sequence_end"], seq2)
                            is not None
                        )
                        match_begin_rev1 = (
                            self.match(
                                self.adapters["adapter_sequence_begin_reverse"], seq1
                            )
                            is not None
                        )
                        match_begin_rev2 = (
                            self.match(
                                self.adapters["adapter_sequence_begin_reverse"], seq2
                            )
                            is not None
                        )
                        match_end_rev1 = (
                            self.match(
                                self.adapters["adapter_sequence_end_reverse"], seq1
                            )
                            is not None
                        )
                        match_end_rev2 = (
                            self.match(
                                self.adapters["adapter_sequence_end_reverse"], seq2
                            )
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

    def count_most_common_sequences(
        self, r1, r2, output_file, max=100000, dependencies=[],
    ):
        if isinstance(output_file, str):
            outfile = Path(output_file)
        outfile.parent.mkdir(parents=True, exist_ok=True)
        our_iter = read_fastq_iterator

        def __dump():
            with gzip.open(r1, "r") as op1:
                with gzip.open(r2, "r") as op2:
                    counter = collections.Counter()
                    examples = {}
                    count = 0
                    for tup in zip(our_iter(op1), our_iter(op2)):
                        seq1, name1, _ = tup[0]
                        seq2, name2, _ = tup[1]
                        seq1 = seq1.decode()
                        seq2 = seq2.decode()
                        key = (seq1, seq2)
                        counter[key] += 1
                        if key not in examples:
                            examples[key] = (name1, name2)
                        count += 1
                        if count >= max:
                            break
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
            df = pd.DataFrame(to_df)
            df = df.sort_values("Count", ascending=False)
            df.to_csv(output_file, sep="\t", index=False)

        return ppg.FileGeneratingJob(outfile, __dump).depends_on(dependencies)
