import gzip
import subprocess

# import pypipegraph as ppg
import pandas as pd
from gzip import GzipFile
from typing import BinaryIO, Union
from pathlib import Path
from typing import List, Callable
from pypipegraph import Job, FileGeneratingJob
from mbf.align import Sample
from pandas import DataFrame

try:
    import string

    maketrans = string.maketrans
except (ImportError, NameError, AttributeError):
    maketrans = bytes.maketrans


rev_comp_table = maketrans(
    b"ACBDGHKMNSRUTWVYacbdghkmnsrutwvy", b"TGVHCDMKNSYAAWBRTGVHCDMKNSYAAWBR"
)


def reverse_complement(seq):
    return seq[::-1].translate(rev_comp_table)


def get_fastq_iterator(filepath: Path):
    if filepath.suffix == ".gz":
        fileobj = gzip.open(filepath, "r")
    else:
        fileobj = filepath.open("rb")
    return read_fastq_iterator(fileobj)


def read_fastq_iterator(file_object: Union[BinaryIO, GzipFile]):
    """A very dump and simple fastq reader, mostly for testing the other more sophisticated variants

    Yield (seq, name, quality)
    """
    row1 = file_object.readline().decode()
    row2 = file_object.readline().decode()
    row3 = file_object.readline().decode()
    row4 = file_object.readline().decode()
    while row1:
        seq = row2[:-1]
        quality = row4[:-1]
        name = row1[1:-1]
        yield (seq, name, quality)
        row1 = file_object.readline().decode()
        row2 = file_object.readline().decode()
        _ = file_object.readline().decode()
        row4 = file_object.readline().decode()


def count_raw_input_reads(gz_filename1):
    if gz_filename1.endswith(".gz"):
        p1 = subprocess.Popen(
            ["gunzip", "-c", gz_filename1], stdout=subprocess.PIPE
        )
        p2 = subprocess.Popen(
            ["wc", "-l"], stdin=p1.stdout, stdout=subprocess.PIPE
        )
        x = int(p2.communicate()[0][:-1]) / 4
    else:
        p2 = subprocess.Popen(
            ["wc", "-l", gz_filename1], stdout=subprocess.PIPE
        )
        t = p2.communicate()[0]
        try:
            x = int(t.split()[0])
            x = x / 4
        except ValueError:
            print("Error ...")
            print(t)
            raise
    return x


def get_reads_for_lanes_df(
    outfile: Path, lanes: List[Sample], dependencies: List[Job] = []
) -> Job:
    #    print_read_loss checks for each lane the number of reads present and writes that to a file.
    outfile.parent.mkdir(exist_ok=True, parents=True)
    df_calc_func = get_reads_for_lanes_callable(lanes, dependencies)

    def __dump():
        df = df_calc_func()
        df.to_csv(outfile, index=False, sep="\t")

    return FileGeneratingJob(outfile, __dump).depends_on(dependencies)


def get_reads_for_lanes_callable(
    lanes: List[Sample], dependencies: List[Job] = []
) -> Callable:
    def calc():
        reads = {"Sample": [], "Reads": []}
        for sample_name in lanes:
            r1 = list(lanes[sample_name].get_aligner_input_filenames())[0]
            cmd = ["wc", "-l", str(r1)]
            try:
                outp = subprocess.check_output(cmd)
                count = int(outp.decode().split()[0]) / 4
            except subprocess.CalledProcessError:
                print(" ".join(cmd))
                raise
            reads["Sample"].append(sample_name)
            reads["Reads"].append(count)
        df = pd.DataFrame(reads)
        return df

    return calc


def get_sequence_df_from_excel_file(
    sequence_file_path: Path, sheet_name: str
) -> DataFrame:
    """Reads a table from excel file with a given sheet name"""
    return pd.read_excel(sequence_file_path, sheet_name=sheet_name)


def identity(something):
    return something
