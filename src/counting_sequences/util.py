import gzip
import subprocess
from gzip import GzipFile
from typing import BinaryIO, Union
from pathlib import Path
from typing import List, Dict
from pypipegraph import Job
import pypipegraph as ppg


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


def generate_stitched_fastq(output_file: Path, r1: Path, r2: Path, dependencies: List[Job] = [], options: Dict[str, str] = {}):
    """
    generate_stitched_fastq wrapper for ngmerge.

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

    def __dump():
        if not output_file.exists():
            cmd = [
                "/project/code/NGmerge/NGmerge",
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
            print(" ".join(cmd))
            subprocess.check_call(cmd)

    job = ppg.FileGeneratingJob(output_file, __dump).depends_on(deps)
    return job
