def read_fastq_iterator(file_object):
    """A very dump and simple fastq reader, mostly for testing the other more sophisticated variants

    Yield (seq, name, quality)
    """
    row1 = file_object.readline()
    row2 = file_object.readline()
    row3 = file_object.readline()
    row4 = file_object.readline()
    while row1:
        seq = row2[:-1]
        quality = row4[:-1]
        name = row1[1:-1]
        yield (seq, name, quality)
        row1 = file_object.readline()
        row2 = file_object.readline()
        _ = file_object.readline()
        row4 = file_object.readline()

