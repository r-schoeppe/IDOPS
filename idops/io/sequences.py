from Bio import SeqIO

_map_formats = {
    "faa": "fasta",
    "fasta": "fasta",
    "gb": "genbank",
    "gbk": "genbank",
}


def get_input_sequences(input_file):
    file_format = input_file.rsplit(".", 1)[-1]
    file_format = _map_formats[file_format]

    if file_format == "fasta":
        for record in SeqIO.parse(input_file, file_format):
            yield record.format("fasta")

    if file_format == "genbank":
        with open(input_file, "r+") as gb_file:
            gb_cds = SeqIO.InsdcIO.GenBankCdsFeatureIterator(gb_file)
            for cds in gb_cds:
                if cds.seq is not None:
                    cds.id = cds.name
                    cds.description = ''
                    yield cds.format("fasta")


def get_selected_sequences(dict_files_seqids):
    for input_file, seq_ids in dict_files_seqids.items():
        file_format = input_file.rsplit(".", 1)[-1]
        file_format = _map_formats[file_format]

        n = 0

        if file_format == "fasta":
            for record in SeqIO.parse(input_file, file_format):
                if record.id in seq_ids:
                    n += 1
                    yield record.format("fasta")

        if file_format == "genbank":
            with open(input_file, "r+") as gb_file:
                gb_cds = SeqIO.InsdcIO.GenBankCdsFeatureIterator(gb_file)
                for cds in gb_cds:
                    if cds.seq is not None and cds.name in seq_ids:
                        n += 1
                        cds.id = cds.name
                        cds.description = ''
                        yield cds.format("fasta")
        if n < len(seq_ids):
            raise RuntimeError("Missed a Sequence")
