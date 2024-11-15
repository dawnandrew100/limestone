#built-in
import typing

def fasta_file_parser(fasta: str)->dict[str, str]:
    fasta_array = []
    with open(fasta, "r") as file:
        for line in file:
            line = line.strip(" \n")
            fasta_array.append(line)
    tag_line = fasta_array[0].strip(">").split()
    accession = tag_line[0]
    fasta_dict = {}
    fasta_dict[accession] = "".join(fasta_array[1:])
    return fasta_dict

def fasta_parser(fasta: typing.TextIO)->dict[str, str]:
    fasta_array = []
    for line in fasta:
      line = line.strip(" \n")
      fasta_array.append(line)
    tag_line = fasta_array[0].strip(">").split()
    accession = tag_line[0]
    fasta_dict = {}
    fasta_dict[accession] = "".join(fasta_array[1:])
    return fasta_dict
