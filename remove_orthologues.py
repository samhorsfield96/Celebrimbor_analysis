import Bio.SeqIO
from BCBio import GFF
import pyrodigal
import argparse
import pandas as pd



def parse_genes(infile):
    gene_dict = {}

    with open(infile, "r") as f:
        for rec in GFF.parse(f):
            for feat in rec.features:
                if feat.type == "gene":
                    gene_dict[feat.id] = (int(feat.location.start), int(feat.location.end))
                    test = 1

    return gene_dict
    # orf_finder = pyrodigal.OrfFinder()
    # orf_finder.train(bytes(record.seq))
    # for pred in orf_finder.find_genes(bytes(record.seq)):


def main():

    #inmat = ""
    infile_list = ["/mnt/c/Users/sth19/PycharmProjects/PhD_project/Celebrimbor_analysis/data/SP_ATCC700669.gff3"]

    gene_pa = pd.read_csv()

    isolate_genes = {}
    for index, entry in infile_list:
        isolate_genes[index] = parse_genes(entry)

if __name__ == "__main__":
    main()