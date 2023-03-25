import Bio.SeqIO
from BCBio import GFF
import pyrodigal
import argparse
import pandas as pd
import json


def get_options():
    description = 'Generates ORFs from a Bifrost graph.'
    parser = argparse.ArgumentParser(description=description,
                                     prog='ggcaller')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--infile',
                    default=None,
                    required=True,
                    help='List of input gff files, one per line. ')
    IO.add_argument('--outpref',
                    default="gene_pa",
                    help='Output prefix. Default = "gene_pa"')

    return parser.parse_args()

def parse_genes(infile):
    member_genome = None
    gene_dict = {}

    with open(infile, "r") as f:
        for rec in GFF.parse(f):
            for feat in rec.features:
                if member_genome == None:
                    member_genome = feat.id.split(".")[0]
                if feat.type == "gene":
                    gene_dict[feat.id] = (int(feat.location.start), int(feat.location.end))

    return member_genome, gene_dict


def main():
    options = get_options()
    infile = options.infile
    outpref = options.outpref

    infile_list = []

    with open(infile, "r") as f:
        for line in f:
            infile_list.append(line.rstrip())

    isolate_genes = {}
    for entry in infile_list:
        member_genome, gene_dict = parse_genes(entry)
        isolate_genes[member_genome] = gene_dict

    # save cluster assignments
    with open(outpref + '.json', 'w') as fp:
        json.dump(isolate_genes, fp)

if __name__ == "__main__":
    main()