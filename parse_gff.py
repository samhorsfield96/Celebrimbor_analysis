from BCBio import GFF
import argparse
import json


def get_options():
    description = 'Generates position json based on Prodigal gffs.'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python parse_gff.py')

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

    contig_pos = 0

    with open(infile, "r") as f:
        for rec in GFF.parse(f):
            record_id = rec.id
            if member_genome == None:
                member_genome = record_id.split(".")[0]
            for feat in rec.features:
                if feat.type == "CDS":
                    CDS_id = feat.id.split("_")[-1]
                    gene_dict[record_id + "_" + CDS_id] = (int(feat.location.start) + contig_pos, int(feat.location.end) + contig_pos)
            contig_pos += len(rec.seq)


    return member_genome, contig_pos, gene_dict


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
        member_genome, total_length, gene_dict = parse_genes(entry)
        isolate_genes[member_genome] = (total_length, gene_dict)

    # save cluster assignments
    with open(outpref + '.json', 'w') as fp:
        json.dump(isolate_genes, fp)

if __name__ == "__main__":
    main()