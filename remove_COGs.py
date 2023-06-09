import argparse
import json
import pandas as pd
import numpy as np

def get_options():
    description = 'Removes COGs from presence/absence matrix based on epirical distribution of completeness.'
    parser = argparse.ArgumentParser(description=description,
                                     prog='python remove_COGs')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--pa-file',
                    default=None,
                    required=True,
                    help='Gene presence/absence .csv file generated by generate_pa_matrix.py. ')
    IO.add_argument('--clusters',
                    default=None,
                    required=True,
                    help='Gene cluster .json file generated by generate_pa_matrix.py. ')
    IO.add_argument('--loc',
                    default=None,
                    required=True,
                    help='Gene location .json file generated by parse_gff.py. ')
    IO.add_argument('--dist',
                    default=None,
                    required=True,
                    help='Genome completeness distribution for sampling. ')
    IO.add_argument('--prop-complete',
                    default=0,
                    type=float,
                    help='Proportion of complete genomes between 0-1. Default = 0. ')
    IO.add_argument('--prop-CDS',
                    default=0.8,
                    type=float,
                    help='Proportion of a gene to be present after slicing to be considered able to cluster. Default = 0.8')
    IO.add_argument('--avg-breaks',
                    type=int,
                    default=0,
                    help='Average number of breaks in section to be removed from genomes. Default = 0')
    IO.add_argument('--outpref',
                    default="gene_pa",
                    help='Output prefix. Default = "gene_pa"')

    return parser.parse_args()

def main():
    options = get_options()
    pa_file = options.pa_file
    clusters = options.clusters
    locations = options.loc
    prop_complete = options.prop_complete
    prop_CDS = options.prop_CDS
    dist_file = options.dist
    outpref = options.outpref
    avg_breaks = options.avg_breaks

    # read in prop_complete, generate distribution
    dist_list = []
    with open(dist_file, "r") as f:
        for line in f:
            dist_list.append(float(line.rstrip()) / 100)

    dist_list = np.array(dist_list)

    #read pa_file
    pa_mat = pd.read_csv(pa_file, index_col=0)
    genome_index = {v: k for v, k in enumerate(pa_mat.index)}

    # read in gene coords file
    with open(locations, "r") as f:
        gene_coords = json.load(f)

    num_genomes = len(gene_coords)

    # determine complete genomes
    number_complete = round(num_genomes * prop_complete)

    # sample from distribution genomes at random
    incomplete_sample = np.random.choice(dist_list, size=(num_genomes - number_complete))
    completeness_array = np.full(num_genomes, 1.0)

    # # set random isolates to non-complete
    incomplete_genomes = np.random.choice(len(gene_coords), size=incomplete_sample.size, replace=False)
    completeness_array[incomplete_genomes] = incomplete_sample

    # read cluster file
    with open(clusters, "r") as f:
        gene_clusters = json.load(f)

    to_remove = {}

    # iterate of genomes to augment
    for index, completeness in enumerate(completeness_array):
        # skip if complete
        if completeness == 1:
            to_remove[genome_index[index]] = (completeness, None, None)
            continue

        # pull out current genome data
        genome_size, gene_pos = gene_coords[genome_index[index]]
        to_remove_list = []

        # generate blocks to remove
        size_removal = round(genome_size * (1 - completeness))
        num_breaks = np.random.poisson(lam=avg_breaks, size=1)
        break_sites = np.random.choice(size_removal, size=num_breaks, replace=False).tolist()
        break_sites.extend([0, size_removal])
        break_sizes = np.diff(np.sort(break_sites))

        # shuffle array to ensure block sizes are randomised
        np.random.shuffle(break_sizes)
        sum_breaks = sum(break_sizes)

        # search between start and furthest possible position to remove all blocks
        start_pos = 1
        end_pos = (genome_size + 1) - sum_breaks
        for block in break_sizes:
            break_start = np.random.choice(range(start_pos, end_pos), size=1)[0]
            break_end = break_start + block

            #identify genes overlapping in range
            for CDS_id, CDS_coords in gene_pos.items():
                CDS_start, CDS_end = CDS_coords
                gene_length = CDS_end - CDS_start
                overlap = len(range(max(break_start, CDS_start), min(break_end, CDS_end) + 1))
                if overlap > (gene_length * (1 - prop_CDS)):
                    to_remove_list.append(gene_clusters[CDS_id])

            # update site values
            start_pos += break_end
            sum_breaks -= block
            end_pos = (genome_size + 1) - sum_breaks

        to_remove[genome_index[index]] = (completeness, to_remove_list, break_sizes.tolist())

    # update pa_mat
    for index, genome_id in genome_index.items():
        completeness, to_remove_list, break_sizes = to_remove[genome_id]
        # move on if complete
        if completeness == 1:
            continue

        row_slice = pa_mat.loc[genome_id, :]
        for element in to_remove_list:
            row_slice[element] = 0

    # write removed CDS json
    with open(outpref + '_removed.json', 'w') as fp:
        json.dump(to_remove, fp)

    # write updated array
    pa_mat.to_csv(outpref + "_removed_gene_pa.csv")

    #write completeness statistics
    completeness_df = pd.DataFrame({'Genome': pa_mat.index, 'Completeness' : completeness_array.tolist()})
    completeness_df.to_csv(outpref + "_completeness_sim.csv")

if __name__ == "__main__":
    main()