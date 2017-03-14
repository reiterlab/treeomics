#!/usr/bin/python
"""Create mutation matrix for benchmarking"""
import logging
import csv
import os

__author__ = 'Johannes REITER'
__date__ = 'January, 2016'


# get logger for application
logger = logging.getLogger('treeomics')


def write_mutation_matrix(phylogeny, filepath):
    """
    Generate ancestor/descendant mutation matrix of the Treeomics solution for efficient benchmarking
    across tools and approaches
    :param phylogeny: data structure around phylogeny
    :param filepath: path to output file
    """

    # for Python 2: newline ", newline=''" argument needs to be removed
    with open(filepath, 'w', newline='') as mm_file:

        logger.debug('Write mutation matrix file {}'.format(filepath))
        mm_writer = csv.writer(mm_file)

        # write header
        header = ['Ancestor/Descendant']
        for mut_idx, (chrom, start_pos, _) in enumerate(phylogeny.patient.mut_positions):
            header.append('{}_{}'.format(chrom, start_pos))
        mm_writer.writerow(header)

        founder_mp = set([sa_idx for sa_idx in range(
            len(phylogeny.patient.sample_names) if phylogeny.patient.sc_names is None
            else len(phylogeny.patient.sc_names))])

        opt_sol = phylogeny.solutions[0]    # optimal solution

        for mut_idx, (chrom, start_pos, _) in enumerate(phylogeny.patient.mut_positions):

            row = ['{}_{}'.format(chrom, start_pos)]
            descendants = set()
            same_subclone = set()

            # check if Treeomics generated a solution for this mutation
            if mut_idx in opt_sol.max_lh_mutations.keys():

                # get inferred mutation pattern of this variant
                mp = opt_sol.max_lh_mutations[mut_idx]

                # check for descending mutations
                for desc_mp in opt_sol.shared_mlh_mps.keys():
                    if desc_mp.issubset(mp) and len(desc_mp) < len(mp):
                        for desc_idx in opt_sol.shared_mlh_mps[desc_mp]:
                            descendants.add(desc_idx)

                    elif desc_mp == mp:     # mutations are acquired on the same edge
                        for desc_idx in opt_sol.shared_mlh_mps[desc_mp]:
                            same_subclone.add(desc_idx)

                # check for other founding mutations if this mutation is also a founder
                if mp == founder_mp:
                    for founder_idx in opt_sol.mlh_founders:
                        same_subclone.add(founder_idx)

                # unique mutation in samples of this pattern are always descendants
                if len(mp) > 1:
                    for sa_idx in mp:
                        for desc_idx in opt_sol.mlh_unique_mutations[sa_idx]:
                            descendants.add(desc_idx)

                # unique mutations are acquired in the same "subclone"
                else:
                    for sa_idx in mp:
                        for desc_idx in opt_sol.mlh_unique_mutations[sa_idx]:
                            same_subclone.add(desc_idx)

            elif opt_sol.conflicting_mutations is None or mut_idx not in opt_sol.conflicting_mutations:
                raise RuntimeError('Mutation {} {} was not assigned during inference procedure!'.format(
                    mut_idx, phylogeny.patient.mut_keys[mut_idx]))

            for desc_idx in range(len(phylogeny.patient.mut_positions)):
                if desc_idx in descendants and desc_idx != mut_idx:
                    row.append(str(1))
                elif desc_idx in same_subclone and desc_idx != mut_idx:
                    row.append(str(0.5))
                else:
                    row.append(str(0))

            mm_writer.writerow(row)

        logger.info('Wrote mutation matrix file {}'.format(os.path.abspath(filepath)))
