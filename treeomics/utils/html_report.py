#!/usr/bin/python
"""Helper functions to generate report in HTML"""

import logging
import os
import numpy as np
import datetime
import math
from collections import Counter

from treeomics.version import __version__
import treeomics.utils.int_settings as def_sets
from treeomics.utils.driver import Driver
import treeomics.settings as settings


__author__ = 'Johannes REITER'
__date__ = 'March 6, 2015'


# get logger
logger = logging.getLogger(__name__)


class HTMLReport(object):
    """
    Generates an HTML formatted file providing an overview about the data and the results
    """

    def __init__(self, file_path, patient_name):
        """
        Opens the output file
        :param file_path: path to the output file of the created HTML report
        :param patient_name: name of the subject
        """

        self.filepath = file_path
        self.file = open(file_path, 'w')
        self.patient_name = patient_name

        self._ind = 0      # current indentation
        self._inds = ['\t'*i for i in range(10)]    # indentation

    def start_report(self):
        """
        Add HTML header and other basics to the report
        """

        # write HTML page header
        self.file.write('<!DOCTYPE html>\n')
        self.file.write('<html lang="en">\n')
        self._ind += 1      # indentation level increases by 1

        self.file.write(self._inds[self._ind]+'<head>\n')
        self._ind += 1      # indentation level increases by 1

        self.file.write(self._inds[self._ind]+'<meta charset="utf-8">\n')
        self.file.write(self._inds[self._ind]+'<meta http-equiv="X-UA-Compatible" content="IE=edge">\n')
        self.file.write(self._inds[self._ind]+'<meta name="viewport" content="width=device-width, initial-scale=1">\n')
        self.file.write(self._inds[self._ind]
                        + '<meta name="description" content="Treeomics analysis report of patient {}">\n'
                        .format(self.patient_name))
        self.file.write(self._inds[self._ind]+'<meta name="author" content="Johannes G. Reiter">\n\n')
        self.file.write(self._inds[self._ind]+'<title>Treeomics analysis report of patient {}</title>\n'
                        .format(self.patient_name))
        self.file.write(self._inds[self._ind]+'<link rel="stylesheet" '
                        'href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.2/css/bootstrap.min.css">\n')
        self.file.write(self._inds[self._ind]+'<style>body{ margin:0 100; background:white; }</style>\n')
        # self.file.write(self._inds[self._ind]+'<style>body{ margin:0 100; background:whitesmoke; }</style>\n')

        self._ind -= 1      # indentation level decreases by 1
        self.file.write(self._inds[self._ind]+'</head>\n\n')

        self.file.write(self._inds[self._ind]+'<body>\n')
        self._ind += 1      # indentation level increases by 1
        self.file.write(self._inds[self._ind]+'<div class="container" style="max-width:900px">\n')
        self._ind += 1      # indentation level increases by 1
        self.file.write(self._inds[self._ind]+'<div class="page-header">\n')
        self._ind += 1      # indentation level increases by 1
        self.file.write(self._inds[self._ind]+'<h2>Treeomics analysis report of patient {}</h2>\n'
                        .format(self.patient_name))
        self._ind -= 1      # indentation level decreases by 1
        self.file.write(self._inds[self._ind]+'</div>\n')

    def add_sequencing_information(self, patient, mut_table_path=None, put_driver_vars=None, put_driver_genes=None,
                                   unlikely_driver_mut_effects=None):
        """
        Adds basic information about the provided sequencing data and if provided add a mutation table plot
        :param patient: instance of class patient
        :param mut_table_path: path to the generated mutation table plot
        :param put_driver_vars: defaultdict with mutation IDs and and instance of driver class
        :param put_driver_genes: set of genes where putative driver gene mutations where identified
        :param unlikely_driver_mut_effects: mutation effects of variants in drivers, likely without effect
        """

        self.file.write(self._inds[self._ind]+'<h4>Input data</h4>\n')

        self.file.write(self._inds[self._ind]
                        + '<table class="table table-striped" '
                        + 'style="text-align: center;width:98%;max-width:800px;font-size:9pt">\n')
        self._ind += 1      # indentation level increases by 1
        col_names = ['Sample', 'Median coverage (mean)', 'Median MAF (mean)', 'Purity',
                     'No variants present', 'absent']
        header = ''.join('<th class="text-center">{}</th>'.format(col_name) for col_name in col_names)
        self.file.write(self._inds[self._ind]+header+'\n')

        # build up output data sequentially
        for sa_idx, sample_name in enumerate(patient.sample_names):

            row = list()
            row.append(sample_name.replace('_', ' '))

            # median coverage
            row.append('{:.1f} ({:.1f})'.format(np.nanmedian(patient.sample_coverages[sample_name]),
                                                np.nanmean(patient.sample_coverages[sample_name])))
            # if patient.sample_dis_phred_coverages is not None:
            #     row.append(np.median(patient.sample_dis_phred_coverages[sample_name]))
            # else:
            #     row.append('n/a')

            # Meadian and mean Variant Allele Frequency (VAF)
            row.append('{:.1%} ({:.1%})'.format(np.median(patient.sample_mafs[sample_name]),
                                                np.mean(patient.sample_mafs[sample_name])))

            # estimated purity
            if sample_name in patient.purities:
                row.append('{:.1%}'.format(patient.purities[sample_name]))
            else:
                row.append('-')

            # Bayesian inference model classification
            # present if probability to be present is greater than 50%
            # row.append(sum(1 for ps in patient.log_p01.values() if ps[sa_idx][1] > pres_lp))
            # row.append(sum(1 for ps in patient.log_p01.values() if ps[sa_idx][1] <= pres_lp))
            row.append(patient.no_present_vars[sa_idx])
            row.append(patient.no_absent_vars[sa_idx])

            # # Previous conventional classification
            # row.append(patient.positives[sample_name])
            # row.append(patient.negatives[sample_name])
            # row.append(patient.unknowns[0][sample_name]+patient.unknowns[1][sample_name])

            # write row to file
            self.file.write(self._inds[self._ind]+'<tr>'+''.join('<td>{}</td>'.format(c) for c in row)+'</tr>\n')

        self._ind -= 1      # indentation level decreases by 1
        self.file.write(self._inds[self._ind]+'</table>\n')

        # provide general summary about the samples
        self.file.write(self._inds[self._ind]+'<p>\n')
        self._ind += 1      # indentation level increases by 1
        self.file.write(self._inds[self._ind]+'Samples that passed the filtering: {}/{}</br>\n'.format(
            len(patient.sample_names), len(patient.sample_names) + len(patient.discarded_samples)))
        # median and mean coverage
        coverages = []
        for mut_key in patient.mut_reads.keys():
            # for sample_name in patient.mut_reads[mut_key].keys():     # all samples
            for sample_name in patient.sample_names:
                if patient.coverage[mut_key][sample_name] >= 0:
                    coverages.append(patient.coverage[mut_key][sample_name])
        self.file.write(self._inds[self._ind]+'Median coverage in the passed samples: {} (mean: {:.2f})'
                        .format(np.median(coverages), np.mean(coverages))+'\n')
        self._ind -= 1      # indentation level decreases by 1
        self.file.write(self._inds[self._ind]+'</p>\n')

        # provide general summary about the variants
        self.file.write(self._inds[self._ind]+'<p>\n')
        self._ind += 1      # indentation level increases by 1
        self.file.write(self._inds[self._ind]+'Total number of passed somatic variants: {} </br>\n'.format(
            len(patient.mutations)))

        if patient.variant_stats is not None:
            if patient.variant_stats[-5] > 0:
                self.file.write(self._inds[self._ind]+'Variants removed due to common variants filter: {} </br>\n'
                                .format(patient.variant_stats[-5]))
            if patient.variant_stats[-2] + patient.variant_stats[-3] + patient.variant_stats[-4] > 0:
                self.file.write(
                    self._inds[self._ind]+'Removed intronic/intergenic variants due to WES filter: {} </br>\n'.format(
                        patient.variant_stats[-2] + patient.variant_stats[-3] + patient.variant_stats[-4]))
            if patient.variant_stats[-1] > 0:
                self.file.write(
                    self._inds[self._ind]+'Variants removed due to never reaching significant level: {} </br>\n'.format(
                        patient.variant_stats[-1]))
            if patient.variant_stats[-6] > 0:
                self.file.write(self._inds[self._ind]+'Variants removed due to detection in normal sample: {} </br>\n'
                                .format(patient.variant_stats[-6]))

        self.file.write(self._inds[self._ind]+'Variants classified as present in at least one of the samples '
                                              'that passed the filtering: {} </br>\n'.format(
                                              len(patient.present_mutations)))

        # self.file.write(self._inds[self._ind]+'Mean number of variants per sample: {:.1f} </br>\n'.format(
        #                 (float(sum(len(muts) for sa_idx, muts in patient.samples.items()))
        #                  / len(patient.sample_names))))

        self.file.write(self._inds[self._ind]+"Founders (variants present in all samples): {} ({:.1%}) </br>\n".format(
            len(patient.founders), float(len(patient.founders))/len(patient.present_mutations) if
            len(patient.present_mutations) > 0 else float('nan')))
        self.file.write(
            self._inds[self._ind] + 'Mean number of unique (private) variants per sample: {:.1f} ({:.1%}) </br>\n'
            .format(float(len(patient.shared_muts[1])) / len(patient.sample_names),
                    (float(len(patient.shared_muts[1])) / len(patient.sample_names)) /
                    (sum(len(muts) for sa_idx, muts in patient.samples.items()) / len(patient.sample_names))))

        def _clamp(x):
            return int(max(0, min(x, 255)))

        def _driver_tag(gene_name, driver, suffix=None):
            if driver is None or driver.mutation_effect == 'unknown':
                return gene_name

            # number of supporting sources
            if Driver.MaxSourceSupport > 0:
                r, g, b, _ = Driver.colors()[len(driver.sources)]
            else:
                r, g, b = (0.8, 0.2, 0.2)
            c = 'color:#{0:02x}{1:02x}{2:02x};'.format(_clamp(r*255), _clamp(g*255), _clamp(b*255))
            tag = '<font style="{}">{}{}</font>'.format(c, gene_name, '' if suffix is None else suffix)
            if driver.cgc_driver:   # is driver in CGC list
                tag = '<strong>'+tag+'</strong>'

            return tag

        # add basic information about driver gene mutations
        if put_driver_vars is not None:
            dr_gene_cnts = Counter([d.gene_name + ('_CGC' if d.cgc_driver else '') for d in put_driver_vars.values()
                                    if d.gene_name is not None])
            driver_strs = []
            added_drivers = set()
            for mut_idx, dr in put_driver_vars.items():
                if dr.cgc_driver:
                    if (dr.gene_name + '_CGC') in added_drivers:
                        continue
                elif dr.gene_name in added_drivers:
                    continue

                driver_strs.append(_driver_tag(
                    patient.gene_names[mut_idx], dr,
                    suffix=None if dr_gene_cnts[patient.gene_names[mut_idx] + ('_CGC' if dr.cgc_driver else '')] == 1
                    else '({})'.format(dr_gene_cnts[patient.gene_names[mut_idx] + ('_CGC' if dr.cgc_driver else '')])))

                added_drivers.add(dr.gene_name + ('_CGC' if dr.cgc_driver else ''))

            self.file.write(
                self._inds[self._ind] + 'Total number of likely driver gene mutations: {} in {} genes ({}; '.format(
                    len(put_driver_vars), len(put_driver_genes),
                    ', '.join(d for d in sorted(driver_strs))) +
                ('more red colored gene names correspond to better supported driver genes. '
                 if Driver.MaxSourceSupport > 0 else '') +
                'Bold names are also present in the Cancer Gene Census list)</br>\n')

        if unlikely_driver_mut_effects is not None and sum(unlikely_driver_mut_effects.values()) > 0:
            self.file.write(self._inds[self._ind] + 'Other mutations in likely driver genes: {}</br>\n'.format(
                ', '.join('{} ({})'.format(d, e) for d, e in unlikely_driver_mut_effects.items())))

        self._ind -= 1      # indentation level decreases by 1
        self.file.write(self._inds[self._ind]+'</p>\n')

        if mut_table_path is not None:

            self.file.write(self._inds[self._ind]+'<div align="center">\n')
            self.file.write(self._inds[self._ind]+'<div style="width:98%;max-width:800px">\n')
            self._ind += 1      # indentation level increases by 1

            self.file.write(self._inds[self._ind]+'<figure>\n')
            self._ind += 1      # indentation level increases by 1
            self.file.write(self._inds[self._ind]+'<img class="img-responsive" src="'+mut_table_path +
                                                  '" alt="Variant classification table" width="800"/>'+'\n')
            self.file.write(self._inds[self._ind]+'<div align="left">\n')
            self.file.write(self._inds[self._ind]+'<figcaption>\n')
            self.file.write(
                self._inds[self._ind]+'<b>Probabilistic variant classification across '
                + '{} samples of patient {}.</b>\n'.format(len(patient.sample_names), patient.name))
            self.file.write(self._inds[self._ind]+'Blue rectangles correspond to present variants, '
                            + 'red to absent variants, and white to unknown mutation status. '
                              'Brighter colors denote higher probability. Reddish colored gene names correspond to '
                              'likely cancer drivers with increasing support. '
                              'Bold names are also present in the Cancer Gene Census list. \n')
            # self.file.write(
            #     self._inds[self._ind] +
            #     'In total {} distinct variants were classified as present in at least one sample'
            #         .format(len(patient.present_mutations)) + ', {} ({:.1%}) of those were founders, \n'.format(
            #         len(patient.founders), float(len(patient.founders))/len(patient.present_mutations))
            #     + 'and {} mutations were unique to single samples.'.format(len(patient.shared_muts[1])))

            self.file.write(self._inds[self._ind]+'</figcaption>\n')
            self.file.write(self._inds[self._ind]+'</div>\n')
            self._ind -= 1      # indentation level decreases by 1
            self.file.write(self._inds[self._ind]+'</figure>\n')

            self._ind -= 1      # indentation level decreases by 1
            self.file.write(self._inds[self._ind]+'</div>\n')
            self.file.write(self._inds[self._ind]+'</div>\n')

        self.file.write(self._inds[self._ind]+'</br>\n\n')

    def add_similarity_information(self, patient):
        """
        Add tables providing the Jaccard similarity coefficient and the genetic distance between all pairs of
        samples to the HTML report.
        :param patient: instance of class patient
        """

        self.file.write(self._inds[self._ind]+'<h4>Genetic similarity</h4>\n')

        # add table with the Jaccard similarity coefficient between all pairs of samples
        self.file.write(self._inds[self._ind] + '<table class="table table-striped" '
                        + 'style="text-align: center;width:98%;max-width:800px;font-size:9pt">\n')
        self._ind += 1      # indentation level increases by 1

        jscs = [patient.bi_sim_coeff[sa1_idx][sa2_idx] for sa1_idx in range(patient.n) for sa2_idx in range(sa1_idx)]

        # add table caption
        self.file.write(
            self._inds[self._ind] + '<caption>Probabilistic Jaccard similarity coefficient between all pairs of ' +
            'samples (median: {:.2f}; mean: {:.2f}; classification threshold: {:.0%}).</caption>\n'.format(
                np.median(jscs), np.mean(jscs), def_sets.CLA_CONFID_TH))

        col_names = ['Sample'] + patient.sample_names
        header = ''.join('<th class="text-center">{}</th>'.format(col_name.replace('_', ' ')) for col_name in col_names)
        self.file.write(self._inds[self._ind]+header+'\n')

        # build up output data sequentially
        for sa1_idx, sample_name in enumerate(patient.sample_names):

            row = list()
            row.append(sample_name.replace('_', ' '))

            for sa2_idx in range(sa1_idx+1):
                row.append('{:.2f}'.format(patient.bi_sim_coeff[sa1_idx][sa2_idx]))

            # write row to file
            self.file.write(self._inds[self._ind]+'<tr>'+''.join('<td>{}</td>'.format(c) for c in row)+'</tr>\n')

        self._ind -= 1      # indentation level decreases by 1
        self.file.write(self._inds[self._ind]+'</table>\n\n')

        # self.file.write(self._inds[self._ind]+'</br>\n\n')

        # add table with the genetic distance between all pairs of samples
        self.file.write(self._inds[self._ind] + '<table class="table table-striped" '
                        + 'style="text-align: center;width:98%;max-width:800px;font-size:9pt">\n')
        self._ind += 1      # indentation level increases by 1

        gds = [patient.bi_gen_dis[sa1_idx][sa2_idx] for sa1_idx in range(patient.n) for sa2_idx in range(sa1_idx)]

        # add table caption
        self.file.write(
            self._inds[self._ind] + '<caption>Genetic distance between all pairs of samples ' +
            '(median: {:.0f}; mean: {:.1f}; classification threshold: {:.0%}).</caption>\n'.format(
                np.median(gds), np.mean(gds), def_sets.CLA_CONFID_TH))

        col_names = ['Sample'] + patient.sample_names
        header = ''.join('<th class="text-center">{}</th>'.format(col_name.replace('_', ' ')) for col_name in col_names)
        self.file.write(self._inds[self._ind]+header+'\n')

        # build up output data sequentially
        for sa1_idx, sample_name in enumerate(patient.sample_names):

            row = list()
            row.append(sample_name.replace('_', ' '))

            for sa2_idx in range(sa1_idx+1):
                row.append('{}'.format(patient.bi_gen_dis[sa1_idx][sa2_idx]))

            # write row to file
            self.file.write(self._inds[self._ind]+'<tr>'+''.join('<td>{}</td>'.format(c) for c in row)+'</tr>\n')

        self._ind -= 1      # indentation level decreases by 1
        self.file.write(self._inds[self._ind]+'</table>\n')

        self.file.write(self._inds[self._ind]+'</br>\n\n')

        # write temporary results to file
        self.file.flush()
        os.fsync(self.file.fileno())

    def add_tree_plot(self, patient, phylogeny):
        """
        Add create ete3 tree of the inferred phylogeny to the HTML report
        :param patient: instance of class patient
        :param phylogeny: data structure around the inferred phylogenetic tree
        """

        self.file.write(self._inds[self._ind]+'<h4>Inferred cancer phylogeny</h4>\n')

        self.file.write(self._inds[self._ind]+'<div align="center">\n')
        self.file.write(self._inds[self._ind]+'<div style="width:98%;max-width:700px">\n')
        self._ind += 1      # indentation level increases by 1

        self.file.write(self._inds[self._ind]+'<figure>\n')
        self._ind += 1      # indentation level increases by 1
        self.file.write(self._inds[self._ind]+'<img class="img-responsive" src="'
                        + os.path.basename(phylogeny.tree_plot) + '" alt="Evolutionary tree" width="750"/>'+'\n')
        self.file.write(self._inds[self._ind]+'<div align="left">\n')
        self.file.write(self._inds[self._ind]+'<figcaption>\n')

        self.file.write(
            self._inds[self._ind]+'<b>Phylogenetic tree illustrating the evolutionary history of '
                                  'the cancer in {}.</b>\n'.format(patient.name))
        self.file.write(
            self._inds[self._ind] + 'Numbers on top of each branch indicate the number of acquired variants '
                                    '(including likely driver gene mutations reported separately in orange).\n')
        self.file.write(
            self._inds[self._ind] + 'Numbers on bottom of each branch indicate the estimated support values.\n')

        if logger.isEnabledFor(logging.DEBUG):
            self.file.write(
                self._inds[self._ind] + 'Frequencies in brackets denote the median VAF of the mutations on a branch.\n')

        if phylogeny.max_no_mps is not None:
            self.file.write(self._inds[self._ind] +
                            'Only the {} most likely mutation patterns were considered for each variant.'.format(
                            phylogeny.max_no_mps))

        self.file.write(self._inds[self._ind]+'</figcaption>\n')
        self.file.write(self._inds[self._ind]+'</div>\n')
        self._ind -= 1      # indentation level decreases by 1
        self.file.write(self._inds[self._ind]+'</figure>\n')

        self._ind -= 1      # indentation level decreases by 1
        self.file.write(self._inds[self._ind]+'</div>\n')
        self.file.write(self._inds[self._ind]+'</div>\n')

        self.file.write(self._inds[self._ind]+'</br>\n\n')

    def add_conflict_graph(self, patient, mp_graph_name, phylogeny=None):
        """
        Add evolutionary conflict plot graph to the HTML report
        :param patient: instance of class patient
        :param mp_graph_name: path to the evolutionary conflict graph plot file
        :param phylogeny: data structure around the inferred phylogenetic tree
        """

        self.file.write(self._inds[self._ind]+'<h4>Evolutionary conflict graph</h4>\n')

        self.file.write(self._inds[self._ind]+'<div align="center">\n')
        self.file.write(self._inds[self._ind]+'<div style="width:98%;max-width:700px">\n')
        self._ind += 1      # indentation level increases by 1

        self.file.write(self._inds[self._ind]+'<figure>\n')
        self._ind += 1      # indentation level increases by 1
        self.file.write(self._inds[self._ind]+'<img class="img-responsive" src="'+mp_graph_name +
                                              '" alt="Evolutionary conflict graph" width="500"/>'+'\n')
        self.file.write(self._inds[self._ind]+'<div align="left">\n')
        self.file.write(self._inds[self._ind]+'<figcaption>\n')
        self.file.write(
            self._inds[self._ind]+'<b>Evolutionary conflict graph of {} samples in patient {}.</b>\n'.format(
                len(patient.sample_names), patient.name))
        self.file.write(
            self._inds[self._ind] + 'Treeomics considered {:.0f} distinct mutation patterns (MPs).\n'.format(
                math.pow(2, len(patient.sample_names))+1))
        self.file.write(self._inds[self._ind]+'Each circular line represents a distinct sample. '
                        + 'Inner to outer lines denote: '
                        + ', '.join(sa_name.replace('_', ' ') for sa_name in patient.sample_names)
                        + '. Marks on these lines denote present variants. \n')
        self.file.write(
            self._inds[self._ind]+'Labels denote the MP reliability scores. ' +
            'Only nodes with the highest reliability score are depicted. ' +
            'Blue colored nodes (MPs) are evolutionarily compatible ' +
            'and red colored nodes are evolutionarily incompatible indicated by edges among the nodes.' +
            (' Minimum reliability score value to be considered as a potential subclone: {:.3f}.'.format(
             phylogeny.min_score) if phylogeny is not None else '') +
            (' Only the {} most likely mutation patterns were considered for each variant.'.format(
                phylogeny.max_no_mps) if phylogeny is not None and phylogeny.max_no_mps is not None else '') + '\n')

        self.file.write(self._inds[self._ind]+'</figcaption>\n')
        self.file.write(self._inds[self._ind]+'</div>\n')
        self._ind -= 1      # indentation level decreases by 1
        self.file.write(self._inds[self._ind]+'</figure>\n')

        self._ind -= 1      # indentation level decreases by 1
        self.file.write(self._inds[self._ind]+'</div>\n')
        self.file.write(self._inds[self._ind]+'</div>\n')

        self.file.write(self._inds[self._ind]+'</br>\n\n')

    def add_inc_mp_information(self, phylogeny, incomp_mps_plot_filepath=None, plot_width=None):
        """
        Add table with sequencing data information about the incompatible mutation patterns and add illustrative plot
        about the VAF and coverage in these incompatible patterns if a path is given.
        :param phylogeny: data structure around the inferred phylogenetic tree
        :param incomp_mps_plot_filepath: path to illustrative sequencing data plot of incompatible mutation patterns
        :param plot_width: width of the plot
        """

        self.file.write(self._inds[self._ind]+'<h4>Evolutionarily incompatible mutation patterns</h4>\n')

        pat = phylogeny.patient

        # add table with the genetic distance between all pairs of samples
        self.file.write(self._inds[self._ind] + '<table class="table table-striped" '
                        + 'style="text-align: center;width:98%;max-width:800px;font-size:9pt">\n')
        self._ind += 1      # indentation level increases by 1

        col_names = ['Variant'] + pat.sample_names
        header = ''.join('<th class="text-center">{}</th>'.format(col_name.replace('_', ' ')) for col_name in col_names)
        self.file.write(self._inds[self._ind]+header+'\n')

        # build up output data sequentially
        if pat.gene_names is not None:
            for mut_idx in sorted(phylogeny.conflicting_mutations, key=lambda k: pat.gene_names[k].lower()):

                self.file.write(
                    self._inds[self._ind] + '<tr><td><em>{}</em> ({})</td> {} </tr>\n'.format(
                        pat.gene_names[mut_idx], pat.mut_keys[mut_idx],
                        ' '.join('<td> {}/{} </td>'.format(
                            pat.mut_reads[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]],
                            pat.coverage[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]]) for sa_idx, sa_name in
                            enumerate(sorted(pat.sample_names)))))
        else:
            for mut_idx in sorted(phylogeny.conflicting_mutations, key=lambda k: pat.mut_keys[k].lower()):
                self.file.write(
                    self._inds[self._ind] + '<tr><td>{}</td> {} </tr>\n'.format(
                        pat.mut_keys[mut_idx], ' '.join('<td> {}/{} </td>'.format(
                            pat.mut_reads[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]],
                            pat.coverage[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]]) for sa_idx, sa_name in
                            enumerate(sorted(pat.sample_names)))))

        self._ind -= 1      # indentation level decreases by 1
        self.file.write(self._inds[self._ind]+'</table>\n')

        # show illustrative plot of the sequencing data of the incompatible mutation patterns
        if incomp_mps_plot_filepath is not None:
            self.file.write(self._inds[self._ind]+'<div align="center">\n')
            self.file.write(self._inds[self._ind]+'<div style="width:98%;max-width:700px">\n')
            self._ind += 1      # indentation level increases by 1

            self.file.write(self._inds[self._ind]+'<figure>\n')
            self._ind += 1      # indentation level increases by 1
            self.file.write(self._inds[self._ind]+'<img class="img-responsive" src="'+incomp_mps_plot_filepath +
                            '" alt="Sequencing data plot of incompatible mutation patterns" width="{}"/>'.format(
                                plot_width)+'\n')
            self.file.write(self._inds[self._ind]+'<div align="left">\n')
            self.file.write(self._inds[self._ind]+'<figcaption>\n')
            self.file.write(self._inds[self._ind]
                            + '<b>Incompatible mutation patterns in patient {}.</b>\n'.format(pat.name))

            present_mutations = pat.present_mutations
            self.file.write(self._inds[self._ind]
                            + 'Treeomics identified {} incompatible mutation patterns '.format(
                len(phylogeny.conflicting_mutations)) + '(out of {} investigated variants; {:.1%}).'.format(
                len(present_mutations), float(len(phylogeny.conflicting_mutations)) / len(present_mutations))+'\n')

            self.file.write(self._inds[self._ind]
                            + ' The color of the border of each rectangle representing a variant illustrates the '
                              'original classification, '
                            + 'the color of the left bar within each rectangle illustrates the VAF, '
                            + 'and the color of the right bar illustrates the coverage.  \n')
            self.file.write(self._inds[self._ind]
                            + 'Blue borders correspond to variants classified as present, red absent variants, '
                              'and light red unknown mutation status.\n')

            self.file.write(self._inds[self._ind]+'</figcaption>\n')
            self.file.write(self._inds[self._ind]+'</div>\n')
            self._ind -= 1      # indentation level decreases by 1
            self.file.write(self._inds[self._ind]+'</figure>\n')

            self._ind -= 1      # indentation level decreases by 1
            self.file.write(self._inds[self._ind]+'</div>\n')
            self.file.write(self._inds[self._ind]+'</div>\n')

    def add_artifacts_information(self, phylogeny, artifacts_plot_filepath=None, plot_width=None):
        """
        Add table with sequencing data information of the putative artifacts and add illustrative plot
        about the VAF and coverage in these artifacts if a path is given.
        :param phylogeny: data structure around the inferred phylogenetic tree
        :param artifacts_plot_filepath: illustrative sequencing data plot for the putative artifacts
        :param plot_width: width of the plot
        """

        self.file.write(self._inds[self._ind]+'<h4>Data artifacts</h4>\n')

        pat = phylogeny.patient
        opt_sol = phylogeny.solutions[0]    # optimal solution

        # add information about putative sequencing errors
        no_fps = sum(len(fps) for mut_idx, fps in opt_sol.false_positives.items())
        no_fns = sum(len(fns) for mut_idx, fns in opt_sol.false_negatives.items())
        no_putative_artifacts = no_fps+no_fns

        # any evolutionarily incompatible variants due to limited search space
        if opt_sol.conflicting_mutations is not None and len(opt_sol.conflicting_mutations) > 0:
            self.file.write(self._inds[self._ind]+'<h5>Unclassified artifacts due to limited solution space:</h5>\n')
            self.file.write(self._inds[self._ind]+'<ul>\n')
            self._ind += 1      # indentation level increases by 1

            for mut_idx in sorted(opt_sol.conflicting_mutations,
                                  key=lambda x: pat.gene_names[x].lower() if pat.gene_names is not None
                                  else pat.mut_keys[x]):
                self.file.write(
                    self._inds[self._ind] + '<li><em>{}</em> {} in samples: {} </li>\n'.format(
                        pat.gene_names[mut_idx] if pat.gene_names is not None else pat.mut_keys[mut_idx],
                        '({})'.format(pat.mut_keys[mut_idx]) if pat.gene_names is not None else '',
                        ', '.join('{} (reads: {}/{})'.format(
                            pat.sample_names[sa_idx].replace('_', ' '),
                            pat.mut_reads[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]],
                            pat.coverage[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]]) for sa_idx in
                            sorted(range(len(pat.sample_names)), key=lambda x: pat.sample_names[x]))))
                # if math.exp(pat.log_p01[mut_idx][sa_idx][1]) > 0.5)))

            self._ind -= 1      # indentation level decreases by 1
            self.file.write(self._inds[self._ind]+'</ul></br>\n')

        # any putative false-positives (sequencing errors)?
        if len(opt_sol.false_positives.keys()) > 0:
            self.file.write(self._inds[self._ind]+'<h5>Putative false-positives:</h5>\n')
            self.file.write(self._inds[self._ind]+'<ul>\n')
            self._ind += 1      # indentation level increases by 1

            for mut_idx, samples in sorted(opt_sol.false_positives.items(),
                                           key=lambda x: pat.gene_names[x[0]].lower() if pat.gene_names is not None
                                           else pat.mut_keys[x[0]]):
                self.file.write(
                    self._inds[self._ind] + '<li><em>{}</em> {} in samples: {} </li>\n'.format(
                        pat.gene_names[mut_idx] if pat.gene_names is not None else pat.mut_keys[mut_idx],
                        '({})'.format(pat.mut_keys[mut_idx]) if pat.gene_names is not None else '',
                        ', '.join('{} (reads: {}/{})'.format(
                            pat.sample_names[sa_idx].replace('_', ' '),
                            pat.mut_reads[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]],
                            pat.coverage[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]]) for sa_idx in
                            sorted(samples, key=lambda x: pat.sample_names[x])
                            if sa_idx in pat.mutations[mut_idx])))

            self._ind -= 1      # indentation level decreases by 1
            self.file.write(self._inds[self._ind]+'</ul></br>\n')

        # any information about putative lost variants?
        if len(opt_sol.false_negatives.keys()) > 0:
            self.file.write(self._inds[self._ind]+'<h5>Putative lost variants:</h5>\n')
            self.file.write(self._inds[self._ind]+'<ul>\n')
            self._ind += 1      # indentation level increases by 1

            for mut_idx, samples in sorted(opt_sol.false_negatives.items(),
                                           key=lambda x: pat.gene_names[x[0]].lower() if pat.gene_names is not None
                                           else pat.mut_keys[x[0]]):
                self.file.write(
                    self._inds[self._ind]+'<li><em>{}</em> {} in samples: {} </li>\n'.format(
                        pat.gene_names[mut_idx] if pat.gene_names is not None else pat.mut_keys[mut_idx],
                        '({})'.format(pat.mut_keys[mut_idx]) if pat.gene_names is not None else '',
                        ', '.join('{} (reads: {}/{})'.format(
                            pat.sample_names[sa_idx].replace('_', ' '),
                            pat.mut_reads[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]],
                            pat.coverage[pat.mut_keys[mut_idx]][pat.sample_names[sa_idx]]) for sa_idx
                            in sorted(samples, key=lambda x: pat.sample_names[x]))))

            self._ind -= 1      # indentation level decreases by 1
            self.file.write(self._inds[self._ind]+'</ul>\n')

        if artifacts_plot_filepath is not None:
            present_mutations = pat.present_mutations
            if len(opt_sol.false_positives.keys()) + len(opt_sol.false_negatives.keys()) > 0:
                self.file.write(self._inds[self._ind]+'<div align="center">\n')
                self.file.write(self._inds[self._ind]+'<div style="width:98%;max-width:700px">\n')
                self._ind += 1      # indentation level increases by 1

                self.file.write(self._inds[self._ind]+'<figure>\n')
                self._ind += 1      # indentation level increases by 1
                self.file.write(self._inds[self._ind]+'<img class="img-responsive" src="'+artifacts_plot_filepath +
                                '" alt="Mutation patterns of putative artifacts" width="{}"/>'.format(plot_width)+'\n')
                self.file.write(self._inds[self._ind]+'<div align="left">\n')
                self.file.write(self._inds[self._ind]+'<figcaption>\n')
                self.file.write(self._inds[self._ind]
                                + '<b>Mutation patterns of putative artifacts in patient {}.</b>\n'.format(pat.name))

                self.file.write(self._inds[self._ind]
                                + 'Treeomics identified {} putative artifacts '.format(no_putative_artifacts)
                                + '(out of {} investigated variants; {:.1%}).'.format(
                                len(pat.sample_names) * len(present_mutations),
                                float(no_putative_artifacts) / (len(pat.sample_names) * len(present_mutations)))+'\n')
                self.file.write(
                    self._inds[self._ind] +
                    'Additionally there were {} putative false-negatives due to insufficient '.format(
                        sum(len(fns) for mut_idx, fns in opt_sol.false_negative_unknowns.items()))
                    + 'coverage (unknowns; data not shown). \n')

                self.file.write(self._inds[self._ind]
                                + ' The color of the border of each rectangle representing a variant illustrates the '
                                  'original classification, '
                                + 'the color of the left bar within each rectangle illustrates the VAF, '
                                + 'and the color of the right bar illustrates the coverage. '
                                + 'If a variant was identified as a putative artifact, a smaller rectangle '
                                  'with the changed classification color is added on top of the bars. \n')
                self.file.write(self._inds[self._ind]
                                + 'Blue borders correspond to variants classified as present, red absent variants, '
                                  'and light red unknown mutation status.\n')

                self.file.write(self._inds[self._ind]+'</figcaption>\n')
                self.file.write(self._inds[self._ind]+'</div>\n')
                self._ind -= 1      # indentation level decreases by 1
                self.file.write(self._inds[self._ind]+'</figure>\n')

                self._ind -= 1      # indentation level decreases by 1
                self.file.write(self._inds[self._ind]+'</div>\n')
                self.file.write(self._inds[self._ind]+'</div>\n')

            else:       # there are no putative well-powered artifacts
                self.file.write(self._inds[self._ind]+'<p>\n')
                self._ind += 1      # indentation level increases by 1
                self.file.write(self._inds[self._ind] + 'Treeomics identified no well-powered artifacts '
                                + '(out of {} investigated variants; {:.1%}).'.format(
                                len(pat.sample_names) * len(present_mutations),
                                float(no_putative_artifacts) / (len(pat.sample_names) * len(present_mutations)))+'\n')
                self.file.write(
                    self._inds[self._ind] +
                    'There were {} putative false-negatives due to insufficient '.format(
                        sum(len(fns) for mut_idx, fns in opt_sol.false_negative_unknowns.items()))
                    + 'coverage (unknowns; data not shown). \n')

                self._ind -= 1      # indentation level decreases by 1
                self.file.write(self._inds[self._ind]+'</p>\n')

    def end_report(self, e, c0, max_absent_vaf, loh_frequency, fpr, fdr,
                   min_absent_cov, min_median_cov, min_median_maf, max_no_mps=None):
        """
        Add used parameter values, HTML file footer, close the file, and create PDF file
        :param e: sequencing error rate for bayesian inference
        :param c0: prior mixture parameter of delta function and uniform distribution for bayesian inference
        :param max_absent_vaf: maximal absent VAF before considering estimated purity
        :param loh_frequency: probability that a SNV along a lineage is lost due loss of heterozygosity
        :param fpr: false-positive rate
        :param fdr: false-discovery rate
        :param min_absent_cov: Minimum coverage for a variant to be called absent
        :param min_median_cov: minimum median coverage of a sample to pass filtering
        :param min_median_maf: minimum median variant allele frequency of a sample to pass filtering
        :param max_no_mps: maximal number of MPs per variant that are considered in the MILP; by default the full
                           solution space is considered and hence 2^(#samples) of MPs are generated
        """

        if self.file is not None:

            self.file.write(self._inds[self._ind]+'</br><p><em>Treeomics settings:</em>\n')
            self.file.write(
                self._inds[self._ind] + ' sequencing error rate e: {}, '.format(e)
                + 'prior absent probability c0: {}, '.format(c0)
                + 'max absent VAF: {}, '.format(max_absent_vaf)
                + 'LOH frequency: {}, '.format(loh_frequency)
                + 'false discovery rate: {}, '.format(fdr)
                + 'false-positive rate: {}. '.format(fpr)
                + ('absent classification minimum coverage: {}. '.format(min_absent_cov) if min_absent_cov > 0 else '')
                + ('sample minimal median coverage: {}. '.format(min_median_cov) if min_median_cov > 0 else '')
                + ('sample minimal median VAF: {}. '.format(min_median_maf) if min_median_maf > 0 else '')
                + ('explored mut patterns per variant: {}.'.format(max_no_mps) if max_no_mps is not None else '')
                + ('filter: minimum VAF per variant in at least one sample {}. '.format(settings.MIN_VAF)
                    if (settings.MIN_VAF is not None and settings.MIN_VAF > 0) else '')
                + ('filter: minimum number of variant reads per variant in at least one sample {}. '.format(
                    settings.MIN_VAR_READS) if settings.MIN_VAR_READS is not None and settings.MIN_VAR_READS > 0
                   else ''))

            self.file.write(self._inds[self._ind]+'</p>\n')

            self._ind -= 1      # indentation level decreases by 1
            self.file.write(self._inds[self._ind]+'</div> <!-- /container -->\n')

            # add footer
            today = datetime.date.today()
            self.file.write(self._inds[self._ind]+'<footer class="footer">\n')
            self._ind += 1      # indentation level increases by 1
            self.file.write(self._inds[self._ind]+'<div class="container">\n')
            self._ind += 1      # indentation level increases by 1
            self.file.write(self._inds[self._ind]+'<p class="text-muted">&copy; Treeomics {}, '.format(
                __version__) + today.strftime('%b %d, %Y')+'</p>\n')
            self._ind -= 1      # indentation level decreases by 1
            self.file.write(self._inds[self._ind]+'</div>\n')
            self._ind -= 1      # indentation level decreases by 1
            self.file.write(self._inds[self._ind]+'</footer>\n')

            self._ind -= 1      # indentation level decreases by 1
            self.file.write(self._inds[self._ind]+'</body>\n')
            self._ind -= 1      # indentation level decreases by 1
            self.file.write(self._inds[self._ind]+'</html>\n')

            self.file.flush()
            self.file.close()

            self.file = None
            logger.info('Successfully created HTML report at {}'.format(self.filepath))

            # try to create a PDF report from the HTML document
            self.create_pdf()

    def create_pdf(self):
        """
        Convert HTML report to PDF
        """

        try:  # check if varcode and pyensembl is available (necessary for Windows)
            import pdfkit

        except ImportError:
            logger.warning('PDFKIT not available and hence no PDF of the HTML report was created.')
            return

        pdf_filepath = self.filepath.replace('.html', '.pdf')

        options = {
            'zoom': settings.ZOOM,
            'dpi': 400,
            # 'print-media-type': '',
            'page-size': 'Letter',
            'margin-top': '1in',
            'margin-right': '0.75in',
            'margin-bottom': '1in',
            'margin-left': '0.75in',
            'disable-smart-shrinking': '',
            'enable-local-file-access': '',
            'quiet': ''

            # 'encoding': "UTF-8",
            # 'custom-header': [
            #     ('Accept-Encoding', 'gzip')
            # ]
            # 'cookie': [
            #     ('cookie-name1', 'cookie-value1'),
            #     ('cookie-name2', 'cookie-value2'),
            # ],
            # 'no-outline': None
        }

        pdfkit.from_file(self.filepath, pdf_filepath, options=options)
        logger.info('Successfully created PDF report at {}'.format(os.path.abspath(pdf_filepath)))
