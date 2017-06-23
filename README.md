## Treeomics: Reconstructing metastatic seeding patterns of human cancers
Developed by: JG Reiter, AP Makohon-Moore, JM Gerold, I Bozic, K Chatterjee, C Iacobuzio-Donahue, B Vogelstein, MA Nowak.

 
========

### What is Treeomics?
Treeomics is a computational tool to reconstruct the phylogeny of metastases with commonly available sequencing technologies.
The tool detects putative artifacts in noisy sequencing data and infers robust evolutionary trees across a variety of evaluated scenarios.
For more details, see our publication *Reconstructing metastatic seeding patterns of human cancers* (Nature Communications, 8, 14114, [http://dx.doi.org/10.1038/ncomms14114](http://dx.doi.org/10.1038/ncomms14114)).

<img align="middle" src="repository_illustration.png">

* <a href="#installation">Installation</a>
* <a href="#getting">Getting started with Treeomics</a>
* <a href="#examples">Examples</a>

#### <a name="releases"> Releases
* Treeomics 1.5.2 2016-10-18: Initial release with acceptance of the manuscript.
* Treeomics 1.6.0 2016-12-09: Improves visualization of generated evolutionary trees by integrating ETE3. ILP solver explores a pool of the best solutions to more efficiently assess the support of the inferred branches.
* Treeomics 1.7.0 2017-02-09: Uses Bayesian inference model for similarity and artifact analyses.
* Treeomics 1.7.1 2017-02-23: Integrated python packages ```pyensembl``` and ```varcode``` to infer the gene names where variants occurred as well as their mutation effect.
* Treeomics 1.7.2 2017-03-02: Improved visualization of predicted driver genes in HTML report and the mutation table.
* Treeomics 1.7.3 2017-03-13: Visualize the 5 most likely evolutionary trees. Improve solution pool usage to better estimate confidence values.
* Treeomics 1.7.4 2017-03-15: Make mutation effect prediction by VarCode optional to reduce dependencies for users.
* Treeomics 1.7.5 2017-04-11: Improved putative driver gene analysis and HTML report. Allow multiple normal samples. Implemented optional filter of common normal variants.
* Treeomics 1.7.6 2017-05-12: Generate new out put file ```<subject>_variants.csv``` with information about the individual variants and how they were classified in the inferred phylogeny. Solved issues with subclone detection and solution pool.
* Treeomics 1.7.7 2017-06-21: Made Treeomics ready for ultra deep targeted sequencing data. Fixed bug in calculation of branch confidence values in partial solution space. Use wkhtmltopdf to create a PDF from the HTML report.

### <a name="installation"> Installation
1. Open a terminal and clone the repository from GitHub with ```git clone https://github.com/johannesreiter/treeomics.git```
2. Install required packages:
  - Install Python 3.4 ([https://www.python.org/downloads](https://www.python.org/downloads); CPLEX explicitly requires Python 3.4)
  - Install NumPy ([http://www.numpy.org](http://www.numpy.org)), 
    SciPy ([http://www.numpy.org](http://www.numpy.org))
  - Install networkx ([https://networkx.github.io/](https://networkx.github.io/))
  - Install matplotlib 1.4 or 1.5 (matplotlib 2 can cause various problems with the layout; [http://matplotlib.org](http://matplotlib.org/))
  - Install pandas ([http://pandas.pydata.org/](http://pandas.pydata.org/))
  - Install seaborn ([https://stanford.edu/~mwaskom/software/seaborn/](https://stanford.edu/~mwaskom/software/seaborn/))
  - Install the IBM ILOG CPLEX Optimization Studio ([http://www-01.ibm.com/support/docview.wss?uid=swg21444285](http://www-01.ibm.com/support/docview.wss?uid=swg21444285))
    and then setup the Python API ([https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.6.3/ilog.odms.cplex.help/CPLEX/GettingStarted/topics/set_up/Python_setup.html](https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.6.3/ilog.odms.cplex.help/CPLEX/GettingStarted/topics/set_up/Python_setup.html));
    An IBM Academic License to freely download CPLEX can be obtained here: [http://www-304.ibm.com/ibm/university/academic/pub/page/academic_initiative](http://www-304.ibm.com/ibm/university/academic/pub/page/academic_initiative)

3. Install optional packages:
  - To create a PDF from the HTML report, install wkhtmltopdf ([https://wkhtmltopdf.org](https://wkhtmltopdf.org)) and pdfkit ([https://github.com/JazzCore/python-pdfkit](https://github.com/JazzCore/python-pdfkit))
  - To automatically generate evolutionary conflict graphs, install circos ((with ```circos``` in your ```PATH``` environment variable; [http://circos.ca/software/installation](http://circos.ca/software/installation))
  - For automatically generating evolutionary tree plots, install LaTeX/TikZ (with ```pdflatex``` in your ```PATH``` environment variable;
    [https://www.tug.org/texlive/quickinstall.html](https://www.tug.org/texlive/quickinstall.html)) and/or ETE3 [https://github.com/etetoolkit/ete](https://github.com/etetoolkit/ete) (installing ETE3 can be tricky; we recommend using Anaconda [https://www.continuum.io](https://www.continuum.io)
  - For annotating only non-synonymous variants in driver genes, install pyensembl ([https://github.com/hammerlab/pyensembl](https://github.com/hammerlab/pyensembl)) and varcode ([https://github.com/hammerlab/varcode](https://github.com/hammerlab/varcode))
    
### <a name="getting"> Getting started with Treeomics
1. Input files: The input to ```__main__.py``` is either
  - two tab-delimited text files -- one for variant read data and one for coverage data. Please see the files ```input/Makohon2017/Pam03_mutant_reads.txt``` and ```input/Makohon2017/Pam03_phredcoverage.txt``` included in this repository for examples.
  - VCF-files of all samples
2. Go into the new folder with ```cd treeomics/src```
3. Type the following command to run the simulation: ```python treeomics -r <mut-reads table> -s <coverage table> -O``` 
where ```<mut-reads table>``` is the path to a tab-separated-value file with the number of 
reads reporting a variant (row) in each sample (column) and ```<coverage table>``` is the path to a tab-separated-value 
file with the sequencing depth at the position of this variant in each sample.

##### Usage: 
```shell
$ python treeomics -r <mut-reads table> -s <coverage table> | -v <vcf file> | -d <vcf file directory> -O
```

##### Optional parameters:
- *-e <sequencing error rate>:* Sequencing error rate *e* in the Bayesian inference model (default 1.0%)
- *-a <max absent VAF>:* Maximum VAF for an absent variant *f<sub>absent</sub>* before considering the estimated purity (default 5%)
- *-z <prior absent probability>:* Prior probability for a variant being absent *c<sub>0</sub> (default 0.5).
- *-o <output directory>:* Provide different output directory (default src/output)
- *-n <normal sample names>:* If a normal sample is provided, variants significantly present in the normal are removed. Additional normal samples (or other samples that should be ignored) can be provided via a space-separated enumeration. E.g. ```-n FIRSTNORMALSAMPLE SECONDNORMALSAMPLE ```
- ```--pool_size``` *<Pool size of ILP solver>:* Number of best solutions explored by ILP solver to assess the support of the inferred branches (default 1000)
- *-b <No bootstrapping samples>:* Number of bootstrapping samples (default 0); Generally using the solution pool instead of bootstrapping seems to be the more efficient way to assess confidence.
- *-u:* Enables subclone detection (default ```False```)
<!--- - ```--no_subclone_detection``` Disables subclone detection) -->
- *-c <min median coverage>:* Minimum median coverage of a sample to be considered (default 0)
- *-f <min median vaf>:* Minimum median mutant allele frequency of a sample to be considered (default 0)
- *-p <false positive rate>:* False-positive rate of conventional binary classification (only relevant for artifact comparison)
- *-i <false discovery rate>:* Targeted false-discovery rate of conventional binary classification  (only relevant for artifact comparison)
- *-y <min absent coverage>:* Minimum coverage for a powered absent variant  (only relevant for artifact comparison)
- *-t <time limit>:* Maximum running time for CPLEX to solve the MILP (in seconds, default ```None```). If not ```None```, the obtained solution is no longer guaranteed to be optimal
- *-l <max no MPS>:* Maximum number of considered mutation patterns per variant (default ```None```). If not ```None```, the obtained solution is no longer guaranteed to be optimal
- ```--wes_filtering``` Removes intronic and intergenic variants in WES data; default ```False```)
- ```--common_vars_file``` Path to file with common variants in normal samples and therefore removed from analysis (default ```None```)
- ```--no_plots``` Disables generation of plots (useful for benchmarking; default ```True```)
- ```--benchmarking``` Generates mutation matrix and mutation pattern files that can be used for automatic benchmarking of silico data (default ```False```)

Default parameter values as well as output directory can be changed in ```treeomics/src/treeomics/settings.py```.
Moreover, the ```settings.py``` provides more options an annotation of driver genes and configuration of plot output names. 
All plots, analysis and logging files, and the HTML report will be in this output directory.

##### Optional input:
- *Driver gene annotation:* Treeomics highlights any non-synonymous or splice-site variants (if VarCode is available, otherwise all) in putative driver genes given in a CSV-file under ```DRIVER_PATH``` in ```treeomics/src/treeomics/settings.py```. As default list, the union of reported driver genes by 20/20+, TUSON, and MutsigCV from Tokheim et al. (PNAS, 2016) is used (see ```treeomics/src/input/Tokheim_drivers_union.csv```). Any CSV-file can be used as long as there is column named 'Gene_Symbol'. Variants in these provided genes will be highlighted in the HTML report as well as in the inferred phylogeny.
- *Cancer Gene Census (CGC) annotation:* Variants that have been identified as likely drivers in the provided genes (under ```DRIVER_PATH```) will be check if they occurred in the reported region in the given CSV-file (default ```treeomics/src/input/cancer_gene_census_grch37_v80.csv```; CGC version 80, reference genome hg19).

### <a name="examples"> Examples
Example 1:
```shell
$ python treeomics -r input/Makohon2017/Pam03_mutant_reads.txt -s input/Makohon2017/Pam03_phredcoverage.txt -e 0.005 -O
```
Reconstructs the phylogeny of pancreatic cancer patient Pam03 based on targeted sequencing data 
of 5 distinct liver metastases, 3 distinct lung metastases, and 2 samples of the primary tumor.

Example 2:
```shell
$ python treeomics -r input/Bashashati2013/Case5_mutant_reads.txt -s input/Bashashati2013/Case5_coverage.txt -e 0.005 -O
```
Reconstructs the phylogeny of the high-grade serous ovarian cancer of Case 5 in Bashashati et al. (2013).

========

### Problems?
If you have any questions, you can contact us ([https://github.com/johannesreiter](https://github.com/johannesreiter)) and we will try to help.


### License
Copyright (C) 2016 Johannes Reiter

Treeomics is licensed under the GNU General Public License, Version 3.
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, 
version 3 of the License.
There is no warranty for this free software.

========

Author: Johannes Reiter, Harvard University, [http://www.people.fas.harvard.edu/~reiter](http://www.people.fas.harvard.edu/~reiter)  
