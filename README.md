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
* Treeomics 1.7.8 2017-10-06: Fixed problem with ete3 visualization of detected subclones. Added additional command line parameters: path to CSV file to highlight given genes in inferred phylogeny and set the maximal number of used threads by CPLEX.
* Treeomics 1.7.9 2017-10-10: Configure the number of top ranked solution trees that are plotted.
* Treeomics 1.7.10 2018-05-15: Improved PDF-report generation. Added support for structural variants. Added support for providing externally estimated sample purities via ```--purities <SAMPLE NAMES>```. Added ```--verbose``` option to run Treeomics in DEBUG logging level. Fixed VCF parsing error thanks to Frank's bug report.
* Treeomics 1.7.11 2018-10-26: Added TCGA consensus driver gene list from Bailey et al., Cell 2018. Added Zoom parameter to ```settings.py``` to better configure PDF report appearance. 
* Treeomics 1.7.12 2018-11-26: Replaced 'nodes_iter()' and 'edges_iter()' calls with 'nodes()' and 'edges()' calls because networkx 2.0+ has no backward compatibility: https://stackoverflow.com/questions/33734836/graph-object-has-no-attribute-nodes-iter-in-networkx-module-python


### <a name="installation"> Installation
1. Open a terminal and clone the repository from GitHub with ```git clone https://github.com/johannesreiter/treeomics.git```
2. Install required packages:
  - Install Python 3.4 ([https://www.python.org/downloads](https://www.python.org/downloads))
  - Install NumPy ([http://www.numpy.org](http://www.numpy.org)), 
    SciPy ([http://www.numpy.org](http://www.numpy.org))
  - Install networkx 2.0 or above ([https://networkx.github.io/](https://networkx.github.io/))
  - Install matplotlib 1.4 or 1.5 (matplotlib 2 can cause various problems with the layout; [http://matplotlib.org](http://matplotlib.org/))
  - Install pandas ([http://pandas.pydata.org/](http://pandas.pydata.org/))
  - Install seaborn ([http://seaborn.pydata.org/](http://seaborn.pydata.org/))
  - Install the IBM ILOG CPLEX Optimization Studio 12.6.3 (or higher) ([http://www-01.ibm.com/support/docview.wss?uid=swg21444285](http://www-01.ibm.com/support/docview.wss?uid=swg21444285))
    and then setup the Python API ([https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.6.3/ilog.odms.cplex.help/CPLEX/GettingStarted/topics/set_up/Python_setup.html](https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.6.3/ilog.odms.cplex.help/CPLEX/GettingStarted/topics/set_up/Python_setup.html));
    An IBM Academic License to freely download CPLEX can be obtained here: [http://www-304.ibm.com/ibm/university/academic/pub/page/academic_initiative](http://www-304.ibm.com/ibm/university/academic/pub/page/academic_initiative). 
    In the new version 12.7.1 apparently Python 3.4 is no longer officially supported, however, cplex seems to work nicely after updating two files and changing the version check from ```(3, 5, 0)``` to ```(3, 4, 0)``` in both the ```setup.py``` (in MacOS here: ```Applications/IBM/ILOG/CPLEX_Studio1271/cplex/python/3.5/x86-64_osx/```) to install cplex in Python 3.4 as well as in your miniconda installation ```miniconda3/lib/python3.4/site-packages/cplex/_internal/_pycplex_platform.py```. You will see where an exception is thrown if you test your installation with ```python3 -c 'import cplex'```.
    You may also need to add cplex to your ```PYTHONPATH``` with: ```export PYTHONPATH="~/Applications/IBM/ILOG/CPLEX_Studio1271/cplex/python/3.5/x86-64_osx/:$PYTHONPATH"```

3. Install optional packages:
  - To create a PDF from the HTML report, install wkhtmltopdf ([https://wkhtmltopdf.org](https://wkhtmltopdf.org)) and pdfkit ([https://github.com/JazzCore/python-pdfkit](https://github.com/JazzCore/python-pdfkit))
  - To automatically generate evolutionary conflict graphs, install circos (with ```circos``` in your ```PATH``` environment variable; [http://circos.ca/software/installation](http://circos.ca/software/installation))
  - For automatically generating evolutionary tree plots, install LaTeX/TikZ (with ```pdflatex``` in your ```PATH``` environment variable;
    [https://www.tug.org/texlive/quickinstall.html](https://www.tug.org/texlive/quickinstall.html)) and/or ETE3 [https://github.com/etetoolkit/ete](https://github.com/etetoolkit/ete) (installing ETE3 can be very frustrating, in particular in Python 3.5+ as it  requires Qt4; we recommend using Python 3.4 and Miniconda [https://www.continuum.io](https://www.continuum.io): ```conda install python=3.4 qt=4``` and then install ete3 ```conda install -c etetoolkit ete3 ete3_external_apps```. You can test your installation with ```python3 -c 'from ete3 import TreeStyle'```.
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
- *-n <normal sample names>:* If a normal sample is provided, variants significantly present in the normal are removed. Additional normal samples can be provided via a space-separated enumeration. E.g. ```-n FIRSTNORMALSAMPLE SECONDNORMALSAMPLE ```
- *-x <samples to exclude>:* Space-separated enumeration of sample names to exclude from the analysis. E.g. ```-x FIRSTEXCLUDEDSAMPLE SECONDEXCLUDEDSAMPLE ```
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
- ```--threads=<N>``` Maximal number of parallel threads that will be invoked by CPLEX (```0```: default, let CPLEX decide; ```1```: single threaded; ```N```: uses up to N threads)

- *-l <max no MPS>:* Maximum number of considered mutation patterns per variant (default ```None```). If not ```None```, the obtained solution is no longer guaranteed to be optimal
- ```--driver_genes=<path to file>``` Path to CSV file with names of putative driver genes highlighted in inferred phylogeny (default ```--driver_genes=../input/Tokheim_drivers_union.csv```)
- ```--wes_filtering``` Removes intronic and intergenic variants in WES data (default ```False```)
- ```--common_vars_file``` Path to file with common variants in normal samples and therefore removed from analysis (default ```None```)
- ```--no_plots``` Disables generation of X11 depending plots (useful for benchmarking; default plots are generated ```plots```)
- ```--no_tikztrees``` Disables generation of latex trees which do not depend on X11 (default latex trees are generated ```tikztrees```)
- ```--benchmarking``` Generates mutation matrix and mutation pattern files that can be used for automatic benchmarking of silico data (default ```False```)
- ```--include``` Provide a list of sample names that should be analyzed (e.g., ```--include PT1 PT2 PT3 PT4```)
- ```--purities``` Provide a list of externally estimated sample purities (e.g., ```--purities 0.7 0.3 0.9 0.8```). Requires ```--include``` argument with the same ordering of samples.


Default parameter values as well as output directory can be changed in ```treeomics/src/treeomics/settings.py```.
Moreover, the ```settings.py``` provides more options an annotation of driver genes and configuration of plot output names. 
All plots, analysis and logging files, and the HTML report will be in this output directory.

##### Optional input:
- *Driver gene annotation:* Treeomics highlights any non-synonymous or splice-site variants (if VarCode is available, otherwise all) in putative driver genes given in a CSV-file under ```DRIVER_PATH``` in ```treeomics/src/treeomics/settings.py```. As default list, the union of reported driver genes by 20/20+, TUSON, and MutsigCV from Tokheim et al. (PNAS, 2016) is used (see ```treeomics/src/input/Tokheim_drivers_union.csv```). Any CSV-file can be used as long as there is column named 'Gene_Symbol'. Variants in these provided genes will be highlighted in the HTML report as well as in the inferred phylogeny.
- *Cancer Gene Census (CGC) annotation:* Variants that have been identified as likely drivers in the provided genes (under ```DRIVER_PATH```) will be check if they occurred in the reported region in the given CSV-file (default ```treeomics/src/input/cancer_gene_census_grch37_v80.csv```; CGC version 80, reference genome hg19).

### <a name="examples"> Examples
Example 1:
```shell
$ python treeomics -r input/Makohon2017/Pam03_1-10_mutant_reads.txt -s input/Makohon2017/Pam03_1-10_phredcoverage.txt -n Pam03N3 -e 0.005 -O
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
Copyright (C) 2017 Johannes Reiter

Treeomics is licensed under the GNU General Public License, Version 3.
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, 
version 3 of the License.
There is no warranty for this free software.

========

Author: Johannes Reiter, Harvard University, [http://www.people.fas.harvard.edu/~reiter](http://www.people.fas.harvard.edu/~reiter)  
