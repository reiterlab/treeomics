## Treeomics: Decrypting somatic mutation patterns to reveal the evolution of metastatic cancer
Developed by: JG Reiter<sup>1,2,3,4</sup>, A Makohon-Moore<sup>5,6</sup>, J Gerold<sup>1</sup>, I Bozic<sup>1,7</sup>, K Chatterjee<sup>2</sup>, C Iacobuzio-Donahue<sup>5,6,8</sup>, B Vogelstein<sup>9,10</sup>, MA Nowak<sup>1,7,11</sup>.

<sup>1</sup> Program for Evolutionary Dynamics, Harvard University, Cambridge, MA, USA.
<sup>2</sup> IST (Institute of Science and Technology) Austria, Klosterneuburg, Austria.
<sup>3</sup> Dana-Farber Cancer Institute, Boston, MA, USA.
<sup>4</sup> Broad Institute of MIT and Harvard, Cambridge, MA, USA.
<sup>5</sup> The David M. Rubenstein Center for Pancreatic Cancer Research, Memorial Sloan Kettering Cancer Center, New York, USA.
<sup>6</sup> Human Oncology and Pathogenesis Program, Memorial Sloan Kettering Cancer Center, New York, USA.
<sup>7</sup> Department of Mathematics, Harvard University, Cambridge, MA, USA.
<sup>8</sup> Department of Pathology, Memorial Sloan Kettering Cancer Center, New York, USA.
<sup>9</sup> The Sol Goldman Pancreatic Cancer Research Center, Johns Hopkins University School of Medicine, Baltimore, MD, USA. 
<sup>10</sup> The Ludwig Center, Johns Hopkins University School of Medicine, Baltimore, MD, USA.
<sup>11</sup> Department of Organismic and Evolutionary Biology, Harvard University, Cambridge, MA, USA.
 
========

### What is Treeomics?
Treeomics is a computational tool to reconstruct the phylogeny of a cancer with commonly available sequencing technologies.
The tool detects putative artifacts in noisy sequencing data and can therefore infer more robust evolutionary trees.

![Observed tumor heterogeneity and inferred](repository_illustration.png "Observed tumor heterogeneity and inferred evolution among spatially-distinct DNA sequencing samples")    


#### Installation
1. Open a terminal and clone the repository from GitHub with ```git clone https://github.com/johannesreiter/treeomics.git```
2. Install required packages:
  - Install Python 3.4 ([https://www.python.org/downloads](https://www.python.org/downloads))
  - Install NumPy ([http://www.numpy.org](http://www.numpy.org)), 
    SciPy ([http://www.numpy.org](http://www.numpy.org))
  - Install the IBM ILOG CPLEX Optimization Studio ([http://www-01.ibm.com/support/docview.wss?uid=swg21444285](http://www-01.ibm.com/support/docview.wss?uid=swg21444285))
    and then setup the Python API ([http://www-01.ibm.com/support/knowledgecenter/SSSA5P_12.2.0/ilog.odms.cplex.help/Content/Optimization/Documentation/CPLEX/_pubskel/CPLEX20.html](http://www-01.ibm.com/support/knowledgecenter/SSSA5P_12.2.0/ilog.odms.cplex.help/Content/Optimization/Documentation/CPLEX/_pubskel/CPLEX20.html));
    An IBM Academic License can be obtained here: [http://www-304.ibm.com/ibm/university/academic/pub/page/academic_initiative](http://www-304.ibm.com/ibm/university/academic/pub/page/academic_initiative)
  - If you want any plots to be automatically generated, install also
    matplotlib ([http://matplotlib.org](http://matplotlib.org/)), LaTeX/TikZ (with ```pdflatex``` in your ```PATH``` environment variable; 
    [https://www.tug.org/texlive/quickinstall.html](https://www.tug.org/texlive/quickinstall.html)), circos ((with ```circos``` in your ```PATH``` environment variable; [http://circos.ca/software/installation](http://circos.ca/software/installation))
    
#### Getting started with Treeomics
1. Input files: The input to ```__main__.py``` is either
  - two tab-delimited text files -- one for variant read data and one for coverage data. Please see the files ```input/Makohon2015/Pam03_mutant_reads.txt``` and ```input/Makohon2015/Pam03_phredcoverage.txt``` included in this repository for examples.
  - VCF-files of all samples
2. Go into the new folder with ```cd treeomics\src```
3. Type the following command to run the simulation: ```python treeomics -r <mut-reads table> -s <coverage table> -O``` 
where ```<mut-reads table>``` is the path to a tab-separated-value file with the number of 
reads reporting a variant (row) in each sample (column) and ```<coverage table>``` is the path to a tab-separated-value 
file with the sequencing depth at the position of this variant in each sample.

##### Usage: 
```python treeomics -r <mut-reads table> -s <coverage table> | -v <vcf file> | -d <vcf file directory> [-n <normal sample name>] [-e <sequencing error rate] [-z <prior absent probability>] [-p false positive rate] [-i false discovery rate] -O```

##### Optional parameters:
- *-e <sequencing error rate>:* Error rate of the sequencing machine.
- *-z <prior absent probability>:* Prior probability for a variant being absent.
- *-p <false positive rate>:* False-positive rate of conventional binary classification.
- *-i <false discovery rate>:* Targeted false-discovery rate of conventional binary classification.

Default parameter values as well as output directory can be changed in ```treeomics\src\settings.py```.
All plots, analysis and logging files, and the HTML report will be in this output directory.

#### Examples
Example 1: ```python treeomics -r input/Makohon2015/Pam03_mutant_reads.txt -s input/Makohon2015/Pam03_phredcoverage.txt```  
Reconstructs the phylogeny of pancreatic cancer patient Pam03 based on targeted sequencing data 
of 5 distinct liver metastases, 3 distinct lung metastases, and 2 samples of the primary tumor.

Example 2: ```python treeomics -r input/Bashashati2013/Case5_mutant_reads.txt -s input/Bashashati2013/Case5_coverage.txt```
Reconstructs the phylogeny of the high-grade serous ovarian cancer of Case 5 in Bashashati et al. (2013).

========

### Problems?
If you have any questions, you can contact us ([https://github.com/johannesreiter](https://github.com/johannesreiter)) and we will try to help.


### License
Copyright (C) 2015 Johannes Reiter

Treeomics is licensed under the GNU General Public License, Version 3.
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, 
version 3 of the License.
There is no warranty for this free software.

========

Author: Johannes Reiter, IST Austria, [http://pub.ist.ac.at/~jreiter](http://pub.ist.ac.at/~jreiter)  
