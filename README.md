## Treeomics: Decrypting somatic mutation patterns to reveal the evolution of cancer
Treeomics reconstructs the phylogeny of a cancer with commonly available sequencing technologies.
 
========

### The easiest way to run Treeomics:
1. Open a terminal and clone the repository from GitHub with ```git clone https://github.com/johannesreiter/treeomics.git```
2. Go into the new folder with ```cd treeomics\src```
3. Type the following command to run the simulation: ```python treeomics -r <mut-reads table> -s <coverage table> -O``` 
where ```<mut-reads table>``` is the path to a tab-separated-value file with the number of 
reads reporting a variant (row) in each sample (column) and ```<coverage table>``` is the path to a tab-separated-value 
file with the sequencing depth at the position of this variant in each sample.

Usage: ```python treeomics -r <mut-reads table> -s <coverage table> | -v vcf_file | "
                "-d vcf_file_directory  [-n normal_sample_name] -O```

##### Examples:
Example 1: ```python treeomics -r input/Makohon2015/Pam03_mutant_reads.txt -s input/Makohon2015/Pam03_phredcoverage.txt```  
Reconstructs the phylogeny of pancreatic cancer patient Pam03 based on targeted sequencing data 
of 5 distinct liver metastases, 3 distinct lung metastases, and 2 samples of the primary tumor.

Example 2: ```-r input/Bashashati2013/Case5_mutant_reads.txt -s input/Bashashati2013/Case5_coverage.txt```
Reconstructs the phylogeny of the high-grade serous ovarian cancer of Case 5 in Bashashati et al. (2013).

========


Author: Johannes Reiter, IST Austria, [http://pub.ist.ac.at/~jreiter](http://pub.ist.ac.at/~jreiter)  
Latest version can be found here: [https://github.com/johannesreiter/treeomics](https://github.com/johannesreiter/treeomics)

=========
