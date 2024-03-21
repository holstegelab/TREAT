# TREAT (Tandem REpeat Annotation Toolkit)

![TREAT](/doc/Figure_TREAT_original.png)

## TREAT in a nutshell
**TREAT** is a command line tool written in **Python** and **R** (for plotting) that can be used to work with tandem repeats and structural variants from long-read sequencing data. **TREAT** was developed specifically for long-read sequencing data. However, it can potentially be used with any sequencing data, including PacBio, Oxford Nanopore, and Illumina.  
**TREAT** integrates a novel targeted local assembler, [**otter**](https://github.com/holstegelab/otter)


## How do you install TREAT
The easiest way to install **TREAT** in your system is to clone the repository in your system, and use the `INSTALL.sh` script. This will install a fresh version of Python along with the required packages, in a separate environment. This script assumes you have `conda` correctly installed in your system. If you are not familiar with Conda, please see [here](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html).
You can run the script by typing:  
`git clone https://github.com/holstegelab/treat.git`  
`cd treat/bin`  
`source INSTALL.sh`  
This script will install:
- `treat environment through conda`, which contains all required packages to execute treat (including Python 3.6 and R)
- `htslib-1.19.1`, which contains `htslib`, required to run `otter` 
- `otter`, which is a novel targeted local assembler. Additional information about `otter` is available at [here](https://github.com/holstegelab/otter).  
Note that `R` needs to be installed in your system. When running the `plot` or `analysis` modules of `TREAT`, it will automatically install the required packages, if not present (`data.table`, `stringr`, `argparse`, `ggplot2`, `dplyr`, `dendextend`, `berryFunctions`). 

You can independently install the required packages by `TREAT` without using the provided `INSTALL.sh`. These are listed in the `treat.yml`, and consists of `samtools`, `python 3.6`, `trf`, and the python packages `pysam`, `pandas`, `scikit-learn`, `numpy`, `biopython`. Please see `otter` [installation](https://github.com/holstegelab/otter) for more information about how to properly install `otter`.

At the end of the installation, you may need to re-run the following lines to enable `TREAT` and `otter` system-wide (otherwise, you need to run the tools using the full path):  
`export PATH=$PWD/:$PATH`  
`export PATH=$PWD/otter/build/:$PATH`  

## What do you need to run TREAT
To run **TREAT**, you need:
- **TREAT** correctly installed and working in your system
- the target regions of interest in the form of a BED file. Please see BED file specifications [here](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
- the target genomes in the form of aligned BAM file(s). Please see BAM file specifications [here](https://genome.ucsc.edu/goldenPath/help/bam.html)
- the reference genomes that was used to align genomes, in the form of a FASTA file. Please see FASTA file specifications [here](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/). GRCh38 is available [here](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/)

## What do you get as output
**TREAT** output consists of:
- `treat_run.log`: a recapitulation of the job along with a replicable command
- `sample.vcf.gz`: the main VCF file summarizing the genotype of the target regions in the target genomes

## How do you run TREAT
To run **TREAT**, you can follow these steps:
- `conda activate treat         # activate the right conda environment`  
- `TREAT.py -h         # run TREAT`  

## Toolkit
**TREAT** contains several tools to manipulate and analyze sequencing data. Three main analysis strategies are available in **TREAT**:
- `reads`: genotype the target region in the target genomes using single reads and a clustering framework. Genotypes are always given as the size of the region of interest
- `assembly`: genotype the target region in the target genomes with targeted local assembly of the single reads. Genotypes are always given as the size of the region of interest
- `merge`: can be used to combine multiple VCF. Genotypes are always given as the size of the region of interest  
- `analysis`: can be used to do downstream analysis of tandem repeats. Currently, two analyses are implemented: outlier analysis and case-control analysis
Typing `/path/to/treat/bin/TREAT.py [reads/assembly/merge] -h` will show the help message specific to the analysis of interest, along with all available run parameters.

## Assembly analysis
The `assembly` analysis take advantage of all sequencing reads aligning to the target region to perform local assembly of the target regions. Local assembly is done with [**otter**](https://github.com/holstegelab/otter). The procedure goes as it follows:
1. extract the reads and relative sequences encompassing the target regions
2. perform haplotype aware local assembly of the target regions in each sample
3. performs motif finding at the individual assembly level using [tandem repeat finder](https://tandem.bu.edu/trf/trf.html)
4. performs haplotype calling

### Required parameters
Same as for the `reads` analysis.

### Optional parameters
- `-w / --window`: the target regions defined in the BED file will be extended by this value upstream and downstream. Default value is 20. Must be an integer.
- `-t / --cpu`: number of parallel threads to be used. Default value is 2.
- `-d / --HaploDev`: during haplotype calling analysis, the median absolute deviation value to assign reads to the same allele. Defaul value is 0.10, that correspond to 10% median absolute deviation.
- `-minSup / --minimumSupport`: during haplotype calling, the minimum number of reads supporting each haplotyping. Default is 2.
- `-minCov / --minimumCoverage`: during haplotype calling, the minimum number of total reads necessary for calling. Default is 5.
- `-wAss / --windowAssembly`: the target regions defined in the BED file by this value upstream and downstream to take reads for assembly. Default value is 20. Must be an integer.
- `-p / --ploidy`: estimated ploidy of the sample. Default value is 2 for autosomal regions. For sex-specific regions, the ploidy is either 1 (for males with chrX and chrY present in the BAM file), or 2 (for females with 2 chrX).

## Reads analysis
The `reads` analysis take advantage of all sequencing reads aligning to the target regions to estimate genotypes. The procedure goes as it follows:
1. extract the reads and relative sequences encompassing the target regions
2. extract the corresponding sequence from the reference genome
3. performs motif finding at the individual read level using [tandem repeat finder](https://tandem.bu.edu/trf/trf.html)
4. performs phasing using SNPs (optional, skipped by default)
5. performs haplotype calling

### Required parameters
- `-b / --bed`: the target regions encoded in a BED file
- `-i / --inBam`: the targer genomes encoded in a BAM file. Multiple comma-separated BAM files can be used. If a folder is provided, all BAM files in the folder will be used.
- `-o / --outDir`: directory where output files will be placed. The output directory must NOT be present. TREAT will automatically create it.
- `-r / --ref`: the reference genome encoded in a FASTA file.

### Optional parameters
- `-w / --window`: the target regions defined in the BED file will be extended by this value upstream and downstream. Default value is 20. Must be an integer.
- `-t / --cpu`: number of parallel threads to be used. Default value is 2.
- `-minSup / --minimumSupport`: during haplotype calling, the minimum number of reads supporting each haplotyping. Default is 2.
- `-minCov / --minimumCoverage`: during haplotype calling, the minimum number of total reads necessary for calling. Default is 5.

## TREAT analysis module
`TREAT` includes a module for downstream analysis of tandem repeats. This takes as input the `VCF` file generated by `TREAT`, and performs either a outlier analysis or a case-control analysis:
- `outlier`: can be invoked with `TREAT.py analysis -a outlier -v input_vcf -r all`. This analysis will identify individuals with extreme deviation in the size of the tandem repeats. Outlier analysis will use the shorter allele, the longer allele, and the joint allele size. Outliers are first identified using Mahalanobis distance, and then a p-value is assigned to each outlier based on the Chi-Squared distribution. The output table reports each outlier along with the corresponding p-value.  
- `case-control`: can be invoked with `TREAT.py analysis -a case-control -v input_vcf -r all -l case_control_labels.txt`. This analysis will compare the allele size of tandem repeats between 2 groups using a logistic regression framework. Both the shorter, the longer, and the joint sum of alleles is compared between cases and controls. Because `TREAT` needs to know the case-control labels of the samples included in the `VCF` file, you must attach a tab-separated file with 2 columns and no header. The first column should include the sample name as found in the `VCF` file, while the second column should contain a binary phenotype. `TREAT` will check the type of phenotype provided, and will binarize it.

### Optional parameters
- `-o`: output directory (by default, the current working directory)  
- `-n`: output name (by default, treat_analysis_output.txt)  
- `-r`: region, which can be `all` (all regions will be analysed), or a specific region (based on the `ID` column of the `VCF` file)  
- `-t`: for outlier analysis, the threshold (number of standard deviations) for calling outliers based on median absolute deviation (by default, 3)  
- `-c`: number of cpus to use (by default, 1)  

## TREAT plot module
`TREAT` includes a module for plotting tandem repeats across samples. This module can be invoked with `TREAT.py plot -v input_vcf -r all`. The plotting module will produce 2 plots: 
- a plot showing tandem repeats allele sizes across individuals, along with hierarchical clustering to order the individuals based on their alleles, and the representative motif, for each allele of each individual  
- a plot showing tandem repeats allele size frequency, in which allele sizes are binned (bins with a size of 20bp)  

### Optional parameters
- `-o`: output directory (by default, the current working directory)  
- `-n`: output name (by default, treat_analysis_output.txt)  
- `-r`: region, which can be `all` (all regions will be analysed), or a specific region (based on the `ID` column of the `VCF` file)  
- `-p`: plot format (`png` or `pdf`)  

## Additional folders in the repository

### test_data folder
The `test_data` folder contains test data that can be use to assess the correct functioning of **TREAT**. A basic test can be the following for a `reads` analysis:
`TREAT.py reads -b test_data/example.bed -i test_data/example.bam -o test_output -r /path/to/reference_genome_hg38.fa`  
While for a `assembly` analysis, the following can be used:  
`TREAT.py assembly -b test_data/example.bed -i test_data/example.bam -o test_output_asm -r /path/to/reference_genome_hg38.fa -s otter`  

### treat_application folder
The `treat_application` folder contains several sub-folders related to the different projects from our group in which **TREAT** was used. Each sub-folder contains additional scripts and downstream analysis scripts that were used. Please look at the README in the specific sub-folders for additional information.