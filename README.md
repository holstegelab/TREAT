# TREAT (Tandem REpeat Annotation Toolkit)

![TREAT](/Figure_TREAT_original.png)

## TREAT in a nutshell
**TREAT** is a command line tool written in **Python** and **R** (for plotting) that can be used to work with tandem repeats and structural variants from long-read sequencing data. **TREAT** was developed specifically for long-read sequencing data. However, it can potentially be used with any sequencing data, including PacBio, Oxford Nanopore, and Illumina.  
**TREAT** integrates a novel targeted local assembler, [**otter**](https://github.com/holstegelab/otter)


## How do you install TREAT
The easiest way to install **TREAT** in your system is to clone the repository in your system, and use the `INSTALL.sh` script. You can do so by typing:  
`git clone https://github.com/holstegelab/treat.git`  
`cd treat/bin`  
`sh INSTALL.sh`  
This script assumes you have `conda` correctly installed in your system. If you are not familiar with Conda, please see [here](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html).
This script will install:
- `treat environment through conda`, which contains all required packages to execute treat (including Python 3.6 and R)
- `htslib-1.19.1`, which contains `htslib`, required to run `otter` 
- `otter`, which is a novel targeted local assembler. Additional information about `otter` is available at [here](https://github.com/holstegelab/otter).

Alternatively, you can independently install the required packages by `TREAT`. These are listed in the `treat.yml`, and main consists of `samtools`, `python 3.6`, `R`, `trf`, and the python packages `pysam`, `pandas`, `scikit-learn`, `numpy`, `biopython`. Please see `otter` [installation](https://github.com/holstegelab/otter) for more information about how to properly install `otter`.

## What do you need to run TREAT
To run **TREAT**, you need:
- **TREAT** correctly installed and working in your system
- the target regions of interest in the form of a BED file. Please see BED file specifications [here](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
- the target genomes in the form of aligned BAM file(s). Please see BAM file specifications [here](https://genome.ucsc.edu/goldenPath/help/bam.html)
- the reference genomes that was used to align genomes, in the form of a FASTA file. Please see FASTA file specifications [here](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/)

## What do you get as output
**TREAT** output consists of:
- `treat_run.log`: a recapitulation of the job along with a replicable command
- `sample.vcf.gz`: the main VCF file summarizing the genotype of the target regions in the target genomes
- `sample.seq.txt.gz`: this file contains the same information as the VCF file, but in a tab-delimited format
- `sample.raw.txt.gz`: this file contains the raw read- or contig-specific information

## How do you run TREAT
To run **TREAT**, you can follow these steps:
- `conda activate treat         # activate the right conda environment`
- `TREAT.py -h         # run TREAT`

## Toolkit
**TREAT** contains several tools to manipulate and analyze sequencing data. Three main analysis strategies are available in **TREAT**:
- `reads`: genotype the target region in the target genomes using single reads and a clustering framework. Genotypes are always given as the size of the region of interest
- `assembly`: genotype the target region in the target genomes with targeted local assembly of the single reads. Genotypes are always given as the size of the region of interest
- `merge`: can be used to combine multiple VCF. Genotypes are always given as the size of the region of interest  
Typing `/path/to/treat/bin/TREAT.py [reads/assembly/merge] -h` will show the help message specific to the analysis of interest, along with all available run parameters.

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
- `-p / --phasingData`: if available, the path to the SNP dataset that will be used for phasing. The SNP dataset must be in PLINK2 format. If you are not familiar with PLINK2, please visit https://www.cog-genomics.org/plink/2.0.
- `-m / --mappingSNP`: path to a 2-column file with SNP IDs and sequencing IDs. If not provided, will assume the IDs between SNP dataset and BAM files are the same.
- `-d / --HaploDev`: during haplotype calling analysis, the median absolute deviation value to assign reads to the same allele. Defaul value is 0.10, that correspond to 10% median absolute deviation.
- `-minSup / --minimumSupport`: during haplotype calling, the minimum number of reads supporting each haplotyping. Default is 2.
- `-minCov / --minimumCoverage`: during haplotype calling, the minimum number of total reads necessary for calling. Default is 5.

## Assembly analysis
The `assembly` analysis take advantage of all sequencing reads aligning to the target region to perform local assembly of the target regions. The procedure goes as it follows:
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
- `-s / --software`: software to be used for assembly [otter/hifiasm]. Default value is otter. 

## Additional folders in the repository

### test_data folder
The `test_data` folder contains test data that can be use to assess the correct functioning of **TREAT**. A basic test can be the following for a `reads` analysis:
`TREAT.py reads -b test_data/example.bed -i test_data/example.bam -o test_output -r /path/to/reference_genome_hg38.fa`  
While for a `assembly` analysis, the following can be used:  
`TREAT.py assembly -b test_data/example.bed -i test_data/example.bam -o test_output_asm -r /path/to/reference_genome_hg38.fa -s otter`  

### treat_application folder
The `treat_application` folder contains several sub-folders related to the different projects from our group in which **TREAT** was used. Each sub-folder contains additional scripts and downstream analysis scripts that were used. Please look at the README in the specific sub-folders for additional information.