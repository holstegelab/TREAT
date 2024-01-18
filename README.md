# TREAT (Tandem REpeat Annotation Toolkit)

## TREAT in a nutshell
TREAT is a command line tool written in Python and R (for plotting) that can be used to work with tandem repeats and structural variants from long-read sequencing data.

## What data can be used with TREAT
Although TREAT was developed for PacBio long-read sequencing data, it can be used with other long-read sequencing technologies (for example, Oxford Nanopore) and potentially short-read sequencing data.

## What do you need to run TREAT
To run TREAT, you need:
- TREAT correctly installed and working in your system
- the target regions of interest in the form of a BED file. Please see BED file specifications at https://genome.ucsc.edu/FAQ/FAQformat.html#format1
- the target genomes in the form of aligned BAM file(s). Please see BAM file specifications at https://genome.ucsc.edu/goldenPath/help/bam.html 
- the reference genomes that was used to align genomes, in the form of a FASTA file. Please see FASTA file specifications at https://www.ncbi.nlm.nih.gov/genbank/fastaformat/ 

## What do you get as output
TREAT output consists of:
- `treat_run.log`: a recapitulation of the job along with a replicable command
- `sample.vcf.gz`: the main VCF file summarizing the genotype of the target regions in the target genomes
- `sample.seq.txt.gz`: this file contains the same information as the VCF file, but in a tab-delimited format
- `sample.raw.txt.gz`: this file contains the raw read- or contig-specific information

## How do you install TREAT
The easiest way to install TREAT in your system is to clone the repository in your system. You can do so by typing:
`git clone https://github.com/holstegelab/treat.git`
All you need to run TREAT is inside the repository.
It is fundamental to install all the tools TREAT uses: this include, among others, aligners (minimap2), assembler (otter, hifiasm, flye). Although it is possible to install all software packages separately, we highly encourage to use the provided Conda environment. If you are not familiar with Conda, please see https://conda.io/projects/conda/en/latest/user-guide/getting-started.html. The easiest way to install all softwares required by TREAT is to replicate the provided Conda environment. You can do so by typing:
`conda env create -f treat_environment.yml`. To make TREAT executable, you may need to type `chmod +x path/to/TREAT.py`. In addition, you may want to add TREAT directory to your bash_profile by typing `export PATH="/path/to/treat/bin/:$PATH"`. Alternatively, you can run TREAT directly from the folder.

## How do you run TREAT
To run TREAT, you can follow these steps:
- activate conda environment by typing `conda activate treat`
- you can run TREAT with `/path/to/treat/bin/TREAT.py -h`

## Toolkit
TREAT contains several tools to manipulate and analyze sequencing data. Three analysis strategies are available in TREAT:
- `reads`: genotype the target region in the target genomes using single reads and a clustering framework. Genotypes are always given as the size of the region of interest.
- `assembly`: genotype the target region in the target genomes with targeted local assembly of the single reads. Genotypes are always given as the size of the region of interest.
- `merge`: can be used to combine multiple VCF. Genotypes are always given as the size of the region of interest.
Typing `/path/to/treat/bin/TREAT.py [reads/assembly/merge] -h` will show the help message specific to the analysis of interest, along with all available run parameters.

## Reads analysis
The `reads` analysis take advantage of all sequencing reads aligning to the target regions to estimate genotypes. The procedure goes as it follows:
1. extract the reads and relative sequences encompassing the region of interest
2. extract the corresponding sequence from the reference genome
3. performs motif finding using tandem repeat finder (https://tandem.bu.edu/trf/trf.html)
4. performs phasing using SNPs (optional, skipped by default)
5. performs haplotype calling

### Required parameters
- `-b / --bed`
### Optional parameters

## General information
### Single and multiple .bam files
.bam files of interest are specified with `--bam-dir path/to/bam` flag. Every analysis type supports a single .bam files, multiple .bam files (that is, a comma-separated list of .bam files) or all .bam files in a directory. For the latter, in case a directory is the input, all .bam files in the directory will be used.

### .bed file
Any tab-separated file that includes at least chromosome (chr1, chr2, ..., chrN), start position and end position is valid. Header is not required (it is assumed that the fields are chromosome, start and end, respectively), but if header is present, it should start with #. If the analysis-type is `--analysis-type trf`, then you should add a fourth field with the expected motif. This can also be a comma-separated list of motifs. The `trf` analysis will compare the expected motif(s) with the observed and report whether the motif is the same, the reverse-complement, any permutation of these, or different. The motif is used mostly to match TRF findings by size. This means that if the motif in the .bed file is 5 nucleotide long, the analysis will report TRF findings with any motif of length 5. If no motifs are found, this is increased to 10 and all motifs are reported.

### Mapping file for phasing
Often the sample identifiers for GWAS SNP array and sequencing are different. To account for this, you can input a mapping file: this should be a 2-column file with header (GWAS ID and PACBIO ID), which should contain, for each row, the mapping information between the two platform. At the moment, the pacbio identifier should report `_step1.bam` instead of `.bam`. This will change.

### Reference genome
In a number of options, the reference genome in fasta format should be provided. 