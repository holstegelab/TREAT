## TREAT (Tandem REpeat Annotation Toolkit)

# TREAT in a nutshell
TREAT is a command line tool written in Python and R (for plotting) that can be used to work with tandem repeats and structural variants from long-read sequencing data.

# What data can be used with TREAT
Although TREAT was developed for PacBio long-read sequencing data, it can be used with other long-read sequencing technologies (for example, Oxford Nanopore) and potentially short-read sequencing data.

# What do you need to run TREAT
To run TREAT, you need:
- TREAT correctly installed and working in your system
- the target regions of interest in the form of a BED file. Please see BED file specifications at https://genome.ucsc.edu/FAQ/FAQformat.html#format1
- the target genomes in the form of aligned BAM file(s). Please see BAM file specifications at https://genome.ucsc.edu/goldenPath/help/bam.html 
- the reference genomes that was used to align genomes, in the form of a FASTA file. Please see FASTA file specifications at https://www.ncbi.nlm.nih.gov/genbank/fastaformat/ 

# How do you install TREAT
The easiest way to install TREAT in your system is to clone the repository in your system. You can do so by typing:
`git clone https://github.com/holstegelab/treat.git`
All you need to run TREAT is inside the repository.
It is fundamental to install all the tools TREAT uses: this include, among others, aligners (minimap2), assembler (otter, hifiasm, flye). Although it is possible to install all software packages separately, we highly encourage to use the provided Conda environment. If you are not familiar with Conda, please see https://conda.io/projects/conda/en/latest/user-guide/getting-started.html. The easiest way to install all softwares required by TREAT is to replicate the provided Conda environment. You can do so by typing:
`conda env create -f treat_environment.yml`. To make TREAT executable, you may need to type `chmod +x path/to/TREAT.py`. In addition, you may want to add TREAT directory to your bash_profile by typing `export PATH="/path/to/treat/bin/:$PATH"`. Alternatively, you can run TREAT directly from the folder.

# How do you run TREAT
To run TREAT, you can follow these steps:
- activate conda environment by typing `conda activate treat`
- you can run TREAT with `/path/to/treat/bin/TREAT.py -h`

