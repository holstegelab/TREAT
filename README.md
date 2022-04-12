# Tandem Repeat Haplotyping Toolkit (TRHT)
A set of tools to manipulate and analyze Pacbio and other sequencing data.

## Setup environment and installation
1. Clone the repository anywhere in your system
2. Type `conda activate py37` to load the conda environment
3. Either add the directory in your bash_profile (`export PATH="/path/to/TRHT/directory/:$PATH"`) or run directly from the folder
4. You may need to make the program executable (`chmod +x path/to/TRHT.py`)

## Dependencies
TRHT will call several external tools. All these tools should be readily available in the cluster. The tools are:
- whatshap (conda environment whatshap)
- tabix (conda environment whatshap)
- bcftools (conda environment py37)
- pbmm2 (conda environment py37)
- plink2 (/project/holstegelab/Software/bin/plink2)
- urchin (running from my environment at the moment)
- samtools (running from my environment at the moment)
- hifiasm (running from my environment at the moment)
- gfatools (running from my environment at the moment)
- samtools (running from my environment at the moment)

## Toolkit
TRHT contains several tools to manipulate and analyze sequencing data. You can choose the type of analysis you want to run with the --analysis-type parameter.
Typing `TRTH.py -h` will show the help message which contains information about the analysis types and the optional parameters. 
Below is presented the list of analyses type along with their optional parameters

1. `extract_reads`: given a .bed file (`--bed file.bed`) and one or multiple .bam file(s) (`--bam-dir path/to/bam`), this analysis will extract reads mapping to the region specified in the bed file, and will output reads in .bam and .fasta format in the output directory (`--out-dir path/to/output/directory`). Optional parameter allows to define the window around your .bed file (`--window`). An example run should look like `TRHT.py --analysis-type extract_reads --bed test_data/example.bed --bam-dir test_data --out-dir test_output`.

2. `measure`: given a .bed file and one or multiple .bam file(s), this analysis will measure the distance between start-position and end-position for the reads spanning the region of interest. An example run should look like `TRHT.py --analysis-type measure --bed test_data/example.bed --bam-dir test_data --out-dir test_output`. By default, temporary files (that is, temporary .bam and .fasta files) are deleted, but you can change this with `--store-temp True`. You can also add a polishing step (using Alex's script) with `--polish True`. The output is a plain-text file with single-read information.

3. `trf`: given a .bed file (containing chromosome, start, end positions, and motif information) and one or multiple .bam file(s), this analysis will extract reads spanning the region in the .bed file and will run Tandem Repeat Finder (TRF) on each read. By default, temporary files (that is, temporary .bam and .fasta files) are deleted, but you can change this with `--store-temp True`. You can also add a polishing step (using Alex's script) with `--polish True`. The output is a plain-text file with single-read information and TRF information. An example command should look like `TRHT.py --analysis-type trf --bed test_data/example.bed --bam-dir test_data --out-dir test_output --polish True`.

4. `assembly`: given a .bed file and one or multiple .bam file(s), this analysis will extract reads mapping to the region of interest and will perform local assembly with `hifiasm`. This procedure will output .bam files of the primary and haplotype-aware contigs. By default, temporary files (that is, temporary .bam and .fasta files) are deleted, but you can change this with `--store-temp True`. The use of a larger window is suggested with `--window 50000`. For the assembly, by default, each .bam will result in an assembly. If you prefer to combine multiple .bam files in the same assembly, please use `--assembly-type` and submit a file with no header and 2 columns: the first column should report, for each line, a comma-separated list of .bam files to  combine in the assembly. The second column, for each line, should report the output prefix of the final assembly for each group. You can define the ploidy with `--assembly-ploidy` (default value is 2). Multiprocessing is supported with option `--threads` followed by the number of cpus you want to use. An example command should look like `TRHT.py --analysis-type assembly --bed test_data/example.bed --bam-dir test_data --out-dir test_output`.

5. `phase_reads`: given a .bed file, one or multiple .bam file(s), and a snp-dataset (`--snp-data`), this analysis will extract reads mapping to the region of interest and will perform phasing (with `whatshap` tool). This is a two-step procedure: in the first step, SNPs will be phased based on the sequencing reads. In the second step, reads will be tagged by their haplotype. The output is plain-text file that reports each individual read ID along with the phase (1, 2 or NA). The SNP dataset (speficied with `--snp-data`) should be a plink2 file. Please submite the .pvar file (see example below). Often, sample identifiers change between sequencing and SNP array datasets: in these cases, you can specify a mapping file (`--snp-data-ids`): this should be a 2-column file with GWAS ID and ID in sequencing data (see example data). If not provided, will assume the IDs are the same. An example run should look like `TRHT.py --analysis-type phase_reads --bed test_data/example.bed --bam-dir test_data --out-dir test_output --snp-data test_data/example.pvar --window 25000 --snp-data-ids test_data/example.map`.

6. `coverage_profile`: given a .bed file and one or multiple .bam file(s), this analysis will partition the region of interest in bins and count the number of reads mapping to each bin. The output is a plain-text file reporting the interval analyzed along with the number of mapping reads observed. Optional parameters are the size of the window (`--window`) to add to the region of interest (it is suggested to use a large window for this analysis), and the size of the bins to partition the interval into (`--coverage-step`, default is 500 nt). An example run should look like `TRHT.py --analysis-type coverage_profile --bed test_data/example.bed --bam-dir test_data --out-dir test_output --window 1000000 --coverage-step 500`.

## General information
### Single and multiple .bam files
.bam files of interest are specified with `--bam-dir path/to/bam` flag. Every analysis type supports a single .bam files, multiple .bam files (that is, a comma-separated list of .bam files) or all .bam files in a directory. For the latter, in case a directory is the input, all .bam files in the directory will be used.

### .bed file
Any tab-separated file that includes at least chromosome (chr1, chr2, ..., chrN), start position and end position is valid. Header is not required (it is assumed that the fields are chromosome, start and end, respectively), but if header is present, it should start with #. If the analysis-type is `--analysis-type trf`, then you should add a fourth field with the expected motif. This can also be a comma-separated list of motifs. The `trf` analysis will compare the expected motif(s) with the observed and report whether the motif is the same, the reverse-complement, any permutation of these, or different. The motif is used mostly to match TRF findings by size. This means that if the motif in the .bed file is 5 nucleotide long, the analysis will report TRF findings with any motif of length 5. If no motifs are found, this is increased to 10 and all motifs are reported.

### 