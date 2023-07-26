# treat (Tandem REpeat Annotation Toolkit)
A set of tools to manipulate and analyze Pacbio and other sequencing data.

## Setup environment and installation
1. Clone the repository anywhere in your system
2. If not present, install the right conda environment using the provided yml file: `conda env create -f treat_environment.yml`
3. Activate the environment typing `conda activate treat`
4. Either add the directory in your bash_profile (`export PATH="/path/to/treat/directory/:$PATH"`) or run directly from the folder
5. You may need to make the program executable (`chmod +x path/to/treat.py`)

## Dependencies
`treat` will call several external tools. If you create the specific conda environment (following point 2 above), then all tools should be correctly installed and `treat` should be ready to use. If you don't use the `treat` conda environment, you can still manually install all required packages. Please refer to `treat_environment.yml` for the complete list of softwares that need installation.

## Toolkit
`treat` contains several tools to manipulate and analyze sequencing data. Although `treat` was specifically designed for PacBio HiFi data, it can work with any type of sequencing data, both long-read and short-read. A summary of all options in `treat` can be invoked with `treat.py --help`. You can choose the type of analysis you want to run with the `--analysis-type` parameter.
Typing `treat.py --help` will show the help message which contains information about the analysis types and the optional parameters. 
Below is presented the list of analyses type along with their optional parameters

1. `extract_reads`: given a `.bed` file (`--bed file.bed`) and one or multiple `.bam` file(s) (`--bam-dir path/to/bam`), this analysis will extract reads mapping to the region specified in the bed file, and will output reads in `.bam` and `.fasta` format in the output directory (`--out-dir path/to/output/directory`). Optional parameter allows to define the window around your `.bed` file (e.g., `--window 10000` for a 10000 bp window). An example run should look like `treat.py --analysis-type extract_reads --bed test_data/example.bed --bam-dir test_data --out-dir test_output`.
Output is one `.bam` file and one `.fasta` file for each sample.

2. `measure`: given a `.bed` file and one or multiple `.bam` file(s), this analysis will measure the distance between start-position and end-position for the reads spanning the region of interest. An example run should look like `treat.py --analysis-type measure --bed test_data/example.bed --bam-dir test_data --out-dir test_output`. By default, temporary files (that is, temporary `.bam` and `.fasta` files) are deleted, but you can change this with `--store-temp True`.
Output is a text file (`measures_spanning_reads.txt`) with the following columns `REGION, SAMPLE_NAME, READ_NAME, PASSES, READ_QUALITY, MAPPING_CONSENSUS, SEQUENCE_WITH_WINDOW, LENGTH_SEQUENCE, PADDING_SIZE`, which codes for the region of interest, the name of the `.bam` file, the read identifier, the number of passes of the read, its quality, its mapping consensus, the actual sequence with 10bp padding window on each side, the length of the sequence and the padding size.

3. `trf`: given a `.bed` file (containing chromosome, start, end positions, and optionally motif information) and one or multiple `.bam` file(s), this analysis will extract reads spanning the region in the `.bed` file and will run Tandem Repeat Finder (TRF) on each read. By default, temporary files (that is, temporary `.bam` and `.fasta` files) are deleted, but you can change this with `--store-temp True`. An example command should look like `treat.py --analysis-type trf --bed test_data/example.bed --bam-dir test_data --out-dir test_output`.
Outputs are text files (`measures_spanning_reads.txt` and `measures_spanning_reads_and_trf.txt`). The TRF specific outputs, in addition to the same columns as in the `measure` analysis, are `START_TRF, END_TRF, LENGTH_MOTIF_TRF, COPIES_TRF, PC_MATCH_TRF, PC_INDEL_TRF, MOTIF_TRF, PADDING_BEFORE, SEQUENCE_TRF,PADDING_AFTER, MATCH_TYPE`, which codes for start position of TRF, end position of TRF, length of the repetitive region, number of copies, percentage of match, percentage of indels, motif found by TRF, padding sequence before repetitive region, repetitive region, padding sequence after repetitive region, and type of match between observed and expected motif (as present in the `.bed` file).

4. `assembly`: given a `.bed` file and one or multiple `.bam` file(s), this analysis will extract reads mapping to the region of interest and will perform local assembly with `hifiasm`. This procedure will output aligned `.bam` files of the primary and haplotype-aware contigs. By default, temporary files (that is, temporary `.bam` and `.fasta` files) are deleted, but you can change this with `--store-temp True`. The use of a larger window is suggested with `--window 50000`. For the assembly, by default, each `.bam` will result in an assembly. If you prefer to combine multiple `.bam` files in the same assembly, please use `--assembly-type` and submit a file with no header and 2 columns: the first column should report, for each line, a comma-separated list of `.bam` files to combine in the assembly. The second column, for each line, should report the output prefix of the final assembly for each group. You can define the ploidy with `--assembly-ploidy` (default value is 2). An example command should look like `treat.py --analysis-type assembly --bed test_data/example.bed --bam-dir test_data --out-dir test_output`.
Outputs: aligned `.bam` file of primary and haplotype-aware contigs, `.fasta` file of primary and haplotype-aware contigs, `.gfa` file of primary contig. Assembly contigs are by default aligned to GRCh38, however, this can be changed to CHM13 with `--reference-genome chm13`.
Strictly speaking, `hifiasm` was developed specifically for Pacbio hifi reads: this doesn't mean it will not produce an output when used with other read types, but the accuracy may be lower. We are working to implement additional methods to deal with different types of reads, such as Nanopore and/or Illumina.

5. `phase_reads`: given a `.bed` file, one or multiple `.bam` file(s), and a snp-dataset (`--snp-data`), this analysis will extract reads mapping to the region of interest and will perform phasing (with `whatshap` tool). This is a two-step procedure: in the first step, SNPs will be phased based on the sequencing reads. In the second step, reads will be tagged by their haplotype. The output is plain-text file that reports each individual read ID along with the phase (1, 2 or NA). The SNP dataset (speficied with `--snp-data`) should be a plink2 file (https://www.cog-genomics.org/plink/2.0/). Please submite the `.pvar` file (see example below). Often, sample identifiers change between sequencing and SNP array datasets: in these cases, you can specify a mapping file (`--snp-data-ids`): this should be a 2-column file with GWAS_ID and ID in sequencing data (see example data). If not provided, will assume the IDs are the same. An example run should look like `treat.py --analysis-type phase_reads --bed test_data/example.bed --bam-dir test_data --out-dir test_output --snp-data test_data/example.pvar --window 25000 --snp-data-ids test_data/example.map`.
Outputs: phased `.VCF` file of the sample(s) of interest and text file (`haplotags.txt`) which shows a summary of the haplotags with the following columns `SAMPLE	READ_ID	HAPLOTYPE`, which code for sample name, read identifier and haplotype (1, 2, NA).

6. `coverage_profile`: given a `.bed` file and one or multiple `.bam` file(s), this analysis will partition the region of interest in bins and count the number of reads mapping to each bin. The output is a plain-text file reporting the interval analyzed along with the number of mapping reads observed. Optional parameters are the size of the window (`--window`) to add to the region of interest (it is suggested to use a large window for this analysis), and the size of the bins to partition the interval into (`--coverage-step`, default is 500 bp). An example run should look like `treat.py --analysis-type coverage_profile --bed test_data/example.bed --bam-dir test_data --out-dir test_output --window 1000000 --coverage-step 500`.
Outputs: `.bed` file without header containing chromosome, start position of bin, end position of bin, sample name and number of reads mapping aligned to the bin.

7. `haplotyping`: this analysis attempts to generalize the output of tandem repeat finder by calling haplotypes for each sample in each specific region. As such, it will display the copy number of each tandem repeat in each sample and each haplotype. This analysis can be run using:
- only trf output of reads-spanning the intervals (`--trf fname_trf.txt`) 
- trf output + phasing output (`--trf fname_trf.txt --phase fname_phase.txt`)
- trf output + phasing + trf output from assembly (`--trf fname_trf.txt --phase fname_phase.txt --asm fname_trf_asm.txt`).
The output, independently from the analysis type, are two files reporting haplotype information at the single-read level as well as at the sample level.

8. `complete`: given a `.bed` file and one or multiple `.bam` file(s), this analysis will perform all previous analyses in a unique workflow. As such, this analysis will:
- extract reads mapping to the regions defined in the `.bed` file
- measure the intervals defined in the `.bed` file in the reads spanning
- tandem repeat finder on the reads spanning the intervals defined in the `.bed` file
- phasing and haplotagging
- assembly and alignment of assembled contigs
- measure of intervals defined in `.bed` file in the assembled contigs
- tandem repeat finder on the assembled contigs
- coverage profile
- haplotyping
This analysis type produces several (sub)folders each containing the relative outputs of a specific sub-analysis.

## General information
### Single and multiple .bam files
.bam files of interest are specified with `--bam-dir path/to/bam` flag. Every analysis type supports a single .bam files, multiple .bam files (that is, a comma-separated list of .bam files) or all .bam files in a directory. For the latter, in case a directory is the input, all .bam files in the directory will be used.

### .bed file
Any tab-separated file that includes at least chromosome (chr1, chr2, ..., chrN), start position and end position is valid. Header is not required (it is assumed that the fields are chromosome, start and end, respectively), but if header is present, it should start with #. If the analysis-type is `--analysis-type trf`, then you should add a fourth field with the expected motif. This can also be a comma-separated list of motifs. The `trf` analysis will compare the expected motif(s) with the observed and report whether the motif is the same, the reverse-complement, any permutation of these, or different. The motif is used mostly to match TRF findings by size. This means that if the motif in the .bed file is 5 nucleotide long, the analysis will report TRF findings with any motif of length 5. If no motifs are found, this is increased to 10 and all motifs are reported.

### Mapping file for phasing
Often the sample identifiers for GWAS SNP array and sequencing are different. To account for this, you can input a mapping file: this should be a 2-column file with header (GWAS ID and PACBIO ID), which should contain, for each row, the mapping information between the two platform. At the moment, the pacbio identifier should report `_step1.bam` instead of `.bam`. This will change.

### Reference genome
In a number of options, the reference genome in fasta format should be provided. 