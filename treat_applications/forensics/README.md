# TREAT used in a forensic project

## Idea of the project
**TREAT** was used with the aim of genotyping short-tandem-repeats (STR) that are normally used in forensics. 

## What STR were used
A complete description of the STR that were used for this project are available in the original publication as well as on the [Netherlands Forensics Institute (NFI) database](https://www.fdstools.nl/strnaming/updates/nl-hym-pyg-v100.html).  
The publication linked with this project will be provided as soon as it gets published.

## What data was used
For this project, long-read sequencing data from Oxford Nanopore technology was used. This confirms that **TREAT** is able to deal with sequencing data from other technologies than PacBio.

## What analysis mode was used
Although we tried to use both the `reads` and the `assembly` analysis types, we used an adapted version of the `reads` analysis. The `assembly` analysis was excluded as in forensics even SNP or small INDELs need to be considered when genotyping. This, couple with some noise in the data, resulted in assemblies deviating too much from the original underlying sequences. Similarly, in the `reads` analysis, we used a different approach for identifying the correct alleles in the data. The reason for this is that *tandem repeat finder*, which is normally implemented in **TREAT** for tandem repeat motif characterization, typically does not take SNP and INDELs into account when generating a consensus representation of the motif.

## Adapted version of TREAT
We therefore implemented an adapted version of **TREAT**. The adapted version takes the output of a typical **TREAT** run, but uses only the raw individual read data information (the `sample.raw.txt.gz` output).  
Briefly, the adapted version of **TREAT** extracts the sequences spanning the regions of interest, and then identifies haplotypes based on the size of the regions using a clustering framework (typical **TREAT** output).  
Sequences resembling each haplotype (`sample.raw.txt.gz`) are then queried for the presence of forensics alleles based on the most recent release from the [Netherlands Forensics Institute database](https://www.fdstools.nl/strnaming/updates/nl-hym-pyg-v100.html).  
We repeated the analysis using different thresholds for the minimum allelic coverage (6, 8, 10, and 12, respectively).  
STR were excluded in case of low coverage or when no matching alleles were found in the database of forensic alleles. Candidate alleles underwent further curation by prioritizing longest alleles and interrogating nearby SNPs. Curated alleles were validated using standard forensics procedures based on deep short read sequencing. Forensics STR were regarded as matching only when all predicted alleles, for a given sample, matched the curated alleles. Finally, we evaluated the effect of GC content, STR size and coverage on allele match fraction, using a linear regression model.

## Reproducibility
To reproduce our analysis, please first run **TREAT**:
`TREAT.py reads -b bed_file.bed -i samples_genome.bam -o output_directory -w 20 -r path/to/reference_genome_hg38.fa -d 0.005`
And then you can run the supplementary script for the adapted version. As the script does multiple operations, we recommend to run it interactively line-by-line through a python console.
`python Adapted_version.py`

## Additional files
- `database_allele.txt`: this file represents the forensics allele catalog as derived from [NFI website](https://www.fdstools.nl/strnaming/updates/nl-hym-pyg-v100.html)
- `forensic_hg38.extended.bed`: this file represents the input BED file for **TREAT**. It contains the location of all STR considered for forensics purposes and included in the project.