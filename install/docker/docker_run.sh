################################################
################################################
# This script contains example scripts to run TREAT from the docker environment.
# You can comment/uncomment lines depending on the script of interest.
################################################
################################################

# Set up variables for running
# Assembly-based analysis
INPUT_BAM=/root/TREAT/data/test_data/HPC_example/
INPUT_BED=/root/TREAT/data/test_data/example.bed
REFERENCE=/root/TREAT/data/GRCh38_full_analysis_set_plus_decoy_hla.fa
OUTPUT_ASM=/root/TREAT/data/TREAT_assembly

# Read-based analysis
OUTPUT_READS=/root/TREAT/data/TREAT_reads

# Outlier analysis
VCF_ASM=/root/TREAT/data/TREAT_assembly/sample.vcf.gz
VCF_READS=/root/TREAT/data/TREAT_reads/sample.vcf.gz
OUTPUT_ANALYSIS_ASM=/root/TREAT/data/outlierAnalysis_Assembly
OUTPUT_ANALYSIS_READS=/root/TREAT/data/outlierAnalysis_Reads
KNOWN_BOUNDARIES=/root/TREAT/data/test_data/HPC_example/known_boundaries.txt

# Case-control analysis
OUTPUT_ANALYSIS_CC_ASM=/root/TREAT/data/CaseControlAnalysis_Assembly
OUTPUT_ANALYSIS_CC_READS=/root/TREAT/data/CaseControlAnalysis_Reads
CC_LABELS=/root/TREAT/data/test_data/HPC_example/case_control_labels.txt

# Plot
OUTPUT_PLOT_ASM=/root/TREAT/data/Plot_Assembly
OUTPUT_PLOT_READS=/root/TREAT/data/Plot_Reads

###############################################
###############################################
# Commands for running

# TREAT assembly based analysis command
docker run -it --rm \
        -v ${INPUT_BAM}:/run_input.bam \
        -v ${INPUT_BAM}.bai:/run_input.bam.bai \
        -v ${INPUT_BED}:/run_input.bed \
        -v ${REFERENCE}:/run_input_ref.fa \
        -v ${REFERENCE}.fai:/run_input_ref.fa.fai \
        -v ${OUTPUT_ASM}:/run_outfolder \
        treat \
        assembly \
        -i /run_input.bam \
        -b /run_input.bed \
        -r /run_input_ref.fa \
        -o /run_outfolder

# TREAT read based analysis command
docker run -it --rm \
        -v ${INPUT_BAM}:/run_input.bam \
        -v ${INPUT_BAM}.bai:/run_input.bam.bai \
        -v ${INPUT_BED}:/run_input.bed \
        -v ${REFERENCE}:/run_input_ref.fa \
        -v ${REFERENCE}.fai:/run_input_ref.fa.fai \
        -v ${OUTPUT_READS}:/run_outfolder \
        treat \
        reads \
        -i /run_input.bam \
        -b /run_input.bed \
        -r /run_input_ref.fa \
        -o /run_outfolder

# TREAT outlier analysis command based on Assembly output
docker run -it --rm \
        -v ${INPUT_BED}:/run_input.bed \
        -v ${OUTPUT_ANALYSIS_ASM}:/run_outfolder \
        -v ${VCF_ASM}:/run_input.vcf.gz \
        -v ${KNOWN_BOUNDARIES}:/run_known_boundaries.txt \
        treat \
        analysis \
        -a outlier \
        -v /run_input.vcf.gz \
        -o /run_outfolder \
        -k /run_known_boundaries.txt

# TREAT outlier analysis command based on Reads output
docker run -it --rm \
        -v ${INPUT_BED}:/run_input.bed \
        -v ${OUTPUT_ANALYSIS_READS}:/run_outfolder \
        -v ${VCF_READS}:/run_input.vcf.gz \
        -v ${KNOWN_BOUNDARIES}:/run_known_boundaries.txt \
        treat \
        analysis \
        -a outlier \
        -v /run_input.vcf.gz \
        -o /run_outfolder \
        -k /run_known_boundaries.txt

# TREAT case-control analysis command based on Assembly output
docker run -it --rm \
        -v ${INPUT_BED}:/run_input.bed \
        -v ${OUTPUT_ANALYSIS_CC_ASM}:/run_outfolder_cc \
        -v ${VCF_ASM}:/run_input.vcf.gz \
        -v ${CC_LABELS}:/run_input_labels.txt \
        treat \
        analysis \
        -a case-control \
        -v /run_input.vcf.gz \
        -l /run_input_labels.txt \
        -o /run_outfolder_cc

# TREAT case-control analysis command based on Reads output
docker run -it --rm \
        -v ${INPUT_BED}:/run_input.bed \
        -v ${OUTPUT_ANALYSIS_CC_READS}:/run_outfolder_cc \
        -v ${VCF_READS}:/run_input.vcf.gz \
        -v ${CC_LABELS}:/run_input_labels.txt \
        treat \
        analysis \
        -a case-control \
        -v /run_input.vcf.gz \
        -l /run_input_labels.txt \
        -o /run_outfolder_cc

# TREAT plot command based on Assembly output
docker run -it --rm \
        -v ${OUTPUT_PLOT_ASM}:/run_outfolder_plot \
        -v ${VCF_ASM}:/run_input.vcf.gz \
        treat \
        plot \
        -v /run_input.vcf.gz \
        -r all \
        -o /run_outfolder_plot

# TREAT plot command based on Reads output
docker run -it --rm \
        -v ${OUTPUT_PLOT_READS}:/run_outfolder_plot \
        -v ${VCF_READS}:/run_input.vcf.gz \
        treat \
        plot \
        -v /run_input.vcf.gz \
        -r all \
        -o /run_outfolder_plot