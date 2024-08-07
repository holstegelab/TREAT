# Example commands to run TREAT/Otter with Docker

# To build the Docker image, use the following command:
docker build --network host -t treat .

# When the image is built, use the following script to run TREAT/Otter docker image

# Set up variables and data to use
INPUT_BAM=/root/TREAT/data/test_data/HPC_example/
INPUT_BED=/root/TREAT/data/test_data/example.bed
REFERENCE=/root/TREAT/data/GRCh38_full_analysis_set_plus_decoy_hla.fa
OUTPUT=/root/TREAT/data/outfolder
VCF=/root/TREAT/data/outfolder/sample.vcf.gz
OUTPUT_ANALYSIS=/root/TREAT/data/outfolder/outlierAnalysis
OUTPUT_ANALYSIS_CC=/root/TREAT/data/outfolder/CaseControlAnalysis
CC_LABELS=/root/TREAT/data/test_data/HPC_example/case_control_labels.txt

# TREAT assembly based analysis command
docker run -it --rm \
       -v ${INPUT_BAM}:/run_input.bam \
       -v ${INPUT_BAM}.bai:/run_input.bam.bai \
       -v ${INPUT_BED}:/run_input.bed \
       -v ${REFERENCE}:/run_input_ref.fa \
       -v ${REFERENCE}.fai:/run_input_ref.fa.fai \
       -v ${OUTPUT}:/run_outfolder \
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
       -v ${OUTPUT}:/run_outfolder \
       treat \
       reads \
       -i /run_input.bam \
       -b /run_input.bed \
       -r /run_input_ref.fa \
       -o /run_outfolder

# TREAT outlier analysis command
docker run -it --rm \
       -v ${INPUT_BED}:/run_input.bed \
       -v ${OUTPUT_ANALYSIS}:/run_outfolder \
       -v ${VCF}:/run_input.vcf.gz \
       treat \
       analysis \
       -a outlier \
       -v /run_input.vcf.gz \
       -o /run_outfolder

# TREAT case-control analysis command
docker run -it --rm \
        -v ${INPUT_BED}:/run_input.bed \
        -v ${OUTPUT_ANALYSIS_CC}:/run_outfolder_cc \
        -v ${VCF}:/run_input.vcf.gz \
        -v ${CC_LABELS}:/run_input_labels.txt \
        treat \
        analysis \
        -a case-control \
        -v /run_input.vcf.gz \
        -l /run_input_labels.txt \
        -o /run_outfolder_cc