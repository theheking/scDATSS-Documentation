#!/bin/bash
#PBS -q normal
#PBS -l ncpus=4
#PBS -l mem=150GB
#PBS -l jobfs=30GB
#PBS -l walltime=48:00:00
#PBS -m be
#PBS -j oe




###---------------
#Brief
###----------------
#Take input csv containing cellbarcodes assigned to pseudobulk population 
#and seperate into individual bam files


###---------------
#LOAD MODULE
###----------------
module load samtools/1.9

 

###---------------
#DIRECTORIES
###----------------
LOC="alt-prom-crispr-fiveprime/"
IN_LOC_LIST="alt-prom-crispr-fiveprime/pseudobulk_bam/cellbarcode_list_ideal/"
OUT_LOC="alt-prom-crispr-fiveprime/pseudobulk_bam/cellbarcode_bam_ideal/"



###---------------
#FILE PATHS
###----------------
BAM_FILE="${LOC}/AP_CRISPR_WHOLE_ENSEMBL/outs/possorted_genome_merged.bam"
FILTERED_SAM_BODY="${SAMPLE}filtered_SAM_body"
FILTERED_SAM="${SAMPLE}filtered.sam"
SAM_HEADER="${OUT_LOC}possorted_genome_bam.header.sam"
SAMPLE_SAM_HEADER="${OUT_LOC}${SAMPLE}possorted_genome_bam.header.sam"
OUTPUT_FILE_BAM="${OUT_LOC}${SAMPLE}.sortedByCoord.out.bam"
CELL_FILE="${IN_LOC_LIST}${SAMPLE}.csv"

# Print file paths for debugging
echo ${CELL_FILE} ${OUTPUT_FILE_BAM} ${BAM_FILE}

# Create output directory if it doesn't exist

mkdir -p ${OUT_LOC}
cd ${OUT_LOC}

# Check if SAM header is present, if not, create it
if [ ! -f "${SAM_HEADER}" ]; then
    # Extract the header of the BAM file
    samtools view -H ${BAM_FILE} > ${SAM_HEADER}
fi

#create a new line in the SAM header file to add then infomation about pseudobulk to another file
NEW_LINE="@CO\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:10X"
scp ${SAM_HEADER} ${SAMPLE_SAM_HEADER}
#append the new line to the file add the end
sed -i -e '$a\'$'\n'"${NEW_LINE}" ${SAMPLE_SAM_HEADER}

#loop through every bam file
# Filter alignments using filter.txt. Use LC_ALL=C to set C locale instead of UTF-8
#search in every bam file for the cell barcodes in the cell file
samtools view -@48 ${BAM_FILE}  | LC_ALL=C grep -F -f ${CELL_FILE} > ${FILTERED_SAM_BODY}

# Combine header and body
cat ${SAM_HEADER} ${FILTERED_SAM_BODY} > ${FILTERED_SAM}

# Convert filtered.sam to BAM format
samtools view -@14 -b ${FILTERED_SAM} > ${OUTPUT_FILE_BAM}

# Remove intermediate files
rm ${FILTERED_SAM_BODY} ${FILTERED_SAM} ${SAMPLE_SAM_HEADER}

# Index the BAM file
samtools index -@14 ${OUTPUT_FILE_BAM}



