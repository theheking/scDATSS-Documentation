#!/bin/bash

#PBS -q normal
#PBS -l ncpus=12
#PBS -l mem=120GB
#PBS -l jobfs=30GB
#PBS -l walltime=48:00:00



###---------------
#Brief
###----------------
#Cellranger running
###---------------
#JOB MODE
###----------------
STAR_LOC="refdata-gex-GRCh38-2024-A/star/"
DIR_LOCATION="alt-prom-crispr-fiveprime/"
LIBRARIES="${DIR_LOCATION}cellranger_files/library.csv"
FEATURE_REF="${DIR_LOCATION}cellranger_files/guides.ensembl.csv"
JOBMODE="local"
SAMPLE_NAME="AP_CRISPR_WHOLE_ENSEMBL"

###---------------
#REFERENCE FILES
###----------------
REF_LOC="reference/homo_sapiens_ensembl/"
REF_LOCATION="reference/homo_sapiens_ensembl/singlecell"
INPUT_GTF="${REF_LOC}Homo_sapiens.GRCh38.77.gtf"
INPUT_FILTERED_GTF="${REF_LOC}Homo_sapiens.filtered.GRCh38.77.gtf"
INPUT_FILTERED_GTF_FINAL="${REF_LOC}Homo_sapiens.filtered.final.GRCh38.77.gtf"
INPUT_FASTA="${REF_LOC}Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
REPEAT_MASKER="weissman_crispr_pbody/files/hg38_rmsk.gtf"
REPEAT_MASKER_NOCHR="weissman_crispr_pbody/files/hg38_rmsk.nochr.gtf"

###---------------
#LOAD MODULES
###---------------
module load samtools

###---------------
#MAKE FOLDER AND CD 
###---------------
mkdir -p ${DIR_LOCATION}
cd ${REF_LOC}
echo ${SAMPLE_NAME} ${REF_LOCATION} ${LIBRARIES} ${FEATURE_REF} 

###---------------
#MAKE REFERENCE 
###---------------
#create if statement that if exists then skip
if [ ! -d "${REF_LOCATION}" ]; then
    gunzip Homo_sapiens.GRCh38.77.gtf.gz Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
    cellranger mkgtf ${INPUT_GTF} ${INPUT_FILTERED_GTF} --attribute=gene_biotype:protein_coding --attribute=gene_biotype:lncRNA  || exit 1
    grep -E '^(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|MT|X|Y)' ${INPUT_FILTERED_GTF} > ${INPUT_FILTERED_GTF_FINAL}
    cellranger mkref --genome="singlecell" --fasta=${INPUT_FASTA} --genes=${INPUT_FILTERED_GTF_FINAL} || exit 1 
else
    echo "Reference folder already exists. Skipping reference creation."
fi

###---------------
#RUN CELLRANGER COUNT
###---------------
cd ${DIR_LOCATION}
cellranger count --id=${SAMPLE_NAME} \
                    --transcriptome=${REF_LOCATION} \
                    --libraries=${LIBRARIES} \
                    --feature-ref=${FEATURE_REF} \
                    --jobmode=${JOBMODE} \
                    --include-introns false \
                    --create-bam true \
                    --chemistry SC5P-R2



###---------------
#RUN SCVELO
###---------------

INPUT_VELO="alt-prom-crispr-fiveprime/${SAMPLE_NAME}"
velocyto run10x -m ${REPEAT_MASKER_NOCHR} ${INPUT_VELO} ${INPUT_FILTERED_GTF_FINAL}