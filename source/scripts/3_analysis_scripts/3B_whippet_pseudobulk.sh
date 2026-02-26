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
###---------------
#Take input csv containing cellbarcodes assigned to pseudobulk population

###---------------
#LOAD MODULES
###---------------
module load java samtools

### Define Paths
WHIPPET_IND="bin/Whippet.jl/bin/whippet-index.jl"
WHIPPET_QUANT="bin/Whippet.jl/bin/whippet-quant.jl"
WHIPPET_DELT="bin/Whippet.jl/bin/whippet-delta.jl"

bam_loc="alt-prom-crispr-fiveprime/pseudobulk_bam/cellbarcode_bam_ideal/"
bam_loc_dedup="alt-prom-crispr-fiveprime/pseudobulk_bam/cellbarcode_bam_ideal_dedup/"
fastq_loc="alt-prom-crispr-fiveprime/pseudobulk_bam/cellbarcode_fastq_ideal_dedup/"
whippet_loc="alt-prom-crispr-fiveprime/pseudobulk_bam/cellbarcode_bam_ideal_whippet_dedup/"
WHIPPET_INDEX_LOC="/g/data/oo78/rw9400/genomes/Hsa38/Homo_sapiens.GRCh38.102_whippet.jls"

fastq_loc_subsample="alt-prom-crispr-fiveprime/pseudobulk_bam/cellbarcode_fastq_ideal_dedup_subsamp/"
whippet_loc_subsample="alt-prom-crispr-fiveprime/pseudobulk_bam/cellbarcode_bam_ideal_whippet_dedup_subsamp/"

mkdir -p ${fastq_loc_subsample} ${whippet_loc_subsample}

cellbarcode_tag="CR"
umi_tag="UB"
remove_duplicates="false"
picard_path="bin/picard/build/libs/picard.jar"

### Process BAM Files
echo "Processing BAM files..."
for i in ${bam_loc}*.bam; do
    SAMPLE=$(basename $i .bam)
    outfile_path="${bam_loc_dedup}${SAMPLE}.deduplicated.bam"
    tmp_bam="${bam_loc_dedup}${SAMPLE}.cbcorrected.tmp.bam"
    metrics_path="${bam_loc_dedup}${SAMPLE}.deduplicated.txt"
    
    if [ ! -f ${tmp_bam} ]; then
        samtools view -h ${i} | sed 's@\(CB:Z:[ACGT]\{9\}\)_\([ACGT]\{9\}\)_\([ACGT]\{9\}\)@\1-\2-\3@' | \
        samtools view -b -o ${tmp_bam} -
    fi
    
    if [ ! -f ${outfile_path} ]; then
        java -jar ${picard_path} MarkDuplicates \
            -I ${tmp_bam} --BARCODE_TAG ${cellbarcode_tag} --CREATE_INDEX true \
            --MOLECULAR_IDENTIFIER_TAG ${umi_tag} --REMOVE_DUPLICATES ${remove_duplicates} \
            -O ${outfile_path} -M ${metrics_path}
    fi

done

### Convert BAM to FASTQ & Run Whippet Quantification
echo "Converting BAM to FASTQ and running Whippet..."
for i in ${bam_loc_dedup}*.bam; do
    SAMPLE=$(basename $i .bam)
    whippet_loc_final="${whippet_loc}${SAMPLE}"
    gene_file="${whippet_loc_final}.gene.tpm.gz"
    FASTQ="${fastq_loc}${SAMPLE}.fastq.gz"
    
    if [ ! -f ${gene_file} ]; then
        java -jar ${picard_path} SamToFastq -I $i -F $FASTQ
        bin/julia-1.6.7/bin/julia ${WHIPPET_QUANT} --biascorrect -x ${WHIPPET_INDEX_LOC} -o ${whippet_loc_final} ${FASTQ}
    fi

done

### Split FASTQ into Three Files
echo "Splitting FASTQ files..."
for i in ${fastq_loc}*deduplicated*.fastq.gz; do
    SAMPLE=$(basename $i .fastq.gz)
    FASTQ1="${fastq_loc_subsample}${SAMPLE}_1.fastq.gz"
    FASTQ2="${fastq_loc_subsample}${SAMPLE}_2.fastq.gz"
    FASTQ3="${fastq_loc_subsample}${SAMPLE}_3.fastq.gz"
    
    if [ ! -f ${FASTQ1} ]; then
        fastqsplitter -i ${i} -o ${FASTQ1} -o ${FASTQ2} -o ${FASTQ3}
    fi

done

### Run Whippet on Subsampled FASTQ Files
echo "Running Whippet on subsampled FASTQ files..."
for i in ${fastq_loc_subsample}*deduplicated*.fastq.gz; do
    SAMPLE=$(basename $i .fastq.gz)
    whippet_loc_final="${whippet_loc_subsample}${SAMPLE}"
    gene_file="${whippet_loc_final}.gene.tpm.gz"
    
    if [ ! -f ${gene_file} ]; then
        bin/julia-1.6.7/bin/julia ${WHIPPET_QUANT} --biascorrect -x ${WHIPPET_INDEX_LOC} -o ${whippet_loc_final} ${i}
    fi

done

echo "Pipeline completed!"
