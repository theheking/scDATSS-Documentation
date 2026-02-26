

# scDATSS Documentation - Single-Cell Discovering Active Transcriptional Start Sites

## Abstract
Transcription start site (TSS) usage is a fundamental driver of transcript diversity and cellular identity. Alternative promoter selection reshapes 5′ untranslated regions, coding potential, and downstream regulatory behaviour, often without proportional changes in total gene-level abundance. While bulk methods such as CAGE-seq provide high-resolution maps of transcription initiation, they average across heterogeneous populations and obscure cell-type–specific regulatory programs.

Here, we present scDATSS, a computational pipeline designed to quantify TSS usage at single-cell resolution from 5′-biased scRNA-seq data. scDATSS models relative TSS usage within genes using a Dirichlet–Multinomial framework, allowing robust inference under sparse, overdispersed count data typical of single-cell experiments. By explicitly modelling multi-TSS competition and background expression variability, scDATSS isolates promoter-specific regulation from global gene-level effects.

![Visual Figure](./source/fig/VisualAbstract_vNAR_v2.png)


## Links
Please access the website: 

- https://scDATSS-Documentation.readthedocs.io/en/latest/index.html

Please access the raw data: 
- **Pending**

Please access the paper:
- **Pending**

## Pipeline
1. TSS x Cell Quantificaion

    **Required inputs:**
    
    - 5′ scRNA-seq aligned BAM files or equivalent fragment-level alignments
    
    - Annotated TSS reference (e.g., refTSS or GENCODE TSS annotations)
    
    - Gene annotation reference (GRCh38 recommended)
    
    **Optional metadata:**
    
    - CRISPR guide assignments (for Perturb-seq analyses), Cell-type annotations or clustering labels, Condition labels (disease vs control, treatment vs baseline)

2. Nextflow pipeline
Reads mapping within a configurable window (e.g., ±25 bp) around annotated TSS are assigned to transcription start sites. UMI counts are aggregated per TSS per cell to construct a cell × TSS count matrix.

Low-confidence TSS (e.g., expressed in <10% of relevant cells) may be filtered to stabilise downstream modelling.

3. Statistical Modelling (R) 





















4. 
