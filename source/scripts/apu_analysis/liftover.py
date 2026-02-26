from pyensembl import EnsemblRelease
# from pyliftover import LiftOver
import pandas as pd
import numpy as np 

#get the chromosome from the gene name using pyensembl
release = EnsemblRelease(75)
def gene_ids_of_gene_name_try(x):
    try:
        return release.gene_ids_of_gene_name(x)[0]
    except:
        return "NA"
    
def gene_by_id_try(x):
    try:
        return release.gene_by_id(x).contig
    except:
        return "NA"
    
def read_promoter_psi(file_path):
    promoter_psi = pd.read_csv(file_path, sep="\t", header=None)
    promoter_psi.columns = ["GENE_NAME", "START_COORD", "END_COORD", "FANTOM_TSS_ANNOTATION", "PROM_PSI"]
    return promoter_psi

def get_gene_ids_and_chromosomes(promoter_psi, release):
    promoter_psi["ensembl_id"] = promoter_psi["GENE_NAME"].apply(lambda x: release.gene_ids_of_gene_name(x)[0] if release.gene_ids_of_gene_name(x) else "NA")
    promoter_psi["chr"] = promoter_psi["ensembl_id"].apply(lambda x: release.gene_by_id(x).contig if release.gene_by_id(x) else "NA")
    return promoter_psi

def convert_coordinates(promoter_psi, converter):
    promoter_psi["hg38_start_full"] = promoter_psi.apply(lambda x: converter.convert_coordinate(x.chr, x.START_COORD), axis=1)
    promoter_psi["hg38_end_full"] = promoter_psi.apply(lambda x: converter.convert_coordinate(x.chr, x.END_COORD), axis=1)
    return promoter_psi

def process_coordinates(promoter_psi):
    promoter_psi["hg38_start"] = promoter_psi["hg38_start_full"].str[0].str[1]
    promoter_psi["hg38_end"] = promoter_psi["hg38_end_full"].str[0].str[1]
    promoter_psi["hg38_strand"] = promoter_psi["hg38_end_full"].str[0].str[2]
    return promoter_psi

def adjust_coordinates(promoter_psi):
    promoter_psi["hg38_start"] = np.where(promoter_psi["hg38_strand"]=="-", promoter_psi["hg38_end"], promoter_psi["hg38_start"])
    promoter_psi["hg38_end"] = np.where(promoter_psi["hg38_strand"]=="-", promoter_psi["hg38_start"], promoter_psi["hg38_end"])
    promoter_psi["hg38_start"] = promoter_psi["hg38_start"] - 500
    promoter_psi["hg38_end"] = promoter_psi["hg38_end"] + 500
    return promoter_psi

def create_bed_file(promoter_psi, output_path):
    promoter_psi_bed = promoter_psi[["chr", "hg38_start", "hg38_end", "FANTOM_TSS_ANNOTATION", "hg38_strand"]]
    promoter_psi_bed.to_csv(output_path, sep="\t", index=False, header=False)
