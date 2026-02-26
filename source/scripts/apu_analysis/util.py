# Perturbseq library for loading and manipulating single-cell experiments
# Copyright (C) 2019  Thomas Norman

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import stats
import statsmodels.stats.weightstats as sm 
from random import sample 
import random
random.seed(10)

###heatmap correlation
def get_cell_barcodes(adata_input, goi, promoter_type, sgRNA_pos): 
    return adata_input.obs.index[(adata_input.obs['promoter_type'] == promoter_type) & (adata_input.obs['perturbation'] == goi) & (adata_input.obs["protospacer_number"]==sgRNA_pos)]

def get_cell_barcodes_minimum(adata_input, goi, promoter_type, sgRNA_pos, length): 
    adata_output=adata_input.obs.index[(adata_input.obs['promoter_type'] == promoter_type) & (adata_input.obs['perturbation'] == goi) & (adata_input.obs["protospacer_number"]==sgRNA_pos)]
    #return only a certain length 
    return adata_output[:length]

def get_random_cell_barcodes(adata_input, goi, length): 
    return sample(list(adata_input.obs.index[(adata_input.obs['perturbation'] != goi)]), length)

def get_cells(adata_input, cell_barcodes, goi): 
    return pd.DataFrame(adata_input[cell_barcodes, goi].X, columns=[goi], index=cell_barcodes)

def get_cells_transcriptome(adata_input, cell_barcodes):
    adata_no_colnames= pd.DataFrame(adata_input[cell_barcodes].X,  index=cell_barcodes)
    adata_no_colnames=adata_no_colnames.rename(columns={x:y for x,y in zip(adata_no_colnames.columns,range(0,len(adata_no_colnames.columns)))})
    return adata_no_colnames

def get_cells_hvg(adata_input, cell_barcodes): 
    hvg = adata_input.var_names[adata_input.var['highly_variable']]
    return pd.DataFrame(adata_input[cell_barcodes,hvg].X, columns=hvg , index=cell_barcodes)

def get_ztest(A_cells, control_neighborhood): 
    return sm.ztest(A_cells.values, x2=control_neighborhood.values, value=0)

def get_rho_test(A_cells, control_neighborhood): 
    #input the z-test between A_cells and control_neighborhood
    # print(control_neighborhood.head())
    #calculate the max value per row for A_cells
    A_cells=A_cells.mean(axis=0)
    control_neighborhood=control_neighborhood.mean(axis=0)
    A_cells = A_cells.sort_values(ascending=False)
    control_neighborhood = control_neighborhood.sort_values(ascending=False)
    # print(control_neighborhood.head(n=1))
    # print(A_cells.index[0], control_neighborhood.index[0])
    return stats.spearmanr(A_cells.index[:500], control_neighborhood.index[:500])
    # return A_cells.corrwith(control_neighborhood,numeric_only=True, drop=True, axis=0, method='pearson')
    # return sm.ztest(A_cells.mean(axis=0).values, x2=control_neighborhood.mean(axis=0).values, value=0)

def process_gene(adata_input, goi, geneset=["goi","transcriptome","hvg"]):
    """"Subset the adata frame with gene of interest 
    """
    if goi in ["ATF5", "CSNK1E", "ARNT", "UPF3B", "SNRPD2", "STAG2", "RPL31", "GINS1"] or adata_input[adata_input.obs["perturbation"]==goi].shape==[0,0]:
        print(goi, "not found")
        return None
    print(goi, "found")
    control_cell_barcodes = adata_input.obs.index[adata_input.obs['perturbation'] == "non-targeting"]
    results = [goi]
    for phase in ["MP", "AP"]:
        sgRNA_pos=["1", "2"]
        for group in sgRNA_pos:
            #get the minimum len of each cellbarcode
            minimum_length= min(len(get_cell_barcodes(adata_input, goi, phase, "1")), len(get_cell_barcodes(adata_input, goi, phase, "2")))
            cell_barcodes = get_cell_barcodes_minimum(adata_input, goi, phase, group, minimum_length)
            cell_barcodes_random = get_random_cell_barcodes(adata_input, goi, minimum_length)
            control_cell_barcodes_subset= sample(list(control_cell_barcodes), minimum_length)
            if len(cell_barcodes) < 15:
                print("skipping", phase, group)
                break
            for cells in [cell_barcodes, cell_barcodes_random]:
                # control_cell_barcodes_subset= sample(list(control_cell_barcodes), 50)
                if geneset == "goi":
                    gene_input = get_cells(adata_input, cells, goi)
                    control_gene_input= get_cells(adata_input, sample(list(control_cell_barcodes), len(cells)), goi)
                    gene_input = get_ztest(gene_input, control_gene_input)
                    results.extend(gene_input)

                elif geneset == "transcriptome":
                    #bootstrap for
                    transcriptome = get_cells_transcriptome(adata_input, cells)
                    control_transcriptome = get_cells_transcriptome(adata_input, control_cell_barcodes_subset)
                    gene_input = get_rho_test(transcriptome, control_transcriptome)
                    results.extend([gene_input.statistic, gene_input.pvalue])

                elif geneset =="hvg":
                    hvg_cells = get_cells_hvg(adata_input, cells)
                    control_hvg = get_cells_hvg(adata_input, control_cell_barcodes_subset)
                    gene_input = get_rho_test(hvg_cells, control_hvg)       
                    results.extend([gene_input.statistic, gene_input.pvalue])
             
                #assign the phase group and test type
                # results.extend([phase, group[0], group[1],cells])
                #output the column name for the z-test
    
    return results
            
### neighbouring genes analysis         
def get_subset(gene_list, goi):
    return gene_list[gene_list["gene_name"]==goi]

def check_promoters(adata_input, goi):
    return "MP" in adata_input.obs["promoter_type"][adata_input.obs["perturbation"]==goi].unique() and "AP" in adata_input.obs["promoter_type"][adata_input.obs["perturbation"]==goi].unique()

def get_cell_barcodes_goi(adata_input, goi, promoter):
    return adata_input.obs.index[(adata_input.obs['perturbation']==goi)&(adata_input.obs['promoter_type']==promoter)]

def get_neighboring_genes(protein_gene_list, goi, n, subset):
    idx = protein_gene_list[protein_gene_list['gene_name'] == goi].index[0]
    neighboring_genes = protein_gene_list.iloc[np.arange(idx - n, idx + n + 1)]
    neighboring_genes["same_chr"]=np.where(neighboring_genes["chr"]==subset.iloc[0]["chr"], True, False)
    neighboring_genes["closest_distance"]=np.where(neighboring_genes["orientation"]=="+", neighboring_genes["start"]-subset.iloc[0]["start"] ,neighboring_genes["end"]-subset.iloc[0]["end"] )
    return neighboring_genes

def get_neighborhood(adata_input, cell_barcodes, neighboring_genes):
    return pd.DataFrame(adata_input[cell_barcodes,neighboring_genes["gene_name"]].X, columns=neighboring_genes["gene_name"], index=cell_barcodes)

def get_control_neighborhood(adata_input, neighboring_genes):
    control_cell_barcodes=adata_input.obs.index[adata_input.obs["perturbation"]=="non-targeting"]
    return pd.DataFrame(adata_input[control_cell_barcodes,neighboring_genes["gene_name"]].X, columns=neighboring_genes["gene_name"], index=control_cell_barcodes)

def perform_ztest(neighborhood, control_neighborhood, neighboring_genes,goi,promoter,n):
    z_test_dict = dict()
    output = []
    i = 0
    for index,row in neighboring_genes.iterrows():
        gene=row["gene_name"]
        distance=row["closest_distance"]
        samechr=row["same_chr"]
        z_test_dict[gene] = sm.ztest(neighborhood[gene].values, x2=control_neighborhood[gene].values, value=0)
        list_output= [goi, promoter,gene, n-i ,distance,samechr, z_test_dict[gene][0], z_test_dict[gene][1]]
        i+=1
        output.append(list_output)
    return output