import pandas as pd 
import numpy as np
import scanpy as sc
from pandas.api.types import is_numeric_dtype
import random
import hdbscan
import gc
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score


def import_label(csv_interest):
    labels_1=pd.read_csv(csv_interest,skiprows=2)
    labels_1=labels_1[labels_1["Population"].isin(["G1","S-ph","G2M"])].sort_values(by="Mean").groupby("Well").head(n=1).dropna()
    labels_1["Sample"]=csv_interest.split("-")[-1].split("_")[0]
    labels_1["cell_name"]=labels_1["Well"]+"_"+labels_1["Sample"]
    labels_1=labels_1[["cell_name","Population","Sample"]]
    return labels_1

    



def strip_low_expression(pop, threshold=0):
    """Remove genes with low or zero expression to reduce memory usage. Modifies the
    target CellPopulation in place.
    
    Args:
        pop: CellPopulation instance
        threshold: all genes with expression <= threshold will be removed
    """
    retain = pop.genes.query('mean > @threshold').index
    remove = pop.genes.query('mean <= @threshold').index
    if len(remove) == 0:
        print('No genes have expression below threshold.')
        return
    pop.matrix = pop.matrix[retain]
    pop.genes.loc[remove, 'in_matrix'] = False
    # set all numeric properties to nan for genes that have been removed 
    for col in np.setdiff1d(pop.genes.columns, ['gene_name', 'in_matrix']):
        if is_numeric_dtype(pop.genes[col]):
            pop.genes.loc[remove, col] = np.nan
    # force garbage collection
    gc.collect()


def import_label(csv_interest):
    labels_1=pd.read_csv(csv_interest,skiprows=2)
    labels_1=labels_1[labels_1["Population"].isin(["G1","S-ph","G2M"])].sort_values(by="Mean").groupby("Well").head(n=1).dropna()
    labels_1["Sample"]=csv_interest.split("-")[-1].split("_")[0]
    labels_1["cell_name"]=labels_1["Well"]+"_"+labels_1["Sample"]
    labels_1=labels_1[["cell_name","Population","Sample"]]
    return labels_1



def add_filter_columns(input_dataframe):
    '''Function that adds T/F columns for same gene, same promoter type sgRNA pos for input into upset plot'''
    # Count the number of different guides for each cell barcode
    num_guides = input_dataframe.groupby('cell_barcode')['guide_identity'].nunique().rename('num_guides')
    # Check if all guides in a cell barcode target the same gene
    same_gene = input_dataframe.groupby('cell_barcode')['guide_target'].nunique() == 1
    # Check if all guides in a cell barcode have the same promoter type
    same_promotertype = input_dataframe.groupby('cell_barcode')['promoter_type'].nunique() == 1
    # Check if the guides in a cell barcode are A and B or C and D
    input_dataframe['sgRNA_pos'] = input_dataframe['guide_identity'].str[-1]
    AorB = input_dataframe.groupby('cell_barcode')['sgRNA_pos'].apply(lambda x: set(x) == {'A', 'B'})
    CorD = input_dataframe.groupby('cell_barcode')['sgRNA_pos'].apply(lambda x: set(x) == {'C', 'D'})
    # Check if the cell barcode satisfies all conditions
    ideal_guide = same_gene & same_promotertype & (num_guides == 2)
    correct_num_guides = num_guides == 2
    correct_partner = CorD | AorB
    # Create the final dataframe with the necessary columns
    result = pd.DataFrame({
        'num_guides': num_guides,
        'same_gene': same_gene,
        'same_promotertype': same_promotertype,
        'correct_partner': correct_partner,
        'ideal_guide': ideal_guide,
        'correct_num_guides': correct_num_guides
    })
    result=input_dataframe.merge(result,left_index=True, right_index=True)
    result.reset_index(inplace=True)
    result["DualGuide"]=np.where(((result["sgRNA_pos"]=="C")|(result["sgRNA_pos"]=="D")), "CD", "Nontargeting") #assigment of which dual guide of interest
    result["DualGuide"]=np.where(((result["sgRNA_pos"]=="A")|(result["sgRNA_pos"]=="B")), "AB", result["DualGuide"] ) #againt repeat 
    no_dup=result["guide_identity"][result["perturbation"]=="non-targeting"].drop_duplicates().reset_index(drop=True)
    no_dup_num = [*range(1, len(no_dup)+1)]
    res = {no_dup[i]: no_dup_num[i] for i in range(len(no_dup_num))}
    #adding in a column with guide assigned to it 
    result["guide_no"]=1
    result["guide_no"]=np.where((result["sgRNA_pos"]=="C") | (result["sgRNA_pos"]=="D"),2,result["guide_no"])
    #np.where for negative controls assign a new number for every guide
    result["guide_no"]=np.where(result["perturbation"]=="non-targeting",result["guide_no"].astype(int)+1,result["guide_no"])
    result["guide_no"]=result['guide_identity'].map(res).fillna(result["guide_no"]) #exhaustive mapping of the guide identity to the guide number
    result["guide_id"]=result["perturbation"]+"_"+result["promoter_type"]
    result["guide_assignment"]=result["perturbation"]+"_"+result["promoter_type"]+"_"+result["guide_no"].astype(int).astype(str)
    #'perturbation'  'promoter_type' 'guide_assignment'  'guide_assignment'
    result=result.sort_values(by="UMI_count",ascending=False).groupby("cell_barcode").head(1)
    result=result[result["perturbation"].isna()==False]
    return result
    
def extract_highly_variable_genes(counts_training_data):
    #extract the highly variable genes
    #read in mitochondrial genes from scanpy
    mito_gene_names = sc.queries.mitochondrial_genes("hsapiens",attrname="ensembl_gene_id")
    #read in training data 
    anndata_cc=sc.read_csv(counts_training_data)
    anndata_cc.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
    #filter the 
    sc.pp.filter_genes(anndata_cc, min_cells=3)
    #mitochondrial genes 
    anndata_cc.var['mt'] = anndata_cc.var_names.isin(mito_gene_names["ensembl_gene_id"]) # annotate the group of mitochondrial genes as 'mt'
    #return a list of the genes to filter for highly varaible only
    sc.pp.calculate_qc_metrics(anndata_cc, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    sc.pl.violin(anndata_cc, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                jitter=0.4, multi_panel=True)
    anndata_cc = anndata_cc[anndata_cc.obs.n_genes_by_counts > 1500, :]
    print(anndata_cc.obs.shape)
    anndata_cc = anndata_cc[anndata_cc.obs.pct_counts_mt > 0.5, :]
    print(anndata_cc.obs.shape)
    sc.pp.normalize_total(anndata_cc, target_sum=1e6)
    sc.pp.log1p(anndata_cc)
    sc.pp.highly_variable_genes(anndata_cc, min_mean=0.0125, max_mean=3, min_disp=0.5)
    highly_variable_genes=anndata_cc.var.index[anndata_cc.var['highly_variable']]
    return highly_variable_genes, anndata_cc



def best_validity(source):
    cols = ['metric','min_cluster_size', 'method','min_samples', 'validity_score', 'n_clusters','randscore','mutualscore']
    df =  pd.DataFrame(source, columns = cols)
    best_validity = df.loc[df['validity_score'].idxmax()]
    return best_validity

def best_mutscore(source):
    cols = ['metric','min_cluster_size', 'method','min_samples', 'validity_score', 'n_clusters','randscore','mutualscore']
    df =  pd.DataFrame(source, columns = cols)
    best_mutualscore = df.loc[df['mutualscore'].idxmax()]
    return best_mutualscore


def hyperparam(df,list_rep):
    import warnings
    warnings.filterwarnings('ignore')  
    n = df.shape[0]
    list_data_list=list()
    for metric in ['euclidean','manhattan']:
        for method in ['eom', 'leaf']:
            for gamma in range (1, int(np.log(n))):
                for ms in range(1, int(2 * np.log(n))):
                    cluster_size=int(gamma * np.sqrt(n))
                    clust_alg = hdbscan.HDBSCAN(algorithm='best',  alpha=1.0,
                                            approx_min_span_tree=True,
                                            gen_min_span_tree=True, 
                                            cluster_selection_method=method,
                                            metric=metric, 
                                            min_cluster_size=cluster_size, 
                                            min_samples=ms,
                                            allow_single_cluster=False).fit(df)
                    min_cluster_size = clust_alg.min_cluster_size 
                    min_samples = clust_alg.min_samples
                    randscore=adjusted_rand_score(np.where(list_rep=="MP", 0, 1), clust_alg.labels_)
                    mutualscore=adjusted_mutual_info_score(np.where(list_rep=="MP", 0, 1), clust_alg.labels_)
                    validity_score = clust_alg.relative_validity_
                    n_clusters = np.max(clust_alg.labels_) 
                    data_list=[metric,min_cluster_size,method,  min_samples,validity_score, n_clusters,randscore,mutualscore]
                    list_data_list.append(data_list)
                    if validity_score >= .5:
                        print (f'metric = {metric} ,min_cluster_size = {min_cluster_size},method = {method} , min_samples = {min_samples}, validity_score = {validity_score} n_clusters = {n_clusters}, randscore = {randscore},mutualscore = {mutualscore}')
    return list_data_list
