import pandas as pd 
import seaborn as sns 
import matplotlib.pyplot as plt
import glob
import numpy as np
import seaborn as sns
from matplotlib.lines import Line2D

import scanpy as sc
import sys
from scipy import stats
from statsmodels.stats.weightstats import ztest
import anndata as ad
import statsmodels.stats.weightstats as ws
from sklearn.model_selection import train_test_split
from sklearn.dummy import DummyRegressor
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error
from sklearn import neighbors, datasets, preprocessing
from sklearn.metrics import accuracy_score
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import roc_curve
from sklearn.metrics import auc

pd.options.display.float_format = '{:.4f}'.format

from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from sklearn.manifold import TSNE
from numpy import reshape
from numpy import array

def import_label(csv_interest):
    labels_1=pd.read_csv(csv_interest,skiprows=2)
    labels_1=labels_1[labels_1["Population"].isin(["G1","S-ph","G2M"])].sort_values(by="Mean").groupby("Well").head(n=1).dropna()
    labels_1["Sample"]=csv_interest.split("-")[-1].split("_")[0]
    labels_1["cell_name"]=labels_1["Well"]+"_"+labels_1["Sample"]
    labels_1=labels_1[["cell_name","Population","Sample"]]
    return labels_1


#ttest between non targeting of the same phase
def ztest_phase(gene_target, promoter,non_targeting, phase_assignment, phase, guide_num):
    """read in negative control and the phase assignment dataframe and return the ttest between the two groups
    """
    targeting = phase_assignment[(phase_assignment['guide_num']==guide_num)& (phase_assignment['phase']==phase) & (phase_assignment['target_gene']==gene_target) & (phase_assignment['promoter_type']==promoter)]
    non_targeting_phase=non_targeting[non_targeting['phase']==phase]
    ztest = ws.ztest(targeting['count'], non_targeting_phase['count'], value=0)
    return ztest

#ttest between non targeting of the same phase
def ttest_phase(gene_target, promoter,non_targeting, phase_assignment, phase, guide_num):
    """read in negative control and the phase assignment dataframe and return the ttest between the two groups
    """
    targeting = phase_assignment[(phase_assignment['guide_num']==guide_num)& (phase_assignment['phase']==phase) & (phase_assignment['target_gene']==gene_target) & (phase_assignment['promoter_type']==promoter)]
    non_targeting_phase=non_targeting[non_targeting['phase']==phase]
    ttest = stats.ttest_ind(targeting['count'], non_targeting_phase['count'], equal_var=False)
    return ttest

def percentage_check(gene_target, promoter, phase_assignment, phase, guide_num):
    targeting=phase_assignment[(phase_assignment['guide_num']==guide_num)& (phase_assignment['phase']==phase) & (phase_assignment['target_gene']==gene_target) & (phase_assignment['promoter_type']==promoter)]
    # print(targeting)
    targeting_all=phase_assignment[(phase_assignment['guide_num']==guide_num) & (phase_assignment['target_gene']==gene_target)& (phase_assignment['promoter_type']==promoter)] 
    # print(targeting_all)
    percentage = (targeting['count'].sum()/targeting_all['count'].sum())*100
    return(percentage)

def num_sim(n1, n2):
  """ calculates a similarity score between 2 numbers """
  return 1 - abs(n1 - n2) / (n1 + n2)

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


def plot_roc_curve(axis, final_df, phase, roc_auc_df, roc_auc_df_seurat, color1, color2):
    palette = sns.color_palette([color1, color2])
    #set whiute style and font size
    sns.set_style("whitegrid")
    sns.set(font_scale=1.5)
    sns.lineplot(data=final_df, hue="Model", x="FPR", y="TPR", palette=palette, ax=axis)
    axis.plot([0, 1], [0, 1],'r--')
    axis.set_xlim([0, 1])
    axis.set_ylim([0, 1])
    axis.set_ylabel('True Positive Rate')
    axis.set_xlabel('False Positive Rate')
    axis.set_title('ROC Curve of '+phase)
    myHandle = [Line2D([], [], color=color1, lw = 3), Line2D([], [], color=color2, lw = 3)]

    # red_patch = mpatches.Patch(color=color1, label=)
    # blue_patch = mpatches.Patch(color=, label=)
    axis.legend(handles=myHandle, labels=['KNN AUC = %0.2f' % roc_auc_df.loc[0,"AUC"], 'Scanpy AUC = %0.2f' % roc_auc_df_seurat.loc[0,"AUC"]])

    # new_labels = ['label 1', 'label 2'] 
    #[', ], 
    # axis.legend(handles=palette,labels=['First line', 'Second line'],loc='center left', bbox_to_anchor=(1, 0.5))

def compare_seurat(X_test, y_test, labels_training_data, knn, color1, color2):
    
    f, axis = plt.subplots(3, 1)
    plt.rcParams['figure.figsize'] = [8,22]
    sns.set(font_scale=2)
    sns.set_style("whitegrid")


    phases = ["G1","G2M","S"]
    pos_labels = [ "G1","G2M","S"]
    #only select
    y_scores = knn.predict_proba(X_test)

    for i in range(3):
        phase = phases[i]
        column_select = phase + "_score"
        pos_label = pos_labels[i]
        fpr, tpr, threshold = roc_curve(y_test, y_scores[:, i], pos_label=pos_label)
        roc_auc = auc(fpr, tpr)
        roc_auc_df = pd.DataFrame({'FPR':fpr, "TPR":tpr, "AUC":roc_auc, "Phase":phase, "Model":"KNN"})
        fpr, tpr, threshold = roc_curve(labels_training_data["Population"], labels_training_data[column_select], pos_label=pos_label)
        roc_auc = auc(fpr, tpr)
        roc_auc_df_seurat = pd.DataFrame({'FPR':fpr, "TPR":tpr, "AUC":roc_auc, "Phase":phase, "Model":"Scanpy"})
        final_df = pd.concat([roc_auc_df, roc_auc_df_seurat]).reset_index()
        plot_roc_curve(axis[i], final_df, phase, roc_auc_df, roc_auc_df_seurat, color1, color2)
    
    # cell_cycle_scores=pd.DataFrame(knn.predict_proba(X_test), columns=["G1.score","G2M.score","S.score","NA.score"])
    cell_cycle_scores=pd.DataFrame(knn.predict_proba(X_test), columns=["G1.score","G2M.score","S.score"])

    cell_cycle_scores["predicted_score"]=knn.predict(X_test)
    cell_cycle_scores.index=X_test.index
    plt.savefig("/Users/helenking/Desktop/Weatheritt_Lab_Y2/alt-prom-crispr-fiveprime/figures/ROC.pdf")
    return cell_cycle_scores
        
        
def extract_genes_training(human_genes_loc, adata,genes_regev,genelist):
    if genelist =="cyclebase":
        # Get biomart annotations
        adata_var_names = pd.DataFrame(adata.var)
        annot = sc.queries.biomart_annotations(
            "hsapiens",
            ["ensembl_gene_id", "hgnc_symbol", "ensembl_peptide_id"],
        )
        
        annot = adata_var_names.merge(annot, on="ensembl_gene_id",  how="left")
        human_genes = pd.read_csv(human_genes_loc, sep="\t")
        annot = annot.merge(human_genes, right_on="gene",left_on="ensembl_peptide_id")
        annot = annot.groupby("ensembl_gene_id").head(n=1)
        annot_sub_100 = annot.sort_values(by="rank", ascending=True).dropna()
    elif genelist=="regev":
        adata_var_names = pd.DataFrame(adata.var)
        annot_sub_100=adata_var_names.merge(genes_regev, on="ensembl_gene_id",  how="inner")
    elif genelist=="both":
        adata_var_names = pd.DataFrame(adata.var)
        annot = sc.queries.biomart_annotations(
            "hsapiens",
            ["ensembl_gene_id", "hgnc_symbol", "ensembl_peptide_id"],
        )
        annot = adata_var_names.merge(annot, on="ensembl_gene_id",  how="left")
        human_genes = pd.read_csv(human_genes_loc, sep="\t")
        #filter for G1 specific genes 
        human_genes=human_genes[human_genes["peaktime"]<45]
        annot = annot.merge(human_genes, right_on="gene",left_on="ensembl_peptide_id")
        annot_sub_100 = annot.sort_values(by="rank", ascending=True).dropna()
        annot_sub_100_add=adata_var_names.merge(genes_regev, on="ensembl_gene_id",  how="inner")
        annot_sub_100=pd.concat([annot_sub_100,annot_sub_100_add]).drop_duplicates()
    return annot_sub_100

def label_from(labels_training_data):
    ##read in all the csv for the labelled data
    list_df=[]
    for csv_interest in glob.glob(labels_training_data):
        list_df.append(import_label(csv_interest))
    label_final=pd.concat(list_df).reset_index(drop=True)
    return label_final

def train_model(counts_training_data, labels_training_data, human_genes_loc,adata,genes_regev,genelist=["regev","cyclebase"]):
    # Load the training data
    GSE146773_counts = sc.read_csv(counts_training_data, first_column_names=True)
    label_final = label_from(labels_training_data)
    label_final["Population"]=label_final["Population"].replace("S-ph","S") #the s-phase
    #add the cell label
    #remove all cells with na 
    GSE146773_counts.obs=GSE146773_counts.obs.merge(label_final, left_index=True, right_on="cell_name", how="left").reset_index(drop=True).set_index('cell_name')
    #drop all cell with na in "population" column in obs
    GSE146773_counts=GSE146773_counts[GSE146773_counts.obs["Population"].isin(["G1","G2M","S"])]
    annot_sub_100 = extract_genes_training(human_genes_loc,adata,genes_regev,genelist)
    # Subset the genes
    gene_subset = GSE146773_counts.var_names[GSE146773_counts.var_names.isin(annot_sub_100["ensembl_gene_id"])]
    GSE146773_counts=GSE146773_counts[:,gene_subset]
    sc.pp.log1p(GSE146773_counts)
    sc.pp.pca(GSE146773_counts, n_comps=30)
    sc.pp.neighbors(GSE146773_counts, n_pcs=30)
    sc.tl.umap(GSE146773_counts)
    cellcycle_palette=sns.color_palette(["#36a047", "#7671b3","#d76127"])
    sc.pl.umap(GSE146773_counts, color="Population", palette=cellcycle_palette, save="cellcycle_assignment.pdf")
    # GSE146773_counts=GSE146773_counts[~GSE146773_counts.obs["Population"].isna()]    # Prepare the features and target
    # X = GSE146773_counts.to_df()
    X = pd.DataFrame(GSE146773_counts.obsm['X_umap'],index=GSE146773_counts.obs.index)
    y = GSE146773_counts.obs["Population"].to_list()
    print( GSE146773_counts.obs["Population"].value_counts())
    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=42)
    # Standardize the features
    scaler = preprocessing.StandardScaler().fit(X_train)
    # X_train = pd.DataFrame(scaler.transform(X_train), columns=X_train.columns, index=X_train.index)
    # X_test = pd.DataFrame(scaler.transform(X_test), columns=X_test.columns, index=X_test.index)
    # y_train=[x if x in ["G1","G2M","S"] else "S" for x in y_train]
    # Initialize the KNN classifier
    knn = neighbors.KNeighborsClassifier(n_neighbors=10)

    # Define the hyperparameters for grid search
    params = { 
        'n_neighbors': np.arange(3,12),
        'metric': ['euclidean','cityblock'],
        "weights":["uniform","distance"]
    } 

    # Perform grid search
    grid = GridSearchCV(estimator=knn, param_grid=params)
    grid.fit(X_train, y_train)
    # Fit the model and make predictions
    knn.fit(X_train, y_train)
    # Print the results
    print(f"Best score: {grid.best_score_}")
    print(f"Best number of neighbors: {grid.best_estimator_.n_neighbors}")
    print(f"Best parameters: {grid.best_params_}")

    return knn, X_train, X_test, y_train, y_test, gene_subset, scaler


