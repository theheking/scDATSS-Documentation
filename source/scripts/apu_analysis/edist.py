
import scperturb
def test_edist_perprom(adata, gene, grouping_variable, negative_control):
    #subset adata for ESR1 and nontargeting
    print(gene)
    negative_control=gene+"_AP"
    grouping_variable="guide_id"
    adata_sub=adata[adata.obs["perturbation"]==gene,:]
    #compare P1 and P2
    estats = scperturb.edist(adata_sub,obs_key=grouping_variable, obsm_key='X_pca', dist='sqeuclidean')
    df_edist = scperturb.etest(adata_sub, obs_key=grouping_variable, obsm_key='X_pca', dist='sqeuclidean', control=negative_control, alpha=0.05, runs=200, n_jobs=-1)
    #change the estats to a dataframe
    del(estats)
    # df_edist.loc[:,"edist"]=estats.iloc[:,1]
    return df_edist