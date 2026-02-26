
import pandas as pd
import numpy as np
from scipy.io import mmread
import os
import six
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster, leaves_list
from time import time
from tqdm import tqdm_notebook as progress
from collections import defaultdict
import warnings
# from .util import _strip_cat_cols, gini

from sklearn import preprocessing as pre
from pandas.api.types import is_numeric_dtype
from six.moves import zip as izip
from time import time
import gc
from tqdm import tqdm_notebook


# CellPopulation class
class CellPopulation:
    
    def __init__(self, matrix, cell_list, gene_list, source='arrays', normalized_matrix=None, calculate_statistics=True):
        """A class for holding single-cell data with separate tables for expression, normalized expression, and data about cells and genes
        
        Args: 
            matrix: expression matrix of UMIs per cells (rows = cell barcodes, columns = Ensembl gene ids)
            cell_list: table of properties of cells in the population (indexed by cell barcode)        
            gene_list: table of properties of genes in the population (indexed by Ensembl gene id)
            source: keeps track of how this population was derived. If it is from raw data, the source will 
                    be 'arrays'. If it is a subpopulation derived from another population, this will contain
                    the list of criteria that were used to create that subpopulation
            normalized_matrix: Add the supplied normalized expression matrix
        """
        self.matrix = matrix
        self.normalized_matrix = normalized_matrix
        
        if calculate_statistics:
            # fill out the list of gene properties
            print("Generating summary statistics...")
            gene_list['mean'] = matrix.mean()
            gene_list['std'] = matrix.std()
            gene_list['cv'] = gene_list['std']/gene_list['mean']
            gene_list['fano'] = gene_list['std']**2/gene_list['mean']
            gene_list['in_matrix'] = True
        self.genes = gene_list
        self.cells = cell_list

        if calculate_statistics:
            self.cells['gem_group'] = self.cells.index.map(lambda x: int(x.split('-')[-1]))
            self.guides = sorted(cell_list[cell_list['single_cell']]['guide_identity'].unique())
        
        self.source = source
        print("Done.")
  
    @classmethod
    def from_file(cls, directory, genome=None, filtered=True, raw_umi_threshold=2000):
        """Load a Perturb-seq data set from cellranger's matrix market exchange format
        
        Args:
            directory: location of base path of 10x experiment. The filtered gene-barcode matrices will automatically be found
            genome: name of genome reference you aligned to (e.g. GRCh38)
            
        Example:
            >>>pop = CellPopulation.from_file('~/sequencing/perturbseq_expt/', genome='GRCh38')
        """ 
        if genome == None:
            genome = 'GRCh38'
        
        # output of cellranger count and cellranger aggr has slightly different structure
        if os.path.isdir(os.path.join(directory, os.path.normpath("outs/filtered_gene_bc_matrices"), genome)):
            if filtered:
                matrix_directory = os.path.join(directory, os.path.normpath("outs/filtered_gene_bc_matrices"), genome)
            else:
                matrix_directory = os.path.join(directory, os.path.normpath("outs/raw_gene_bc_matrices"), genome)
        else:
            if filtered:
                matrix_directory = os.path.join(directory, os.path.normpath("filtered_feature_bc_matrix"), genome)
            else:
                matrix_directory = os.path.join(directory, os.path.normpath("outs/raw_gene_bc_matrices_mex"), genome)
        
        print('Loading digital expression data: {0}...'.format(os.path.join(matrix_directory, "matrix.mtx")))
        genes_path = os.path.join(matrix_directory, "features.tsv.gz")
        gene_list = pd.read_csv(genes_path, sep='\t', header=None, names=['gene_id', 'gene_name'])
        
        barcodes_path = os.path.join(matrix_directory, "barcodes.tsv.gz")
        cell_barcodes = pd.read_csv(barcodes_path, sep='\t', header=None, names=['cell_barcode'])
        # for the moment we keep matrix as dense... should be sparse
        matrix = mmread(os.path.join(matrix_directory, "matrix.mtx.gz"))
        
        print('Densifying matrix...')
        if filtered:
            matrix = pd.DataFrame(matrix.transpose().todense(),
                                  columns=gene_list['gene_id'],
                                  index=cell_barcodes['cell_barcode'],
                                  dtype='int32')
        else:
            # raw gene-barcode matrices generally too big to fit in memory so must do some thresholding
            m = pd.Series(np.asarray(np.sum(matrix, axis=0)).flatten())
            ind = m[m >= raw_umi_threshold].index.values
            print('Filtering cell barcodes with fewer than {0} UMIs...'.format(raw_umi_threshold))
            matrix = matrix.tocsc()[:, ind]
            matrix = pd.DataFrame(matrix.transpose().todense(),
                                  columns=gene_list['gene_id'],
                                  index=cell_barcodes['cell_barcode'].iloc[ind],
                                  dtype='int32')
            
        # table of gene properties
        gene_list = gene_list.set_index('gene_id')
        
        if filtered:
            identity_filename = 'cell_identities.csv'
        else:
            identity_filename = 'raw_cell_identities.csv'
        
        print("Loading guide identities:" + os.path.join(directory,  identity_filename) + '...')
        # adjust for historical differences in column names...
        guide_identities = pd.read_csv(os.path.join(directory, identity_filename)) \
            .rename(columns={'cell BC': 'cell_barcode',
                             'read count': 'guide_read_count',
                             'UMI count': 'guide_UMI_count',
                             'coverage': 'guide_coverage',
                             'cell_BC': 'cell_barcode',
                             'read_count': 'guide_read_count',
                             'UMI_count': 'guide_UMI_count'}).set_index('cell_barcode')
        cols = guide_identities.columns
        cols = cols.map(lambda x: x.replace(' ', '_') if isinstance(x, six.string_types) else x)
        guide_identities.columns = cols
                
        # table of cell properties
        if filtered:
            cell_list = pd.merge(cell_barcodes, guide_identities, left_on='cell_barcode', right_index=True, how='left').set_index('cell_barcode')
        else:
            cell_list = pd.merge(cell_barcodes.iloc[ind], guide_identities, left_on='cell_barcode', right_index=True, how='left').set_index('cell_barcode')
        guide_targets = cell_list['guide_identity'].map(lambda x: str(x).split('_')[0])
        guide_targets.name = 'guide_target'
        cell_list['guide_target'] = guide_targets
        
        cell_list['single_cell'] = (cell_list['number_of_cells'] == 1) & (cell_list['good_coverage']) & (~(cell_list['guide_identity'] == '*'))
        cell_list['UMI_count'] = matrix.sum(axis=1)             
        
        return cls(matrix, cell_list, gene_list, source=directory)
        
        def densify_matrix(self, fill_value=0):
            print('Densifying matrix...')
            self.matrix = self.matrix.to_dense()


    def where(self, cells=None, genes=None, normalized=False, gene_names=False, densify=True, return_query_str=False, dropna=False, **kwargs):
        """Return expression matrix or normalized expression matrix sliced based on properties of cells or genes

        Args:
            cells: list of gene names or query string for cell properties (i.e. executes pop.cells.query(cells))
            genes: likewise for genes
            normalized: if True, return values from the normalized expression matrix    
            gene_names: if True, return a table with gene names as columns instead of Ensembl ids (default: False)
            densify: If True, densify matrix while subsetting (which can be faster)
            return_query_str: If true, return the conditions that were used internally for subsetting
            dropna: If true, drop any genes from output that have NaN or infinite values (only relevant to 
                normalized matrices)
            **kwargs: dict of variables to inject into name space accessed by Pandas' query, necessary if you
                want to for example test membership against a predefined list
        
        Returns:
            A dataframe containing expression subsetted based on conditions and (optionally) the final queries that
            were used to create it
        
        Examples:
            >>>pop.where(cells='perturbation == "control"', genes='mean > 0.25')
            
            would return the expression matrix genes of mean > 0.25 within control cells
            
            >>>pop.where(cells=interesting_cells,
                      genes=interesting_genes)
        
            would return an expression matrix with cell barcodes taken from the list interesting_cells and genes 
            whose names are in the list interesting_genes.
        
            >>>pop.where(genes='mean > 0.5 or index in @interesting_genes',
                         interesting_genes=interesting_genes)
        
            This is an example leveraging Pandas' query syntax and will return genes that either have mean
            expression > 0.25 (a query executed on the genes table) or in the list interesting_genes. This
            list has to be provided as a keyword argument which is then injected into the namespace searched by 
            query.

        """
        # which matrix to use
        if (normalized):
            matrix = self.normalized_matrix
        else:
            matrix = self.matrix
            
        if densify and isinstance(matrix, pd.SparseDataFrame):
            if not normalized:
                print('(where) Densifying matrix...')
            else:
                print('(where) Densifying normalized matrix...')
            matrix = matrix.to_dense()
                
        # if no queries, just return the expression matrix
        if (genes == None) & (cells == None):
            if gene_names:
                out_matrix = matrix.copy()
                out_matrix.columns = self.gene_names(out_matrix.columns)
            else:
                return matrix
        
        # are we performing a query based on traits?
        complex_gene_indexing = isinstance(genes, str)
        # complex_gene_indexing = isinstance(genes, list)
        complex_cell_indexing = isinstance(cells, str)
        cell_query_str = ''
        gene_query_str = ''
        
        # construct cell list
        if complex_cell_indexing: # query on metadata
            cell_index = self.cells.query(cells, global_dict=kwargs, engine='numexpr').index
            cell_query_str = '| ' + cells + ' '
        else: # list of cell barcodes or no condition
            cell_index = cells
            if cells !=  None:
                cell_query_str = '| cells in list '
            else: # no condition
                cell_query_str = ''
        
        # construct gene list
        if complex_gene_indexing: # query on metadata
            gene_index = self.genes.query('(' + genes + ') and (in_matrix)', global_dict=kwargs, engine='numexpr').index
            gene_query_str = '| '+ genes + ' '
            # print(gene_index.head())
            # print(gene_query_str[0])
        else: # list of genes or no condition
            if (genes != None): # list
                test_gene = genes[0]
                # print(genes[0])
                if test_gene[0:4] != 'ENSG': # we already have a list of ensembl ids
                    gene_names = True # return gene names when passed a list of gene names
                    # genes = self.gene_ids(genes)
                    genes = genes

                    if isinstance(genes, six.string_types):
                        genes = [genes,]
                gene_index = self.genes.index[self.genes.index.isin(genes) & self.genes['in_matrix']]
                gene_query_str = '| genes in list ' 
            else: # no condition
                gene_query_str = ''
                gene_index = genes
        
        # construct output, separate cases because faster
        if genes == None:
            out_matrix = matrix.loc[cell_index]
        elif cells == None:
            out_matrix = matrix[gene_index]
        else: # querying on both
            out_matrix = matrix.loc[cell_index, gene_index]  
        
        # if supplied either a list of gene names (rather than Ensembl ids) or requested, return a
        # table that has human readable gene names
        if gene_names:
            out_matrix.columns = out_matrix.columns
            # print(gene_names)
            # out_matrix.columns = self.gene_names(out_matrix.columns)

        # drop any columns that have NaN or infinite values
        if dropna:
            out_matrix = out_matrix.replace({np.inf: np.nan, -np.inf: np.nan}).dropna(axis=1)
        
        if return_query_str:
            return out_matrix, (cell_query_str, gene_query_str)
        else:
            return out_matrix

    def subpopulation(self, cells=None, genes=None, normalized_matrix=None, **kwargs):
        """Return a new CellPopulation instance based on properties of cells or genes
        
        Note: This function internally uses pop.where and enables the same slicing.
        
        Args:
            cells: list of gene names or query string for cell properties (i.e. executes pop.cells.query(cells))
            genes: likewise for genes
            normalized_matrix: if set to 'inherit', subset the parent normalized matrix. If a dataframe, this
                will be used as the normalized matrix.
            **kwargs: dict of variables to inject into name space accessed by pop.where, see 
                documentation for pop.where
       
        Example: 
        
            pop = pop.subpopulation(cells='(single_cell) and \
                                   (guide_identity in @perturbation_list)',
                                    normalized_matrix='inherit',
                                    perturbation_list=perturbation_list)
            
            would return a new CellPopulation instance consisting of cell barcodes that are called as singletons (i.e. no
            doublets), whose guide_identity is in the list perturbation_list. The normalized matrix will be subsetted
            from the parent population.
        """
        new_matrix, (cell_query_str, gene_query_str) = self.where(cells=cells,
                          genes=genes,
                          return_query_str=True,
                          **kwargs)
        
        new_pop = self.__class__(new_matrix,
                                 self.cells.loc[new_matrix.index],
                                 self.genes.loc[new_matrix.columns],
                                 source=self.source + '  ||' + gene_query_str + cell_query_str + ' || ')
        
        if normalized_matrix == 'inherit':
            print('Inheriting from parent normalized matrix...')
            new_pop.normalized_matrix = self.normalized_matrix.loc[new_matrix.index, new_matrix.columns].copy()
        else:
            new_pop.normalized_matrix = normalized_matrix

        # if some genes have been removed for memory purposes, add the names back to the genes table
        if not self.genes['in_matrix'].all():
            missing = self.genes.query('~in_matrix').copy()
            new_pop.genes = pd.concat([missing, new_pop.genes])
            new_pop.genes = new_pop.genes.loc[self.genes.index]
            
        return new_pop
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
        
