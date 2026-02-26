# Analysis for APU  library for loading and manipulating single-cell experiments
# Copyright (C) Helen king

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

from __future__ import absolute_import

__version__ = 0.1

# from reference_prep import transcriptomic_start_end, find_genomic_coord
from .perturbseq_matrix_qc import strip_low_expression, add_filter_columns, import_label,extract_highly_variable_genes 
from .guide_calling_prep import clean_cellranger_guidecalling
from .reference_prep import transcriptomic_start_end, find_genomic_coord
from .cell_import import CellPopulation
from .sctransform import SCTransform
from .perturbseq_matrix_qc import hyperparam, best_validity
from .util import process_gene, get_subset, check_promoters, get_cell_barcodes_goi, get_neighboring_genes, get_neighborhood, get_control_neighborhood, perform_ztest
from .liftover import gene_ids_of_gene_name_try, gene_by_id_try, read_promoter_psi, get_gene_ids_and_chromosomes, convert_coordinates, process_coordinates, adjust_coordinates, create_bed_file
from .edist import test_edist_perprom
from .cellcycle_assignment import num_sim,extract_highly_variable_genes, percentage_check,train_model,label_from,extract_genes_training,compare_seurat,extract_highly_variable_genes,percentage_check, ztest_phase, ttest_phase, train_model, extract_genes_training, compare_seurat 
# from .expression_normalization import z_normalize_expression, normalize_to_control, normalize_matrix_to_control, log_normalize_expression, equalize_UMI_counts, normalize_to_gemgroup_control, inherit_normalized_matrix, normalize_matrix_by_key