.. _analysis-scripts:

Analysis Scripts & Pipelines
============================

This section contains the core computational notebooks and scripts used for data processing, 
statistical analysis, and visualization.

.. _preprocessing-scripts:

Preprocessing & Guide Calling
-----------------------------

.. toctree::
   :maxdepth: 1

   scripts/3_analysis_scripts/2A_cellranger_guidecalling.ipynb
   scripts/3_analysis_scripts/2B_guide_calling_cellrangerfilter4_idealguides.ipynb
   scripts/3_analysis_scripts/12_negativecontrol.ipynb

**Additional Files:**

* **Shell Pipeline:** :download:`1_cellranger.sh <scripts/3_analysis_scripts/1_cellranger.sh>`

.. _transcript-scripts:

Transcript Quantification
-------------------------

.. toctree::
   :maxdepth: 1

   scripts/3_analysis_scripts/3C_whippet_pseudobulk

**Additional Files:**

* **Shell Scripts:** :download:`3A_seperate_bam_file_ideal.sh <scripts/3_analysis_scripts/3A_seperate_bam_file_ideal.sh>` | :download:`3B_whippet_pseudobulk.sh <scripts/3_analysis_scripts/3B_whippet_pseudobulk.sh>`

.. _estat-scripts:

Transcriptional Divergence (E-Stats & Differential Expression)
---------------------------------------------------------------

.. toctree::
   :maxdepth: 1

   scripts/3_analysis_scripts/4_estatistic_gene_guide.ipynb
   scripts/3_analysis_scripts/5_genekd_neighbouring_gene_expression.ipynb
   scripts/3_analysis_scripts/6_perturbseq_tsne_kldivergence_umap.ipynb
   scripts/3_analysis_scripts/8_perturbseq_differentialexpression.ipynb

.. _cellphase-pathway-scripts:

Cell Phase & Pathway Modeling
------------------------------

.. toctree::
   :maxdepth: 1

   scripts/3_analysis_scripts/7_cellphase_model.ipynb
   scripts/3_analysis_scripts/9A_cnv_score.ipynb
   scripts/3_analysis_scripts/9B_velocity_loom.ipynb
   scripts/3_analysis_scripts/11_spectra.ipynb

.. _clinical-de-scripts:

Clinical & Promoter Prevalence
-------------------------------

.. toctree::
   :maxdepth: 1

   scripts/3_analysis_scripts/10_survival_curve_isoform