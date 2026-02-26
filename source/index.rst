scDATSS Documentation - Single-Cell Discovering Active Transcriptional Start Sites
===================================================================================

.. attention::
   ***Under construction***

Abstract 
-------------------------

Transcription start site (TSS) usage is a fundamental driver of transcript diversity and cellular identity. Alternative promoter selection reshapes 5′ untranslated regions, coding potential, and downstream regulatory behaviour, often without proportional changes in total gene-level abundance. While bulk methods such as CAGE-seq provide high-resolution maps of transcription initiation, they average across heterogeneous populations and obscure cell-type–specific regulatory programs.

Here, we present scDATSS, a computational pipeline designed to quantify TSS usage at single-cell resolution from 5′-biased scRNA-seq data. scDATSS models relative TSS usage within genes using a Dirichlet–Multinomial framework, allowing robust inference under sparse, overdispersed count data typical of single-cell experiments. By explicitly modelling multi-TSS competition and background expression variability, scDATSS isolates promoter-specific regulation from global gene-level effects.

.. figure:: ../fig/figure-01.png
   :class: with-border
   :alt: Visual abstract summarising single-cel scDATSS pipeline
   :align: center
   :width: 90%


.. toctree::
   :maxdepth: 1
   :caption: Analysis Pipeline:

   promoter_identification
   guide_design
   analysis_scripts
