.. _guide-design:

Guide Design
====================================

Procedures for selecting ideal promoter-specific guides using FlashFry. We prioritize 
protospacers within the first 100 bp of the TSS to maximize CRISPRi efficiency. 

.. _flashfry-analysis:

FlashFry Analysis
-----------------

.. toctree::
   :maxdepth: 1

   scripts/2_guide_design/2_FlashFry_dualguide

Design Strategy
---------------

**Negative Controls**

With regards to negative controls: in previous papers the % varies from ~2% to ~5% of total sgRNA with a minimum number of 5 negative controls. For the controls here, we are using the negative controls from Weissman hCRISPRi library which are scrambled. To make these guides, the frequency of each DNA base at each position along the sgRNA protospacer sequence was calculated. Random sgRNA protospacer sequences weighted by these base frequencies were then generated to mirror the composition of the targeting sgRNAs. These were then filtered for sgRNAs based on mismatch scores proximal to TSS.

**Guide Selection Workflow**

Using the P1 and P2 regions that are found through the promoter identification pipeline, the guides chosen were decided through these following steps:

1. **Quality Filtering:**
   * Removed enzyme sites, polyT (length >= 4bp) and GC content between 80% and 20%
   * Removed sites close to one another 

2. **Database Priority:**
   * First chose guides that are found in hCRISPRiv2.1 database
   * Then chose guides in CRISPRDO with greater than 80% efficacy score 

3. **Final Selection:**
   * Then chose the top 4 guides by rank (the rank is an aggregated score of all the different on-target and off-target metrics)
   * E.g. guide with rank 1 would be chosen first

**Quick Downloads:**

* :download:`FlashFry Analysis Report (HTML) <scripts/2_guide_design/2_FlashFry_dualguide.html>`

.. note::
   Refer to the internal ``README.md`` in the guide design folder for specific library construction details.

