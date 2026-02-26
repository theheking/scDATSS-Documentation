# Dual guides are designed by Flashfry 

Negative Control With regards to negative controls: in previous papers the % varies from ~2% to ~5% of total sgRNA with a minimum number of 5 negative controls. For the controls here, we are using the negative controls from Weissman hCRISPRi library which are scrambled. To make these guides, the frequency of each DNA base at each position along the sgRNA protospacer sequence was calculated. Random sgRNA protospacer sequences weighted by these base frequencies were then generated to mirror the composition of the targeting sgRNAs. These were then filtered for sgRNAs based on mismatch scores proximal to TSS

Both P1 and P2
Using the P1 and P2 regions that are found through html step (2). The guides chosen were decided through these following steps:
•	Removed enzyme sites, polyT (length >= 4bp) and GC content between 80% and 20%
•	Removed sites close to one another 
•	First chose guides that are found in hCRIPSRiv2.1 database
•	Then chose guides in CRISPRDO with greater than 80% efficacy score 
•	Then chose the top 4 guides by rank (the rank is an aggregated score of all the different on-target and off-target metrics). E.g. guide with rank 1 would be chosen first.
