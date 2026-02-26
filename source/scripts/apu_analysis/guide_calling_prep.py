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

import pandas as pd
import numpy as np
import seaborn as sns
import anndata as ad

#read in cell identites
def clean_cellranger_guidecalling(fileloc, oligosloc):
    """Read in the sgRNA sgRNA_candidates.xlsx file for the oligos and the cell assignment for very
    """
    oligos=pd.read_excel(oligosloc)
    cell_identites=pd.read_csv(fileloc)
    cell_identites["guide_target"]=cell_identites["guide_identity"].str.split("_").str.get(0).str.split("sg").str.get(1)
    cell_identites["guide_target"]=np.where(cell_identites["guide_target"].str.startswith("non"), "NegCtrl3b", cell_identites["guide_target"]  )
    cell_identites["guide_target"][cell_identites["guide_target"]=="NegCtrl3b"]="Negative_Control"
    #form table of promoter type nicely formatted
    cell_identites["promoter_type"]=cell_identites["guide_identity"].str.split("_").str.get(1).str.split("sg").str.get(0).str.slice(0,2)
    cell_identites["promoter_type"]=np.where((cell_identites["promoter_type"]=="P1")|(cell_identites["promoter_type"]=="P2")|(cell_identites["promoter_type"]=="MP")|(cell_identites["promoter_type"]=="AP"),cell_identites["promoter_type"], "Control"  )
    cell_identites["promoter_type"]=np.where(cell_identites["promoter_type"]=="MP", "P1", cell_identites["promoter_type"]  )
    cell_identites["promoter_type"]=np.where(cell_identites["promoter_type"]=="AP", "P2", cell_identites["promoter_type"]  )
    cell_identites["guide_target"]=np.where((cell_identites["promoter_type"]=="Control")&(cell_identites["guide_target"]!="Negative_Control"),cell_identites["guide_target"].str.slice(0,-1),cell_identites["guide_target"])
    num_guides=cell_identites.reset_index()[["cell_barcode","guide_identity"]].drop_duplicates()
    #for each barcode how many differetnt guides are called 
    num_guides=pd.DataFrame(num_guides["cell_barcode"].value_counts())
    num_guides.columns=["num_guides"]
    #merge num_guides and cell 
    cell_identites=num_guides.merge(cell_identites, left_index=True, right_on="cell_barcode")
    #cell barcode as a column
    cell_identites["cell_bc"]=cell_identites["cell_barcode"]
    cell_identites["cell_barcode"]=cell_identites.index
    #check if same gene present 
    cell_identites["same_gene"]=cell_identites.groupby('cell_barcode').guide_target.nunique() == 1
    #check if A and B or C and D and MP/AP
    #form another column with A, B, C, D 
    cell_identites["same_promotertype"]=cell_identites.groupby('cell_barcode').promoter_type.nunique() == 1
    #get the position letter
    cell_identites["sgRNA_pos"]=cell_identites["guide_identity"].str.get(-1)
    L = ['A','B']
    s = frozenset(L)
    cell_identites["AorB"] = cell_identites.groupby('cell_barcode')['sgRNA_pos'].apply(s.issubset)
    L = ['C','D']
    c = frozenset(L)
    cell_identites["CorD"] = cell_identites.groupby('cell_barcode')['sgRNA_pos'].apply(c.issubset)
    #check the combo
    #provide correct partner colum with negative control section 
    #check if A an C or C and D depending if they are within the origignal guides
    cell_identites["AorD"]=np.where(cell_identites["guide_identity"].isin(oligos["target_A"]),'A',False)
    cell_identites["AorD"]=np.where(cell_identites["guide_identity"].isin(oligos["target_B"]),'D',cell_identites["AorD"])
    cell_identites['neg_correct_part'] = cell_identites['cell_barcode'].isin(cell_identites.loc[cell_identites['AorD']=="A", 'cell_barcode'])
    cell_identites['neg_correct_part'] = np.where(cell_identites['cell_barcode'].isin(cell_identites.loc[cell_identites['AorD']=="D", 'cell_barcode']), True, False)
    #####gets t
    cell_identites["correct_partner"]=np.where(((cell_identites["CorD"]==True) | (cell_identites["AorB"]==True) )   ,True,False)
    cell_identites["correct_partner"][cell_identites['guide_target']=="Negative_Control"]= np.where( (cell_identites["neg_correct_part"]==True)&(cell_identites['guide_target']=="Negative_Control"),True,False)
    cell_identites["ideal_guide"]=np.where(((cell_identites["correct_partner"]==True) & (cell_identites["same_promotertype"]==True) & (cell_identites["same_gene"]==True) & (cell_identites["num_guides"]==2)) | (cell_identites['guide_target']=="Negative_Control") ,True,False)
    cell_identites["correct_num_guides"]=np.where((cell_identites["num_guides"]==2) ,True,False)
    return cell_identites





