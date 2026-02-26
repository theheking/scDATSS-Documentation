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
from collections import OrderedDict
import seaborn as sns
import pyensembl

# Initialize Ensembl database
ensembl = pyensembl.EnsemblRelease(76)
protein_coding=pd.read_table("~/Desktop/Weatheritt_Lab_Y2/alt-prom-crispr-fiveprime/files/reference/protein_coding.txt",low_memory=False)
# Read in the information
file_location = "/Users/helenking/Desktop/Weatheritt_Lab_Y2/alt-prom-crispr-fiveprime/cellranger_input/candidate_AP_MP_additionalinfo.xlsx"
candidate_MP = pd.read_excel(file_location)



def transcriptomic_start_end(gene_id,genomic_coord1,genomic_coord2):
    """Return the transcripomic start and end of the genomic coordinates given
    """
    # Specify the gene ID and the genomic coordinates
    # gene_id = 'ENSG00000040341'  # Replace with the actual gene ID
    # genomic_coord1 = 73747442	  # Replace with the first genomic coordinate
    # genomic_coord2 = 73582830  # Replace with the second genomic coordinate

    # Get the gene object
    gene = ensembl.gene_by_id(gene_id)

    # Initialize variables
    transcriptomic_distance = 0
    cumulative_length = 0
    transcriptomic_coord1 = 0

    # Check the gene strand
    if gene.strand == "+":
        # Iterate over the exons of the gene
        for exon in gene.exons:
            exon_id=exon.exon_id
            if exon_id not in list(protein_coding["Exon stable ID"]):
                print("notproteincoding")
                continue
            else:
                # Check if the exon end coordinate is before the first genomic coordinate
                if exon.end < genomic_coord1:
                    cumulative_length += exon.end - exon.start + 1
                # Check if the exon start coordinate is after the second genomic coordinate
                elif exon.start > genomic_coord2:
                    break
                # Check if the exon spans the first genomic coordinate
                elif exon.start <= genomic_coord1 <= exon.end:
                    transcriptomic_coord1 = cumulative_length + genomic_coord1 - exon.start
                    cumulative_length += exon.end - exon.start + 1
                # Check if the exon spans the second genomic coordinate
                elif exon.start <= genomic_coord2 <= exon.end:
                    transcriptomic_coord2 = cumulative_length + genomic_coord2 - exon.start
                    transcriptomic_distance = transcriptomic_coord2 - transcriptomic_coord1
                    break
                # Check if the exon is between the two genomic coordinates
                elif genomic_coord1 < exon.start < genomic_coord2:
                    transcriptomic_distance += exon.end - exon.start + 1

    else:
        # Iterate over the exons of the gene in reverse order (for negative strand)
        for exon in reversed(gene.exons):
            exon_id=exon.exon_id
            if exon_id not in list(protein_coding["Exon stable ID"]):
                print("notproteincoding")
                continue
            else:

                # Check if the exon start coordinate is after the first genomic coordinate
                if exon.start > genomic_coord1:
                    cumulative_length += exon.end - exon.start + 1
                # Check if the exon end coordinate is before the second genomic coordinate
                elif exon.end < genomic_coord2:
                    break
                # Check if the exon spans the first genomic coordinate
                elif exon.start <= genomic_coord1 <= exon.end:
                    transcriptomic_coord1 = cumulative_length - (genomic_coord1 - exon.start)
                    cumulative_length += exon.end - exon.start + 1
                # Check if the exon spans the second genomic coordinate
                elif exon.start <= genomic_coord2 <= exon.end:
                    transcriptomic_coord2 = cumulative_length - (genomic_coord2 - exon.start)
                    transcriptomic_distance = transcriptomic_coord1 - transcriptomic_coord2
                    break
                # Check if the exon is between the two genomic coordinates
                elif genomic_coord2 < exon.start < genomic_coord1:
                    transcriptomic_distance += exon.end - exon.start + 1
    print(f"The transcriptomic distance between the two coordinates is: {abs(transcriptomic_distance)} bp")
    return transcriptomic_distance

def find_genomic_coord(gene_id,start_coord,distance=500):
    """returns the genomic coord
    """
    # Get the transcript object
    #needs to become more positive
    genomic_coord = None
    exon_length=0
    gene = ensembl.gene_by_id(gene_id)
    if gene.strand=="+":
        # Iterate over the exons of the gene
        for exon in gene.exons:
            exon_id=exon.exon_id
            if exon_id not in list(protein_coding["Exon stable ID"]):
                print("notproteincoding")
                continue
            else:
                #filter for exons after the start coord
                if exon.start >= start_coord:
                    length= exon.end - exon.start + 1 #calculate length of exon
                    exon_length = exon_length + exon.end - exon.start + 1  # calculate culmative length of the exon
                    print(exon.exon_id,exon.start,exon.end,exon_length)
                    if exon_length>=distance: #check if culmative length is greater than or equal to 500 then find the internal transcript
                        # If the transcriptomic coordinate falls within the exon
                        genomic_coord = exon.start + (exon_length-distance)
                        print(exon.start , genomic_coord ,exon_length)
                        break
    elif gene.strand=="-":
        # Iterate over the exons of the gene in reverse order (for negative strand)
        for exon in reversed(gene.exons):
            exon_id=exon.exon_id
            if exon_id not in list(protein_coding["Exon stable ID"]):
                print("notproteincoding")
                continue
            else:
                # Filter for exons before the start coordinate
                if start_coord >= exon.start:
                    length = exon.end - exon.start + 1  # Calculate the length of the exon
                    exon_length += length  # Calculate cumulative length of the exon
                    print(length, exon_length)
                    if exon_length >= distance:  # Check if cumulative length is greater than or equal to 500
                        # If the transcriptomic coordinate falls within the exon
                        genomic_coord = exon.end - (exon_length - distance)
                        print(exon.end, genomic_coord, exon_length)
                        break




    # Print the result
    print("Genomic coordinates ("+str(distance)+"bp downstream/upstream):", start_coord, genomic_coord)
    return genomic_coord


def find_final_end(gene_id,start_coord,  promoter_type):
    """find the end coordinate
    """
    #if a negative strand gene and AP, find the end_coordinate 
    gene = ensembl.gene_by_id(gene_id)
    if gene.strand=="+":
        #check if MP or AP promoter
        #if MP then find the final coordinate by extracting the AP coordinate from the candidate_MP file 
        if promoter_type=="MP":
            final_coord=candidate_MP.loc[candidate_MP["Ensembl_ID"]==gene_id,"TSSend_AP"].values[0]
        #if AP then find the final coordinate by find the end coordinate of the final exon of the gene
        elif promoter_type=="AP":
            final_coord=gene.exons[-1].end
    elif gene.strand=="-":
        # Apply same logic but in the reverse order for the negative strand
        if promoter_type=="MP":
            final_coord=candidate_MP.loc[candidate_MP["Ensembl_ID"]==gene_id,"TSSend_AP"].values[0]
        elif promoter_type=="AP":
            final_coord=gene.exons[0].start
    # Print the result
    print("Genomic coordinates", start_coord,promoter_type , final_coord)
    return final_coord


# df=pd.DataFrame({"gene_name":gene.gene_name, "transcript_support_level":transcript.support_level, "transcript_biotype":transcript.biotype ,"transcript_name":transcript.transcript_name ,"exon_number":range(1,len(transcript.exons)+1) ,"transcript_id":transcript.transcript_id,"gene_id":[i.gene_id for i in transcript.exons], "chr":[i.contig for i in transcript.exons] ,"start":[i.start for i in transcript.exons],"end":[i.end for i in transcript.exons],"strand":[i.strand for i in transcript.exons]})
# transcript_support_level transcript_biotype transcript_name   
def format_to_gtf(dataframe):
    dataframe["dot"] = "."
    dataframe["dot_2"] = "."
    dataframe["exon"] = "exon"
    dataframe["source"] = "havana"
    dataframe["describe"] = 'gene_id "' + dataframe["gene_id"] +'"; transcript_id "' + dataframe["transcript_id"]+  '"; gene_name "' + dataframe["gene_name"] + '"; transcript_name "' + dataframe["transcript_name"] +'"; exon_number '+dataframe["exon_number"].astype(str)+';'
    gtf_columns=['chr','source','exon','start', 'end', "dot", "strand", "dot_2", "describe"]
    dataframe=dataframe[gtf_columns].dropna()
    return dataframe    


def create_df(gene,transcript,gene_id,promoter_type="_mp"):
    """Create a dataframe of the exons and the name of the transcript
    """
    df=pd.DataFrame({"gene_name":gene.gene_name,"transcript_name":transcript.transcript_name ,"exon_number":range(1,len(transcript.exons)+1) ,"transcript_id":transcript.transcript_id,"gene_id":[i.gene_id for i in transcript.exons], "chr":[i.contig for i in transcript.exons] ,"start":[i.start for i in transcript.exons],"end":[i.end for i in transcript.exons],"strand":[i.strand for i in transcript.exons]})
    df["gene_id"]=df["gene_id"]+promoter_type
    df["transcript_id"]=df["transcript_id"]+promoter_type
    df=format_to_gtf(df)
    #add final line of the dataframe of transcript
    df.loc[len(df.index)] = [transcript.contig, "havana", "transcript", transcript.start, transcript.end, ".", transcript.strand, ".", 'gene_id "' + gene_id + promoter_type+'"; transcript_id "' + transcript.transcript_id +promoter_type+'"; gene_name "' + gene.gene_name + '"; transcript_name '+transcript.transcript_name+ '"; transcript_support_level "' + str(transcript.support_level) +'"; transcript_biotype "' + str(transcript.biotype)  +'";'] 
    return df


def related_exons(start, end,  gene_id, promoter_type,strand,bp=500):
    """"Use gtf coordinated to extract all the transcripts that are related to the promoter
    """
    distance_from_prom=10000
    if promoter_type=="AP":
        # Get the gene object
        gene = ensembl.gene_by_id(gene_id)
        # Loop through the transcripts and exons
        output_df_list=[]
        if strand == "+":
            beginning=start
            for transcript in gene.transcripts:
                distance=(transcript.start - (beginning-bp))
                distance_bool=(transcript.start > (beginning-bp))
                #remove very small transcripts
                if (transcript.length > np.quantile([i.length for i in gene.transcripts],0.1)):
                    if distance < distance_from_prom:
                        if distance_bool:
                            df=create_df(promoter_type="_ap",gene=gene,gene_id=gene_id,transcript=transcript)
                            output_df_list.append(df)
                        else:
                            df=create_df(promoter_type="_mp",gene=gene,gene_id=gene_id,transcript=transcript)
                            output_df_list.append(df)
        if strand == "-":
            beginning=end
            num=0
            for transcript in gene.transcripts:
                distance_bool=(transcript.end < (beginning+bp ))
                distance=((beginning+bp)-transcript.end )
                # if transcript.transcript_id=="ENST00000398671":
                #     print(distance)
                #     print(distance_from_prom)
                #     print(transcript)
                #     print(transcript.end)
                #     print((beginning+bp ))
                #     print(transcript.length)
                #     print(np.quantile([i.length for i in gene.transcripts],0.25))
                #     break
                #remove any transcripts that are below the median length of the transcript length
                if (transcript.length > np.quantile([i.length for i in gene.transcripts],0.1)): 
                    #remove extreme distance
                    if distance < distance_from_prom:
                        if distance_bool:
                                df=create_df(promoter_type="_ap",gene=gene,gene_id=gene_id,transcript=transcript)
                                output_df_list.append(df)
                        else:
                            df=create_df(promoter_type="_mp",gene=gene,gene_id=gene_id,transcript=transcript)
                            output_df_list.append(df)
        #concatenate the two
        # print( len(ap), len(mp))
        output_df=pd.concat(output_df_list)
        output_df.reset_index(drop=True, inplace=True)
        return output_df                   


def num_related_transcript(gtf, start_col, end_col,  gene_id, promoter_type,strand,bp=500):
    """"Use gtf coordinated to extract all the transcripts that are related to the promoter
    """
    gene = ensembl.gene_by_id(gene_id)
    #input the melted dataframe 
    subset=gtf[ (gtf["Ensembl_ID"]==gene_id) & (gtf["strand"]==strand)]
    if 1>=subset.shape[0]:
        print("skip")
    else:
        #make dictionary of the coordinate of the MP promoter and AP promoter
        prom_dict={"start": subset[start_col][subset["Prom"]=="AP"].values[0],"end":subset[end_col][subset["Prom"]=="AP"].values[0]}
        if promoter_type=="AP":
            mp=[]
            ap=[]
            distance_from_prom=10000
            # Get the gene object
            # Loop through the transcripts and exons
            if strand == "+":
                beginning=prom_dict["start"]
                for transcript in gene.transcripts:
                    distance=(transcript.start - (beginning-bp))
                    distance_bool=(transcript.start > (beginning-bp))
                    #remove very small transcripts
                    if (transcript.length > np.quantile([i.length for i in gene.transcripts],0.1)):
                        if distance < distance_from_prom:
                            if distance_bool:
                                ap.append(transcript)
                            else:
                                mp.append(transcript)
            if strand == "-":
                beginning=prom_dict["end"]
                num=0
                for transcript in gene.transcripts:
                    distance_bool=(transcript.end < (beginning+bp ))
                    distance=((beginning+bp)-transcript.end )
                    #remove any transcripts that are below the median length of the transcript length
                    if (transcript.length > np.quantile([i.length for i in gene.transcripts],0.1)): 
                        #remove extreme distance
                        if distance < distance_from_prom:
                            if distance_bool:
                                ap.append(transcript)
                            else:
                                mp.append(transcript)
            #concatenate the two
            # print( len(ap), len(mp))
            final_dim=str(len(mp))+"_"+str(len(ap))+"_"+str(len(set(mp+ap))-len(ap)-len(mp))
            return final_dim 



def num_related_transcripts_old(gtf, gene_id, strand,start_col="TSS_start",  proportion_up=0.1,proportion_down=0.2):
    """"Use gtf coordinated to extract all the transcripts that are related to the promoter and export the number
    """
    gene = ensembl.gene_by_id(gene_id)

    upstream_bp=gene.length*proportion_up
    downstream_bp=gene.length*proportion_down
    bp=500

    #input the melted dataframe 
    subset=gtf[ (gtf["Ensembl_ID"]==gene_id) & (gtf["strand"]==strand)]
    #make dictionary of the coordinate of the MP promoter and AP promoter
    prom_dict={"MP": subset[start_col][subset["Prom"]=="MP"].values[0],"AP":subset[start_col][subset["Prom"]=="AP"].values[0]}
    # print(prom_dict)
    # Get the gene object
    # Loop through the transcripts and exons
    mp=[]
    ap=[]
    for transcript in gene.transcripts:
        diff_mp=abs(transcript.start-prom_dict["MP"])
        diff_ap=abs(transcript.start-prom_dict["AP"])
        if strand == "+":
            if (transcript.start > (prom_dict["AP"]-bp)):
                ap.append(transcript)
            else:
                mp.append(transcript)
        elif strand == "-":
            if (transcript.end < (prom_dict["AP"]+bp )):
                ap.append(transcript)
            else:
                mp.append(transcript)

            
    final_dim=str(len(mp))+"_"+str(len(ap))+"_"+str(len(set(mp+ap))-len(ap)-len(mp))
    return final_dim 

