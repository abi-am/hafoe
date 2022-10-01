import pandas as pd
import numpy as np
#import seaborn as sns
import os
from Bio import SeqIO

from header_html import write_header
from bokeh_plots import bokeh_composite, bokeh_histogram, view_alignment, view_cluster_size_distribution

import sys

out_dir = sys.argv[1] 

def generate_header(config_df, output_paths_df):
    chimeric_orf_summary_df = pd.read_csv(output_paths_df.loc[output_paths_df[0] == "chimeric_orf_summary"].iloc[0,1])
    enriched1_orf_summary_df = pd.read_csv(output_paths_df.loc[output_paths_df[0] == "enriched1_orf_summary"].iloc[0,1])

    write_header(config_df, chimeric_orf_summary_df, enriched1_orf_summary_df)

def get_df_from_fasta(path):
    with open(path) as fasta_file:  
        identifiers = []
        #sequences = []
        lengths = []
        for seq_record in SeqIO.parse(fasta_file, 'fasta'): 
            identifiers.append(seq_record.id)
            #sequences.append(seq_record.seq)
            lengths.append(len(seq_record.seq))

    df = pd.DataFrame({'Header':identifiers,'Length':lengths})
    return df

def get_vd_dictionary(csv_path):
    chimeric_library_df = pd.read_csv(csv_path, header = None)
    ### Create a dictionary with chimeric_name:[serotype ids]
    variant_description_dictionary = {}
    for chimeric_library_item in chimeric_library_df[0]:
        
        chimeric_library_item_name = chimeric_library_item.split()[0]
        chimeric_library_item_parents = np.array(chimeric_library_item.split()[1:]).astype("int")
        
        variant_description_dictionary[chimeric_library_item_name] = chimeric_library_item_parents
    return variant_description_dictionary

def get_clustering_df(clustering_data_path):
    clustering_data = pd.read_csv(clustering_data_path)
    ### Introduce an additional column with member count
    clustering_data["member_count"] = clustering_data["Members_char"].str.count('AAV')
    ### Sort the entries by the member count
    clustering_data = clustering_data.sort_values(by=['member_count'], ascending=False)  
    return clustering_data

def generate_main_report(config_df):
    chimeric_orf_df = get_df_from_fasta(output_paths_df.loc[output_paths_df[0] == "final_chimeric_orf_fasta"].iloc[0,1]) 
    enriched1_orf_df = get_df_from_fasta(output_paths_df.loc[output_paths_df[0] == "final_enriched_orf_fasta"].iloc[0,1]) 
    chimeric_lib_vd_dictionary = get_vd_dictionary(output_paths_df.loc[output_paths_df[0] == "representatives_variant_description"].iloc[0,1])
    chimeric_lib_top20_vd_dictionary = get_vd_dictionary(output_paths_df.loc[output_paths_df[0] == "representatives_variant_description_top20_msa"].iloc[0,1])
    chimeric_lib_clustering_info = get_clustering_df(output_paths_df.loc[output_paths_df[0] == "chimeric_lib_clustering_info"].iloc[0,1])

    

    figure_list = []
    figure_list += [bokeh_histogram("Chimeric library ORF length distribution", chimeric_orf_df, 50, plot_width=600),
                    bokeh_histogram("Enriched library ORF length distribution", enriched1_orf_df, 50, plot_width=600)]
    figure_list += [view_alignment("Compositions of chimeric library representatives", chimeric_lib_vd_dictionary, serotype_colors, serotype_names), None]
    figure_list += [view_alignment("Compositions of the 20 most abundant representatives in chimeric library (with multiple sequence alignment)", chimeric_lib_top20_vd_dictionary, serotype_colors, serotype_names), None]
    figure_list += [view_cluster_size_distribution(chimeric_lib_clustering_info, lower_cut=0), None]

    title = config_df.loc[config_df[0] == "title_of_the_run"].iloc[0, 1]
    bokeh_composite(title + "_main",  
                    figure_list, 
                    os.path.join(config_df.loc[config_df[0] == "output.dir"].iloc[0,1], "reports", title + "_main.html"), 
                    ncols=2)

if __name__ == '__main__':
    #os.chdir("c:/Users/Tatevik/Desktop/ABI/hafoe")
    config_df = pd.read_table(os.path.join(out_dir,"config", "config_file"), sep=" ", header=None)
    output_paths_df = pd.read_table(os.path.join(out_dir, "config", "output_paths"), sep=" ", header=None)

    serotype_names = ["No alignment", "AAV1", "AAV2", "AAV3", "AAV4", "AAV5", "AAV6", 
                  "AAV7", "AAV8", "AAV9", "AAV10", "AAV11", "AAV12", "AAV13", 
                  "AAVrh8", "AAVrh10", "AAVrh32", "Multiple alignments", "gap"]
    serotype_colors = ['#124984', '#2065ab', '#327cb7', '#4393c3', '#6bacd1', '#90c4dd', '#b1d5e7', '#d1e5f0', '#e4eef4', '#f7f6f6', '#fae9df', '#fddbc7', '#f8bfa4', '#f3a481', '#e48066', '#d6604d', '#c43b3c', '#b1182b', '#8a0b25'] #sns.color_palette("RdBu_r", len(serotype_names)).as_hex()

    generate_header(config_df, output_paths_df)
    generate_main_report(config_df)