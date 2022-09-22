import pandas as pd
import os
from Bio import SeqIO

from header_html import write_header
from bokeh_plots import bokeh_composite, bokeh_histogram

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

def generate_main_report(config_df):
    chimeric_orf_df = get_df_from_fasta(os.path.join(out_dir,"files", "Chimeric_ORF.fasta"))
    enriched1_orf_df = get_df_from_fasta(os.path.join(out_dir,"files", "Enriched1_ORF.fasta"))


    figure_list = []
    figure_list += [bokeh_histogram("Chimeric library ORF length distribution", chimeric_orf_df, 50),
                    bokeh_histogram("Enriched library ORF length distribution", enriched1_orf_df, 50)]

    title = config_df.loc[config_df[0] == "title_of_the_run"].iloc[0, 1]
    bokeh_composite(title, 
                    figure_list, 
                    os.path.join(config_df.loc[config_df[0] == "output.dir"].iloc[0,1], "reports", title + "_main.html"), 
                    ncols=2)

if __name__ == '__main__':
    #os.chdir("c:/Users/Tatevik/Desktop/ABI/hafoe")
    config_df = pd.read_table(os.path.join(out_dir,"config", "config_file"), sep=" ", header=None)
    output_paths_df = pd.read_table(os.path.join(out_dir, "config", "output_paths"), sep=" ", header=None)

    generate_header(config_df, output_paths_df)
    generate_main_report(config_df)