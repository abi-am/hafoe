import pandas as pd
import numpy as np
import os
from Bio import SeqIO
import bokeh as bk

from header_html import write_header
from bokeh_plots import view_radar, bokeh_composite, bokeh_histogram, view_composition, view_cluster_size_distribution, view_serotype_abundance, view_positional_serotype_abundance

import sys

out_dir = sys.argv[1] 

def generate_header(config_df, output_paths_df, enriched1_names):
    chimeric_orf_summary_df_path = output_paths_df.loc[output_paths_df[0] == "chimeric_orf_summary"].iloc[0,1]
    #enriched1_orf_summary_df_path = output_paths_df.loc[output_paths_df[0] == "enriched1_orf_summary"].iloc[0,1]

    enriched1_orf_summary_path_keys = ["Enriched1_"+ n + "_orf_summary" for n in enriched1_names]    
    enriched1_orf_summary_path_paths = output_paths_df.loc[output_paths_df[0].isin(enriched1_orf_summary_path_keys)].iloc[:,1]

    os.makedirs(os.path.join(config_df.loc[config_df[0] == "output.dir"].iloc[0,1], "reports", "main"), exist_ok=True)

    write_header(config_df, chimeric_orf_summary_df_path, enriched1_orf_summary_path_paths)

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

def get_clustering_df(clustering_data_path, rep_counts_path):
    clustering_data = pd.read_csv(clustering_data_path)
    rep_counts = pd.read_csv(rep_counts_path)
    ### Introduce an additional column with member count
    clustering_data["member_count"] = clustering_data["Members_char"].str.count('AAV')
    ### Sort the entries by the member count
    clustering_data = clustering_data.sort_values(by=['member_count'], ascending=False)
    clustering_data = pd.merge(clustering_data, rep_counts, on = ['Representative'])

    return clustering_data

def get_serotype_abundance_df(csv_path):
    chimeric_library_df = pd.read_csv(csv_path, header = None)
    ### Create a dictionary with chimeric_name:[serotype ids]
    variant_description_dictionary = {}
    for chimeric_library_item in chimeric_library_df[0]:
        
        chimeric_library_item_name = chimeric_library_item.split()[0]
        chimeric_library_item_parents = np.array(chimeric_library_item.split()[1:]).astype("int")
        
        variant_description_dictionary[chimeric_library_item_name] = chimeric_library_item_parents

    ### We compute the distribution of serotype abundance per representative
    ### It's a bit contrived but works fine
    chimeric_composition_array = np.array(list(variant_description_dictionary.values()))
    chimeric_representative_number, chimeric_length = np.shape(chimeric_composition_array)

    distribution_serotypes = []
    distribution_representatives = []

    for serotype_indx in range(len(serotype_names)):
        
        ### This is not very correct for several reasons. 1) it's normlized by chimeric_length, which is the same
        ### for all the reprepsentatives because they are padded by gaps (=> many are artificially enlarged),
        ### 2) we assume all the members behave exactly as the reprepsentatives,
        ### 3) we don't account for combined counts provided by Tatev
        
        composition_tmp = 100.*np.sum(chimeric_composition_array == serotype_indx, axis = 1)/chimeric_length
        composition_tmp = np.around(composition_tmp, decimals = 2)
        
        ### This is a required for the sns violin plot
        distribution_serotypes += [serotype_names[serotype_indx]]*chimeric_representative_number
        distribution_representatives.append(list(composition_tmp))
        
    distribution_serotypes = np.array(distribution_serotypes)
    distribution_representatives = np.array(distribution_representatives).flatten()
        
    distribution_serotypes.astype('str');
    distribution_representatives.astype('float64');
    
    ### Create a data-frame
    serotype_dist_data_frame = pd.DataFrame(np.array([distribution_serotypes, distribution_representatives]).T, 
                                                    columns = ["Serotype", "Frequency"])
    serotype_dist_data_frame['Frequency'] = serotype_dist_data_frame['Frequency'].astype('float64')

    #Remove columns "gap","Multiple alignments","No alignment"
    serotype_dist_data_frame_del = serotype_dist_data_frame[~serotype_dist_data_frame.Serotype.isin(["gap","Multiple alignments","No alignment"])]
    return serotype_dist_data_frame_del

def get_positional_serotype_abundance_matrix(csv_path):
    chimeric_library_df = pd.read_csv(csv_path, header = None)
    ### Create a dictionary with chimeric_name:[serotype ids]
    variant_description_dictionary = {}
    for chimeric_library_item in chimeric_library_df[0]:
        
        chimeric_library_item_name = chimeric_library_item.split()[0]
        chimeric_library_item_parents = np.array(chimeric_library_item.split()[1:]).astype("int")
        
        variant_description_dictionary[chimeric_library_item_name] = chimeric_library_item_parents

    ### 2-D array of chimeric composition
    chimeric_composition_array = np.array(list(variant_description_dictionary.values()))

    ### Initialize the abundance matrix ---> rows:serotypes, columns:positional abundance 
    positional_serotype_abundance = np.zeros((len(serotype_names), np.shape(chimeric_composition_array)[1]))

    ### Populate the abundance matrix
    for serotype_indx in np.arange(0, len(serotype_names), 1):
        positional_serotype_abundance[serotype_indx] = np.sum(chimeric_composition_array == serotype_indx, axis = 0)
    
    # Remove gap, no alignment, multiple alignment
    positional_serotype_abundance_del = np.delete(positional_serotype_abundance, idx_del, 0)

    return positional_serotype_abundance_del

def get_radar_df(csv_path, fc_threshold=1):
    # quantile_prob , count_threshold = ? 

    counts_mor_df = pd.read_csv(csv_path)
    numeric_cols = counts_mor_df.columns[counts_mor_df.columns.str.contains('Count|FC')]
    counts_mor_df[numeric_cols] = counts_mor_df[numeric_cols].astype('float32')

    count_cols = counts_mor_df.columns[counts_mor_df.columns.str.contains('Count')][1:]
    fc_cols = counts_mor_df.columns[counts_mor_df.columns.str.contains('FC')]

    ### Filter enriched variants 
    # !!! maybe add filtering by normalized counts
    counts_mor_df = counts_mor_df[(counts_mor_df[fc_cols] > fc_threshold).any(axis=1)].reset_index()

    return counts_mor_df   


def generate_main_report(config_df):
    os.makedirs(os.path.join(config_df.loc[config_df[0] == "output.dir"].iloc[0,1], "reports", "supplement"), exist_ok=True)

    chimeric_orf_df = get_df_from_fasta(output_paths_df.loc[output_paths_df[0] == "final_chimeric_orf_fasta"].iloc[0,1]) 
    chimeric_lib_vd_dictionary = get_vd_dictionary(output_paths_df.loc[output_paths_df[0] == "representatives_variant_description"].iloc[0,1])
    chimeric_lib_top20_vd_dictionary = get_vd_dictionary(output_paths_df.loc[output_paths_df[0] == "representatives_variant_description_top20"].iloc[0,1])
    chimeric_lib_clustering_info = get_clustering_df(output_paths_df.loc[output_paths_df[0] == "chimeric_lib_clustering_info"].iloc[0,1], output_paths_df.loc[output_paths_df[0] == "chimeric_lib_representatives_counts"].iloc[0,1])
    serotype_dist_data_frame_del = get_serotype_abundance_df(output_paths_df.loc[output_paths_df[0] == "representatives_variant_description"].iloc[0,1])
    positional_serotype_abundance_del = get_positional_serotype_abundance_matrix(output_paths_df.loc[output_paths_df[0] == "representatives_variant_description_top20"].iloc[0,1])
    counts_mor_df = get_radar_df(output_paths_df.loc[output_paths_df[0] == "counts_mor_df"].iloc[0,1], fc_threshold=2)

    enriched1_hist_list = []
    for n in enriched1_names:
        final_enriched1_orf_path_key = "final_Enriched1_"+ n + "_orf_fasta" 
        final_enriched1_orf_path_path = output_paths_df.loc[output_paths_df[0] == final_enriched1_orf_path_key].iloc[0,1]
        enriched1_orf_df = get_df_from_fasta(final_enriched1_orf_path_path) 
        enriched1_hist_list += [[bokeh_histogram(n + " enriched library ORF length distribution", enriched1_orf_df, 20, plot_width=600)]]
        # enriched1_hist_list += [[bokeh_histogram(n + " enriched library ORF length distribution", enriched1_orf_df, 20, svg_dir=os.path.join(config_df.loc[config_df[0] == "output.dir"].iloc[0,1], "reports", "supplement"), plot_width=600)]]
    
    figure_layout = []
    figure_layout += [[view_cluster_size_distribution(chimeric_lib_clustering_info, lower_cut_size_ = 100, lower_cut_abundance_ = 100)]]
    figure_layout += [[view_composition("Compositions of the 20 most abundant representatives in chimeric library (with multiple sequence alignment)", chimeric_lib_top20_vd_dictionary, serotype_colors, serotype_names)]]
    #figure_layout += [[view_composition("Compositions of chimeric library representatives", chimeric_lib_vd_dictionary, serotype_colors, serotype_names)]]
    figure_layout += [[bk.layouts.row(view_serotype_abundance("Serotype abundance in chimeric library",
                                            data=serotype_dist_data_frame_del, 
                                            serotype_names=serotype_names_del, 
                                            serotype_colors=serotype_colors_del),
                                    view_positional_serotype_abundance("Positional serotype abundance in 20 most abundant representatives of chimeric library",
                                                                    positional_serotype_abundance_del, 
                                                                    serotype_names_del, 
                                                                    serotype_colors_del), 
                                                                    width = 1500)]]
    figure_layout += [[bk.layouts.row(view_radar(counts_mor_df, type="Count"), 
                                      view_radar(counts_mor_df, type="FC"), 
                                      width = 1500)]]
    #figure_layout += [[[bokeh_histogram("Chimeric library ORF length distribution", chimeric_orf_df, 20, plot_width=600)]] + enriched1_hist_list]

    title = config_df.loc[config_df[0] == "title_of_the_run"].iloc[0, 1]
    bokeh_composite(title + "_main",  
                    figure_layout = figure_layout, 
                    filename = os.path.join(config_df.loc[config_df[0] == "output.dir"].iloc[0,1], "reports", "main", title + "_main.html"))

    #supplement
    bokeh_composite(title + "_chimeric_library_composition",  
                    figure_layout = [view_composition("Compositions of chimeric library representatives", chimeric_lib_vd_dictionary, serotype_colors, serotype_names)], 
                    filename = os.path.join(config_df.loc[config_df[0] == "output.dir"].iloc[0,1], "reports", "supplement", title + "_chimeric_library_composition.html"))

    bokeh_composite(title + "_orf_quality_control",  
                    figure_layout = [[[bokeh_histogram("Chimeric library ORF length distribution", chimeric_orf_df, 20, plot_width=600)]] + enriched1_hist_list], 
                    filename = os.path.join(config_df.loc[config_df[0] == "output.dir"].iloc[0,1], "reports", "supplement", title + "_orf_quality_control.html"))

if __name__ == '__main__':
    config_df = pd.read_table(os.path.join(out_dir,"config", "config_file"), sep=" ", header=None)
    output_paths_df = pd.read_table(os.path.join(out_dir, "config", "output_paths"), sep=" ", header=None)

    # Enriched1 names
    enriched1_path = os.path.abspath(config_df.loc[config_df[0] == "enriched"].iloc[0, 1])
    def list_full_paths(directory):
        return [os.path.join(directory, file) for file in os.listdir(directory)]
    enriched1_files = list_full_paths(os.path.abspath(enriched1_path))
    enriched1_files = [os.path.abspath(s) for s in enriched1_files]
    enriched1_names = [os.path.basename(s) for s in enriched1_files]
    enriched1_names = [s.split(".", 1)[0] for s in enriched1_names]

    # Serotypes
    serotype_names = np.array(["No alignment", "AAV1", "AAV2", "AAV3", "AAV4", "AAV5", "AAV6", 
                  "AAV7", "AAV8", "AAV9", "AAV10", "AAV11", "AAV12", "AAV13", 
                  "AAVrh8", "AAVrh10", "AAVrh32", "Multiple alignments", "gap"])
    serotype_colors = np.array(["gray", "#AA4488", "#4477AA", "#AAAA44", 
                 "#AA7744", "#AA4455", "#44AAAA", "#771155", 
                 "#114477", "#777711", "#774411", "#771122", 
                 "#117777", "#A6CEE3", "#B2DF8A", "#F1B6DA", 
                 "#B2ABD2", "black", "white"])

    ### remove gap, multiple, no alignments from serotype_names, serotype_colprs
    #find index location of first occurrence of each value of interest
    sorter = np.argsort(serotype_names)
    idx_del = sorter[np.searchsorted(serotype_names, ["gap", "Multiple alignments", "No alignment"], sorter=sorter)]
    serotype_names_del = np.delete(serotype_names, idx_del, 0)
    serotype_colors_del = np.delete(serotype_colors, idx_del, 0)

    generate_header(config_df, output_paths_df, enriched1_names)
    generate_main_report(config_df)