import numpy as np
import pandas as pd

from Bio import AlignIO, SeqIO

import bokeh as bk

import panel as pn
import panel.widgets as pnw
pn.extension()




### load the data

clustering_data = pd.read_csv("Data/clustering_info.csv")

### Introduce an additional column with member count
clustering_data["member_count"] = clustering_data["Members_char"].str.count('AAV')
### Sort the entries by the member count
clustering_data = clustering_data.sort_values(by=['member_count'])

# set output to static HTML file
bk.plotting.output_file(filename = "orf_distribution.html", title="Static HTML file")

### Lower and upper cluster sizes to be included in the bar plot
lower_cut = 100
upper_cut = max(clustering_data["member_count"])

### Filter the data
cut_tmp = clustering_data[(clustering_data.member_count > lower_cut) & 
                          (clustering_data.member_count < upper_cut)]
### Prepare the ColumnDataSource
### reps is the representative names, counts is the size of the corresponding cluster
source = bk.models.ColumnDataSource(data = dict(reps = cut_tmp.Representative, counts = cut_tmp.member_count))

### Make a range slider for changing the lower and upper cuts
range_slider = bk.models.RangeSlider(start = 1, 
                                     end = upper_cut, 
                                     value = (lower_cut, upper_cut), 
                                     step=1, 
                                     title="Range of the cluster sizes")


### Make the log and linear bar plots
### VV: Ideally, we should have changed the y_axis_type property, but I could now access it
p_linear = bk.plotting.figure(width = 1800, 
                              height = 1300, 
                              y_range = [1e-2, 5e3], 
                              x_range = source.data["reps"].values, 
                              title="Cluster sizes", 
                              sizing_mode="scale_both", 
                              y_axis_type="linear")
    
    
p_linear.vbar(x = 'reps', 
              bottom = 1e-2, 
              top = 'counts', 
              width = 0.7, 
              source = source,
              line_color = 'white', 
              fill_color = 'lightblue')

p_log = bk.plotting.figure(width = 1800, 
                              height = 1300, 
                              y_range = [1e-2, 5e3], 
                              x_range = source.data["reps"].values, 
                              title="Cluster sizes", 
                              sizing_mode="scale_both", 
                              y_axis_type="log")
    
    
p_log.vbar(x = 'reps', 
              bottom = 1e-2, 
              top = 'counts', 
              width = 0.7, 
              source = source,
              line_color = 'white', 
              fill_color = 'lightblue')

### Rotate the x-axis labels
p_linear.xaxis.major_label_orientation = np.pi/4
p_log.xaxis.major_label_orientation = np.pi/4

### Combine the log and linear plots into separate tabs of a single panel
panels = [bk.models.widgets.Panel(child = p_linear, title="Linear scale"), 
          bk.models.widgets.Panel(child = p_log, title="Log scale")]

tabs = bk.models.widgets.Tabs(tabs = panels)

### Also show the data in a table format
    
columns = [
    bk.models.TableColumn(field = 'reps', title = 'Representative'),
    bk.models.TableColumn(field = 'counts', title = 'Cluster size')
    ]

representative_table = bk.models.DataTable(source=source, columns=columns, fit_columns=False)



def update_data(attrname, old, new):

    ### Get the current slider values
    lower_cut = range_slider.value[0]
    upper_cut = range_slider.value[1]

    ### And introduce the cut to the data
    cut_tmp = clustering_data[(clustering_data.member_count>lower_cut) & 
                              (clustering_data.member_count<upper_cut)]
    
    ### There is a bug when the filtered data is empty
    ### Should figure out how to correct this efficiently
    #if not len(cut_tmp):
    #    lower_cut = 100
    #    upper_cut = 1000
        
    #    range_slider.value[0] = lower_cut
    #    range_slider.value[1] = upper_cut
        
    #    cut_tmp = clustering_data[(clustering_data.member_count>lower_cut) & 
    #                         (clustering_data.member_count<upper_cut)]
        
    source.data = dict(reps=cut_tmp.Representative, counts=cut_tmp.member_count)

    ### Prepare new list, and update the figure
    new_list = list(source.data["reps"].values)
    p_log.x_range.factors = new_list
    p_linear.x_range.factors = new_list

    
    



range_slider.on_change('value', update_data)

bk.io.curdoc().add_root(bk.layouts.column(range_slider, tabs, representative_table, width=800))


