import os
from turtle import width
#import seaborn as sns
import bokeh as bk

import panel as pn
import panel.widgets as pnw
import colorlover as cl
import matplotlib as mpl
import matplotlib.cm as cm
import numpy as np
import pandas as pd
from bokeh.colors import RGB
from bokeh.io import output_file, save, export_svgs, export_png
from bokeh.io import show
from bokeh.layouts import gridplot, row, widgetbox
from bokeh.models import HoverTool, Arrow, NormalHead, ColumnDataSource, LinearColorMapper, FactorRange, DataTable, \
    PanTool, BoxZoomTool, WheelZoomTool, SaveTool, ResetTool, Div, Legend, Panel, Tabs, LabelSet
from bokeh.plotting import figure, show
from bokeh.transform import transform
from bokeh.core.properties import value, field
from header_html import get_div_logo, get_div_footer, get_div_title

tools = [PanTool(), BoxZoomTool(), WheelZoomTool(), SaveTool(), ResetTool()]

LEGEND_POSITION = 'above'
LEGEND_ALIGNMENT = 'bottom_left'
COMBINED = "combined"


def bokeh_composite(title, figure_list, filename, ncols=2):

    output_file(title + ".html", title=title)

    p = gridplot(figure_list, ncols=ncols, toolbar_location="left")
    div_logo = Div(text=get_div_logo())
    div_title = Div(text=get_div_title(title))
    save([div_logo, div_title, p], filename=filename)

def bokeh_histogram(title, df, n_bins, plot_width=800):
    hist, edges = np.histogram(df['Length'], density=False, bins=n_bins)
    p = figure(title=title, tools='') #, width=plot_width
    p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],
           fill_color="lightblue", line_color="white")

    p.y_range.start = 0
    #p.legend.location = "center_right"
    #p.legend.background_fill_color = "#fefefe"
    p.xaxis.axis_label = 'Length'
    p.yaxis.axis_label = 'Frequency'
    #p.grid.grid_line_color="white"

    return p

def view_alignment(title, serotype_dictionary, serotype_colors, serotype_names, include_text = False, plot_width = 800, max_num_x = 3000):
    
    ### seqs: Serotype list
    ### keys: Library name
    seqs = np.array(list(serotype_dictionary.values()))
    keys = np.array(list(serotype_dictionary.keys()))
    
    N = len(seqs[0])
    S = len(seqs)    
    width = 1.
    
    tmp = [i for s in seqs for i in s]
    colors = [serotype_colors[i] for i in tmp]
    serotypes = [serotype_names[i] for i in tmp]
    
    x = np.arange(1, N + 1)
    y = np.arange(0, S, 1)    
    xx, yy = np.meshgrid(x, y)
    
    gx = xx.flatten()
    gy = yy.flatten()
    recty = gy + .1
    h = 1/S
    
    plot_height = len(seqs)*15 + 50
    x_range = bk.models.Range1d(0, N + 1, bounds = 'auto')
    if N > max_num_x:
        viewlen = max_num_x
    else:
        viewlen = N
    view_range = (0, viewlen) 
    
    
    
    source = bk.models.ColumnDataSource(dict(x=gx, y=gy, recty=recty, colors = colors, serotypes = serotypes))
    
    p = bk.plotting.figure(title = title, 
               plot_width = plot_width, plot_height = plot_height, 
               x_range = view_range, y_range = keys, tools = "xwheel_zoom, xpan, reset", 
               min_border = 0, toolbar_location = 'below')
    
    if include_text == True:
        glyph = bk.models.glyphs.Text(x="x", y="y", text="text", 
                                   text_align='center',
                                   text_color="black",
                                   text_font=value("monospace"))
                                   #text_font_size=fontsize)
    elif include_text == False:
        glyph = bk.models.glyphs.Text(x="x", y="y", text="text", 
                                   text_align='center',
                                   text_color="black",
                                   text_font=value("monospace"))
        
    rects = bk.models.glyphs.Rect(x = "x", y = "recty",  
                 width = 1., height = 1., 
                 fill_color = "colors", 
                 line_color = None, fill_alpha = 1.)
    
    p.add_glyph(source, glyph)
    p.add_glyph(source, rects)
  
    p.grid.visible = False
    p.xaxis.axis_label = 'Nucleotide position'
    p.xaxis.major_label_text_font_style = "bold"
    p.yaxis.minor_tick_line_width = 0
    p.yaxis.major_tick_line_width = 0
    
    
    tooltips = [("serotype", "@serotypes"),]
    p.add_tools(bk.models.HoverTool(tooltips=tooltips))
    
    return p

def view_cluster_size_distribution(clustering_data, lower_cut = 100):
    ### Lower and upper cluster sizes to be included in the bar plot
    lower_cut = lower_cut
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
                                #y_range = [1e-2, 5e3], 
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
                                #y_range = [1e-2, 5e3], 
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
    
    return bk.layouts.column(range_slider, tabs, representative_table, width=800)
    #return p_linear