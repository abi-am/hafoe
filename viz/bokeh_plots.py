import bokeh as bk
import numpy as np
import pandas as pd
import seaborn as sns
import math
from bokeh.io import output_file, save, export_svgs, export_png
from bokeh.io import show
from bokeh.layouts import gridplot, row, widgetbox, layout, grid
from bokeh.models import HoverTool, Arrow, NormalHead, ColumnDataSource, LinearColorMapper, FactorRange, DataTable, \
    PanTool, BoxZoomTool, WheelZoomTool, SaveTool, ResetTool, Div, Legend, Panel, Tabs, LabelSet
from bokeh.plotting import figure, show
from bokeh.transform import transform, jitter, factor_cmap
from bokeh.core.properties import value, field
from header_html import get_div_logo, get_div_footer, get_div_title

tools = [PanTool(), BoxZoomTool(), WheelZoomTool(), SaveTool(), ResetTool()]

LEGEND_POSITION = 'above'
LEGEND_ALIGNMENT = 'bottom_left'
COMBINED = "combined"


def bokeh_composite(title, figure_layout, filename):

    output_file(title + ".html", title=title)

    #p = gridplot(figure_list, ncols=ncols, toolbar_location="left")
    p = layout(figure_layout, sizing_mode="scale_width")
    # ([
    #             [bollinger],
    #             [sliders, plot],
    #             [p1, p2, p3],
    #         ])
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

def view_composition(title, serotype_dictionary, serotype_colors, serotype_names, include_text = False, plot_width = 1500, max_num_x = 3000):
    
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

    # layout = bk.layouts.column(p, height=600, width=1500)
    return p


def view_serotype_abundance(title, data, serotype_names, serotype_colors):
    cats = serotype_names
    # find the quartiles and IQR for each category
    groups = data.groupby('Serotype')
    q1 = groups.quantile(q=0.25)
    q2 = groups.quantile(q=0.5)
    q3 = groups.quantile(q=0.75)
    iqr = q3 - q1
    upper = q3 + 1.5*iqr
    lower = q1 - 1.5*iqr

    q1 = q1.loc[cats]
    q2 = q2.loc[cats]
    q3 = q3.loc[cats]
    iqr = iqr.loc[cats]
    lower = lower.loc[cats]
    upper = upper.loc[cats]

    # find the outliers for each category
    # def outliers(group):
    #     cat = group.name
    #     return group[(group.Frequency > upper.loc[cat][0]) | (group.Frequency < lower.loc[cat][0])]['Frequency']
    # out = groups.apply(outliers).dropna()

    # # prepare outlier data for plotting, we need coordinates for every outlier.
    # if not out.empty:
    #     outx = list(out.index.get_level_values(0))
    #     outy = list(out.values)

    p = figure(title = title, 
               tools="", x_range=cats, width=650, height=600,
               toolbar_location=None,
               x_axis_label = "Serotype", 
               y_axis_label = "Frequency")

    # if no outliers, shrink lengths of stems to be no longer than the minimums or maximums
    qmin = groups.quantile(q=0.00)
    qmax = groups.quantile(q=1.00)

    qmin = qmin.loc[cats]
    qmax = qmax.loc[cats]

    upper.Frequency = [min([x,y]) for (x,y) in zip(list(qmax.loc[:,'Frequency']),upper.Frequency)]
    lower.Frequency = [max([x,y]) for (x,y) in zip(list(qmin.loc[:,'Frequency']),lower.Frequency)]

    # stems
    p.segment(cats, upper.Frequency, cats, q3.Frequency, line_color="black")
    p.segment(cats, lower.Frequency, cats, q1.Frequency, line_color="black")

    # boxes
    p.vbar(cats, 0.7, q2.Frequency, q3.Frequency, fill_color=serotype_colors, line_color="black") #"#E08E79"
    p.vbar(cats, 0.7, q1.Frequency, q2.Frequency, fill_color=serotype_colors, line_color="black") #"#3B8686"

    # whiskers (almost-0 height rects simpler than segments)
    p.rect(cats, lower.Frequency, 0.2, 0.01, line_color="black")
    p.rect(cats, upper.Frequency, 0.2, 0.01, line_color="black")

    # outliers
    # if not out.empty:
    #     p.circle(outx, outy, size=6, color="black", fill_alpha=0.3) #"#F38630"

    p1 = p.scatter(x=jitter('Serotype', width=0.5, range=p.x_range), y='Frequency', source=data, 
                  size=4, color=factor_cmap('Serotype', serotype_colors, serotype_names), line_color="black",
                  alpha=0.4)


    p.xgrid.grid_line_color = None
    p.xaxis.major_label_orientation = math.pi/4

    ### Add a tooltip
    tooltips = [("frequency", "@Frequency"),]
    p.add_tools(bk.models.HoverTool(renderers=[p1], tooltips = tooltips))
    return p

def view_positional_serotype_abundance(title, positional_serotype_abundance, serotype_names, serotype_colors):
    # def update_data(attrname, old, new):

    #     ### Get the current slider value
    #     smoothing_window_size = smoothing_slider.value
            
    #     smoothed_serotype_abundance = np.copy(positional_serotype_abundance)
    #     for serotype_indx in np.arange(0, len(serotype_names), 1):
    #         abundance_tmp = smoothed_serotype_abundance[serotype_indx]
    #         smoothed_serotype_abundance[serotype_indx] = np.convolve(abundance_tmp, 
    #                                                                 np.ones(smoothing_window_size), 
    #                                                                 'same')/smoothing_window_size
            
    #     ### update the data container
    #     source.data["abundance"] = list(smoothed_serotype_abundance)

    ### The smoothing slider
    smoothing_slider = bk.models.Slider(start = 1, end = 500, value = 100, step = 1, title = "Smoothing length")

    ### Prepare the plotting data
    # placeholder array
    positions_base = np.arange(0, np.shape(positional_serotype_abundance)[1])
    # the main data container
    data_dict = {}
    data_dict["positions"] = [positions_base for _ in np.arange(0, np.shape(positional_serotype_abundance)[0])]
    data_dict["abundance"] = list(positional_serotype_abundance)
    data_dict["color"] = serotype_colors
    data_dict["serotypes"] = serotype_names
    source = bk.plotting.ColumnDataSource(data = data_dict)
    
    # smooth the initially displayed data with default value
    smoothed_serotype_abundance = np.copy(positional_serotype_abundance)
    for serotype_indx in np.arange(0, len(serotype_names), 1):
        abundance_tmp = smoothed_serotype_abundance[serotype_indx]
        smoothed_serotype_abundance[serotype_indx] = np.convolve(abundance_tmp, 
                                                                np.ones(100), 
                                                                'same')/100  
    source.data["abundance"] = list(smoothed_serotype_abundance)

    ### Start plotting
    p = bk.plotting.figure(title = title, 
                            width = 850, height=550,
                            x_axis_label = "Nucleotide position", 
                            y_axis_label = "Serotype abundance")

    p.multi_line(xs = 'positions', ys = 'abundance', source = source,
                                                    line_width = 3, 
                                                    color = 'color',
                                                    legend = 'serotypes')

    p.legend.location = "center"
    #p.legend.click_policy="hide" not working

    leg = p.legend[0]
    p.add_layout(leg,'right') 

    ### Add a tooltip
    tooltips = [("serotype", "@serotypes"),]
    p.add_tools(bk.models.HoverTool(tooltips = tooltips))


    update_data = bk.models.CustomJS(args = dict(source=source, external_data = data_dict), code="""
        var data = source.data;
        var smoothing_window_size = cb_obj.value;
        
        const x_external = external_data["positions"];
        const y_external = external_data["abundance"];
        
        var chunk; var chunk_start; var chunk_end;
        
        data["positions"] = x_external;
        
        console.log("entering loops")
        
        for (var i = 0; i < x_external.length; i++) {

            for (var j = 0; j <= y_external[i].length - Math.floor(smoothing_window_size/2); j++){
            
                if (j - Math.floor(smoothing_window_size/2) < 0)
                    {chunk_start = 0;}
                else
                    {chunk_start = j - Math.floor(smoothing_window_size/2);}
                if (j + Math.floor(smoothing_window_size/2) > y_external[i].length)
                    {chunk_end = y_external[i].length;}
                else
                    {chunk_end = j + Math.floor(smoothing_window_size/2);}
                    
                chunk = y_external[i].slice(chunk_start, chunk_end + 1);
                console.log(chunk)
                data["abundance"][i][j] = chunk.reduce((total, current) => total + current)/chunk.length;
                
            }
        }
        
        
        
        source.change.emit();
        
    """)


    #smoothing_slider.on_change('value', update_data)
    smoothing_slider.js_on_change('value', update_data)


    # show the results
    layout = bk.layouts.column(smoothing_slider, p, width = 750)

    return layout

def view_cluster_size_distribution(clustering_data, topn = 20, lower_cut_size_ = 100, lower_cut_abundance_ = 100):
    ### Lower and upper cluster sizes to be included in the bar plot
    # lower_cut_size = lower_cut_size_
    # upper_cut_size = max(clustering_data["member_count"])+10
    # lower_cut_abundance = lower_cut_abundance_
    # upper_cut_abundance = max(clustering_data["Chimeric.Count"])+10
    # ### Filter the data
    # cut_tmp = clustering_data[(clustering_data["member_count"] > lower_cut_size) & 
    #                           (clustering_data["member_count"] < upper_cut_size) &
    #                           (clustering_data["Chimeric.Count"] > lower_cut_abundance) &
    #                           (clustering_data["Chimeric.Count"] < upper_cut_abundance)
    #                          ]

    # take top 20 by cluster size
    clustering_data = clustering_data.sort_values(by=['member_count'], ascending=False)
    cut_tmp = clustering_data.head(topn)
    ### Prepare the ColumnDataSource
    ### reps is the representative names, sizes is the size of the corresponding cluster
    source = bk.models.ColumnDataSource(data = dict(reps = cut_tmp.Representative, sizes = cut_tmp["member_count"], abundances = cut_tmp["Chimeric.Count"]))

    ### Make a range slider for changing the lower and upper cuts
    # range_slider_size = bk.models.RangeSlider(start = 1, 
    #                                      end = upper_cut_size, 
    #                                      value = (lower_cut_size, upper_cut_size), 
    #                                      step=1, 
    #                                      title="Range of the cluster sizes")
    # range_slider_abundance = bk.models.RangeSlider(start = 1, 
    #                                      end = upper_cut_abundance, 
    #                                      value = (lower_cut_abundance, upper_cut_abundance), 
    #                                      step=1, 
    #                                      title="Range of the cluster abundances")


    ### Make the log and linear bar plots
    ### VV: Ideally, we should have changed the y_axis_type property, but I could now access it

    ### Size plots
    p_linear_size = bk.plotting.figure(height=500, width=1100,
                                  y_range = [1e-2, math.ceil(max(source.data["sizes"].values)/100)*100],  
                                  x_range = source.data["reps"].values, 
                                  title="Chimeric library cluster size for top "+str(topn)+" clusters", 
                                  sizing_mode="scale_both", 
                                  y_axis_type="linear")


    p_linear_size.vbar(x = 'reps', 
                  bottom = 1e-2, 
                  top = 'sizes', 
                  width = 0.7, 
                  source = source,
                  line_color = 'white', 
                  fill_color = 'lightblue')

    p_log_size = bk.plotting.figure(height=500, width=1100,
                                  y_range = [1e-2, math.ceil(max(source.data["sizes"].values)*10)], 
                                  x_range = source.data["reps"].values, 
                                  title="Chimeric library cluster size for top "+str(topn)+" clusters", 
                                  sizing_mode="scale_both", 
                                  y_axis_type="log")


    p_log_size.vbar(x = 'reps', 
                  bottom = 1e-2, 
                  top = 'sizes', 
                  width = 0.7, 
                  source = source,
                  line_color = 'white', 
                  fill_color = 'lightblue')

    ### Rotate the x-axis labels
    p_linear_size.xaxis.major_label_orientation = np.pi/4
    p_log_size.xaxis.major_label_orientation = np.pi/4

    ### Abundance plots
    p_linear_abundance = bk.plotting.figure(height=500, width=1100,
                                          y_range = [1e-2, math.ceil(max(source.data["abundances"].values)/100)*100], 
                                          x_range = source.data["reps"].values, 
                                          title="Chimeric library cluster abundance for top "+str(topn)+" clusters", 
                                          sizing_mode="scale_both", 
                                          y_axis_type="linear")


    p_linear_abundance.vbar(x = 'reps', 
                          bottom = 1e-2, 
                          top = 'abundances', 
                          width = 0.7, 
                          source = source,
                          line_color = 'white', 
                          fill_color = 'lightblue')

    p_log_abundance = bk.plotting.figure(height=500, width=1100,
                                  y_range = [1e-2, max(source.data["abundances"].values)*10], 
                                  x_range = source.data["reps"].values, 
                                  title="Chimeric library cluster abundance for top "+str(topn)+" clusters", 
                                  sizing_mode="scale_both", 
                                  y_axis_type="log")


    p_log_abundance_bars = p_log_abundance.vbar(x = 'reps', 
                                                bottom = 1e-2, 
                                                top = 'abundances', 
                                                width = 0.7, 
                                                source = source,
                                                line_color = 'white', 
                                                fill_color = 'lightblue')

    ### Add log transformed values to source for hovering
    log_sizes = np.log10(source.data["sizes"].values)
    source.data['log_sizes'] = log_sizes
    log_abundances = np.log10(source.data["abundances"].values)
    source.data['log_abundances'] = log_abundances

    ### Rotate the x-axis labels
    p_linear_abundance.xaxis.major_label_orientation = np.pi/4
    p_log_abundance.xaxis.major_label_orientation = np.pi/4

    tooltips1 = [("Cluster size", "@sizes"),]
    p_linear_size.add_tools(bk.models.HoverTool(tooltips=tooltips1))

    tooltips2 = [("Cluster abundance", "@abundances"),]
    p_linear_abundance.add_tools(bk.models.HoverTool(tooltips=tooltips2))

    tooltips3 = [("Cluster size", "@log_sizes"),]
    p_log_size.add_tools(bk.models.HoverTool(tooltips=tooltips3))

    tooltips4 = [("Cluster abundance", "@log_abundances"),]
    p_log_abundance.add_tools(bk.models.HoverTool(renderers=[p_log_abundance_bars], tooltips=tooltips4))

    ### Combine the log and linear plots into separate tabs of a single panel
    panels = [bk.models.widgets.Panel(child = p_linear_size, title="Cluster size, linear scale"), 
              bk.models.widgets.Panel(child = p_linear_abundance, title="Cluster abundance, linear scale"), 
              bk.models.widgets.Panel(child = p_log_size, title="Cluster size, log scale"),
              bk.models.widgets.Panel(child = p_log_abundance, title="Cluster abundance, log scale")]

    tabs = bk.models.widgets.Tabs(tabs = panels)

    ### Also show the data in a table format

    columns = [
        bk.models.TableColumn(field = 'reps', title = 'Representative'),
        bk.models.TableColumn(field = 'sizes', title = "Cluster size"),
        bk.models.TableColumn(field = 'abundances', title = "Cluster abundance")
        ]

    representative_table = bk.models.DataTable(source=source, columns=columns, fit_columns=True, height=500, width=400)


    # def update_data(attrname, old, new):

    #     ### Get the current slider values
    #     lower_cut_size = range_slider_size.value[0]
    #     upper_cut_size = range_slider_size.value[1]
    #     lower_cut_abundance = range_slider_abundance.value[0]
    #     upper_cut_abundance = range_slider_abundance.value[1]


    #     ### And introduce the cut to the data
    #     cut_tmp = clustering_data[(clustering_data["member_count"]>=lower_cut_size) & 
    #                               (clustering_data["member_count"]<=upper_cut_size) &
    #                               (clustering_data["Chimeric.Count"]>=lower_cut_abundance) & 
    #                               (clustering_data["Chimeric.Count"]<=upper_cut_abundance)]

    #     ### There is a bug when the filtered data is empty
    #     ### Should figure out how to correct this efficiently
    #     #if not len(cut_tmp):
    #     #    lower_cut = 100
    #     #    upper_cut = 1000

    #     #    range_slider.value[0] = lower_cut
    #     #    range_slider.value[1] = upper_cut

    #     #    cut_tmp = clustering_data[(clustering_data.member_count>lower_cut) & 
    #     #                         (clustering_data.member_count<upper_cut)]

    #     source.data = dict(reps=cut_tmp.Representative, sizes=cut_tmp["member_count"], abundances=cut_tmp["Chimeric.Count"])

    #     ### Prepare new list, and update the figure
    #     new_list = list(source.data["reps"].values)
    #     p_log_size.x_range.factors = new_list
    #     p_linear_size.x_range.factors = new_list
    #     p_log_abundance.x_range.factors = new_list
    #     p_linear_abundance.x_range.factors = new_list

    # range_slider_size.on_change('value', update_data)
    # range_slider_abundance.on_change('value', update_data)

    # layout = bk.layouts.column(bk.layouts.column(range_slider_size, range_slider_abundance, width=1500), bk.layouts.row(tabs, representative_table, width=1500, height=550))
    layout = bk.layouts.row(tabs, representative_table, width=1500, height=550)

    return layout


def view_radar(df, type, fig_width = 720, fig_height = 700):
	sample_names = [i.replace("."+type,"") for i in df.columns[df.columns.str.contains(type)]]

	### Assign colors to samples
	num_samples = len(df.columns[df.columns.str.contains('FC')])
	palette = sns.color_palette("Set1", num_samples+1).as_hex()
	sample_color =  {} 

	if type == "FC":
		cols = df.columns[df.columns.str.contains('FC')]
	elif type == "Count":
		sample_color["Chimeric.Count"] = palette[0]
		cols = df.columns[df.columns.str.contains('Count')][1:]

	for i in range(len(cols)):
		sample_color[cols[i]] = palette[i+1]
	
	alpha = np.pi/25.

	inner_radius = 20
	start_radius = 90
	outer_radius = 300 - 10

	if type == "FC":
		title = "log2 Fold Changes of AAV variants (enriched library/chimeric library)"
	elif type == "Count":
		title = "log2 Normalized Counts of AAV variants"

	# main frame of the fig
	p = bk.plotting.figure(width = fig_width, height = fig_height, title = title,
			x_axis_type = None, y_axis_type = None,
			x_range = (-420, 420), y_range = (-420, 420),
			min_border = 0, match_aspect = True) # outline_line_color = "black", 

	# bace wedge
	p.annular_wedge(0, 0, inner_radius, outer_radius, 
					np.pi/2. + alpha, np.pi/2. - alpha, 
					color = "#eaeaea", name = "Base wedge")

	# circular grids and their labels
	max_value = round(np.log2(df[df.columns[df.columns.str.contains(type)]].replace(0., np.nan)).max().max())
	labels = np.linspace(0, max_value, 5)
	radii = np.linspace(start_radius, outer_radius, 5)
	p.circle(0, 0, radius = radii, fill_color = None, line_color = "black", line_width = 1)
	p.text(0, radii + 15, [str(ax_label) for ax_label in labels],
			text_font_size = "14px", text_align = "center", text_baseline = "middle")


	angle = (2.0*np.pi - 2.*alpha)/(len(df) + 1)
	angles = np.pi/2 - alpha - angle/2 - df.index.to_series()*angle

	counter=0
	# for col_name in df.columns[df.columns.str.contains(type)]:
	for sample_name in sample_names:
		df['log2_'+type+'_'+sample_name] = start_radius + np.log2(df[sample_name+"."+type].replace(0., np.nan))*np.floor((outer_radius-start_radius)/max_value)
		df['real_log_'+type+'_'+sample_name] = np.log2(df[sample_name+"."+type].replace(0., np.nan))

		df['start_angle_'+sample_name] = -angle + angles + (4.5-counter*0.5)*angle/6.
		df['end_angle_'+sample_name] = -angle + angles + (6.+counter*0.5)*angle/6.
		counter+=1

	source = ColumnDataSource(df)

	# small wedges
	counter=0
	for col_name in df.columns[df.columns.str.contains('log2_'+type)]:
		p.annular_wedge(x = 0, y = 0, 
						inner_radius = start_radius, outer_radius = col_name,
						start_angle = 'start_angle_'+col_name.replace('log2_'+type+'_', ""), end_angle = 'end_angle_'+col_name.replace('log2_'+type+'_', ""), 
						source = source, color = sample_color[col_name.replace('log2_'+type+'_', "")+"."+type], name = col_name, alpha = 1-counter*0.2) #"#0d3362"
		counter+=1

	# serotype labels
	xr = 3.8*radii[0]*np.cos(np.array(-angle + angles + 10.5*angle/12.))
	yr = 3.8*radii[0]*np.sin(np.array(-angle + angles + 10.5*angle/12.))

	label_angles = np.array(-angle/2. + angles)
	label_angles[label_angles < -np.pi/2] += np.pi
	p.text(xr, yr, df.Representative, angle = label_angles,
			text_font_size = "11px", text_align = "center", text_baseline = "middle")


	# legend
	x_coord = -380
	y_coord = -400
	for group, color in sample_color.items():
		p.rect(x_coord, y_coord + 6.5, width=30, height=13, color=color)
		p.text(x_coord + 30, y_coord + 6.5, text=[group.replace("."+type, "")], text_font_size="10px", text_align="left", text_baseline="middle")
		# Update coordinates for the next group
		y_coord += 30

	# add a tooltip
	tooltips = [("AAV variant", "@Representative")] + [("log2"+type+"-"+i, "@real_log_"+type+"_"+i) for i in sample_names]
	p.add_tools(bk.models.HoverTool(tooltips = tooltips, mode="mouse", point_policy="follow_mouse", 
									names = list(df.columns[df.columns.str.contains('log2_'+type)]))) 
	
	return p
    
