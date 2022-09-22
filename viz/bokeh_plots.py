import os

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

def bokeh_histogram(title, df, n_bins):
    hist, edges = np.histogram(df['Length'], density=False, bins=n_bins)
    p = figure(title=title, tools='') #, background_fill_color="#fafafa"
    p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],
           fill_color="navy", line_color="white", alpha=0.5)

    p.y_range.start = 0
    #p.legend.location = "center_right"
    #p.legend.background_fill_color = "#fefefe"
    p.xaxis.axis_label = 'Length'
    p.yaxis.axis_label = 'Frequency'
    #p.grid.grid_line_color="white"

    return p
