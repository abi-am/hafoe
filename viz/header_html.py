import os
import pandas as pd

def human_readable_size(size, decimal_places=1):
    #convert byte to higher unit
    for unit in ['B','KB','MB','GB','TB']:
        if size < 1024.0:
            break
        size /= 1024.0
    return f"{round(size, 1):.{decimal_places}f}{unit}"

def get_dir_size(path='.'):
    # in bytes
    total = 0
    with os.scandir(path) as it:
        for entry in it:
            if entry.is_file():
                total += entry.stat().st_size
            elif entry.is_dir():
                total += get_dir_size(entry.path)
    return total

def get_version():
    version_path = os.path.abspath(os.path.join(os.curdir, "version.txt"))
    f = open(version_path, 'r')
    version = f.read()

    return version

def write_header(config_df, chimeric_orf_summary_df, enriched1_orf_summary_df):
    html_file_path = os.path.join(config_df.loc[config_df[0] == "output.dir"].iloc[0,1], "reports", config_df.loc[config_df[0] == "title_of_the_run"].iloc[0,1] + ".html")
    html_file = open(html_file_path, "w")

    html_file.write(get_html_body(config_df, chimeric_orf_summary_df, enriched1_orf_summary_df))

    html_file.close()


def get_div_logo(w=200, h=110):
    div = """<div>
    <img src="logo.png" alt="hafoe logo" width="%d" height="%d">
    </div>
    &nbsp;
    &nbsp;""" % (w, h)
    return div


def get_div_title(title):
    div = """
    <h2 style="color:#404040;" align="center">%s</h2>
    """ % title
    return div


def get_div_footer():
    div = """<hr><div>
    These plots are generated with <it>hafoe</it> version %s. 
    For more information on usage and citation, please visit <a href = "https://github.com/TatevikJ/hafoe">
    our github page</a>.
    </div>""" % get_version()
    return div


def get_libsize_table(orf_summary_df):
    table = "<table>" # \
            # "<tr>" \
            # "<th> </th>" \
            # "<th> </th>" \
            # "</tr>"

    for i in range(len(orf_summary_df)):
        table += "<tr><td>%s:</td><td>%d</td></tr>" \
                    % (orf_summary_df.iloc[i,0], pd.to_numeric(orf_summary_df.iloc[i,1]))
    table += "</table>"

    return table


def get_html_body(config_df, chimeric_orf_summary_df, enriched1_orf_summary_df):

    curr_dir = os.getcwd()
    with open(os.path.join(curr_dir, "viz", "template.html"), "r") as f:
        template = f.read()

        if ('title_of_the_run' in config_df[0].values) and config_df.loc[config_df[0] == "title_of_the_run"].iloc[0, 1]  is not None and config_df.loc[config_df[0] == "title_of_the_run"].iloc[0, 1] != 'none':
            title = config_df.loc[config_df[0] == "title_of_the_run"].iloc[0, 1]
        else:
            title = "None"
        
        if ('parent' in config_df[0].values) and config_df.loc[config_df[0] == "parent"].iloc[0, 1]  is not None and config_df.loc[config_df[0] == "parent"].iloc[0, 1] != 'none':
            parent = os.path.abspath(config_df.loc[config_df[0] == "parent"].iloc[0, 1])
            parent_size = human_readable_size(os.path.getsize(parent))
        else:
            parent = "None"
            parent_size = "None"
        
        if ('chimera' in config_df[0].values) and config_df.loc[config_df[0] == "chimera"].iloc[0, 1]  is not None and config_df.loc[config_df[0] == "chimera"].iloc[0, 1] != 'none':
            chimera = os.path.abspath(config_df.loc[config_df[0] == "chimera"].iloc[0, 1])
            chimera_size = human_readable_size(os.path.getsize(chimera))
        else:
            chimera = "None"
            chimera_size = "None"

        if ('enriched' in config_df[0].values) and config_df.loc[config_df[0] == "enriched"].iloc[0, 1]  is not None and config_df.loc[config_df[0] == "enriched"].iloc[0, 1] != 'none':
            enriched = os.path.abspath(config_df.loc[config_df[0] == "enriched"].iloc[0, 1])
            enriched_size = human_readable_size(os.path.getsize(enriched))
        else:
            enriched = "None"
            enriched_size = "None"

        if ('explore.out' in config_df[0].values) and config_df.loc[config_df[0] == "explore.out"].iloc[0, 1]  is not None and config_df.loc[config_df[0] == "explore.out"].iloc[0, 1] != 'none':
            explore_out = os.path.abspath(config_df.loc[config_df[0] == "explore.out"].iloc[0, 1])
            explore_out_size = human_readable_size(get_dir_size(explore_out))
        else:
            explore_out = "None"
            explore_out_size = "None"

        if ('output.dir' in config_df[0].values) and config_df.loc[config_df[0] == "output.dir"].iloc[0, 1]  is not None and config_df.loc[config_df[0] == "output.dir"].iloc[0, 1] != 'none':
            out_dir = config_df.loc[config_df[0] == "output.dir"].iloc[0, 1]
        else:
            out_dir = "None"

        if ('explore' in config_df[0].values) and config_df.loc[config_df[0] == "explore"].iloc[0, 1]  is not None and config_df.loc[config_df[0] == "explore"].iloc[0, 1] != 'none':
            explore = config_df.loc[config_df[0] == "explore"].iloc[0, 1]
        else:
            explore = "None"
        
        if ('identify' in config_df[0].values) and config_df.loc[config_df[0] == "identify"].iloc[0, 1]  is not None and config_df.loc[config_df[0] == "identify"].iloc[0, 1] != 'none':
            identify = config_df.loc[config_df[0] == "identify"].iloc[0, 1]
        else:
            identify = "None"

        if ('read_length' in config_df[0].values) and config_df.loc[config_df[0] == "read_length"].iloc[0, 1]  is not None and config_df.loc[config_df[0] == "read_length"].iloc[0, 1] != 'none':
            read_length = config_df.loc[config_df[0] == "read_length"].iloc[0, 1]
        else:
            read_length = "None"

        if ('overlap' in config_df[0].values) and config_df.loc[config_df[0] == "overlap"].iloc[0, 1]  is not None and config_df.loc[config_df[0] == "overlap"].iloc[0, 1] != 'none':
            overlap = config_df.loc[config_df[0] == "overlap"].iloc[0, 1]
        else:
            overlap = "None"

        if ('step_size' in config_df[0].values) and config_df.loc[config_df[0] == "step_size"].iloc[0, 1]  is not None and config_df.loc[config_df[0] == "step_size"].iloc[0, 1] != 'none':
            step_size = config_df.loc[config_df[0] == "step_size"].iloc[0, 1]
        else:
            step_size = "None"

        if ('vd_criterion' in config_df[0].values) and config_df.loc[config_df[0] == "vd_criterion"].iloc[0, 1]  is not None and config_df.loc[config_df[0] == "vd_criterion"].iloc[0, 1] != 'none':
            vd_criterion = config_df.loc[config_df[0] == "vd_criterion"].iloc[0, 1]
        else:
            vd_criterion = "None"

        result = template.format(svg=get_div_logo(),
                                 version=get_version(),
                                 title=get_div_title(title),
                                 footer=get_div_footer(),
                                 fa=parent,
                                 fa_size=parent_size,
                                 csv_fq=chimera,
                                 csv_fq_size=chimera_size,
                                 fq=enriched,
                                 fq_size=enriched_size,
                                 explore_out=explore_out,
                                 explore_out_size=explore_out_size,
                                 out_dir=out_dir,
                                 explore=explore,
                                 identify=identify,
                                 read_length=read_length,
                                 step_size=step_size,
                                 overlap=overlap,
                                 vd_criterion=vd_criterion,
                                 chimeric_orf_summary=get_libsize_table(chimeric_orf_summary_df),
                                 enriched_orf_summary=get_libsize_table(enriched1_orf_summary_df),
                                 main_report=os.path.abspath(os.path.join(config_df.loc[config_df[0] == "output.dir"].iloc[0,1], "reports", title + "_main.html")))
        return (result)

