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

def write_header(config_df, chimeric_orf_summary_df_path, enriched1_orf_summary_path_paths):
    html_file_path = os.path.join(config_df.loc[config_df[0] == "output.dir"].iloc[0,1], "reports", "main", config_df.loc[config_df[0] == "title_of_the_run"].iloc[0,1] + ".html")
    html_file = open(html_file_path, "w")

    html_file.write(get_html_body(config_df, chimeric_orf_summary_df_path, enriched1_orf_summary_path_paths))

    html_file.close()


def get_div_logo(w=280, h=130):
    div = """<div>
    <svg  #width="%d" height="%d" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1" id="Layer_1" x="0px" y="0px" viewBox="0 0 337.5 224.2" style="enable-background:new 0 0 337.5 224.2;" xml:space="preserve">
    <style type="text/css">
        .st0{fill:none;stroke:#000000;stroke-width:2;stroke-linecap:round;stroke-linejoin:round;stroke-miterlimit:10;}
        .st1{opacity:0.95;}
        .st2{fill:#FFFFFF;}
        .st3{fill:#9AB8DA;stroke:#034EA2;stroke-miterlimit:10;}
        .st4{fill:none;stroke:#034EA2;stroke-miterlimit:10;}
        .st5{fill:#CEE9D1;stroke:#84C98D;stroke-miterlimit:10;}
        .st6{fill:none;stroke:#84C98D;stroke-miterlimit:10;}
        .st7{fill:#DDE9F2;stroke:#8FB6D5;stroke-miterlimit:10;}
        .st8{fill:none;stroke:#8FB6D5;stroke-miterlimit:10;}
        .st9{fill:none;stroke:#0C1525;stroke-width:0.75;stroke-miterlimit:10;}
    </style>
    <g>
        <g>
            <path class="st0" d="M295.5,135.2c2.2,2.2,4.4,4.5,6.8,6.5c2.9,2.6,4.8,1.9,5.6-1.9c1.1-5.6,2-11.3,3-16.9c1.3-7,5.1-13,7.9-19.5    c3.3-7.6,3.9-14.4-2.2-21.1c-3.6-4-4.9-9.2-4.7-14.9c0.3-8-1.6-15.6-9.3-19.8c-5-2.7-10.5-4.7-16.5-2.7"/>
            <path class="st0" d="M272.6,191.6c-5.4,0.3-10.8,1.2-16.1,2.2c-2.6,0.4-5.3,0.5-7.9,0.7c-8,0.4-16,1.1-24,1.8c-15,1.3-30,1-45,1.7    c-15.6,0.8-31.2-0.5-46.9-0.6c-8.5,0-17-0.5-25.5-1.4c-8.6-1-17.3-1.6-25.9-2.7c-8-1-15.9-2.7-24-3.2"/>
            <path class="st0" d="M298.5,70.5c-1.8-3.3-1.6-7.1-3-10.5c-4.1-9.9-11.4-16.7-20.6-21.8c-10.7-5.8-22.2-8.5-34.1-8.9    c-7.2-0.2-14.4,1.4-21.4,3.3c-1.8,0.5-3.1,0.2-4.9-0.2c-6.8-1.6-13.6-1.7-19.8,2.4c-3.8,2.6-8.2,3.9-12.7,4.1    c-12.4,0.4-18.6,8.3-22.8,18.5c-2.9,6.9-3.7,14.2-3.9,21.7c-0.4,12.2,1.4,24.1,3.8,36c1.4,6.8,6.6,8.3,11.6,5.2    c4.3-2.7,8.5-5.7,12.8-8.6"/>
            <path class="st0" d="M216.4,33.7c-8.1,8.4-14.3,18-17.6,29.3c-1.9,6.5-2.6,13.2-3,19.9c-0.5,8.5-2.8,16.4-8.4,22.8    c-3.4,3.9-3.5,8.2-3.7,12.8c-0.7,19.5,10.8,30.8,28.9,35.6c9.6,2.6,19.4,3.8,29.3,3.3c11.9-0.5,23.7-1.3,34.9-5.9    c8.1-3.3,14.4-8.8,18.1-16.3c2.9-5.7,3.9-12.5,2.7-19c-2.2-11.8-2.6-23.5,0.1-35.3c0.4-1.9,0.1-4,0-6c-0.1-2.1,0.3-3.7,2.2-5.4    c3.1-2.8,4.1-6.8,3-10.7c-1-3.8-4-5.7-7.8-6.9c-1.1-0.3-1.8,0.1-2.7,0.3"/>
            <path class="st0" d="M154.1,84.7c-6.9-0.8-13.7-2.1-20.6-2.1c-15.8-0.2-31.1,1.8-43.5,13.3c-8.6,8-12.7,18.2-15.4,29.3    c-1.2,4.8-1.7,9.8-2.5,15.1c-5.8-0.2-11.5-0.1-17,1.5c-7,2-13.5,4.4-15.4,12.7c-2.7,11.8,2,18.1,12,21.3c7.8,2.5,15.9,3.1,24,3.4    c12,0.4,24-0.4,36-0.8c5.3-0.2,10.8-0.6,15.7-2.4c4.1-1.5,7.6-4.8,9.1-9.5c0-0.1,0-0.2,0-0.4"/>
            <path class="st0" d="M85.5,98.8c-4.7-1.6-9.4-0.4-13.8,0.7c-14.8,3.7-26.7,12.2-36.4,23.9c-3.2,3.8-6,7.8-8.7,12    c-1.5,2.4-1.1,4.2,0.9,5.6c1.8,1.2,3.8,0.6,5.2-1.1c5.8-7.2,13-12.8,21-17.2c6.8-3.7,14-6.6,22.1-5.7"/>
            <path class="st0" d="M129.8,175.5c2.5,1,4.9,0.7,7.5,0.4c4.9-0.6,9.7-1,14.6-1.1c6.1-0.2,12.3,0.1,18.4-0.1c5.6-0.2,11,1,16.5,1.7    c4,0.6,8.1,1.4,12,1.9c10.3,1.4,20.7,2.9,31.1,2.8c4,0,8-0.7,11.6-2.3c5.3-2.4,7.7-10.5,5.3-16.4c-0.5-1.3-1.2-2.5-1.9-3.8"/>
            <path class="st0" d="M288.8,145.1c3.9,4.8,5.3,10.2,3.7,16.1c-1.5,5.5-5.8,8.2-11.2,9.4c-5.6,1.2-11.2,1.1-16.9,0.9    c-5.5-0.2-10.8-1.1-16.1-2"/>
            <path class="st0" d="M66,141c1.7,6.5,3.7,13,3.3,19.9c-0.3,5.9-1.8,11.5-3,17.3"/>
            <path class="st0" d="M259.5,120.7c0,1.2,0,2.5,0,3.8c0,4.3-3.1,8.3-7.1,9c-4.6,0.8-9.2-1.7-10.9-6c-0.4-1.1-1.1-2.1-0.8-3.4"/>
            <path class="st0" d="M214.9,94.8c1.2,2,1.4,4.5,3.7,6.1c5.6,3.6,11.9,0.7,12.8-6.1"/>
            <path class="st0" d="M272.6,94.8c0.5,3.5,3.4,7.1,6.4,7.5c4.1,0.5,7.4-0.8,9-4.8c0.3-0.9,0.7-1.8,1.1-2.6"/>
            <path class="st0" d="M275.6,124.1c-0.7,2-0.3,4.2-1.8,6c-3.6,4.3-7.7,4.6-11.7,0.7c-0.6-0.6-1.3-1.2-1.9-1.9"/>
            <path class="st0" d="M253.7,114.3c-0.1,3.6,1.6,5.4,5.4,5.6c3.4,0.2,6-1.6,6.4-4.5c0.3-2.3-1.7-3.5-3-3.6c-2.5-0.2-5.3-1.2-7.5,1"/>
            <path class="st0" d="M265.9,155.6c0.2,3.7-2,6.9-2.4,10.5c-0.2,1.6-0.5,3.2-0.6,4.9"/>
            <path class="st0" d="M219.4,120.3c-4.2,0.8-8.5,1.4-12.8,0"/>
            <path class="st0" d="M219.4,111.3c-4.1-0.6-8.1-0.5-12,0.8"/>
            <path class="st0" d="M217.9,168.3c-0.1,3.9-1.7,7.7-1.5,11.6"/>
            <path class="st0" d="M292.5,164.2c3.2,1.7,6.8,2.2,10.1,3.4"/>
        </g>
    </g>
    <g class="st1">
        <g>
            <path d="M51.8,77.8c-2,1.4-4,2.9-5.9,4.4c-3-3.7-6.1-7.3-9.3-10.8c-1.3-1.5-2.6-2.3-3.8-2.6c-1.2-0.3-2.4,0.1-3.6,1    s-1.7,2.1-1.7,3.5s0.6,2.7,1.8,3.9c3.3,3.5,6.4,7,9.5,10.6c-1.9,1.5-3.7,3.1-5.6,4.7c-8-8.9-16.5-17.5-25.7-25.6    c1.9-1.8,3.9-3.5,5.9-5.2c3.7,3.4,7.4,7,10.9,10.5c0,0,0.1,0,0.1-0.1c-0.3-1.7-0.2-3.2,0.4-4.6c0.5-1.4,1.5-2.6,2.9-3.7    c2.4-1.9,4.8-2.7,7.2-2.3c2.4,0.4,4.8,1.9,7.1,4.5C45.3,69.8,48.6,73.8,51.8,77.8z"/>
            <path d="M75.5,62.6c-2.1,1.2-4.2,2.4-6.3,3.7c-0.7-1-1.4-2-2.1-3c0,0-0.1,0-0.1,0.1c0.1,1.6-0.2,3-0.9,4.3    c-0.7,1.3-1.7,2.4-3.1,3.2c-1.9,1.2-4,1.8-6.1,1.6c-2.1-0.2-3.9-1.2-5.4-3s-2.1-3.8-1.7-5.8c0.4-2,1.9-4.2,4.5-6.4    c1.7-1.4,3.3-2.8,5.1-4.2c-0.1-0.1-0.2-0.3-0.3-0.4c-0.8-1-1.8-1.6-3-1.6c-1.2-0.1-2.4,0.3-3.7,1.2c-1,0.6-1.9,1.4-2.8,2.5    c-0.9,1-1.6,2.2-2.2,3.6c-1.7-0.6-3.4-1.2-5.2-1.7c0.8-1.8,1.7-3.4,2.9-4.8c1.2-1.4,2.7-2.7,4.5-3.9c3.4-2.2,6.5-3.1,9.3-2.7    c2.8,0.4,5.3,2,7.5,4.9C69.5,54.3,72.6,58.4,75.5,62.6z M63.1,58c-0.4-0.5-0.7-0.9-1.1-1.4c-1.3,1.1-2.6,2.1-3.9,3.2    c-1.1,0.9-1.7,1.8-1.9,2.6c-0.2,0.9,0,1.7,0.6,2.5c0.6,0.7,1.3,1.1,2.2,1.3c0.9,0.1,1.9-0.1,2.8-0.8c1.4-0.9,2.2-2,2.4-3.4    C64.4,60.7,64.1,59.3,63.1,58z"/>
            <path d="M78,25.7c-0.6,0.1-1.2,0.1-1.6,0.3c-0.4,0.1-0.9,0.3-1.4,0.5c-1.1,0.6-1.8,1.3-1.9,2.1c-0.2,0.8,0.1,1.7,0.9,2.7    c0.6,0.7,1.1,1.5,1.7,2.2c1.9-1,3.9-2,5.8-2.9c1.1,1.5,2.2,3,3.3,4.6c-2,0.9-3.9,1.9-5.8,2.9c4,5.5,7.8,11.1,11.4,16.7    c-2.3,1.1-4.7,2.3-7,3.5c-3.7-5.5-7.6-11-11.7-16.4c-1.4,0.8-2.7,1.5-4.1,2.3c-1.1-1.5-2.3-2.9-3.5-4.4c1.4-0.8,2.8-1.6,4.1-2.3    c-0.7-0.8-1.3-1.7-2-2.5c-1.7-2.2-2.3-4.4-1.5-6.7c0.8-2.3,2.6-4.3,5.6-5.8c1.1-0.6,2.1-1,2.9-1.3s1.8-0.5,2.8-0.7    C76.8,22.1,77.4,23.9,78,25.7z"/>
            <path d="M112.7,45.7c-4.3,1.7-8.2,2-11.8,1c-3.7-1-6.6-3.2-8.9-6.7c-2.4-3.6-3.1-7-1.9-10.1c1.1-3.2,4-5.8,8.5-7.6    c4.4-1.8,8.7-2.2,12.6-1c3.9,1.1,6.9,3.5,9,7.1c2.2,3.8,2.5,7.2,1.1,10.3C119.9,41.7,117,44,112.7,45.7z M110,40.5    c2-0.8,3.2-2,3.6-3.5s-0.1-3.5-1.5-5.8c-1.3-2.2-2.9-3.6-4.6-4.3c-1.8-0.7-3.7-0.6-5.8,0.3c-2,0.8-3.2,2.1-3.6,3.7    c-0.4,1.6,0.2,3.5,1.5,5.6c1.4,2.2,3,3.6,4.7,4.3C106,41.5,107.9,41.4,110,40.5z"/>
            <path d="M152.8,21.5c-6.4,1.5-12.8,3.1-19,4.9c1.1,1.8,2.5,2.9,4.3,3.4c1.7,0.5,3.7,0.5,6-0.1c0.8-0.2,1.9-0.7,3.1-1.5    c1.2-0.8,2.4-1.7,3.3-2.9c1.9,1.1,3.8,2.2,5.6,3.3c-1.4,1.7-3,3-4.7,4.1c-1.7,1-3.8,1.9-6.2,2.5c-4.5,1.2-8.4,1.1-11.9-0.3    c-3.5-1.4-6.3-3.9-8.4-7.6s-2.6-7.1-1.3-10.2c1.3-3.1,4-5.4,8.2-6.6c4.1-1.2,7.9-1.1,11.4,0.3s6.3,4,8.2,7.8    C151.9,19.6,152.4,20.6,152.8,21.5z M142.7,19.1c-1-1.9-2.2-3.1-3.5-3.7c-1.4-0.6-2.9-0.6-4.7-0.1c-1.5,0.4-2.6,1.3-3.1,2.6    c-0.6,1.3-0.5,2.7,0.2,4.4C135.2,21.1,138.9,20,142.7,19.1z"/>
        </g>
    </g>
    <g>
        <g>
            <path class="st2" d="M296.1,159.1c-2.4,1.7-3.4,4.4-4.7,6.9c-1.3,2.7-2.6,5.3-4,8c-0.6,1.2-0.4,2.8,0.9,3.4    c1.1,0.6,2.8,0.4,3.4-0.9c1.2-2.5,2.5-5,3.7-7.5c0.6-1.2,1.2-2.4,1.8-3.6c0.2-0.5,0.5-0.9,0.8-1.4c0.2-0.2,0.2-0.2,0,0    c0.1-0.1,0.2-0.2,0.3-0.3c0.1-0.1,0.7-0.6,0.3-0.3c1.1-0.8,1.7-2.2,0.9-3.4C298.9,158.9,297.3,158.3,296.1,159.1L296.1,159.1z"/>
        </g>
    </g>
    <g>
        <g>
            <path class="st2" d="M305.2,171.6c-0.1-0.9-0.2-1.9-0.3-2.8c-0.1-1.1-0.2-2.1-0.8-3.1c-0.6-1.1-1.8-2-3.1-1.9    c-0.9,0-1.6,0.4-2.3,0.9c-0.3,0.2-0.5,0.5-0.7,0.8c0.2-0.1,0.3-0.3,0.5-0.4c-0.1,0.1-0.2,0.1-0.2,0.2c0.2-0.1,0.4-0.2,0.6-0.3    c-0.1,0-0.2,0.1-0.4,0.1c0.2,0,0.4-0.1,0.7-0.1c-0.1,0-0.2,0-0.3,0c0.2,0,0.4,0.1,0.7,0.1c-0.2,0-0.3-0.1-0.5-0.1    c0.2,0.1,0.4,0.2,0.6,0.3c-0.1-0.1-0.3-0.2-0.4-0.3c0.2,0.1,0.3,0.3,0.5,0.4c-0.1-0.1-0.2-0.2-0.3-0.3c0.1,0.2,0.3,0.3,0.4,0.5    c-0.1-0.1-0.2-0.3-0.2-0.4c0.1,0.2,0.2,0.4,0.3,0.6c0-0.1-0.1-0.2-0.1-0.3c0,0.2,0.1,0.4,0.1,0.7c0-0.1,0-0.3,0-0.4    c0,0.2-0.1,0.4-0.1,0.7c0-0.2,0.1-0.4,0.1-0.5c-0.1,0.2-0.2,0.4-0.3,0.6c0.1-0.2,0.2-0.4,0.4-0.6c-0.1,0.2-0.3,0.3-0.4,0.5    c0.2-0.2,0.3-0.3,0.5-0.5c-0.2,0.1-0.3,0.3-0.5,0.4c0.2-0.2,0.5-0.3,0.7-0.4c-0.2,0.1-0.4,0.2-0.6,0.3c0.3-0.1,0.6-0.2,1-0.3    c-0.2,0-0.4,0.1-0.7,0.1c0.8-0.1,1.6,0,2.4,0c0.6,0.1,1.4-0.3,1.8-0.7s0.8-1.1,0.7-1.8c0-0.6-0.2-1.3-0.7-1.8s-1.1-0.7-1.8-0.7    c-1.7-0.1-3.3-0.1-4.8,0.7c-1.6,0.8-2.7,2.5-2.7,4.2c0,1.9,1.3,3.5,3,4c0.9,0.3,2,0.2,2.8-0.3c0.4-0.2,0.7-0.5,1-0.8    c0.2-0.2,0.3-0.4,0.5-0.6c-0.2,0.1-0.3,0.3-0.5,0.4c0.1,0,0.1-0.1,0.2-0.1c-0.2,0.1-0.4,0.2-0.6,0.3c0.1,0,0.1,0,0.2,0    c-0.2,0-0.4,0.1-0.7,0.1c0.1,0,0.1,0,0.2,0c-0.2,0-0.4-0.1-0.7-0.1c0.1,0,0.2,0,0.2,0.1c-0.2-0.1-0.4-0.2-0.6-0.3    c0.1,0,0.1,0.1,0.2,0.1c-0.2-0.1-0.3-0.3-0.5-0.4c0.1,0.1,0.1,0.1,0.2,0.2c-0.1-0.2-0.3-0.3-0.4-0.5c0.1,0.2,0.2,0.3,0.3,0.5    c-0.1-0.2-0.2-0.4-0.3-0.6c0.1,0.3,0.2,0.6,0.3,0.9c0-0.2-0.1-0.4-0.1-0.7c0.1,1.1,0.2,2.3,0.4,3.4c0.1,0.7,0.2,1.3,0.7,1.8    c0.4,0.4,1.1,0.8,1.8,0.7c0.6,0,1.3-0.2,1.8-0.7C304.9,172.9,305.3,172.3,305.2,171.6L305.2,171.6z"/>
        </g>
    </g>
    <g>
        <g>
            <path class="st2" d="M56.8,193.5c14.3-1,27.8,4.8,41.9,5.6c7.4,0.4,14.6-1.1,22-1.7c3.8-0.3,7.5-0.2,11.2,0.2    c3.6,0.4,7.2,0.9,10.9,0.9c7.3-0.1,14.6-1.7,21.7-3.1c3.4-0.7,6.9-1.4,10.3-1.9c3.6-0.5,7.3-0.8,10.8,0c3.4,0.8,6.5,2.3,9.8,3.2    c3.4,0.9,7,1.4,10.6,1.6c7.3,0.4,14.7-0.6,21.8-1.9c8.3-1.5,16.5-3.1,24.7-4.7c4.1-0.8,8.2-1.8,12.4-2.3c2.8-0.3,7-0.1,8,3.1    c0.8-1.1,1.6-2.1,2.4-3.2c-20.8,2.4-41.6,4.5-62.4,6.5c-5.2,0.5-10.3,1-15.5,1.3c-5.1,0.3-10.3,0.1-15.4-0.3    c-10.5-0.8-20.8-1.5-31.3-1.4c-11.8,0.1-23.5,0.9-35.2,1.8c-1.3,0.1-2.5,1.1-2.5,2.5c0,1.3,1.1,2.6,2.5,2.5    c20.9-1.7,41.7-2.5,62.6-0.8c5.1,0.4,10.2,0.8,15.4,0.7s10.3-0.6,15.5-1.1c10.5-1,20.9-2,31.3-3c11.7-1.2,23.3-2.5,35-3.8    c1.5-0.2,2.9-1.5,2.4-3.2c-1-3.2-3.6-5.4-6.8-6.4c-3.8-1.2-7.8-0.3-11.5,0.5c-7.7,1.5-15.5,3-23.2,4.5c-14.1,2.7-29,6-43,1    c-3.5-1.3-6.9-2.5-10.6-2.7c-3.8-0.2-7.6,0.3-11.3,0.9c-7.1,1.1-14.1,3-21.3,4c-3.7,0.5-7.4,0.7-11.1,0.4c-3.8-0.3-7.6-1-11.4-1.2    c-7.7-0.4-15.3,1.5-23,1.9c-16.2,0.8-31.6-6.7-47.9-5.6c-1.3,0.1-2.5,1.1-2.5,2.5C54.3,192.3,55.5,193.6,56.8,193.5L56.8,193.5z"/>
        </g>
    </g>
    <g>
        <g>
            <path class="st2" d="M200.8,193c-6.2,0.6-12.4,1.4-18.7,1.7c-6.1,0.3-12.1-0.4-18-1.9c-3.1-0.8-4.4,4-1.3,4.8    c6.1,1.6,12.5,2.4,18.8,2.1c6.4-0.2,12.8-1.1,19.2-1.7c1.3-0.1,2.5-1.1,2.5-2.5C203.3,194.2,202.1,192.8,200.8,193L200.8,193z"/>
        </g>
    </g>
    <g>
        <g>
            <path class="st2" d="M127.4,199.6c-12.6-5.6-25.9-11.4-40-10.7c-2.2,0.1-3.2,2.6-1.8,4.3c6.3,7.2,16.4,8.8,25.4,8.4    c5.3-0.3,10.6-1,15.9-1.6c5.3-0.6,10.5-1.2,15.8-1.7c5.2-0.5,10.3-0.7,15.4,0.3c5,1,9.8,1.6,14.9,1.2c9.9-0.8,19.3-4.2,29.1-5.6    c11.6-1.7,23.5,0.2,35,2.1c1.3,0.2,2.7-0.3,3.1-1.7c0.3-1.2-0.4-2.9-1.7-3.1c-10.3-1.7-20.8-3.4-31.2-2.8    c-10,0.6-19.4,3.8-29.2,5.5c-4.5,0.8-9.3,1.1-13.9,0.4c-2.5-0.4-5-1-7.5-1.3c-2.6-0.4-5.2-0.5-7.8-0.4c-5.3,0.2-10.5,0.9-15.7,1.5    s-10.4,1.2-15.5,1.7c-9.7,1-21.4,1.8-28.5-6.3c-0.6,1.4-1.2,2.8-1.8,4.3c13.2-0.7,25.7,4.7,37.5,10c1.2,0.5,2.7,0.4,3.4-0.9    C128.9,202,128.7,200.2,127.4,199.6L127.4,199.6z"/>
        </g>
    </g>
    <g id="Layer_2_00000178889659175815906840000006661522956699967401_">
    </g>
    <g>
        <g>
            <path class="st3" d="M127.1,93.6l-8,4.5l-7.4,4.1l0.3,17.6l15.3,8.4l15.3-8.6v-17.5C142.6,102.1,141.7,101.6,127.1,93.6    C126.8,93.4,127.4,93.8,127.1,93.6"/>
            <polygon class="st4" points="111.8,102.4 142.5,102.1 142.5,119.6 127.1,122.5 111.8,119.8   "/>
        </g>
        <polyline class="st4" points="130,96.5 135.6,102.1 127.4,122.5 127.4,128.2 127.4,122.5 135.6,102.1 127.2,93.6 119.3,102.3    127.1,122.4 142.5,119.5 135.3,102  "/>
        <line class="st4" x1="111.8" y1="119.8" x2="119.3" y2="102.3"/>
    </g>
    <g>
        <g>
            <path class="st5" d="M148.7,129.3l-8,4.5l-7.4,4.1l0.3,17.6l15.3,8.4l15.3-8.6v-17.5C164.2,137.8,163.3,137.3,148.7,129.3    C148.4,129.1,149,129.5,148.7,129.3"/>
            <polygon class="st6" points="133.4,138.1 164.1,137.8 164.1,155.3 148.7,158.2 133.4,155.5   "/>
        </g>
        <polyline class="st6" points="151.6,132.2 157.1,137.8 148.9,158.2 148.9,163.9 148.9,158.2 157.1,137.8 148.8,129.3 140.9,138    148.7,158.1 164.1,155.2 156.9,137.7  "/>
        <line class="st6" x1="133.4" y1="155.5" x2="140.9" y2="138"/>
    </g>
    <g>
        <g>
            <path class="st7" d="M105.3,129.3l-8,4.5l-7.4,4.1l0.3,17.6l15.3,8.4l15.3-8.6v-17.5C120.8,137.8,119.9,137.3,105.3,129.3    C105,129.1,105.6,129.5,105.3,129.3"/>
            <polygon class="st8" points="90,138.1 120.7,137.8 120.7,155.3 105.3,158.2 90,155.5   "/>
        </g>
        <polyline class="st8" points="108.2,132.2 113.8,137.8 105.6,158.2 105.6,163.9 105.6,158.2 113.8,137.8 105.4,129.3 97.5,138    105.3,158.1 120.7,155.2 113.5,137.7  "/>
        <line class="st8" x1="90" y1="155.5" x2="97.5" y2="138"/>
    </g>
    <g>
        <path class="st9" d="M118.8,113.8c-0.2,0.4-0.5,1.2-0.8,1.2s-0.6-0.8-0.6-1.3c0-0.6,0.2-1.5,0.5-1.5c0.2,0,0.4,0.1,0.7,0.7"/>
        <path class="st9" d="M135.5,113.9c0.2,0.4,0.5,1.2,0.9,1.2c0.3,0,0.6-0.8,0.6-1.3c0-0.6-0.2-1.5-0.6-1.5c-0.2,0-0.4,0.1-0.8,0.7"/>
        <path class="st9" d="M118.8,113.9c5.6,0,11.1,0,16.7,0"/>
        <path class="st9" d="M118.7,113c0.5,0,1,0,1.5,0"/>
        <path class="st9" d="M134,113c0.5,0,1,0,1.5,0"/>
    </g>
    <g>
        <path class="st9" d="M97.2,149.7c-0.2,0.4-0.5,1.2-0.8,1.2s-0.6-0.8-0.6-1.3c0-0.6,0.2-1.5,0.5-1.5c0.2,0,0.4,0.1,0.7,0.7"/>
        <path class="st9" d="M113.8,149.8c0.2,0.4,0.5,1.2,0.9,1.2c0.3,0,0.6-0.8,0.6-1.3c0-0.6-0.2-1.5-0.6-1.5c-0.2,0-0.4,0.1-0.8,0.7"/>
        <path class="st9" d="M97.2,149.8c5.6,0,11.1,0,16.7,0"/>
        <path class="st9" d="M97.1,148.9c0.5,0,1,0,1.5,0"/>
        <path class="st9" d="M112.4,148.9c0.5,0,1,0,1.5,0"/>
    </g>
    <g>
        <path class="st9" d="M140.2,149.7c-0.2,0.4-0.5,1.2-0.8,1.2s-0.6-0.8-0.6-1.3c0-0.6,0.2-1.5,0.5-1.5c0.2,0,0.4,0.1,0.7,0.7"/>
        <path class="st9" d="M156.8,149.8c0.2,0.4,0.5,1.2,0.9,1.2c0.3,0,0.6-0.8,0.6-1.3c0-0.6-0.2-1.5-0.6-1.5c-0.2,0-0.4,0.1-0.8,0.7"/>
        <path class="st9" d="M140.2,149.8c5.6,0,11.1,0,16.7,0"/>
        <path class="st9" d="M140.1,148.9c0.5,0,1,0,1.5,0"/>
        <path class="st9" d="M155.4,148.9c0.5,0,1,0,1.5,0"/>
    </g>
    </svg>
    </div>""" % (w, h)
    return div

def get_div_title(title):
    div = """
    <h2 style="color:#404040;margin-top: auto; display: block;" align="center">%s</h2>
    """ % title
    return div


def get_div_footer():
    div = """<hr><div>
    These plots are generated with <it>hafoe</it> version %s. 
    For more information on usage and citation, please visit <a href = "https://github.com/abi-am/hafoe">
    our github page</a>.
    </div>""" % get_version()
    return div


def get_chimeric_summary_table(orf_summary_df_path):
    orf_summary_df = pd.read_csv(orf_summary_df_path)

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


def get_enriched_summary_table(enriched1_orf_summary_path_paths):

    tables = ""
    for f in enriched1_orf_summary_path_paths:
        enriched1_name = os.path.abspath(f)
        enriched1_name = os.path.basename(enriched1_name) 
        enriched1_name = enriched1_name.split(".", 1)[0] 
        
        orf_summary_df = pd.read_csv(f)

        table = "<h5>%s</h5> \
                <table>" % (enriched1_name) # \
                # "<tr>" \
                # "<th> </th>" \
                # "<th> </th>" \
                # "</tr>"

        for i in range(len(orf_summary_df)):
            table += "<tr><td>%s:</td><td>%d</td></tr>" \
                        % (orf_summary_df.iloc[i,0], pd.to_numeric(orf_summary_df.iloc[i,1]))
        table += "</table>"
        table += "&nbsp;"
        tables += table
    return tables

def get_html_body(config_df, chimeric_orf_summary_df_path, enriched1_orf_summary_path_paths):

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
            enriched1 = os.path.abspath(config_df.loc[config_df[0] == "enriched"].iloc[0, 1])
            def list_full_paths(directory):
                return [os.path.join(directory, file) for file in os.listdir(directory)]
            
            enriched1_files = list_full_paths(enriched1)

            enriched1_files = [os.path.abspath(s) for s in enriched1_files]
            # enriched_names = [s.split(".", 1)[0] for s in enriched_files]
            enriched1_sizes = [os.path.getsize(s) for s in enriched1_files]
            enriched1_sizes = [human_readable_size(s) for s in enriched1_sizes]
            enriched1_sizes_str = ', '.join([str(elem) for elem in enriched1_sizes])
            # enriched_size = human_readable_size(os.path.getsize(enriched))
        else:
            enriched1 = "None"
            enriched1_sizes = "None"

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
                                 title=title,
                                 title_div=get_div_title(title),
                                 footer=get_div_footer(),
                                 fa=parent,
                                 fa_size=parent_size,
                                 csv_fq=chimera,
                                 csv_fq_size=chimera_size,
                                 fq=enriched1,
                                 fq_size=enriched1_sizes_str,
                                 explore_out=explore_out,
                                 explore_out_size=explore_out_size,
                                 out_dir=out_dir,
                                 explore=explore,
                                 identify=identify,
                                 read_length=read_length,
                                 step_size=step_size,
                                 overlap=overlap,
                                 vd_criterion=vd_criterion,
                                 chimeric_orf_summary=get_chimeric_summary_table(chimeric_orf_summary_df_path),
                                 enriched_orf_summary=get_enriched_summary_table(enriched1_orf_summary_path_paths), 
                                 main_report=os.path.join(title + "_main.html"))
        return (result)

