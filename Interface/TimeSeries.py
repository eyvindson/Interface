import os
import numpy as np
import matplotlib.pyplot as plt
import imageio
import plotly.express as px
from IPython.display import Image

def time_series_map(geo1,limit = [0,0],variable = "V",FPS = 2):
    if limit == [0,0]:
        limit = [geo1["Harvested_V"].min(),geo1[variable].max()]
    filenames = []
    cmap = (px.colors.sequential.Plasma)

    colors = [cmap[k] for k in range(0,10)]
    cc_scale = ([[i/10, colors[i]] for i in range(0,9)] + [[1, colors[9]]])

    im = []
    for i in [2016+i for i in range(5,105,5)]:

        # create file name and append it to a list
        filename = f'{i}.png'
        filenames.append(filename)

        geodf2 = geo1[(geo1['year'] == i) ]
        fig = px.choropleth_mapbox(geodf2.set_index("id"),    geojson=geodf2.geometry,    locations=geodf2.index,    color=variable,    center=dict(lat= 62.82633, lon=21.259906),    mapbox_style="open-street-map",opacity =0.4,    zoom=12,title="Year "+str(i),range_color=limit,)#,color_continuous_scale=cc_scale)

        # save frame
        fig.write_image(filename)


    # build gif
    with imageio.get_writer(variable+'.gif', mode='I',fps=FPS) as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)

    # Remove files
    for filename in set(filenames):
        os.remove(filename)
    display(Image(data=open(variable+".gif",'rb').read(), format='gif'))
