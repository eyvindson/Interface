#DISPLAYING UNOPTIMIZED MAPS & Tables
#MAP
import requests
from pathlib import Path
from zipfile import ZipFile
import geopandas as gpd
import numpy as np
import plotly.express as px
import plotly.graph_objs as go
from ipywidgets import Button, HBox, VBox,interact
from IPython.display import display, HTML
from IPython.display import clear_output
import ipywidgets as widgets
import numpy
import pandas as pd



RCP = "RCP0" # no climate change
filename = "DATA/DATA.h5" # Test data from Central Finland with 3579 forest stands

scenario = "BDS" 
extension = "test" # some additional info to the saved output 

import wget
import os
import pandas as pd
import sys

module_path = os.path.abspath(os.path.join(''))


class Mapping_table:
    path_to_zip_file = module_path + "/DATA.zip"
    directory_to_extract_to = module_path+"/DATA"

    try:
        import zipfile
        with zipfile.ZipFile(path_to_zip_file, 'r') as zip_ref:
            zip_ref.extractall(directory_to_extract_to)
        print("Zip Extracted")
    except:
        print("Zip already existed")
    
    def __init__(self):
        
        self.geodf = gpd.read_file(module_path+"/Region.geojson") ## Data with geometry
        # shape file is a different CRS,  change to lon/lat GPS co-ordinates
        self.geodf = self.geodf.to_crs("WGS84")
        #dat = pd.read_csv("C:/Users/kyey/Documents/DATA.csv",sep=";") #Simulated data
        self.dat = pd.read_hdf(module_path+"/"+filename)
        
        #Merging databases together - linking spatial information with simulated
        self.geometry = pd.DataFrame(list(set(self.dat['id'])),self.geodf['geometry'])
        self.d = {'geometry':self.geodf['geometry'], 'id1': list(set(self.dat['id']))}
        self.geometry = pd.DataFrame(self.d)
        self.geodf1 = self.dat.merge(self.geometry,left_on="id",right_on="id1")
        self.geometry = self.geodf1['geometry']
        self.geodf1.drop('geometry',axis=1)
        self.geodf1 = gpd.GeoDataFrame(self.geodf1,crs="WGS84",geometry=self.geometry)
    
    def selection_fn(self,trace,points,selector):
        #geodf2 = geodf1[(geodf1['year'] == year) & (geodf1['regime'] == regime)]
        bb = [[self.geodf2.iloc[points.point_inds][col] for col in ['id']]]
        self.sel_stand = [bb[0][0].iloc[i] for i in range(0,len(bb[0][0]))]
        
    def on_display_map(self,**kwargs):
        display(kwargs)
        clear_output()
        self.geodf2 = self.geodf1[(self.geodf1['year'] == kwargs['year']) & (self.geodf1['regime'] == kwargs['regime'])]
        self.geodf2 = self.geodf2.set_index("id")
        fig = px.choropleth_mapbox(self.geodf2,    geojson=self.geodf2.geometry,    locations=self.geodf2.index,    color=kwargs['values'],    center=dict(lat= 62.82633, lon=21.259906),    mapbox_style="open-street-map",opacity =0.4,    zoom=13,)
        fig.update_layout(    height=500,    autosize=True,    margin={"r": 0, "t": 0, "l": 0, "b": 0})#
        ff =go.FigureWidget(fig)
        scatter = ff.data[0]
        scatter.on_selection(self.selection_fn)
        display(VBox([ff]))
        
    def on_display_table(self,**kwargs):
        clear_output()
        if 'sel_stand' in globals():
            t = go.FigureWidget([go.Table(
                header=dict(values=['id','year']+[kwargs['1']+kwargs['2']+kwargs['3']],#variables of interest in map -- could be a cross box (?)
                            fill = dict(color='#C2D4FF'),
                            align = ['left'] * 5),
                cells=dict(values=[self.geodf2[col] for col in ['id']+[kwargs['1']+kwargs['2']+kwargs['3']]],#variables of interest in map -- could be a cross box (?)
                           fill = dict(color='#F5F8FF'),
                           align = ['left'] * 5))])
        display(VBox([t]))
    
    def on_display_table(self,**kwargs):
        display(kwargs)
        clear_output()
        #if 'sel_stand' in globals():
        
        if kwargs['Year of interest'] != "ALL" and kwargs['Regime of interest'] != "ALL":
            geodfv = self.geodf1[(self.geodf1['year'] == kwargs['Year of interest']) & (self.geodf1['regime'] == kwargs['Regime of interest'])]
        elif kwargs['Year of interest'] == "ALL" and kwargs['Regime of interest'] == "ALL":
            geodfv = self.geodf1
        elif kwargs['Year of interest'] == "ALL":
            geodfv = self.geodf1[(self.geodf1['regime'] == kwargs['Regime of interest'])]
        elif kwargs['Regime of interest'] == "ALL":
            geodfv = self.geodf1[(self.geodf1['year'] == kwargs['Year of interest'])]
        else:
            geodfv = self.geodf1
        try:
            geodfv = geodfv[geodfv['id'].isin(self.sel_stand)]
        except:
            geodfv = geodfv
        t = go.FigureWidget([go.Table(
            header=dict(values=['id',"regime","year"]+[kwargs['Column 1'],kwargs['Column 2'],kwargs['Column 3']],#variables of interest in map -- could be a cross box (?)
                        fill = dict(color='#C2D4FF'),
                        align = ['left'] * 5),
            cells=dict(values=[geodfv[col] for col in ['id','regime','year']+[kwargs['Column 1'],kwargs['Column 2'],kwargs['Column 3']]],#variables of interest in map -- could be a cross box (?)
                       fill = dict(color='#F5F8FF'),
                       align = ['left'] * 5))])
        display(VBox([t]))
        '''else:
            print("T")
            if kwargs['Year of interest'] != "ALL" and kwargs['Regime of interest'] != "ALL":
                geodfv = self.geodf1[(self.geodf1['year'] == kwargs['Year of interest']) & (self.geodf1['regime'] == kwargs['Regime of interest'])]
            elif kwargs['Year of interest'] == "ALL":
                geodfv = self.geodf1[(self.geodf1['regime'] == kwargs['Regime of interest'])]
            elif kwargs['Regime of interest'] == "ALL":
                geodfv = self.geodf1[(self.geodf1['year'] == kwargs['Year of interest'])]
            else:
                geodfv = self.geodf1
            
            t = go.FigureWidget([go.Table(
                header=dict(values=['id',"regime","year"]+[kwargs['Column 1'],kwargs['Column 2'],kwargs['Column 3']],#variables of interest in map -- could be a cross box (?)
                            fill = dict(color='#C2D4FF'),
                            align = ['left'] * 5),
                cells=dict(values=[geodfv[col] for col in ['id',"regime",'year']+[kwargs['Column 1'],kwargs['Column 2'],kwargs['Column 3']]],#variables of interest in map -- could be a cross box (?)
                           fill = dict(color='#F5F8FF'),
                           align = ['left'] * 5))])
            display(VBox([t]))'''
    
    def showGUI(self):
        colTypeChooser = widgets.interactive(self.on_display_map,{"manual":True,"manual_name": "Display / Update Map"},**{'values':[l for l in list(self.geodf1.columns)[1:] if l not in ['year','regime','geometry','id1']],'year':list(set(self.geodf1.year)),'regime':list(set(self.geodf1.regime))})
        TablecolTypeChooser = widgets.interactive(self.on_display_table,{"manual":True,"manual_name": "Display / Update Table"},**{'Year of interest':["ALL"]+list(set(self.geodf1.year)),'Regime of interest':["ALL"]+list(set(self.geodf1.regime)),'Column 1':list(self.geodf1.columns)[3:],'Column 2':list(self.geodf1.columns)[3:],'Column 3':list(self.geodf1.columns)[3:]})
        display(HTML('''<style>
        .widget-label { min-width: 60% !important; }
        </style>'''))
        display(colTypeChooser)
        display(HTML('''<style>
        .widget-label { min-width: 60% !important; }
        </style>'''))
        display(TablecolTypeChooser)

MT = Mapping_table()
MT.showGUI()
