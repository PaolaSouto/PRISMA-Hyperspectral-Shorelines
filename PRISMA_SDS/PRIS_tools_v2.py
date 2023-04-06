#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 19:45:00 2021

This module contains functions that will be used in the other modules

@author: paolasouto
"""

import os 
import numpy as np
#import nest_asyncio
import time
import sys

import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from pylab import ginput
import matplotlib as mpl
import matplotlib
mpl.rcParams['toolbar'] = 'None' 

from osgeo import gdal,osr, ogr
import pandas as pd
import geopandas as gpd
import fiona
from fiona.crs import from_epsg
from shapely.geometry import mapping, Polygon

#from geojson import Point, Feature, FeatureCollection, dump
#import geojson
import json
sys.path.append('/media/sf_VBox_Shared/PAOLA/PRISMA3')

sys.path.append('/home/paola/PRISMA_FOAM/')

from PRISMA_SDS import PRIS_tools_v2, PRIS_img_v2

#nest_asyncio.apply()


#fiona.drvsupport.supported_drivers['kml'] = 'rw' # enable KML support which is disabled by default
#fiona.drvsupport.supported_drivers['LIBKML'] = 'rw' # enable KML support which is disabled by default


################################################################################
################# GENERATES THE GENERAL FOLDER STRUCTURE #######################
################################################################################


def generate_structure(path0, scene, beach): # To be modify as new results will be produced
    
    """
    
    This function generates the folder structure per beach
    
    Inputs:
        
        path0 = string where the folder PRISMA is stored
        scene = string with the PRISMA SCENE ANALIZED (e.g Oristano)
        beach = string with the name of the scene beach analyzed (e.g Arborea)
        
    Outputs:
        
        directories = dictionary containing the different paths used in the SDS retrival
    
    """
    
    # Initialize the dictionary
    
    directories = {}
    
    
    # General SCENES folder
    
    path_scenes = os.path.join(path0,'SCENES')
    
    if not os.path.exists(path_scenes):
        raise Exception('The SCENES folder does not exist or the path is wrong')
        
    # Scene of interest folder
    
    path_sceneI = os.path.join(path_scenes, scene)
    
    if not os.path.exists(path_sceneI):
        raise Exception('The scene of Interest does not exist or the path is wrong')
        
    directories['scene'] = os.path.join(path_scenes, scene)
        
    # Beach folder
                        
    path_beach = os.path.join(path_scenes, scene, beach)
    
    if not os.path.exists(path_beach):
        
        os.makedirs(path_beach)        
        
    # VHR images for co-registration
    
    path_vhr = os.path.join(path_scenes, scene, 'VHR')
    
    # if not os.path.exists(path_vhr):
    #     raise Exception('The vhr folder does not exists or the path is wrong')
        
    directories['VHR'] = os.path.join(path_scenes, scene, 'VHR')
    
        
        
    ###################### SHAPEFILES PATHS
    
        # Folder with the BB shapefile used to crop the images
    
    path_CropShp = os.path.join(path_beach, 'QGIS', 'Crop_shapefile')
    
    if not os.path.exists(path_CropShp):
        os.makedirs(path_CropShp)
    
    directories['Crop_shp'] = path_CropShp
    
    
        # Folder where the profiles will be stored
        
    path_pr = os.path.join(path_beach, 'QGIS', 'Profiles')
    
    if not os.path.exists(path_pr):
        os.makedirs(path_pr)
        
    directories['Profiles'] = path_pr
    
        # Folder where the profiles pix coordinates will be stored
        
    path_pr_pix = os.path.join(path_beach, 'QGIS', 'Profiles_pix')
    
    if not os.path.exists(path_pr_pix):
        os.makedirs(path_pr_pix)
        
    directories['Profiles_pix'] = path_pr_pix
    
        # Folder where the baseline will be stored
        
    path_bs = os.path.join(path_beach, 'QGIS', 'Baseline')
    
    if not os.path.exists(path_bs):
        os.makedirs(path_bs)
        
    directories['Baseline'] = path_bs
        
        
        
    ###################### IMAGES FOLDERS
    
    
    # ORIGINAL IMAGE IN TIFF FORMAT
    
    image_type_folder = ['PAN', 'HS']
    
    processing_folders = ['CoReg']

    
    # Define the folder names of the manipulated images for each image "type"
    
    
    for folder_img in image_type_folder:
        
        
        path = os.path.join(directories['scene'], folder_img)
        
        if not os.path.exists(path):
            
            os.makedirs(path)
            
        directories[folder_img] = path
            
        path_coreg = os.path.join(directories[folder_img], 'Coreg')
        
        name = folder_img + '_CoReg'
        
        if not os.path.exists(path_coreg):
            
            os.makedirs(path_coreg)
            
        directories[name] = path_coreg
        

    # Define the types of images contained inside the folder
    
    image_type_folders = ['HS', 'PAN','PanSharp']

    processing_folders = ['Square_Cropped', 'Beach_Cropped', 'RGB']
    
    # Generate the folders that will contain the images
    
    for folder_img in image_type_folders:
        
        path_folder_img = os.path.join(path_beach, folder_img)
        
        if (folder_img == 'PanSharp'):
        
            if not os.path.exists(path_folder_img):
                os.makedirs(path_folder_img)
                #directories[folder_img] = os.path.join(path_beach, folder_img)
            
        # Generate the subfolders inside each type of image folder
        
        for processing_folder in processing_folders:
            
            if folder_img == 'PAN' and processing_folder == 'RGB': # The RGB images are only for the HS images
                
                continue
            
            if folder_img == 'HS' and processing_folder == 'Beach_Cropped':
                
                continue
            
            if folder_img == 'PAN' and processing_folder == 'Beach_Cropped':
                
                continue
            

        
            path_folder_processing = os.path.join(path_folder_img, processing_folder)
            
            name_directory = folder_img + '_' + processing_folder
            
            directories[name_directory] = os.path.join(path_folder_img, processing_folder)
            
            if not os.path.exists(path_folder_processing):
                os.makedirs(path_folder_processing)
                
                
    ###################### EuHydro SHORELINE FOLDER
    
    path_EuHydro = os.path.join(path0, 'EuHydro', 'Shapefile', 'EUHYDRO_Coastline_EEA39_v013.shp')
    
    # if not os.path.exists(path_EuHydro):
    #     raise Exception('The EuHydro path does not exit or the path is wrong')
    
    directories['EuHydro'] = path_EuHydro
                
                
    ###################### RESULTS FOLDERS (¡¡¡¡¡¡¡ THIS WILL CHANGE!!!!!!!)
    
    result_types_folders = ['Img_profiles', 'K_means', 'CoReg_Shifts']
    
    for result_types_folder in result_types_folders:
        
        path_result = os.path.join(path_beach, 'Results', result_types_folder)
        
        directories[result_types_folder] = path_result
        
        if not os.path.exists(path_result):
            
            os.makedirs(path_result)
            
            
    
    ###################### RESULTS FOLDERS (¡¡¡¡¡¡¡ THIS WILL CHANGE!!!!!!!)    
    
    subfolders = ['csv', 'plot']
    result_types_folders = ['Interface_pix', 'PR_pix', 'SS_pix']
    
    for result_types_folder in result_types_folders:
        for subfolder in subfolders:
            path = os.path.join(path_beach, 'Results', result_types_folder, subfolder)
            
            if not os.path.exists(path):
                os.makedirs(path)
                
            name_path = result_types_folder + '_' + subfolder
            directories[name_path] = path
            
            
    ##################### SHORELINES
    
    path_shorelines = os.path.join(path_beach, 'Results', 'SDS')
    
    if not os.path.exists(path_shorelines):
        
        os.makedirs(path_shorelines)
        
    directories['SDS'] = path_shorelines
    
       
            
    return directories


################################################################################
################# DISTANCES BETWEEN SS #######################
################################################################################

def Bray_Curtis(fx, fy):
    
    fx = np.array(fx)
    fy = np.array(fy)
    
    fy = fy.astype('float32')
    fx = fx.astype('float32')
    
    #score = 100 - (np.sum(np.abs(fx - fy))/(np.sum(fx) + np.sum(fy)))*100
    score = 100 - (np.sum(np.abs(fx - fy))/(np.sum(fx) + np.sum(fy)))*100
    
    return score


        
        

################################################################################
################# GENERATES THE GENERAL FOLDER STRUCTURE #######################
################################################################################

def coor_to_geojson(x, y):
    
        
    geoson = {'type':'FeatureCollection', 'features':[]}
    
    for xx, yy in zip(x,y):
        feature = {'type':'Feature',
                   'geometry':{'type':'Point',
                               'coordinates':[]}}
        feature['geometry']['coordinates'] = [xx,yy]
    
        geoson['features'].append(feature)
    
    return geoson


def save_geojson(path, filename, geoson):
    
    output_file = os.path.join(path, filename)
    
    with open(output_file, 'w') as f:
        
        json.dump(geoson, f)
        

     
        

################################################################################
############# FUNCTIONS RELATED WITH FEATURES PROJECTION #######################
################################################################################

def get_GeoTransform(geo, img2, res):
    
    """
    This could be changed after compare with PRISMAREAD [Hay que explicarla bien]
    
    inputs:
        
        * geo = image extension
        * img2 = np.array containing the image
        * res = image resolution
        * crop = Boolean indicating if we have to obtain the GeoTransform for the crop
    
    I think that this could be applied to the original image.
    If we are using the crop, this doesn't apply'
    
    https://github.com/ranghetti/prismaread/blob/master/R/pr_create_pan.R
    
    """

    ex = {'xmin' : geo['xmin'] - res/2,
          'xmax': geo['xmin'] - res/2 + img2.shape[1] * res,
          'ymin': geo['ymin'] - res/2,
          'ymax': geo['ymin'] - res/2    + img2.shape[0] * res}



    # Set the resolution
    
    
    GeoT = (ex['xmin'], res, 0, ex['ymax'], 0, -res)
    
    driver = gdal.GetDriverByName(str("GTiff"))
    

    
    Projj = osr.SpatialReference()
    Projj.ImportFromEPSG(int(geo['proj_epsg'])) #4326
    Projj.ExportToPrettyWkt()

    
    return GeoT, driver, Projj

def world_to_pixel(GeoT, x, y):
    """
    Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate
    the pixel location of a geospatial coordinate
    
    inputs:
        
        * GeoT: gdal geomatrix
        * x,y: point geographic coordinates to be transformed
        
    outputs:
        
        col,row: point image coordinates

    """
    
    ul_x= GeoT[0]
    
    ul_y = GeoT[3]
    
    x_dist = GeoT[1]
    
    y_dist = GeoT[5]
    
    col= int((x - ul_x) / x_dist)
    
    row = -int((ul_y - y) / y_dist)
    
    return col, row

def pixel_to_world(ds, col, row):
    
    """
    
    Transform from pixel coordinates to world coordinates using gdal
    
    inputs:
        
        * ds : image open with gdal
        * col : col coorinate
        * row : row coodinate
        
    ouputs:
        
        * x,y: lon/lat world coordinates 
        
    Based on:
        
        https://stackoverflow.com/questions/59052516/find-lat-long-coordinates-from-pixel-point-in-geotiff-using-python-and-gdal
    
    """
    
    # Find the transform inverse
    
    source = osr.SpatialReference(wkt = ds.GetProjection())
    target = osr.SpatialReference()
    
    # Get the target EPSG
    
    proj = osr.SpatialReference(wkt=ds.GetProjection())
    EPSG = int(proj.GetAttrValue('AUTHORITY',1))
    
    target.ImportFromEPSG(EPSG)
    
    transform = osr.CoordinateTransformation(source, target)
    
    GeoT = ds.GetGeoTransform()
    
    
    
    ul_x = GeoT[0] 
    
    ul_y = GeoT[3] 
    
    x_dist = GeoT[1] 
    
    y_dist = GeoT[5] 
    
    world_x = col * x_dist + ul_x #x pixel
    world_y = row * y_dist + ul_y #y pixel
    
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(world_x, world_y)
    point.Transform(transform)
    
    
    return point.GetX(), point.GetY()


def pixel_to_world_resized(ds,GeoT,col,row):
    
    # Find the transform inverse
    
    source = osr.SpatialReference(wkt = ds.GetProjection())
    target = osr.SpatialReference()
    
    # Get the target EPSG
    
    proj = osr.SpatialReference(wkt=ds.GetProjection())
    EPSG = int(proj.GetAttrValue('AUTHORITY',1))
    
    target.ImportFromEPSG(EPSG)
    
    transform = osr.CoordinateTransformation(source, target)
    
    #GeoT = ds.GetGeoTransform()
    
    
    
    ul_x = GeoT[0] 
    
    ul_y = GeoT[3] 
    
    x_dist = GeoT[1] 
    
    y_dist = GeoT[5] 
    
    world_x = col * x_dist + ul_x #x pixel
    world_y = row * y_dist + ul_y #y pixel
    
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(world_x, world_y)
    point.Transform(transform)
    
    
    return point.GetX(), point.GetY()


def get_img_epsgCode(directories):
    
    """
    
    Obtain the epsg code of the scene images
    
    input:
        
        * directories: dictionary containing the different paths used in the SDS retrival
        
    output:
        
        * EpsgCode_img: int epsg code of the scene images

    """
    
    # Initialize the list that will contail all the epsg codes
    
    epsg_img = []
    
    # Otain the co-registrated PAN images 
    
    pan_images = [f for f in os.listdir(directories['PanSharp_RGB']) if (f.endswith('.tiff') and not f.startswith('_.'))]

    # Obtain the epsg codes of all the PAN images available
    
    for pan_img in pan_images:
        
        path_pan_img = os.path.join(directories['PanSharp_RGB'], pan_img)

        img = gdal.Open(path_pan_img)
        Projj = osr.SpatialReference(wkt=img.GetProjection())
        epsg_img.append(Projj.GetAttrValue('AUTHORITY', 1))

    if len(np.unique(epsg_img)) != 1:

        raise Exception('No all the images have the same epsg code')

    else:

        EpsgCode_img = int(np.unique(epsg_img)[0])

    return EpsgCode_img


################################################################################
######################### LOAD OTHER DATASETS ##################################
################################################################################

### OLD FUNCTION
        
def load_kml(directories):
    
    """
    This function load the beach kml
    
    inputs : dictionary with the paths
    
    outputs:
        
        col_beach:
        row_beach:
        
    """
    
    path_kml_file = directories['Beach_kml']
    
    # Load the beach polygon
    
    ff = gpd.read_file(path_kml_file)
    
    # Transform the ff epsg code to the same of the images
    
    EpsgCode_img = get_img_epsgCode(directories)
    
    # Converts the spatial reference using the epsg codes
    
    ff = ff.to_crs(epsg = EpsgCode_img) # ver si da error
    
    # Obtain the kml geographic coordinates
    
    c = ff['geometry'][0]
    
    shell_coords = np.array(c.exterior)
    
    polygon = []
    
    for x,y in zip(shell_coords[:,:1], shell_coords[:,1:2]):
        
        polygon.append([list(x)[0],list(y)[0]])
        
        
    # Obtain the polygon extent
    
    x_pol = []
    y_pol = []
    
    for n in polygon:
        
        x_pol.append(n[0])
        y_pol.append(n[1])
        
    
    minx = np.min(x_pol)
    maxx = np.max(x_pol)
    miny = np.min(y_pol)
    maxy = np.max(y_pol)
    
    extent = (minx, maxx, miny, maxy)
    
        
    return extent
        
    
    
    
def load_EuHydro(directories):
    
    """
    
    This function loads the EuHydro shoreline.
    
    inputs:
        
        * directories: dictionary with the paths
    
    outputs: 
        
        * EuHydro: geopandas  DataFrame with the EuHydro domains
    
    """
    
    EuHydro = gpd.read_file(directories['EuHydro'])
    
    epsg_code_EuHydro = EuHydro.crs.to_epsg()
    
    # Check if the EuHydro shoreline epsg code is the same one of the images
    
    epsg_img = get_img_epsgCode(directories)
    
    if epsg_code_EuHydro != epsg_img:
        
        EuHydro = EuHydro.to_crs(epsg = epsg_img)
        
    return EuHydro

    
    
    
################################################################################
######################### EXTENT CROPPING ##################################
################################################################################   
    

def modify_cropping_extent(extent):
    
    """
    
    Modify the kml extent in order to mantain the original image resolution
    
    inputs:
        
        * extent: tuple with the kml extension in geographic coordinates
        
    output:
        
        * extent_modify: tuple with the kml geographic coordinates

    """
    
    # Y dimension
    
    hiR_y = round(( extent[3] - extent[2])/5)
    LoR_y = round(( extent[3] - extent[2])/30)
    
    # X dimension
    
    hiR_x = round((extent[1] - extent[0])/5)
    LoR_x = round((extent[1] - extent[0])/30)
    
    # Ratio
    
    rrW_y = hiR_y/LoR_y
    rrW_x = hiR_x/LoR_x
    
    
    if (rrW_x < 6) or (rrW_x > 6):
        
        hiR_x = LoR_x * 6
        hiR_y = LoR_y * 6
        
        # Check if now the ratio is correct
        
        # rrW_y = hiR_y/LoR_y
        # rrW_x = hiR_x/LoR_x
        

        
    minx = extent[0]
    miny = extent[2]
                
    # Obtain the new extension
    
    maxx_new = minx + (hiR_x * 5)
    maxy_new = miny + (hiR_y * 5)
    
    extent_modi = (minx, maxx_new, miny, maxy_new)
    
    return extent_modi


############ GUI TO SELECT THE EXTENT

class Labeler2(object):
    
    """
    GUI to select the beach of interest inside the SCENE image
    
    inputs:
        
        * data = paths to the Co-Registered PAN images
        
    outputs:
        
        * self.extension_modi = extension in geographic coordinates of the BB cotaining
                                the beach of interest
    
    """
    
    def __init__(self, data, directories):
        
        self.i = 0 # image counter
        self.dirs = data
        self.directories = directories
        
        # Genrate the initial figure and axes
        
        self.fig = plt.figure(figsize=(9,9))
        self.ax = self.fig.add_subplot(111)
        
        # Boolena close loop
        
        self.boolean = True
        
        # Generate the initial image
        
        self.img = gdal.Open(self.dirs[self.i])
        self.GeoT = self.img.GetGeoTransform()
        self.img2 =  np.array(self.img.GetRasterBand(1).ReadAsArray())
        self.ax.imshow(self.img2, cmap = 'gray')
        self.ax.axis('off')
        self.ax.set_title(self.dirs[self.i].split('/')[-1])
        self.fig.canvas.draw()
        
        plt.ion()
        plt.show()  

        self.i += 1
        
        # Initialize the variables that will store the coordinates
        
        self.x_points = []
        self.y_points = []
        
        
        ## Generate the buttons
        
        self.visible = True
        self.no_visible = False
        
        ## Buttons
        
        ## Position = posx, posy, width, height
        
        # a) NEXT IMAGE button
        
        self.next_ax = self.fig.add_axes([0.9, 0.78, 0.075, 0.1])
        self.button_next = Button(self.next_ax, 'Next', color='#32CD32')
        self.next_ax.set_visible(self.visible)
        self.button_next.on_clicked(self._next_img)
        
        # b) SELECT IMAGE button
        
        self.done_ax = self.fig.add_axes([0.9, 0.65, 0.075, 0.1])
        self.button_done = Button(self.done_ax, 'Select', color='orange')
        self.done_ax.set_visible(self.visible)
        self.button_done.on_clicked(self._select_img)

        # c) ACCEPT BB button

        self.accept_ax = self.fig.add_axes([0.45, 0.01, 0.075, 0.075])
        self.button_accept = Button(self.accept_ax, 'Accept', color='green')
        self.accept_ax.set_visible(self.no_visible)
        self.button_accept.on_clicked(self._accept_bb)
        
        # d) REJECT BB button
        
        self.reject_ax = self.fig.add_axes([0.55, 0.01, 0.075, 0.075])
        self.button_reject = Button(self.reject_ax, 'Reject', color = 'red')
        self.reject_ax.set_visible(self.no_visible)
        self.button_reject.on_clicked(self._reject_bb)
        

        
        # SET THE CALLBACKS
        
        plt.waitforbuttonpress()
        
        
        self.cid_next = self.fig.canvas.mpl_connect('button_press_event', self._next_img)
        self.cid_select = self.fig.canvas.mpl_connect('button_press_event', self._select_img)
        # self.cid_accept = self.fig.canvas.mpl_connect('button_press_event', self._accept_bb)
        # self.cid_reject = self.fig.canvas.mpl_connect('button_press_event', self._reject_bb)
        
        ### FUNCTIONS
        
    def _next_img(self, event):
        
        self.img = gdal.Open(self.dirs[self.i])
        self.GeoT = self.img.GetGeoTransform()
        self.img2 =  np.array(self.img.GetRasterBand(1).ReadAsArray())
        self.ax.imshow(self.img2, cmap = 'gray')
        self.ax.axis('off')
        self.ax.set_title(self.dirs[self.i].split('/')[-1])
        
        self.i += 1
        self.xlim = self.ax.get_xlim()
        self.ylim = self.ax.get_ylim()
        self.fig.canvas.draw()

        
        
        #plt.waitforbuttonpress()
        
        
    def _onclick(self,event):
        
        # OPCIÓN GINPUT
        
        ### GIUNPUT PROBLEMS WITH BACKEND AND JUPYTER NOTEBOOK
        
        
        # Make visible the second pair of buttons

        self.reject_ax.set_visible(self.visible)
        self.accept_ax.set_visible(self.visible)
 
        
        sys.stdout.flush()
        

        self.ix, self.iy = event.xdata, event.ydata
        
        
        
        # Check each time if the point can be added or not

        
        if (self.ix is None) and (self.iy is None):
            
            self.ix = 0
            self.iy = 0
                        
        if (self.ix > 1) and (self.iy > 1):
            
            
            self.x_points.append(self.ix)
            self.y_points.append(self.iy)
            
            self.plot = self.ax.plot(self.ix, self.iy, 'r*')
            self.fig.canvas.draw()


                
            # Check each time which is the length of the list storing the clicked points
      
                
            if len(self.x_points) > 2:
                
                del self.x_points
                del self.y_points
                
                self.x_points = []
                self.y_points = []
                

                    
            if len(self.x_points) == 2:
                
                
    
                #print(self.x_points)
                
                xx = self.x_points
                yy = self.y_points
                
            
                # If we have clicked 2 valid points, plot them
                
                xmin, xmax = np.array(sorted(self.x_points))
                ymin, ymax = np.array(sorted(self.y_points))
                
    
    
                self.ax.set_xlim(xmin, xmax)
                self.ax.set_ylim(ymax, ymin)
                self.ax.set_title('Selected area')
                self.fig.canvas.draw()
                
                
                #Transform extent from pix coor into geo coordinates
                
                xx = []
                yy = []
                
                for col, row in zip(self.x_points, self.y_points): 
                    
                    x,y = pixel_to_world(self.img, col, row)
                    
                    xx.append(x)
                    yy.append(y)
                    
                # Sort them
                
                xmin_geo, xmax_geo = sorted(xx)
                ymin_geo, ymax_geo = sorted(yy)
                
                self.extent_geo = (xmin_geo, xmax_geo, ymin_geo, ymax_geo)
                
                # Obtain the extent modify

                
                self.extent_modi = modify_cropping_extent(self.extent_geo)
                
                
                
                # Activate the buttons 
                
                #self.cid_accept = self.fig.canvas.mpl_connect('button_press_event', self._accept_bb)
                self.cid_reject = self.fig.canvas.mpl_connect('button_press_event', self._reject_bb)
                


        
        # Si acepto despuès tengo que desconectar tb el on_click

            
        
    def _select_img(self, event):
            
        # Hidde next and select buttons
        
        self.next_ax.set_visible(self.no_visible)
        self.done_ax.set_visible(self.no_visible)
        
          
        # Change image title
        
        self.ax.set_title('Select the beach of interest')
        
        #p Stop the _next_img method
        
        self.fig.canvas.mpl_disconnect(self.cid_next)
        
        self.fig.canvas.mpl_disconnect(self.cid_select)
        
        # Launch the on_click method
        
        self.cid = self.fig.canvas.mpl_connect('key_press_event', self._onclick)
        

        # Now activate the event
        
        self._onclick(event)
        

        
    
    def _reject_bb(self, event):
        
        
        self.ax.clear()
        self.ax.imshow(self.img2, cmap = 'gray')
        self.ax.axis('off')
        self.ax.set_title('Select the beach of interest')
        self.fig.canvas.draw()
        
        
        # Delete the x_points and y_points variables and generate them again
        
        del self.x_points
        del self.y_points
        
        self.x_points = []
        self.y_points = []
        
        
       # Launch the on_click method
        
        self.cid = self.fig.canvas.mpl_connect('key_press_event', self._onclick)
        

        # Now activate the event
        
        self._onclick(event)
        
        
        
    def _accept_bb(self, event):
        
        # Hiddde the buttons
        
        self.reject_ax.set_visible(self.no_visible)
        self.accept_ax.set_visible(self.no_visible) 
        

        # Disconnect the key press event
        
        self.fig.canvas.mpl_disconnect(self.cid)
        
        create_cropping_shapefile(self.directories, self.extent_modi)
        
        ## Now retrieve the BB info
        
        plt.close()
        
        
def define_beach(directories):
    
    
    pan_images = [f for f in os.listdir(directories['PAN_CoReg']) if f.endswith('.tiff')]
    data = []
    for i in range(len(pan_images)):
        data.append(os.path.join(directories['PAN_CoReg'], pan_images[i]))
        
    labeler = Labeler2(data, directories)
    
    return labeler
        
    
def create_cropping_shapefile(directories ,extent_modi):
    
    """
    
    Generates a temporary folder with the shapefile used to crop the images
    
    """
    
    ###### OLD METHOD
    
    # Load the original kml file
    
    #extent = load_kml(directories)
    
    # Obtain the modify extent
    
    #extent_modi = modify_cropping_extent(extent)
    
    #### NEW METHOD CALLING THE GUI
    

    
    
    # Define the new polygon
    
    poly = Polygon([(extent_modi[0], extent_modi[2]), 
                (extent_modi[0], extent_modi[3]), 
                (extent_modi[1], extent_modi[3]), 
                (extent_modi[1], extent_modi[2]),
                (extent_modi[0], extent_modi[2])])
    
    # Save in a temporary folder that will be delete one the images will be cropped
    
    polygon_name = 'Square_polygon.shp'

    # Get the images epsg code
    
    EpsgCode_img = get_img_epsgCode(directories)
    
    # Define a polygon feature geometry with one attribute
    
    Schema = {
        'geometry' : 'Polygon',
        'properties': {'id':'int'},
        }
    
    # Write the shapefile
    
    with fiona.open(os.path.join(directories['Crop_shp'], polygon_name), 'w', crs=from_epsg(EpsgCode_img) , driver = 'ESRI Shapefile', schema = Schema ) as c:
    
    ## If there are multiple geometries, put the "for" loop here
        c.write({
            'geometry': mapping(poly),
            'properties': {'id': 1},
        })
        

        
def crop_shp_geo(directories):
    
    path_poly = os.path.join(directories['Crop_shp'], 'Square_polygon.shp')

    ds = ogr.Open(path_poly)
    poly = ds.GetLayer(0)
    
    crs = poly.GetSpatialRef()
        
    geo = dict()

    geo['xmin'] = poly.GetExtent()[0]
    geo['xmax'] = poly.GetExtent()[1]
    geo['ymin'] = poly.GetExtent()[2]
    geo['ymax'] = poly.GetExtent()[3]
    geo['proj_epsg'] = int(crs.GetAttrValue('AUTHORITY',1))
    
    return geo


def beach_crop_shapefile(directories):
    
    """
    
    This function allows to crop the PanSharp images to the profiles extent
    
    inputs:
        
        * directories: dictionary with the paths
    
    """
    
    # Load the profiles
    
    profiles_names = [f for f in os.listdir(directories['Profiles']) if f.endswith('_pr.shp') and '90_0' in f]
    
    for profile_name in profiles_names:
    
    
        dir_profile = os.path.join(directories['Profiles'], profile_name)
        
        profiles = gpd.read_file(dir_profile, encoding="utf-8")
                                 
                                 
        ## Obtain the coordinates of the area defined by the profiles
        
        minx, miny, maxx, maxy = profiles.geometry.total_bounds
        
        pto_x_start = []
        pto_x_end = []
        pto_y_start = []
        pto_y_end = []
        
        for line in profiles.iterrows():
            
            coords = line[1]['geometry'].coords.xy
            pto_x_start.append(coords[0][0])
            pto_x_end.append(coords[0][1]) # coord x final
            pto_y_start.append(coords[1][0])
            pto_y_end.append(coords[1][1])
        
            
            
        ## Define the geometry
        
        polygon_geom = Polygon([(pto_x_start[0], pto_y_start[0]), 
                         (pto_x_start[-1], pto_y_start[-1]), 
                         (pto_x_end[-1], pto_y_end[-1]), 
                         (pto_x_end[0], pto_y_end[0]),
                         (pto_x_start[0], pto_y_start[0])])
        
        
        #polygon_geom = Polygon([(minx, miny),
                               # (minx, maxy),
                                #(maxx, maxy),
                                #(maxx, miny),
                                #(minx, miny)])
        
        ## Save the polygon 
        
        crs = get_img_epsgCode(directories)
        polygon = gpd.GeoDataFrame(index=[0], crs=crs, geometry=[polygon_geom]) 
        
        name_polygon = 'Beach_polygon_' + profile_name
        
        dir_polygon = os.path.join(directories['Crop_shp'], name_polygon)
        
        polygon.to_file(filename= dir_polygon, driver="ESRI Shapefile")
    
    
    
    
#################################################################
################ PRISMA RUNNING ON SAET #########################
#################################################################

def create_Tiff(output, img2, driver, GeoT, Projj, band_num):
    
    """
    
    This function saves the Panchromatic band of the PRISMA he5 images
    
    inputs:
        
        output = string with the path where the image will be stored
        img2 = numpy array; image to be stored as tiff
        driver = type of image; in this case always "Tiff"
        GeoT = parameters defininig the image geolocation
        Projj = projection
        
    outputs:
        
        Tiff image
        
    """
    
    driver = gdal.GetDriverByName(str('GTiff'))
    

    rows = img2.shape[0]
    cols = img2.shape[1]
    
    
    
    DataSet = driver.Create(output, cols, rows, 1, gdal.GDT_Float32)
    DataSet.SetGeoTransform(GeoT)
    DataSet.SetProjection(Projj.ExportToWkt())
    
    band = img2[:,:,band_num]
    DataSet.GetRasterBand(1).WriteArray(band)
    DataSet.FlushCache()



def GenerateSAETBandsCombination(path0_prisma, path_saet_data, scene, beach, index, sat):
    
    path_save = path_saet_data
    
    path_hs = os.path.join(path0_prisma, 'SCENES', scene, beach, 'HS_Georef') 
    file = [f for f in os.listdir(path_hs) if 'hs_cropped_georef' in f and (f.endswith('tiff') or f.endswith('tif'))]
    date = file[0].split('_')[1]
    
        
    # Open the tiff file
    
    dir_file = os.path.join(path_hs, file[0])
    
    ds = gdal.Open(dir_file)
    GeoT = ds.GetGeoTransform()
    driver = gdal.GetDriverByName(str("GTiff"))
    
    proj = osr.SpatialReference(wkt=ds.GetProjection())
    crs =  proj.GetAttrValue('AUTHORITY',1)
    
    Projj = osr.SpatialReference()
    Projj.ImportFromEPSG(int(crs)) #4326
    Projj.ExportToPrettyWkt()
        
    arr = PRIS_img_v2.numpy_from_tiff(dir_file, False)
    
    CW = pd.read_csv(os.path.join(path0_prisma, 'cw.txt')).to_numpy().flatten()
    
    if index == 'MNDWI':
        
        SWIR1_range = np.arange(120,140)
        Green_range = np.arange(16, 26)
        
        for green in Green_range:
            for swir in SWIR1_range:
                
                green_band = arr[:,:,green]
                swir1_band = arr[:,:,swir]
                
                if sat == 'l8':
                
                    name_folder = scene + '_L1TP_192033_' + date + '_G' + str(green) + '_S' + str(swir) + '_02_T1'
                    
                    name_band_green = name_folder +'_B03.TIF'
                    name_band_swir = name_folder +'_B06.TIF'
                    
                if sat == 's2':
                    
                    name_folder = scene + '_MSIL1C_' + date + '_G' + str(green) + '_S' + str(swir) + 'T32UMG' 
                    
                    name_band_green = name_folder +'_B03.TIF'
                    name_band_swir = name_folder +'_B11.TIF'
                
                path_folder = os.path.join(path_save, name_folder)
                
                if not os.path.exists(path_folder):
                    os.makedirs(path_folder)
                    
                

                
                output_green = os.path.join(path_folder, name_band_green)
                output_swir = os.path.join(path_folder, name_band_swir)
        
                create_Tiff(output_green, arr, driver, GeoT, Projj, green)
                create_Tiff(output_swir, arr, driver, GeoT, Projj, swir)
    
    # if index == 'AWEISH':
        
    # if index == 'AWEINSH':
        
    
    
def FromSaetShapefile2Txt(path0_prisma,dir_prisma_sds, dir_saet_sds, scene, beach, sat, index):
    
    dir_store = os.path.join(path0_prisma, 'SCENES_open')
    paths = PRIS_tools_v2.generate_structure(dir_store, scene, beach)
    epsg_code = PRIS_tools_v2.get_img_epsgCode(paths)


    sds_folders = [f for f in os.listdir(dir_saet_sds)]

    for folder in sds_folders:
        
        path_folder = os.path.join(dir_saet_sds, folder)
        
        if len(sds_folders) != 0:
        
        # find the sds shapefile

            sds_file = [f for f in os.listdir(os.path.join(dir_saet_sds, folder)) if 'shp' in f and ('points' in f or 'cp' in f)]
    
            print(beach, folder)
    
         # load the sds shapefile
            
            dir_shp = os.path.join(dir_saet_sds, folder, sds_file[0])
            
            sds = gpd.read_file(dir_shp)
            sds = sds.to_crs(epsg_code)
            
            ss = [[a.x, a.y] for a in sds.geometry.values]
            
            X, Y = [], []
            
            for s in ss:
                X.append(s[0])
                Y.append(s[1])
    
            
            if sat == 'l8':
                
                filename = beach + '_' + '_'.join(folder.split('_')[3:6]) + '_' + sat  + '_' + index + '_saet.txt'
                
            if sat == 's2':
                
                filename = beach + '_' + '_'.join(folder.split('_')[3:6]).split('T')[0] + '_' + sat  + '_' + index + '_saet.txt'
            
            dir_fileout = os.path.join(dir_prisma_sds, filename)
            
    
            fileout = open(dir_fileout, "w")
        
            for lo, la in zip(X, Y):
                
                # sub_col = (lo - upper_left[0])/geoT[1]
                # sub_row = (upper_left[1] - la)/geoT[1]
                fileout.write('{:9.3f} {:9.3f}\n'.format(lo, la))
                
            fileout.close()               
                    
    
    print('**************************************************************')
    
    
    # sds_saet = 
        
    

def create_Tiff(output, img2, driver, GeoT, Projj, band_num):
    
    """
    
    This function saves the Panchromatic band of the PRISMA he5 images
    
    inputs:
        
        output = string with the path where the image will be stored
        img2 = numpy array; image to be stored as tiff
        driver = type of image; in this case always "Tiff"
        GeoT = parameters defininig the image geolocation
        Projj = projection
        
    outputs:
        
        Tiff image
        
    """
    
    driver = gdal.GetDriverByName(str('GTiff'))
    

    rows = img2.shape[0]
    cols = img2.shape[1]
    
    
    
    DataSet = driver.Create(output, cols, rows, 1, gdal.GDT_Float32)
    DataSet.SetGeoTransform(GeoT)
    DataSet.SetProjection(Projj.ExportToWkt())
    
    band = img2[:,:,band_num]
    DataSet.GetRasterBand(1).WriteArray(band)
    DataSet.FlushCache()



def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array-value)).argmin()
    
    return array[idx], idx

def Generate_S2_Img(arr, CW):
    
    red_min, red_max =650, 680
    xx, red_min_idd =  find_nearest(CW, red_min)
    xx, red_max_idd = find_nearest(CW, red_max)
    mean_red = np.mean(arr[:,:,red_min_idd:red_max_idd], axis = 2)
    
    blue_min, blue_max = 458, 523
    xx, blue_min_idd = find_nearest(CW, blue_min)
    xx, blue_max_idd = find_nearest(CW, blue_max)
    mean_blue = np.mean(arr[:,:,blue_min_idd:blue_max_idd], axis=2)
    
    green_min, green_max = 543, 578
    xx, green_min_idd = find_nearest(CW, green_min)
    xx, green_max_idd = find_nearest(CW, green_max)
    mean_green = np.mean(arr[:,:,green_min_idd:green_max_idd], axis = 2)
    
    nir_min, nir_max = 785, 899
    xx, nir_min_idd = find_nearest(CW, nir_min)
    xx, nir_max_idd = find_nearest(CW, nir_max)
    mean_nir = np.mean(arr[:,:,nir_min_idd:nir_max_idd], axis = 2)
    
    swir1_min, swir1_max = 1565, 1655
    xx, swir1_min_idd = find_nearest(CW, swir1_min)
    xx, swir1_max_idd = find_nearest(CW, swir1_max)
    mean_swir1 = np.mean(arr[:,:,swir1_min_idd:swir1_max_idd], axis= 2)
    
    swir2_min, swir2_max = 2100, 2280
    xx, swir2_min_idd = find_nearest(CW, swir2_min)
    xx, swir2_max_idd = find_nearest(CW, swir2_max)
    mean_swir2 = np.mean(arr[:,:,swir2_min_idd:swir2_max_idd], axis= 2)
    
    ## Seleccionamos y promediamos las bandas
    
    PRISMA_S2 = []
    
    PRISMA_S2.append(mean_red)
    PRISMA_S2.append(mean_green)    
    PRISMA_S2.append(mean_blue)
    PRISMA_S2.append(mean_nir)
    PRISMA_S2.append(mean_swir1)
    PRISMA_S2.append(mean_swir2)
    
    
    PRISMA_S2 = np.stack(PRISMA_S2)
    PRISMA_S2 = np.transpose(PRISMA_S2, [1,2,0])
    
    return PRISMA_S2


def Generate_L8_Img(Prism_arr, CW):
    
    red_min, red_max = 636, 673
    xx, red_min_idd = find_nearest(CW, red_min)
    xx, red_max_idd= find_nearest(CW, red_max)
    mean_red = np.mean(Prism_arr[:,:,red_min_idd:red_max_idd], axis=2)
    
    blue_min, blue_max = 452, 512
    xx, blue_min_idd = find_nearest(CW, blue_min)
    xx, blue_max_idd = find_nearest(CW, blue_max)
    mean_blue = np.mean(Prism_arr[:,:,blue_min_idd:blue_max_idd], axis=2)
    
    green_min, green_max = 533, 590
    xx, green_min_idd = find_nearest(CW, green_min)
    xx, green_max_idd = find_nearest(CW, green_max)
    mean_green = np.mean(Prism_arr[:,:,green_min_idd:green_max_idd], axis = 2)
    
    nir_min, nir_max = 851, 880
    xx, nir_min_idd = find_nearest(CW, nir_min)
    xx, nir_max_idd = find_nearest(CW, nir_max)
    mean_nir = np.mean(Prism_arr[:,:,nir_min_idd:nir_max_idd], axis = 2)
    
    swir1_min, swir1_max = 1566, 1651
    xx, swir1_min_idd = find_nearest(CW, swir1_min)
    xx, swir1_max_idd = find_nearest(CW, swir1_max)
    mean_swir1 = np.mean(Prism_arr[:,:,swir1_min_idd:swir1_max_idd], axis= 2)
    
    swir2_min, swir2_max = 2107, 2294
    xx, swir2_min_idd = find_nearest(CW, swir2_min)
    xx, swir2_max_idd = find_nearest(CW, swir2_max)
    mean_swir2 = np.mean(Prism_arr[:,:,swir2_min_idd:swir2_max_idd], axis= 2)
    
    ## Seleccionamos y promediamos las bandas
    
    PRISMA_L8 = []
    
    PRISMA_L8.append(mean_red)
    PRISMA_L8.append(mean_green)    
    PRISMA_L8.append(mean_blue)
    PRISMA_L8.append(mean_nir)
    PRISMA_L8.append(mean_swir1)
    PRISMA_L8.append(mean_swir2)
    
    
    PRISMA_L8 = np.stack(PRISMA_L8)
    PRISMA_L8 = np.transpose(PRISMA_L8, [1,2,0])
    
    return PRISMA_L8

def GenerateSAETBandsPRspDeg(path0_prisma, path_saet_data, scene, beach, index, sat):
    
    path_save = path_saet_data
    
    path_hs = os.path.join(path0_prisma, 'SCENES', scene, beach, 'HS_Georef')
    file = [f for f in os.listdir(path_hs) if 'cropped_georef' in f and (f.endswith('tiff') or f.endswith('tif'))]
    date = file[0].split('_')[1]
    
    # Open the tiff files
    
    dir_file = os.path.join(path_hs, file[0])
    
    ds = gdal.Open(dir_file)
    GeoT = ds.GetGeoTransform()
    driver = gdal.GetDriverByName(str('GTiff'))
    
    proj = osr.SpatialReference(wkt=ds.GetProjection())
    crs = proj.GetAttrValue('AUTHORITY', 1)
    
    Projj = osr.SpatialReference()
    Projj.ImportFromEPSG(int(crs))
    Projj.ExportToPrettyWkt()
    
    arr = PRIS_img_v2.numpy_from_tiff(dir_file, False)
    CW = pd.read_csv(os.path.join(path0_prisma, 'cw.txt')).to_numpy().flatten()
        
    if index == 'AWEINSH':
        
        if sat == 'l8':
            
            arr_mean = Generate_L8_Img(arr, CW)
    
            # Involved bands: B3(green), B5(nir), B6(swir1), B7(swir2)
            name_folder = scene + 'L1TP_192033_' + date + '_X' + 'xxx' + '_' + index + '_02_T1'
            
            green_band = arr_mean[:,:,1]
            name_band_green = name_folder + '_B3.TIF'
            
            nir_band = arr_mean[:,:,3]
            name_band_nir = name_folder + '_B5.TIF'
            
            swir1_band = arr_mean[:,:,4]
            name_band_swir1 = name_folder + '_B6.TIF'
            
            swir2_band = arr_mean[:,:,5]
            name_band_swir2 = name_folder + '_7.TIF'

        if sat == 's2':
            
            arr_mean = Generate_S2_Img(arr, CW)
            
            #Involved bands: B03(green),B08(nir), B11(swir1), B12(swir2).
            
            name_folder = scene + '_MSIL1C_' + date + '_X' + 'xxx' + '_' + index + '_T32UMG'
            
            green_band = arr_mean[:,:,1]
            name_band_green = name_folder + '_B03.TIF'
            
            nir_band = arr_mean[:,:,3]
            name_band_nir = name_folder + '_B08.TIF'
            
            swir1_band = arr_mean[:,:,4]
            name_band_swir1 = name_folder + '_B11.TIF'
            
            swir2_band = arr_mean[:,:,5]
            name_band_swir2 = name_folder + '_B12.TIF'
            
                
        path_folder = os.path.join(path_save, name_folder)
        
        if not os.path.exists(path_folder):
            os.makedirs(path_folder)
            
        output_green = os.path.join(path_folder, name_band_green)
        output_nir = os.path.join(path_folder, name_band_nir)
        output_swir1 = os.path.join(path_folder, name_band_swir1)
        output_swir2 = os.path.join(path_folder, name_band_swir2)
        
        create_Tiff(output_green, arr_mean, driver, GeoT, Projj, 1)
        create_Tiff(output_nir, arr_mean, driver, GeoT, Projj, 3)
        create_Tiff(output_swir1, arr_mean, driver, GeoT, Projj, 4)
        create_Tiff(output_swir2, arr_mean, driver, GeoT, Projj, 5)
            
    if index == 'AWEISH':
        
        if sat == 'l8':
            
            #B2(blue), B3(green), B5(nir), B6(swir1), B7(swir2)
            
            arr_mean = Generate_L8_Img(arr, CW)
            
            name_folder = scene + 'L1TP_192033_' + date + '_X' + 'xxx' + '_' + index + '_02_T1'
            
            blue_band = arr_mean[:,:,2]
            name_band_blue = name_folder + '_B2.TIF'
                        
            green_band = arr_mean[:,:,1]
            name_band_green = name_folder + '_B3.TIF'
            
            nir_band = arr_mean[:,:,3]
            name_band_nir = name_folder + '_B5.TIF'
            
            swir1_band = arr_mean[:,:,4]
            name_band_swir1 = name_folder + '_B6.TIF'
            
            swir2_band = arr_mean[:,:,5]
            name_band_swir2 = name_folder + '_7.TIF'
            
        if sat == 's2':
            
             #Involved bands: B02(blue), B03(green), B08(nir), B11(swir1), B12(swir2).
            arr_mean = Generate_S2_Img(arr, CW)
                        
            name_folder = scene + '_MSIL1C_' + date + '_X' + 'xxx' + '_' + index + '_T32UMG'
            
            blue_band = arr_mean[:,:,2]
            name_band_blue = name_folder + '_B02.TIF'
            
            green_band = arr_mean[:,:,1]
            name_band_green = name_folder + '_B03.TIF'
            
            nir_band = arr_mean[:,:,3]
            name_band_nir = name_folder + '_B08.TIF'
            
            swir1_band = arr_mean[:,:,4]
            name_band_swir1 = name_folder + '_B11.TIF'
            
            swir2_band = arr_mean[:,:,5]
            name_band_swir2 = name_folder + '_B12.TIF'
            
        path_folder = os.path.join(path_save, name_folder)
        
        if not os.path.exists(path_folder):
            os.makedirs(path_folder)
            
        output_blue = os.path.join(path_folder, name_band_blue)
        output_green = os.path.join(path_folder, name_band_green)
        output_nir = os.path.join(path_folder, name_band_nir)
        output_swir1 = os.path.join(path_folder, name_band_swir1)
        output_swir2 = os.path.join(path_folder, name_band_swir2)
        
        create_Tiff(output_blue, arr_mean, driver, GeoT, Projj, 2)
        create_Tiff(output_green, arr_mean, driver, GeoT, Projj, 1)
        create_Tiff(output_nir, arr_mean, driver, GeoT, Projj, 3)
        create_Tiff(output_swir1, arr_mean, driver, GeoT, Projj, 4)
        create_Tiff(output_swir2, arr_mean, driver, GeoT, Projj, 5)
        
    if index == 'MNDWI':
        
        if sat == 'l8':
            
            # B3(green), B6(swir1)
            
            arr_mean = Generate_L8_Img(arr, CW)
            
            name_folder = scene + 'L1TP_192033_' + date + '_X' + 'xxx' + '_' + index + '_02_T1'
            
            green_band = arr_mean[:,:,1]
            name_band_green = name_folder + '_B3.TIF'
            
            swir1_band = arr_mean[:,:,4]
            name_band_swir1 = name_folder + '_B6.TIF'
            
        if sat == 's2':
            
            arr_mean = Generate_S2_Img(arr, CW)
                        
            name_folder = scene + '_MSIL1C_' + date + '_X' + 'xxx' + '_' + index + '_T32UMG'
            
            #Involved bands: B03(green), B11(swir1).
                    
    if index == 'MNDWI':
        
        if sat == 'l8':
            #Involved bands: B3(green), B6(swir1)
            
            arr_mean = Generate_L8_Img(arr, CW)
            
            name_folder = scene + 'L1TP_192033_' + date + '_X' + 'xxx' + '_X' + 'xxx' + '_02_T1'
            
            green_band = arr_mean[:,:,1]
            name_band_green = name_folder + '_B3.TIF'
            
            
            swir1_band = arr_mean[:,:,4]
            name_band_swir1 = name_folder + '_B6.TIF'
            
            
        if sat == 's2':
            
            #Involved bands: B03(green), B11(swir1).
            
            arr_mean = Generate_S2_Img(arr, CW)
                        
            name_folder = scene + '_MSIL1C_' + date + '_X' + 'xxx' + '_X' + 'xxx' + 'T32UMG'

            green_band = arr_mean[:,:,1]
            name_band_green = name_folder + '_B03.TIF'
            
            swir1_band = arr_mean[:,:,4]
            name_band_swir1 = name_folder + '_B11.TIF'
            
            
            
    path_folder = os.path.join(path_save, name_folder)
    
    if not os.path.exists(path_folder):
        os.makedirs(path_folder)
            
    
    output_green = os.path.join(path_folder, name_band_green)
    output_swir1 = os.path.join(path_folder, name_band_swir1)

    create_Tiff(output_green, arr_mean, driver, GeoT, Projj, 1)
    create_Tiff(output_swir1, arr_mean, driver, GeoT, Projj, 4)
           
           
           
            
        

        
