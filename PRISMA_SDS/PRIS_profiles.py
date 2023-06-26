#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 11:05:23 2021

This module generate the profiles and their related information

@author: paolasouto
"""


import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import pandas as pd
import geopandas as gpd
#from PIL import Image
import osgeo.gdalnumeric as gdn
import json
import random
from matplotlib.colors import hex2color
import math


from osgeo import  osr, ogr, gdal
from shapely.affinity import rotate
from shapely.geometry import LineString, Point
from scipy import interpolate


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 21:01:14 2021

This module contains the functions needed to handle the PRISMA images

@author: paolasouto
"""

import shutil
import pyproj

## Used in the panSharpening

from scipy import signal
from scipy import ndimage
import PIL
import h5py


from PRISMA_SDS import PRIS_tools, PRIS_img




################################################################################
####################### PROFILES DEFINITION  ##################################
################################################################################



##### OPTION I




def generate_profiles_baseline(directories, profiles_car):
    
    """
    
    Generate the profiles from the baseline
    
    inputs:
        
        * img: image numpy array
        
        * directories: dict with the path folders
        
        * baseline : points defining the baseline
        
        * profiles_car: dictionary containing the characteristics of the profiles to be generated
    
    outputs:
        
        * geojson/shapefile containing the profiles
        
        
    The method for interpolate any curve based on the distance is taken from:
        
        https://stackoverflow.com/questions/52014197/how-to-interpolate-a-2d-curve-in-python
    
    """
    
    # Hay que cargar la Baseline
    
    # dir_baseline = os.path.join(directories['Baseline'], 'baseline.txt')
    # baseline = np.loadtxt(dir_baseline)

    
    pan_img_cr = [f for f in os.listdir(directories['PanSharp_RGB']) if f.endswith('.tiff')]
    
    #print(pan_img_cr)
    
    for pan_img in pan_img_cr:
        
        # Load the baseline
        
        baseline_name = [f for f in os.listdir(directories['Baseline']) if f.endswith('.txt') and pan_img.split('_')[2] in f]
        
        if len(baseline_name) == 0:
            continue
        else:
            baseline_name = baseline_name[0]
        
        dir_baseline = os.path.join(directories['Baseline'], baseline_name)
        
        baseline = np.loadtxt(dir_baseline)

        print(pan_img)
        
        img = gdal.Open(os.path.join(directories['PanSharp_RGB'], pan_img))
        
        
        # Extract the baseline pix coordinates and convert into geographic coordinates
        
        x = []; y = [] # Points in pixel coordinates
        
        
        for pt in baseline:
        
            # Transform to world coordinates
            
            X, Y = PRIS_tools.pixel_to_world(img, pt[0], pt[1])
            
            # Baseline
            
            x.append(X)
            y.append(Y)
            
        # Parametric interpolation based on the distance
        
        points = np.array([x,y]).T
        
        # Linear length along the line
        
        distance = np.cumsum( np.sqrt(np.sum( np.diff(points, axis=0)**2, axis=1 )))
        
        distance = np.insert(distance, 0, 0)
            
        num_points = int((int(np.max(distance)) - np.min(distance))/5)
        
            
        alpha = np.linspace(np.min(distance), int(np.max(distance)), num_points)
    
        interpolator =  interpolate.interp1d(distance, points, kind='slinear', axis=0)
        interpolated_points = interpolator(alpha)
        
        x_new =  interpolated_points.T[0]
        y_new = interpolated_points.T[1]
        
        
        # Save the Baseline (as geojson, then will be changed to store it as a shapefile)
        
        geoson = PRIS_tools.coor_to_geojson(x_new, y_new)
        
        ff = pan_img.split('_')
        
        baseline_name = 'Baseline_' + ff[1] + '_' + ff[2] + '.geojson'
        
        output_baseline = os.path.join(directories['Baseline'], baseline_name)
        
        with open(output_baseline, 'w') as f:
            
            json.dump(geoson, f)
            
        
        # For each origin point calculate the profile
        
        # Based on StackOverflow
        
        lines = []
        
        for xx, yy in zip(x_new,y_new):
            
            start = Point(xx, yy)
            length = profiles_car['length']
            angle = profiles_car['orientation']
            
            end = Point(start.x + length, start.y)
            line = LineString([start, end])
            line = rotate(line, angle, origin = start, use_radians = False)
            
            lines.append(line)
            
            
            
        # Generate the Geopandas containing the profiles
        
        EpsgCode_img = PRIS_tools.get_img_epsgCode(directories)
        
        Profiles = gpd.GeoDataFrame(crs = EpsgCode_img, geometry = lines)
        
        Profiles['PR'] = ['PR' + str(i+1) for i in range(len(lines))]
        
        # Save the GeoDataFrame as a shapefile
        
        output_pr = os.path.join(directories['Profiles'], directories['PAN'].split('/')[-2]+ '_' + ff[2] + '_pr.shp')
        
        Profiles.to_file(driver = 'ESRI Shapefile', filename= output_pr)
        


class BaselineGenerator(object):
    
    """
    
    This class allows the generation of the baseline where from which the profiles
    will be generated
    
    inputs:
        
        * directories: dict with the all the pahs contained in the structure
        * profiles_car = dict with the profiles caracteristics
        
    """
    

    def __init__(self, directories): #, profiles_car
        
        # Save the profiles_car dict in the object
        
        #self.profiles_car = profiles_car
        
        # Save the directories in the object
        
        self.directories = directories
        
        # Save the scene
        
        self.scene = directories['scene'].split('/')[-1]
        self.beach = directories['PanSharp_Square_Cropped'].split('/')[-3]
        
        # Obtain the cropped PAN images paths
        
        pan_img_cr = [f for f in os.listdir(directories['PanSharp_RGB']) if f.endswith('.tiff')]
        
        data = []
        
        for img in pan_img_cr:
            
            data.append(os.path.join(directories['PanSharp_RGB'], img))
            
        
        self.i = 0 #image counter
        
        self.dirs = data
        self.img_names = pan_img_cr
        
        # Generate the initial figure and axes
        
        self.fig = plt.figure(figsize=(9,9))
        self.ax = self.fig.add_subplot(111)
        
        # Generate the initial image
        
        
        self.img = gdal.Open(self.dirs[self.i])
        
        self.GeoT = self.img.GetGeoTransform()
        
        bands = [self.img.GetRasterBand(o) for o in range(1, self.img.RasterCount + 1)]
        
        arr = np.array([gdn.BandReadAsArray(band) for band in bands]).astype('float32')
        
        self.img2 = np.transpose(arr, [1, 2, 0])
        
        #self.img2 =  np.array(self.img.GetRasterBand(1).ReadAsArray())
        
        self.ax.imshow(self.img2**(1/2))
        self.ax.axis('off')
        self.ax.set_title(self.dirs[self.i].split('/')[-1])
        self.fig.canvas.draw()
        
        plt.ion()
        plt.show()  

        self.i += 1
        
        ## Generate the buttons
        
        self.visible = True
        self.no_visible = False
        
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
        
        # c) FINISH BUTTON
        
        # self.finish_ax = self.fig.add_axes([0.55, 0.01, 0.075, 0.075])
        # self.button_finish = Button(self.finish_ax, 'Finish', color = 'red')
        # self.finish_ax.set_visible(self.no_visible)
        # self.button_finish.on_clicked(self._finish)
        
        
    def _next_img(self, event):
        
        self.img = gdal.Open(self.dirs[self.i])
        self.GeoT = self.img.GetGeoTransform()
        
        bands = [self.img.GetRasterBand(o) for o in range(1, self.img.RasterCount + 1)]
        
        arr = np.array([gdn.BandReadAsArray(band) for band in bands]).astype('float32')
        
        self.img2 = np.transpose(arr, [1, 2, 0])
        
        #self.img2 =  np.array(self.img.GetRasterBand(1).ReadAsArray())
        self.ax.imshow(self.img2**(1/2))
        self.ax.axis('off')
        self.ax.set_title(self.dirs[self.i].split('/')[-1])
        
        self.i += 1
        self.xlim = self.ax.get_xlim()
        self.ylim = self.ax.get_ylim()
        self.fig.canvas.draw()

        
    def _select_img(self, event):
            
        # Hidde next and select buttons
        
        self.next_ax.set_visible(self.visible)
        self.done_ax.set_visible(self.visible) # no_visible
        
        # Stop the _next_img method
        
        # self.fig.canvas.mpl_disconnect(self.cid_next)
        
        # self.fig.canvas.mpl_disconnect(self.cid_select)
        
        # Make visible the finish button
        
        # self.finish_ax.set_visible(self.visible)
        
          
        # Change image title
        
        self.ax.set_title('Click baseline points, and press ENTER')
    
        
        # Select the baseline with the ginput method
        
        
        self.pts = self.fig.ginput(n = 0, timeout = -1, show_clicks = True)
        
        
        # Save the baseline
        
        
        # Once the points have been saved we generate the profiles and plot it on the image
        
        # Estamos cambiando el nombre de lla baseline salvada
        
        name_baseline =  self.beach + '_' + self.img_names[self.i-1].split('_')[1] + '_bas'+ '.txt'
        
        baseline_dir = os.path.join(self.directories['Baseline'], name_baseline)
        
        np.savetxt(baseline_dir, self.pts, fmt='%.18g', delimiter=' ', newline=os.linesep)
        
        
        
        # Save the baseline also as geojson
        
        x = []; y = [] # Points in pixel coordinates
        
        
        for pt in self.pts:
        
            # Transform to world coordinates
            
            X, Y = PRIS_tools.pixel_to_world(self.img, pt[0], pt[1])
            
            # Baseline
            
            x.append(X)
            y.append(Y)
            

        
        # Parametric interpolation based on the distance
        
        points = np.array([x,y]).T
        
        # Linear length along the line
        
        distance = np.cumsum( np.sqrt(np.sum( np.diff(points, axis=0)**2, axis=1 )))
        
        distance = np.insert(distance, 0, 0)
            
        num_points = int((int(np.max(distance)) - np.min(distance))/5)
        
            
        alpha = np.linspace(np.min(distance), int(np.max(distance)), num_points)
    
        interpolator =  interpolate.interp1d(distance, points, kind='slinear', axis=0)
        interpolated_points = interpolator(alpha)
        
        x_new =  interpolated_points.T[0]
        y_new = interpolated_points.T[1]
        
        
        # Save the Baseline (as geojson, then will be changed to store it as a shapefile)
        
        geoson = PRIS_tools.coor_to_geojson(x_new, y_new)
        
        
        baseline_name = name_baseline.split('.')[0] + '.geojson'
        
        output_baseline = os.path.join(self.directories['Baseline'], baseline_name)
        
        with open(output_baseline, 'w') as f:
            
            json.dump(geoson, f)
        
        
        #return self.pts, self.profiles_car # Esto se podrà mejorar, cargar txt

        
        #generate_profiles_baseline(self.directories, self.pts , self.profiles_car)
        
        
        
        
        
    
    # def _finish(self, event):
        
    #     plt.close()
        
        
        # Save the points into the geojson file
        
        
        



################################################################################
################ PROFILES INTERPOLATION AND PIX COORDINATES ####################
################################################################################


def redistribute_vertices(geom, distance):
    
    """
    
    This function interpolates the profiles
    
    """
    
    if geom.geom_type == 'LineString':
        
        num_vert = int(round(geom.length / distance))
        
        if num_vert == 0:
            num_vert = 1
            
        return LineString(
            [geom.interpolate(float(n) / num_vert, normalized=True)
             for n in range(num_vert + 1)])
    
    elif geom.geom_type == 'MultiLineString':
        parts = [redistribute_vertices(part, distance)
                 for part in geom]
        
        return type(geom)([p for p in parts if not p.is_empty])
    
    else:
        
        raise ValueError('unhandled geometry %s', (geom.geom_type,))
        
        
def linestring_to_points(feature,line):
    
    """
    
    This function transform the LineString features into point features
    
    """
    return {feature:line.coords}


def DataFrameProfileCoor(profiles_id, profiles_number, gdf,  ds, transform):

    point = ogr.Geometry(ogr.wkbPoint)
    
    # Genetate lists to the needed parameters
    
    name_pr = []
    row = []
    col = []
    xx = []
    yy = []
    
    for idd, num in zip(profiles_id, profiles_number):

        for idxs in range(len(gdf['points'][num][idd])):
            
            name_pr.append(idd)
            
            pt = gdf['points'][num][idd][idxs]
            
            x = pt[0]
            y = pt[1]
            
            xx.append(x)
            yy.append(y)
    
            point.AddPoint(x, y)
            
            point.Transform(transform)
            
            lon, lat = PRIS_tools.world_to_pixel(ds.GetGeoTransform(), point.GetX(), point.GetY())
            
            row.append(lat)
            col.append(lon)
            
    transects_info = pd.DataFrame({'pr': name_pr, 'row': row, 'col': col,
                              'X': xx, 'Y': yy})
    
    return transects_info


def plot_profiles(directories, prof_coor_pan):
    
    # Load one PanSharp image
    
    PanSharp_imgs = [f for f in os.listdir(directories['PanSharp_Square_Cropped']) if f.endswith('.tiff')]
    
    path_img = os.path.join(directories['PanSharp_Square_Cropped'], PanSharp_imgs[0])

    _, raster_rgb = PRIS_img.numpy_from_tiff(path_img, True)
    
    # Groupby the pix per profile
    
    groups = prof_coor_pan.groupby('pr')
    
    # Color the profiles
    
    colors = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
         for i in range(len(groups))]
    
    colors_rgb = []
    
    
    for coc in colors:
        colors_rgb.append(hex2color(coc))
    
    o = 0
    
    for idd,group in groups:
        
        group.reset_index(inplace = True)
        raster_rgb[group['row'], group['col']] = colors_rgb[o]
        o = o + 1
        
        
        
    plt.imshow(raster_rgb)



def profiles_pix_coor(directories, plot):
    
    """
    
    This function retrieves the profiles pixels coordinates
    
    inputs:
        
        * directories
        * plot: Boolean. True if the pixels ploted on the image are desidered
        
    
    """
    
    # Load the profiles
    
    names_profiles = [f for f in os.listdir(directories['Profiles']) if ('90_0' in f) and  (f.endswith('_pr.shp'))]
    
    for name_pr in names_profiles:
    
    # profiles_path = os.path.join(directories['Profiles'],  directories['PAN'].split('/')[-2] + '_pr.shp')
    
        profile_path = os.path.join(directories['Profiles'], name_pr)
        
        profiles = gpd.read_file(profile_path, encoding="utf-8")
        
        # Interpolate the profiles
        
        ddis = 2 # distance used to interpolate the profiles
        
        Inter_profiles = profiles.copy()
        
        
        Inter_profiles.geometry = profiles.geometry.apply(redistribute_vertices,distance=ddis)
        
        
        Inter_profiles['nverts'] = Inter_profiles.geometry.apply(lambda x: len(x.coords)) # nverts = numero di punti fianli
        
        # Get the points to be interpolated
        
        gdf = pd.DataFrame()
        
        gdf['points'] = Inter_profiles.apply(lambda l: linestring_to_points(l['PR'],l['geometry']),axis=1)
        
        # Get the indexes
        
        profiles_id = Inter_profiles['PR'].values
        
        profiles_number = Inter_profiles.index.values
        
        # Get the pix coordinates, each image wil have different profile pix coordinates? Yes, AROCSIS doesn't work and each image is at one world's place
        
        epsg = PRIS_tools.get_img_epsgCode(directories)
        
        
        panSharp_img = [f for f in os.listdir(directories['PanSharp_Square_Cropped']) if f.endswith('.tiff')]
    
        
        for pans_img in panSharp_img:
            
            if name_pr.split('_')[1] in pans_img:
        
                path_img = os.path.join(directories['PanSharp_Square_Cropped'], pans_img)
                
                # Get the transform parameters needed to transform from geographic coor to pix coor
                
                source = osr.SpatialReference()
                source.ImportFromEPSG(int(epsg))
                source.ExportToPrettyWkt()
                
                # Load the target image ingo
                
                ds_pan = gdal.Open(path_img)
                target_pan = osr.SpatialReference(wkt = ds_pan.GetProjection())
                transform_pan = osr.CoordinateTransformation(source, target_pan)
                
                # Obtain pixels from world coordinates
                
                prof_coor_pan = DataFrameProfileCoor(profiles_id, profiles_number, gdf, ds_pan, transform_pan)
                
                # Delete the duplicate pixels 
                
                prof_coor_pan.drop_duplicates(subset=['pr', 'row', 'col'], keep = 'first', inplace = True)
                
                # Mantain only  the Dataframe rows where pix coordinates are greater than 0 
                
                prof_coor_pan = prof_coor_pan[(prof_coor_pan['row'] >=0) & (prof_coor_pan['col'] >= 0)]
                
                # Save the info
                
                csv_coor = name_pr.split('.')[0] + '_pix_Coor' + '.csv'
                
                path_csv_coor = os.path.join(directories['Profiles_pix'], csv_coor)
            
                prof_coor_pan.to_csv(path_csv_coor, index = False)    
            
                if plot:
                    
                    plot_profiles(directories, prof_coor_pan)


################################################################################
######################### PROFILES ANALYSIS ####################################
################################################################################    


### PLOT THE PROFILE PIXELS


def plotPR_Pixs(raster, groups, colors, buffer, directories, filename):
    
    """
    
    This function colour the image pixels which spectral signatures will be analyzed
    
    inputs:
        
        * raster = RGB image raster
        * buffer = int; Number of pixels visualize around the pixel profiles
        * groups =  pixels grouped by profile
        * colors = colors assigned to the pixels
        * buffer = pixels around the profile to be shown on the image
        * directories = dict with the folders structure
        
    
    
    """
    
    pathStore = os.path.join(directories['PR_pix'], filename.split('_')[1] + '_' + filename.split('_')[2])
    
    if not os.path.exists(pathStore):
        
        os.makedirs(pathStore)
    
   
    
    colors_rgb = []
    
    for coc in colors:
        colors_rgb.append(hex2color(coc))
    
    
    
    for idd, group in groups:
        
        cc = 0 # color counter
        
        # Reset DataFrame index
        
        group.reset_index(inplace = True, drop = True)
        
        img_rgb = raster.copy()
        
        # Change pixel colors
              
        for row, col in zip(group['row'], group['col']):
        
            img_rgb[row,col,:] = colors_rgb[cc]
        
            cc = cc + 1
        
        # Delete the last pix (water) if the number of pux is not pair

            
        if (len(group) % 2) != 0: # si es impar!!!!!
            group = group[0:len(group)-1]
        
        

        # Find the transects limits
        
        max_r = np.max(group['col'])
        min_r = np.min(group['col'])
        max_c = np.max(group['row'])
        min_c = np.min(group['row'])
        
        if (min_r < 0) or (min_c) < 0:
            continue
        

         # Crop the image
         
        
        
        crop = img_rgb[min_c - buffer : max_c + buffer , min_r - buffer : max_r + buffer, :]
        
        print(crop.shape)
        
        
        
        # SGenerate the image
        
        fig2,ax2 = plt.subplots(1,1,figsize = (20,10))
        plt.title(idd + ' pixels', fontsize = 20)
        plt.axis('off')
        ax2.imshow(crop)
        
        # Save the image
        
        name_fig_crop = idd + '_Pix_' + filename[11:-14] + '.png'
    
        plt.savefig(os.path.join(pathStore, name_fig_crop), tight_layout=True )
        
        plt.close()
        
        del img_rgb
        
        
        
def plotSS_PR_Pix(raster, groups, CW, colors, directories, filename):
    
    
    
    pathStore = os.path.join(directories['SS_pix'], filename.split('_')[1] + '_' + filename.split('_')[2])
    
    if not os.path.exists(pathStore):
        os.makedirs(pathStore)
    
    
    # Reset group index 
    
    for idd, group in groups:
    
        group.reset_index(inplace = True, drop = True)
        
        raster = raster.copy()
        
        
        # Delete the last pixel (water) when the number of pixels is not pair
        
        if (len(group) % 2) != 0:
            group = group[0:len(group)-1]
            
            
        o = 0 # figures rows counter
        u = 0 # figures col counter
        cc = 0 # color counter
        
        fig, axs = plt.subplots(len(group)//2, 2 ,figsize = (20,30),
                                sharex = True)
        plt.suptitle('Spectral signatures profile' + idd + 'pixels', fontsize = 24)
        
        for row, col in zip(group['row'], group['col']):
            
            # Obtain SS
            
            ff = raster[row, col, :]
            axs[o,u].plot(CW, ff, color = colors[cc], linewidth = 5)
            axs[o,u].set_ylim([0,0.45])
            axs[o,u].tick_params(axis = 'both', which = 'major',
                                 labelsize = 14 )
            
            o = o + 1
            cc = cc + 1
            
            if o == len(group)//2:
                u = 1
                o = 0
            
        fig.tight_layout()
        fig.subplots_adjust(top=0.96)
        plt.gcf().text(0.5,-0.1,"Central Wave lengths (nm)", ha="center", fontsize = 16)
        
        name_fig_ss = 'SS_' + idd + '_' + filename[11:-14] + '.png'
        
        plt.savefig(os.path.join(pathStore, name_fig_ss) ,tight_layout=True)
    
        plt.close()
        
        
def SaveSS(raster,groups, directories, filename):
    
    PathStore = os.path.join(directories['SS_pix_csv'], filename.split('_')[1] + '_' + filename.split('_')[2])
    
    if not os.path.exists(PathStore):
        os.makedirs(PathStore)
        
    # Generate the lists needed for the Data
        
    PR = []
    ROW = []
    COL = []
    SIGNA = []
    
    for idd, group in groups:
        
        group.reset_index(inplace = True, drop = True)
        
        if len(group) % 2 != 0:
            
            group = group[0 : len(group) -1]
            
        PR += [idd]*len(group)
        
        for i in range(len(group)):
            
            ROW.append(group['row'][i])
            COL.append(group['col'][i])
            
            ffx = list(raster[group['row'][i], group['col'][i], :])
            
            SIGNA.append(ffx)
            
    SS = pd.DataFrame({'pr': PR, 'row':ROW, 'col':COL, 'signature': SIGNA})
    
    csv_filename = 'SS_Profiles' + '_' + filename.split('.')[0] + '.csv'
    
    SS.to_csv(os.path.join(PathStore, csv_filename), index = False)


        
def profiles_analysis(directories, plot_pr_pix, plot_pix_ss):
    
    # Here we will have to load only once the pixel coordinates (remember to change)
    
    pansharp_imgs = [f for f in os.listdir(directories['PanSharp_Square_Cropped']) if f.endswith('.tiff')]
    
    profiles_pix = [f for f in os.listdir(directories['Profiles_pix']) if f.endswith('.csv')]
    
    
    # Import CW information
    
    original_img = [f for f in os.listdir(directories['scene']) if f.endswith('.he5')]
            
    CW = PRIS_img.ObtainCW(directories ,original_img[0])

    
    tot_img = len(pansharp_imgs)
    
    o = 1
    
    for i in range(len(pansharp_imgs)):
        
        pansharp_img = pansharp_imgs[i]
        
         # FIND THE CORRESPONDENT CSV
        
        for pr_pix in profiles_pix:
            
            if pansharp_img.split('_')[1] in pr_pix:
                
                print(pansharp_img, pr_pix)
            
                
                # Import the profiles info
            
                prof_coor_pan = pd.read_csv(os.path.join(directories['Profiles_pix'], pr_pix))
                
                # Group by profile
                
                groups = prof_coor_pan.groupby(['pr']) 
                
                # Get the number of pixels per profile
                
                len_groups = []
                for idds, group in groups:
                    len_groups.append(len(group))
                    
            
                max_points = np.max(len_groups)    
                
                # Generate as many colors as the maximum number of pixels
                
                colors = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
                         for i in range(max_points)]
                
                
                
                # LOAD THE IMAGE
                
                path_img = os.path.join(directories['PanSharp_Square_Cropped'], pansharp_img)
                
        
                # Get the RGB numpy array
                
                raster_panS, raster_panS_rgb = PRIS_img.numpy_from_tiff(path_img , True)
                #print(raster_panS_rgb.shape)
                
                # PLOT THE IMAGES WITH COLORED PIXELS
                
                if plot_pr_pix:
                
                    buffer = 2
                    
                    
                    plotPR_Pixs(raster_panS_rgb, groups, colors, buffer, directories, pansharp_img )
                
                
                # PLOT THE SPECTRAL SIGNATURES
                
                if plot_pix_ss:
                
                    plotSS_PR_Pix(raster_panS, groups, CW, colors, directories, pansharp_img)
        
                 
                # SAVE THE SS
                
                SaveSS(raster_panS, groups, directories, pansharp_img)
        
                print('Obtaining SS: ', (o/tot_img)*100, '%', end='\r')
        
        
        
        
        

def GetProfilesSS(directories):
    
    profiles_pix_coor(directories, False)
    profiles_analysis(directories, False, False)
    
        
def BaselineCoords(baseline):
    
    baseline_x = []
    baseline_y = []

    for line in baseline['geometry']:
        baseline_x.append(line.x)
        baseline_y.append(line.y)

    return baseline_x, baseline_y

def Cal_ang(lineA, lineB):

    vA = [(lineA[0][0] - lineA[1][0]), (lineA[0][1] - lineA[1][1])]
    vB = [(lineB[0][0] - lineB[1][0]), (lineB[0][1] - lineB[1][1])]

    # Get dot profuct
    dot_pro = dot(vA, vB)
    # Get magnitudes
    magA = dot(vA, vA)**0.5
    magB = dot(vB, vB)**0.5
    # Get cosine value
    cos_ = dot_pro/magA/magB

    # Get angle in radians and then convert to degrees
    if dot_pro == 0:

        ang_deg = np.degrees(np.arccos(0))
    else:
        angle = math.acos(dot_pro/magB/magA)
        # Doing angle <- angle mod 360
        ang_deg = math.degrees(angle)%360

    if ang_deg - 180 >= 0:
        return 360 - ang_deg
    else:
        return ang_deg

def GenerateProfiles(paths, date, symbol, length, rotated, additional):

    """

    Generates authomatically the profiles

    paths: dictionary with the paths
    date: string representing the image date to be analysed
    symbol: + rotates clockwise and - anti
    rotate: degrees that the profiles will be rotated
    additional: adding inclination to the perpendicular profiles
    """

    baselineName = [f for f in os.listdir(paths['Baseline']) if f.endswith('geojson') and date in f][0]
    baseline = gpd.read_file(os.path.join(paths['Baseline'], baselineName))
    baseline_x, baseline_y = BaselineCoords(baseline)

    xmin, xmax, ymin, ymax, north = PRIS_img.ImageBB(paths, date)
    beach_line = [(baseline_x[0], baseline_y[0]), (baseline_x[-1], baseline_y[-1])]
    beach_angle = Cal_ang(north, beach_line)

    lines = []

    for xx, yy in zip(baseline_x, baseline_y):
    
        start = Point(xx,yy)

        if symbol == 'negative':

            angle = -(beach_angle + rotated) + additional

        if symbol == 'positive':

            angle = +(beach_angle + rotated) + additional

        end = Point(start.x, start.y + length)

        Line = LineString([start, end])

        line = rotate(Line, angle, origin = [xx, yy], use_radians = False)

        lines.append(line)


    EpsgCode_img = PRIS_tools.get_img_epsgCode(paths)
    Profiles = gpd.GeoDataFrame(crs = EpsgCode_img, geometry = lines)
    Profiles['PR'] = ['PR' + str(i+1) for i in range(len(lines))]

    return Profiles, beach_angle

def SaveProfiles(paths, filename, Profiles):

    output_pr = os.path.join(paths['Profiles'], filename)
    
    Profiles.to_file(driver = 'ESRI Shapefile', filename = output_pr)

def dot(vA, vB):
    return vA[0]*vB[0] + vA[1] * vB[1]

def AzimuthAngle(coor_pr_i, coor_pr_f):
    angle = (180/np.pi) * math.atan2(coor_pr_f[0] - coor_pr_i[0], coor_pr_f[1] -coor_pr_i[1])
    
    if angle < 0:
        angle = (angle + 360)%360
    return(angle)


def AutomaticProfiles(paths, symbol, length, rotate, additional):

    beach = paths['PanSharp_Square_Cropped'].split('/')[-3]

    imgs = [f for f in os.listdir(paths['PanSharp_Square_Cropped']) if f.endswith('tiff')]
    dates = [f.split('_')[1] for f in imgs]

    for date in dates:
        Profiles, beach_angle = GenerateProfiles(paths, date, symbol, length, rotate, additional)
        filename = beach + '_' + date + '_' + str(rotate) + '_' + 'profiles.shp'
        SaveProfiles(paths, filename, Profiles)
        print('The profiles have been generated...')


def ProfileBounds(pr):
    
    """
    Returns the profile BB
    """
    
    xmin = pr.bounds['minx'].values[0]
    xmax = pr.bounds['maxx'].values[0]
    ymin = pr.bounds['miny'].values[0]
    ymax = pr.bounds['maxy'].values[0]
    
    return xmin, xmax, ymin, ymax

def ProfileBeginning(pr):
    
    """
    Gives the profile starting coordinate
    """
    
    pr.reset_index(inplace = True)
    coords = [(coords) for coords in list(pr.geometry[0].coords)]
    x1 = coords[0][0]
    y1 = coords[0][1]
    
    return x1, y1

def MidPointDistances(overlayed, x_prx, y_prx):
    
    """
    Obtain the midpoint of a segment
    
    """
    
    distances = []

    for point in overlayed['midpoint']:

        mid_point = list(point.coords)
        x2, y2 = mid_point[0][0], mid_point[0][1]

        distances.append(PRIS_tools.EuclideanDistance(x_prx, y_prx, x2 ,y2))
        
    overlayed['dist_midpoint'] = distances
    overlayed.sort_values(by = 'dist_midpoint', ascending = True, inplace = True)
    
    return overlayed


def MeanReflectances(overlayed, arr):
    
    """
    Add to the overlated geodatabase the mean reflectance value of the particular pixel
    """
    
    mean_reflectance = []

    for r, c in zip(overlayed['row'], overlayed['col']):
        ss = arr[r,c,:]
        mean_reflectance.append(np.mean(ss))

    overlayed['mean_R'] = mean_reflectance
    
    return overlayed

def Interpolate_R_PR(overlayed, interval):
    
    x_old = np.array(overlayed['dist_midpoint'])
    y_old = np.array(overlayed['mean_R'])
    
    x_new = np.arange(x_old[0], x_old[-1], interval)
    f = interpolate.interp1d(x_old, y_old, kind = 'cubic')
    y_new =  f(x_new)
    
    deriv_new = y_new[1:] - y_new[:-1]
    idd_inter_new = np.argmin(deriv_new)
    inter = y_new[idd_inter_new]
    
    return idd_inter_new, inter, x_new, y_new, deriv_new

def SSInterfacePixe(GeoT, arr, lon, lat):

    # Pixel containing the transition coordinate
    
    col, row = PRIS_tools.world_to_pixel(GeoT, lon, lat)
    
    # Get the spectral signature
    
    ss = arr[row, col, :]
    
    return col, row, ss
    
    

def PlotProfileValues(paths, date, idd_inter_new, inter, x_new, y_new, deriv_new, pr):

        
    pathS0 = os.path.join(paths['Interface_pix_plot'], 'MeanSS_pr')
    
    if not os.path.exists(pathS0):
        os.mkdir(pathS0)
    
    
    name = 'MeanSS_'+ pr + '_' + date + '.png'

    pathS = os.path.join(pathS0, date)

    if not os.path.exists(pathS):
        os.mkdir(pathS)


        
    plt.scatter(idd_inter_new, inter, s=150, c='firebrick', alpha = 0.5, label = 'interface')
    plt.plot(y_new, marker = 'o', color='royalblue', label = 'mean reflectance') # royalblue
    plt.plot(deriv_new, marker = 'o', color = 'slategray', label = 'derived') # #  skyblue

    ax.spines['top'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_color('grey')
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('grey')
    ax.grid(color='gainsboro', linestyle = '--', axis = 'y')

    plt.title('Profile '+ pr, fontsize = 15)
    plt.xlabel("d (m)", fontsize = 15)
    plt.ylabel("Mean Reflectance", fontsize = 15)
    plt.legend()
    plt.savefig(os.path.join(pathS, name))
    plt.show()
    plt.close()
    
