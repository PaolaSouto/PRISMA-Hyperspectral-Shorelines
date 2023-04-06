#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 21:01:14 2021

This module contains the functions needed to handle the PRISMA images

@author: paolasouto
"""

import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import shutil
import pyproj

## Used in the panSharpening

from scipy import signal
from scipy import ndimage
import PIL
import h5py

from osgeo import osr, ogr, gdal
import osgeo.gdalnumeric as gdn

from PRISMA_SDS import PRIS_tools_v2



################################################################################
############################# BASIC FUNCTIONS ##################################
################################################################################


def create_Tiff(output, img2, driver, GeoT, Projj, num_bands):
    
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
    
    if num_bands == 1:
    
        rows = img2.shape[0]
        cols = img2.shape[1]
        
        
        
        DataSet = driver.Create(output, cols, rows, 1, gdal.GDT_Float32)
        DataSet.SetGeoTransform(GeoT)
        DataSet.SetProjection(Projj.ExportToWkt())
        
        DataSet.GetRasterBand(1).WriteArray(img2)
        DataSet.FlushCache()
        
    if (num_bands != 1 and num_bands != 3): 
        
        rows = img2.shape[0]
        cols = img2.shape[1]
        band = img2.shape[2]
        
        
        DataSet = driver.Create(output, cols, rows, band, gdal.GDT_Float32)
        DataSet.SetGeoTransform(GeoT)
        DataSet.SetProjection(Projj.ExportToWkt())
        
        for i in range(band):
            
            image = img2[:,:,i]
            DataSet.GetRasterBand(i+1).WriteArray(image)
            
        DataSet.FlushCache()
        
    if num_bands == 3: # PRISMA RGB IMAGES 
    
        rows = img2.shape[0]
        cols = img2.shape[1]
        
        
        DataSet = driver.Create(output, cols, rows, 3, gdal.GDT_Float32)
        DataSet.SetGeoTransform(GeoT)
        DataSet.SetProjection(Projj.ExportToWkt()) #.ExportToWkt()
  
        for i,o in zip([29,19,9], [0,1,2]): # Bands defining the RGB image [20,19,9]
            
            image = img2[:,:,i]
            DataSet.GetRasterBand(o+1).WriteArray(image)
            
        DataSet.FlushCache()
    
    
    
def numpy_from_tiff(img_path, RGB):

    
    """
    The function gives the image np.array
    
    inputs:
        
        * img_path : image path
        * RGB : Boolean (True/False). True if the desidered output is the RGB np.array
        
    outputs:
        
        * arr: image np.array
        
    Based on : https://gis.stackexchange.com/questions/32995/fully-load-raster-into-a-numpy-array 
        
    """
    
    img = gdal.Open(img_path)
    
    bands = [img.GetRasterBand(i) for i in range(1, img.RasterCount + 1)]
    
    arr = np.array([gdn.BandReadAsArray(band) for band in bands]).astype('float32')
    
    arr = np.transpose(arr, [1, 2, 0])
    
    if RGB:
        
        arr_rgb = []
        
        for n in [29,19,9]:
            
            ff = arr[:,:,n]
            arr_rgb.append(ff)
            
        arr_rgb = np.stack(arr_rgb)
        arr_rgb = np.transpose(arr_rgb, [1,2,0])
        
        return arr, arr_rgb
    
    else:
        return arr
    
    
    
def ObtainCW(directories ,img):
    
    
    dir_image = os.path.join(directories['scene'], img)
    
    # Load the image
    
    f = h5py.File(dir_image, 'r')
    
    # Load the CW images
    
    CW = np.concatenate([f.attrs['List_Cw_Vnir'][::-1], f.attrs['List_Cw_Swir'][::-1]])
    

    # Load the Flags
    
    Flag = np.concatenate([f.attrs['CNM_VNIR_SELECT'][::-1], f.attrs['CNM_SWIR_SELECT'][::-1]])
    
    CW2 = CW[Flag.astype(bool)]
    
    CW3 = CW2[np.argsort(CW2)]
    
    return CW3

################################################################################
################# GENERATE AND SAVE THE TIFFS ##################################
################################################################################


def generate_PAN(directories, error_matrix):
    
    """
    This function reads the panchromatics bands of each one of the PRISMA images stored on the defined scene and store them in 
    the PAN folder in tiff format
    
    Inputs:
        
        directories: dictionary with all the paths in the in the scene of interest folder
        error_matrix: boolean indicating if the image error matrix must be applied or not. Mask the satured pixels
        
    Outputs:
        
        PAN images stored as 'Tiff' on the 'PAN' folder
    
    """
    
    #### READ THE PRISMA IMAGES CONTAINED ON THE SCENE FOLDER
    
    images_names = [f for f in os.listdir(directories['scene']) if (f.endswith('.he5') and not f.startswith('._'))]
    
    tot_img = len(images_names)
    
    print('The total number of images: ', tot_img)
    
    #### LOAD THROUGH THE IMAGES
    
    o = 1# countr
    
    for img_name in images_names:
        
        path_img = os.path.join(directories['scene'], img_name)
        
        if not os.path.exists(path_img):
            raise Exception('The image does not exists')
        
        # Load the image
        
        f = h5py.File(path_img, 'r')
        
        # Load the PAN band
        
        panchro = np.array(f['HDFEOS']['SWATHS']['PRS_L2D_PCO']['Data Fields']['Cube'])
        
        #  Transformation from unit16 DN to reflectance values
        
        ScalePanMin = f.attrs['L2ScalePanMin']
        ScalePanMax = f.attrs['L2ScalePanMax']
        
        # Transform from DN to reflectance
    
        #print('Transforming from DN to reflectance ...')
        
        panchro2 = ScalePanMin + panchro*(ScalePanMax-ScalePanMin)/65535
        
        # Apply the Error matrix
    
        if error_matrix :
        
            #print('Applying the error matrix ...')
            
            ERR_Pan = np.array(f['HDFEOS']['SWATHS']['PRS_L2D_PCO']['Data Fields']['PIXEL_L2_ERR_MATRIX'])
            
            idxs = np.where( ERR_Pan != 0)
        
            panchro2[idxs[0], idxs[1]] = 0
            
            
        # Geo with the attributes given in the image
    
        geo = {'proj_code':f.attrs['Projection_Id'], 
               'proj_name':f.attrs['Projection_Name'],
               'proj_epsg':f.attrs['Epsg_Code'], 
               'xmin':np.min([f.attrs['Product_ULcorner_easting'], f.attrs['Product_LLcorner_easting']]),
               'xmax':np.max([f.attrs['Product_LRcorner_easting'], f.attrs['Product_URcorner_easting']]),
               'ymin':np.min([f.attrs['Product_LLcorner_northing'], f.attrs['Product_LRcorner_northing']]),
               'ymax':np.max([f.attrs['Product_ULcorner_northing'], f.attrs['Product_URcorner_northing']])}
    
        res = 5
        
    
        
        GeoT, driver, Proje = PRIS_tools_v2.get_GeoTransform(geo, panchro2, res)
        
        ## Save the images
        
        
        dir_save = directories['PAN']
        
        if not os.path.exists(dir_save):
            raise Exception('Wrong PAN path')
            
        new_name = img_name[:-4] + 'pan.tiff'
        
        output = os.path.join(dir_save, new_name)
        
        print('Saved PAN images: ', (o/tot_img)*100, '%', end='\r')
        
        #print('Saved images: ', (o/tot_img)*100, '%')
        
        o += 1
        
        # Define the number of bands
        
        num_bands = 1

        create_Tiff(output, panchro2, driver, GeoT, Proje, num_bands)
        
     
        
     
## IMPORTANTE MIRAR SI TENGO QUE GUARDAR LAS BANDAS ORDENADAS O NO
## REMINDER PARA MI: BANDAS DEL SWIR Y DEL VNIR CUBE SE SOBREPONEN, EN EL SCRIPT INICIAL NO ORDENO LAS BANDAS.
## HABRÀ QUE VOLVER A MIRAR ESTO

def generate_HS(directories, error_matrix, save_RGB):
    
    """
    This function reads the HS images of each one of the PRISMA images stored on the defined scene and store them in 
    the HS folder in tiff format
    
    Inputs:
        
        directories: dictionary with all the paths in the in the scene of interest folder
        error_matrix: boolean indicating if the image error matrix must be applied or not
        save_rgb: True if the rgb image must be stored
        
    Outputs:
        
        HS images stored as 'Tiff' on the 'HS' folder
    
    """
    
    #### READ THE PRISMA IMAGES CONTAINED ON THE SCENE FOLDER
    
    images_names = [f for f in os.listdir(directories['scene']) if (f.endswith('.he5') and not f.startswith('._'))]
    
    tot_img = len(images_names)
    
    print('The total number of images: ', tot_img)
    
    #### LOAD THROUGH THE IMAGES
    
    o = 1 # count HS images
    u = 1 # count rgb images
    
    for img_name in images_names:
        
        path_img = os.path.join(directories['scene'], img_name)
        
        if not os.path.exists(path_img):
            raise Exception('The image does not exists')
        
        # Load the image
        
        f = h5py.File(path_img, 'r')
        
        # Geo with the attributes given in the image
    
        geo = {'proj_code':f.attrs['Projection_Id'], 'proj_name':f.attrs['Projection_Name'],
                              'proj_epsg':f.attrs['Epsg_Code'], 'xmin':np.min([f.attrs['Product_ULcorner_easting'], f.attrs['Product_LLcorner_easting']]),
          'xmax':np.max([f.attrs['Product_LRcorner_easting'], f.attrs['Product_URcorner_easting']]),
          'ymin':np.min([f.attrs['Product_LLcorner_northing'], f.attrs['Product_LRcorner_northing']]),
          'ymax':np.max([f.attrs['Product_ULcorner_northing'], f.attrs['Product_URcorner_northing']])}
        
        # Load the CW images
        
        CW = np.concatenate([f.attrs['List_Cw_Vnir'][::-1], f.attrs['List_Cw_Swir'][::-1]])
        
        # Load the band width
        
        #BandWidth = np.concatenate([f.attrs['List_Fwhm_Vnir'][::-1], f.attrs['List_Fwhm_Swir'][::-1]])
        
    
        
        # Load the Flags (add description if we have time)
        
        Flag = np.concatenate([f.attrs['CNM_VNIR_SELECT'][::-1], f.attrs['CNM_SWIR_SELECT'][::-1]])
        
        # Load the Geolocation Information
        
        #Lat = np.array(f['HDFEOS']['SWATHS']['PRS_L2D_HCO']['Geolocation Fields']['Latitude'])
        #Lon = np.array(f['HDFEOS']['SWATHS']['PRS_L2D_HCO']['Geolocation Fields']['Longitude'])
        
        # Load the HS cube
        
        SWIR_bands =np.array(f['HDFEOS']['SWATHS']['PRS_L2D_HCO']['Data Fields']['SWIR_Cube'])
        VNIR_bands = np.array(f['HDFEOS']['SWATHS']['PRS_L2D_HCO']['Data Fields']['VNIR_Cube'])
        SWIR_bands_C = np.swapaxes(SWIR_bands, 1, 2)
        VNIR_bands_C = np.swapaxes(VNIR_bands, 1, 2)
        VNIR_bands_CC = VNIR_bands_C[:, :, ::-1]
        SWIR_bands_CC = SWIR_bands_C[:, :, ::-1]
        
        # Load the parameters for scale the DN
        
        L2ScaleSwirMax = f.attrs['L2ScaleSwirMax']
        L2ScaleSwirMin = f.attrs['L2ScaleSwirMin']
        L2ScaleVnirMax = f.attrs['L2ScaleVnirMax']
        L2ScaleVnirMin = f.attrs['L2ScaleVnirMin']
        
        # Aplly the correction
    
        # On the SWIR Cube
    
        SWIR_bands_R = np.float32(SWIR_bands_CC.copy())
        
        for n in range(SWIR_bands_CC.shape[2]): # Pensar si le queremos añadir los tqdm o no
            
            SWIR_bands_R[:,:,n] = L2ScaleSwirMin + SWIR_bands_CC[:,:,n]*\
                (L2ScaleSwirMax-L2ScaleSwirMin)/65535
                
        # On the VNIR Cube
        
        VNIR_bands_R = np.float32(VNIR_bands_CC.copy())
        
        for n in range(VNIR_bands_CC.shape[2]): # Pensar si le queremos añadir los tqdm o no
            
            VNIR_bands_R[:,:,n] = L2ScaleVnirMin + VNIR_bands_CC[:,:,n]*\
                (L2ScaleVnirMax - L2ScaleVnirMin)/65535
                
                
        # Concatenate the Cubes
        
        img = np.concatenate([VNIR_bands_R,SWIR_bands_R], axis=2)
        
        # Delete "incorrect" bands
        
        img2 = img[:,:,Flag.astype(bool)]
        CW2 = CW[Flag.astype(bool)]
        #BandWidth2 = BandWidth[Flag.astype(bool)] 
        

        
        
        # Apply the error matrix (¡¡¡¡ MAYBE IT IS BETTER TO DO IT BY DEFAULT, THINK!!!!)
        
        if error_matrix:
                    
            ERR_VNIR = np.array(f['HDFEOS']['SWATHS']['PRS_L2D_HCO']['Data Fields']['VNIR_PIXEL_L2_ERR_MATRIX'])  
            ERR_VNIR_C = np.swapaxes(ERR_VNIR, 1, 2)
            ERR_VNIR_CC = ERR_VNIR_C[:, :, ::-1]
            
            
            ERR_SWIR = np.array(f['HDFEOS']['SWATHS']['PRS_L2D_HCO']['Data Fields']['SWIR_PIXEL_L2_ERR_MATRIX'])  
            ERR_SWIR_C = np.swapaxes(ERR_SWIR, 1, 2)
            ERR_SWIR_CC = ERR_SWIR_C[:, :, ::-1]
            
            ERR = np.concatenate([ERR_VNIR_CC,ERR_SWIR_CC], axis=2)
            ERR_C = ERR[:,:,Flag.astype(bool)]
    
            
            for n in range(img2.shape[2]): #(THINK ABOUT TQDM)
                idx = np.where(ERR_C[:,:,n] != 0)
                # Estamos probando con img3
                img2[idx[0], idx[1],n] = -999 # Tb lo podriamos poner a 0 (tengo que explorar)
                #print('He cambiado a nan')
          
        # Re - ordering the bands (PORQUE NO USA ESTA ???)
        
        CW3 = CW2[np.argsort(CW2)]
        img3 = img2[:,:, np.argsort(CW2)]
        
        
        # Resolution of the HS cube
        
        res = 30
        
    
        GeoT, driver, Proje = PRIS_tools_v2.get_GeoTransform(geo, img3, res)  
        
        
        ## Save the images
        
        
        dir_save = directories['HS']
        
        if not os.path.exists(dir_save):
            raise Exception('Wrong HS path')
            
        new_name = img_name[:-4] + 'hs.tiff'
        
        output = os.path.join(dir_save, new_name)

        # Define the number of bands (aleatory, must be higher than 3)

        num_bands = 150    

        create_Tiff(output, img3, driver, GeoT, Proje, num_bands )
        
        print('Saved HS images: ', (o/tot_img)*100, '%', end='\r')
        
        o += 1
        
        if save_RGB:
            
            img3[img3 == -999] = 0

            dir_save_RGB = directories['HS_RGB']
            
            if not os.path.exists(dir_save_RGB):
                
                raise Exception('Wrong HS RGB path')
            
            new_name_RGB = img_name[:-4] + 'rgb.tiff'
            
            output_RGB = os.path.join(dir_save_RGB, new_name_RGB)
            
            # Define the number of bands
            
            num_bands = 3
            
            create_Tiff(output_RGB, img3, driver, GeoT, Proje, num_bands)
            
            print('Saved HS RGB images: ', (u/tot_img)*100, '%', end='\r')
            
            u = u + 1
            
            
            
################################################################################
############################# CROPPING IMAGES ##################################
################################################################################

# Aqui ahora cargamos unas coordenadas no un shapefile, ver como lo hacemos

def crop_images(directories):
    
    """
    
    This function crops the PAN and HS images to the beach selected using the GUI
    
    inputs:
        
        * directories: dict with all the paths 
        * extent_modi: extension of the AOI obtained by clicking in the GUI
        
    outputs:
    
        * Cropped HS and PAN images
        
    LOOK HOW TO DELETE THE FOLDER WITH THE SHAPEFILE
    LOOK TO DELETE THE TEMPORAL FOLDER IN DIRECTORIES DICT
    TRY TO AVOID CALL THE SHELL https://gis.stackexchange.com/questions/262021/how-to-replicate-gdalwarp-cutline-in-python-to-cut-out-tif-image-based-on-the-e
    
    """
    
        
    # Load the polygon
    
    inshape = os.path.join(directories['Crop_shp'], 'Square_polygon.shp')
    
    ds = ogr.Open(inshape)
    lyr = ds.GetLayer(0)
    
    
    ###### HS IMAGES
    
    ## Load the HS co-registered images
    
    hs_images = [f for f in os.listdir(directories['HS_CoReg']) if (f.endswith('.tiff') and not f.startswith('_.'))]
    
    # Loop throught the images
    
    o = 0
    
    tot_img = len(hs_images)
    
    for hs_img in hs_images:
        
        inraster = os.path.join(directories['HS_CoReg'], hs_img)
        
        # Define the image output name
        
        file_out = 'SquareCrop_' + hs_img[12:-4] + 'tiff'
        
        outraster = os.path.join(directories['HS_Square_Cropped'], file_out)
        
        cmd = 'gdalwarp '+ '-cutline ' + inshape  +' -crop_to_cutline' + ' -dstnodata "-999.0" '+ inraster + ' ' + outraster
        
        os.system(cmd)
        
        o += 1
        
        print('Cropped HS images: ', (o/tot_img)*100, '%', end='\r')
        
        
    ###### PAN IMAGES
    
    pan_images = [f for f in os.listdir(directories['PAN_CoReg']) if (f.endswith('.tiff') and not f.startswith('_.'))]
    
    
    o = 0
    
    for pan_img in pan_images:
        
        inraster = os.path.join(directories['PAN_CoReg'], pan_img)
        
        file_out = 'SquareCrop_' + pan_img[12:-4] + 'tiff'
        
        outraster = os.path.join(directories['PAN_Square_Cropped'], file_out)
        
        cmd = 'gdalwarp '+ '-cutline ' + inshape +  ' -crop_to_cutline' + ' -dstnodata "-999.0" '+ inraster + ' ' + outraster 

        os.system(cmd)
        
        o += 1
        
        print('Cropped PAN images: ', (o/tot_img)*100, '%', end='\r')
        

    ## Delete the temporal path
    
    #shutil.rmtree(directories['Temporal'], ignore_errors=True)
    
    
    
################################################################################
############################# PAN-SHARPENING ##################################
################################################################################


def alpha_estimation(a, imageHR0):
    
    IHc = imageHR0.reshape((imageHR0.shape[0]*imageHR0.shape[1],1), order='F')
    ILRc = a.reshape((a.shape[0]*a.shape[1], a.shape[2]),order='F')
    alpha = np.linalg.lstsq(ILRc,IHc)
    alpha = alpha[0]
    
    return alpha



def PanSharpening_GSA(HS, PAN): 
    
    """
    
    Performs the PanSahrpening GSA
    
    inputs:
        
        * HS: np.array with the HS image
        * PAN: np.array with the PAN image
        
    outputs:
        
        * I_Fus_GSA: Pansharpened image with the PAN spatial resolution and the 
        hyperspectral cube spectral resolution
    
    """
    
    ratio1 = PAN.shape[0]/HS.shape[0]
    
    
    if not ratio1.is_integer():
        
        raise Exception('The ratio is not an iteger')
 
        
    ratio1 = int(ratio1)
    
    
    # Aplicamos el filtro
    
    r_im, c_im, b_im = HS.shape
    r_pan, c_pan, b_pan = PAN.shape
    
    L = 45
    
    BaseCoeff = ratio1*signal.firwin(L,1/ratio1)
    
    I1LRU = np.zeros([ratio1*r_im, ratio1*c_im, b_im])
    
    I1LRU[0:r_pan:ratio1, 0:c_pan:ratio1, :] = HS
    
    m = L//2
    
    Filtro = np.zeros([L,L])
    Filtro[m,:] = BaseCoeff
    
    print('Paso1')
    
    for nn in range(I1LRU.shape[2]):
        
        
        t = I1LRU[:,:,nn]
        t = ndimage.filters.convolve(t.T, Filtro, mode = 'wrap')
        I1LRU[:,:,nn] = ndimage.filters.convolve(t.T, Filtro, mode = 'wrap')
        
    
    imageLR = I1LRU
    
    #REMOVE MEANS FROM imageLR
    
    print('Paso 2')
    
    imageLR0 = np.zeros([imageLR.shape[0], imageLR.shape[1], imageLR.shape[2]])
        
    for i in range(imageLR.shape[2]):
        imageLR0[:,:,i] = imageLR[:,:,i] - np.mean(imageLR[:,:,i])
        
        
    # REMOVE MEANS FROM imageLR_LP
    
    print('Paso 3')
    
    imageLR_LP0 = np.zeros([HS.shape[0], HS.shape[1], HS.shape[2]])
    
    for ii in range(HS.shape[2]):
        imageLR_LP0[:,:,ii] = HS[:,:,ii] -np.mean(HS[:,:,ii])
        
    
    # Synthetic intensity
    
    print('Paso 4')
    
    imageHR = PAN[:,:,0]
    
    imageHR0 = imageHR - np.mean(imageHR)
    
    image_pil = PIL.Image.fromarray(imageHR0)
    imageHR0 = image_pil.resize((HS.shape[1], HS.shape[0]), resample = PIL.Image.BICUBIC )
    imageHR0 = np.asarray(imageHR0)
    
    a = np.dstack((imageLR_LP0, np.ones([HS.shape[0], HS.shape[1]])))
    
    alpha = alpha_estimation(a, imageHR0)
    alpha = np.reshape(alpha, (1,1,len(alpha)))
    
    kk2 = np.tile(alpha,(imageLR.shape[0],imageLR.shape[1], 1))
    kk = np.dstack((imageLR0, np.ones([imageLR.shape[0], imageLR.shape[1]])))
    gg2 = kk * kk2
    I = np.sum(gg2, axis=2)
    

    # REMOVE MEAN FROM I
    
    I0 = I - np.mean(I)
    
    # COEFFICIENTS
    
    print('Paso 5')
    
    g = np.ones([1,1, imageLR.shape[2] +1])
    
    for i in range(imageLR.shape[2]):
        h = imageLR0[:,:,i]
        c = np.cov(I0.flatten(), h.flatten())
        g[0,0,i+1] = c[0,1]/np.var(I0.flatten())
        
    imageHR = imageHR - np.mean(imageHR)
    
    # DETAIL EXTRACTION
    
    delta = imageHR - I0
    deltam = np.tile(delta.T.flatten(),(imageLR.shape[2]+1));
    deltam = deltam.reshape([delta.T.flatten().shape[0], imageLR.shape[2]+1], order = 'F')
    
    
    # FUSION
    
    print('Paso 6')
    
    V = I0.T.flatten()

    for ii in range(imageLR.shape[2]):
        h = imageLR0[:,:,ii]
        V = np.concatenate((V,h.T.flatten()))
        
    V = V.reshape([delta.T.flatten().shape[0], imageLR.shape[2]+1], order = 'F')
    
    gm = np.zeros([V.shape[0], V.shape[1]])

    for ii in range(g.shape[2]):
        pp = np.squeeze(g[0,0,ii]) * np.ones([imageLR.shape[0] * imageLR.shape[1],1])
        gm[:,ii] = pp.T
        
    V_hat = V + deltam * gm
    
    
    
    # RESHAPE FUSION RESULT
    

    V_hat_r = V_hat[:,1:V_hat.shape[1]]
    
    I_Fus_GSA = V_hat_r.reshape([imageLR.shape[0], imageLR.shape[1], imageLR.shape[2]],
                               order = 'F');
    
    # FINAL MEAN EQUALIZATION
    
    print('Paso 7')
    
            
    for ii in range(imageLR.shape[2]):
        h = I_Fus_GSA[:,:,ii]
        I_Fus_GSA[:,:,ii] = h - np.mean(h) + np.mean(imageLR[:,:,ii])
        
        
        
    return I_Fus_GSA


def PanSharpening(directories):
    
    """
    
    This is the function that performs the panSharpening
    
    inputs:
        
        * directories: dict with the folders paths folder structure
        * save_RGB : boolean (True/False). True if the RGB imahe has to be saved
    
    """
    
    # Call the images
    
    PAN_cropped_img = [f for f in os.listdir(directories['PAN_Square_Cropped']) if f.endswith('.tiff')]
    
    HS_cropped_img= [f for f in os.listdir(directories['HS_Square_Cropped']) if f.endswith('.tiff')]
    
    tot_img = len(PAN_cropped_img)
    
    # Perform the panSharpening
    
    # Load the images
    
    u = 1
    
    for (PAN_img, HS_img) in zip(PAN_cropped_img, HS_cropped_img):
        
        # Numpy array of the PAN and the HS images
        
        PAN_path = os.path.join(directories['PAN_Square_Cropped'], PAN_img)
        PAN_np = numpy_from_tiff(PAN_path, False)
        
        HS_path = os.path.join(directories['HS_Square_Cropped'], HS_img)
        HS_np = numpy_from_tiff(HS_path, False)
        
        # Perform the PanSharpening
        
        PanSharp_np = PanSharpening_GSA(HS_np, PAN_np)
        
        
        
        
        # Save the PanSharpened image and its RGB version
        
        #geo = PRIS_tools_v2.crop_shp_geo(directories)
        
        
        # Obtain the GeoTransform
        
        res = 5
        
        ds = gdal.Open(PAN_path)
        GeoT = ds.GetGeoTransform()
        #GeoT = (geo['xmin'], res, 0, geo['ymax'], 0, -res)
    
        driver = gdal.GetDriverByName(str("GTiff"))
        
        proj = osr.SpatialReference(wkt=ds.GetProjection())
        crs =  proj.GetAttrValue('AUTHORITY',1)
        
        Projj = osr.SpatialReference()
        Projj.ImportFromEPSG(int(crs)) #4326
        Projj.ExportToPrettyWkt()
        
        #GeoT, driver, Projj = PRIS_tools.get_GeoTransform(geo, PanSharp_np, res)
        
        # Save both PanSharp and PanSharp_rgb
        
        name_panS = 'PanSharpened' + PAN_img[10:-14] + '.tiff'
        out_panS = os.path.join(directories['PanSharp_Square_Cropped'], name_panS)
        
        name_panS_rgb = 'PanSharpened_RGB' + PAN_img[10:-14] + '.tiff'
        out_panS_rgb = os.path.join(directories['PanSharp_RGB'], name_panS_rgb)
        
        # Ahora las salvamos
        
        create_Tiff(out_panS, PanSharp_np, driver, GeoT, Projj, PanSharp_np.shape[2])
        create_Tiff(out_panS_rgb, PanSharp_np, driver, GeoT, Projj, 3)


        print('Saved PanSharp RGB images: ', (u/tot_img)*100, '%', end='\r')
        
        u += 1
        

def crop_PanSharp_to_beach(directories):
    
    """
    
    This function crops the PanSharpened 
    
    """
    
    # Generate the shapefile
    
    PRIS_tools_v2.beach_crop_shapefile(directories)
    
    # Shapefiles beach crop
    
    beaches_polygons = [f for f in os.listdir(directories['Crop_shp']) if f.startswith('Beach_polygon') and f.endswith('.shp')]
    
    
    # Retrive the PanSharp images
    
    PanSharp_imgs = [f for f in os.listdir(directories['PanSharp_Square_Cropped']) if f.endswith('.tiff')]
    
    tot_img = len(PanSharp_imgs)
    
    
    for beach in beaches_polygons:
        
        beach_and_date = beach.split('_')[2] + '_' + beach.split('_')[3]

        # Generate the Beach Cropped PanSharp images
        
        u = 1
        
        for img in PanSharp_imgs:
            
            if  beach_and_date in img:
                
                print(beach, img)
            
                new_name = 'BeachCrop_' + beach_and_date + '.tiff'
        
                
                OutTile = gdal.Warp(os.path.join(directories['PanSharp_Beach_Cropped'], new_name),
                                    os.path.join(directories['PanSharp_Square_Cropped'], img),
                                    cutlineDSName = os.path.join(directories['Crop_shp'],beach),
                                    cropToCutline = True,
                                    dstNodata = 0)
                
            print('Saved Beach PanSharp  images: ', (u/tot_img)*100, '%', end='\r')
                
            u += 1
            
            OutTile = None
        
        
    
def Spectral_Angle_Mapper(si, sj):
    
    """The smaller the spectral angle the more similar pixel isto a given endmember class"""
    
    den = 0
    sptr1 = 0
    sptr2 = 0
    
    for l in range(len(si)):
        
        den += si[l] * sj[l]
        sptr1 += si[l]**2
        sptr2 += sj[l]**2
        
    angle = np.arccos(den/(np.sqrt(sptr1) * np.sqrt(sptr2)))
    
    return angle


        
