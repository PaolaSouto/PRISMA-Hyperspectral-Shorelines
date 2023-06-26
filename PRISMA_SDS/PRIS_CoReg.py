#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 18:35:13 2021

This module performs the co-registration of the PRISMA images

The module uses the package AROSICS (https://git.gfz-potsdam.de/danschef/arosics/-/blob/master/arosics/CoReg.py) written by Daniel Scheffler

This module uses the GLOBAL Co-registration method that performs the shifts maintaining the 
original image values

@author: paolasouto
"""

import os
import numpy as np
from osgeo import gdal, osr
import pandas as pd
import shutil
import subprocess

from arosics import  COREG

from PRISMA_SDS import PRIS_img, PRIS_tools

def coreg_PAN(directories, PAN_ref, reference_img):
    
    """
    
    This function co-registrate the panchromatic to a reference one (the one with higher geoaccuracy or the one similer to a given ortophoto)
    
    inputs:
        
        * directories = dictionary with all the paths in the in the scene of interest folder
        * PAN_ref = boolean
        * reference_img = string with the name of the PAN image used to corregistrate the other ones
        
        
    outputs:
        
        * PAN images co-registered in the /PAN/CoReg path
        * csv file with the shifts applied to each image. The csv will be stored in /Results/CoReg_Shifts
    

    """
    
    # Check the PRISMA EPSG code
    
    EpsgCode = PRIS_tools.get_img_epsgCode(directories)
        
    code = np.unique(EpsgCode)
    
    pan_images = [f for f in os.listdir(directories['PAN']) if (f.endswith('.tiff') and not f.startswith('._')) and f != reference_img]
    
    if PAN_ref:
    
        im_pan = os.path.join(directories['PAN'], reference_img)
        
    else:
        
        
        im_pan = os.path.join(directories['VHR'], reference_img)
        


        # If it 
    
    # List that will contain the results obtained from the co-registration process
    
    date = []
    corrected_shift_x = []
    corrected_shift_y = []
    corrected_shift_row = []
    corrected_shift_col = []
    shift_reliability = []
    img_similarity_bef = []
    img_similarity_aft = []
    
    # Co-registration process
    
    for pan_img in pan_images:
        
        new_name = pan_img[0:-5] + '_coreg.tiff'
        
        output = os.path.join(directories['PAN_CoReg'], new_name)
    
        date.append(pan_img[12:26])
        
        im_target = os.path.join(directories['PAN'], pan_img)
        
        # 1500, 1500
        CR = COREG(im_pan, im_target, ws=(1200, 1200), 
                   max_shift = 7000, max_iter = 200, v = False);
        
        # Save the outputs of the co-registration process
        
        corrected_shift_x.append(CR.coreg_info['corrected_shifts_map']['x'])
        corrected_shift_y.append(CR.coreg_info['corrected_shifts_map']['y'])
        corrected_shift_col.append(CR.coreg_info['corrected_shifts_px']['x'])
        corrected_shift_row.append(CR.coreg_info['corrected_shifts_px']['y'])
        shift_reliability.append(CR.shift_reliability)
        img_similarity_bef.append(CR.ssim_orig)
        img_similarity_aft.append(CR.ssim_deshifted)
        
        # Correct the shifts
        
        GeoArray = CR.correct_shifts();
        
        # Extract the corrected image 
        
        ff = np.asarray(GeoArray['GeoArray_shifted'])
        ff = ff.astype('float32')
        
        # Save the co-registered images
        
        driver = gdal.GetDriverByName(str('GTiff'))
        GeoT = GeoArray['updated geotransform']
        Projj = osr.SpatialReference()
        Projj.ImportFromWkt(GeoArray['updated projection'])
        Projj.ExportToPrettyWkt()
        num_bands = 1
        
        PRIS_img.create_Tiff(output, ff, driver, GeoT, Projj, num_bands)
        
    # Move the reference image to the /PAN/CoReg folder
    
    if PAN_ref:
    
        im_pan_dest = os.path.join(directories['PAN_CoReg'], reference_img)
    
        shutil.move(im_pan, im_pan_dest)
        
        # Save the csv with the calculated shifts
        
        filename = 'PAN' + '_Pix_Shifts_from_' + reference_img.split('/')[-1][0:46] + '.csv'
        
    else:   
        
        filename = 'PAN' + 'Pix_Shifts_from_VHRimage.csv'
    
    shifts = pd.DataFrame({'Date':date, 'x_shift':corrected_shift_x, 'y_shift':corrected_shift_y, 'row_shift':corrected_shift_row,
                      'col_shift':corrected_shift_col, 'shift_real': shift_reliability, 'SSIM bef': img_similarity_bef,
                      'SSIM aft': img_similarity_aft})
    
    shifts.to_csv(os.path.join(directories['CoReg_Shifts'], filename), index = False)
    
    

def coreg_HS(directories): 
    
    """
    
    Co-registrate each HS image to its correspondet pan image that has been already 
    co-registered
    
    inputs:
        
        * directories
    
    outputs:
        
        * Co-registered hs images on path /HS/CoReg
        * csv file with the shifts applied to each image. The csv will be stored in /Results/CoReg_Shifts
    
    """
    
    pan_images = [f for f in os.listdir(directories['PAN_CoReg']) if (f.endswith('.tiff') and not f.startswith('._')) ]
    hs_images = [f for f in os.listdir(directories['HS']) if (f.endswith('.tiff') and not f.startswith('._')) ]
    
    
    
    # List that will contain the results obtained from the co-registration process
    
    date = []; corrected_shift_x = []; corrected_shift_y = []
    corrected_shift_row = []; corrected_shift_col = []
    shift_reliability = []; img_similarity_bef = []
    img_similarity_aft = []   
    
    # Co-registration process
    
    for img_pan in pan_images:
        
        for img_hs in hs_images:
            
            if img_hs[0:46] == img_pan[0:46]:
                
                img_reference = os.path.join(directories['PAN_CoReg'], img_pan )
                img_target = os.path.join(directories['HS'], img_hs)
                
                date.append(img_hs[12:26])
                
                CR = COREG(img_reference, img_target, ws=(600,600), max_shift = 100, max_iter=800)
                
                corrected_shift_x.append(CR.coreg_info['corrected_shifts_map']['x'])
                corrected_shift_y.append(CR.coreg_info['corrected_shifts_map']['y'])
                corrected_shift_col.append(CR.coreg_info['corrected_shifts_px']['x'])
                corrected_shift_row.append(CR.coreg_info['corrected_shifts_px']['y'])
                shift_reliability.append(CR.shift_reliability)
                img_similarity_bef.append(CR.ssim_orig)
                img_similarity_aft.append(CR.ssim_deshifted)
                
                # Correct the shifts
                
                GeoArray = CR.correct_shifts()
                
                # Obtain the corrected image
                
                ff = np.asarray(GeoArray['GeoArray_shifted'])
                ff = ff.astype('float32')
                
                # Save the image
                
                new_name = img_hs[0:-5] + '_coreg.tiff'
                output = os.path.join(directories['HS_CoReg'], new_name)
                
                driver = gdal.GetDriverByName(str('GTiff'))
                GeoT = GeoArray['updated geotransform']
                Projj = osr.SpatialReference()
                Projj.ImportFromWkt(GeoArray['updated projection'])
                Projj.ExportToPrettyWkt()
                
                num_bands = 150 # aleatory number greater than 3
                
                PRIS_img.create_Tiff(output, ff, driver, GeoT, Projj, num_bands)
    
    # Save the calculated shifts between images
    
    filename = 'HS' + 'Pix_Shifts.csv'
    
    shifts = pd.DataFrame({'Date':date, 'x_shift':corrected_shift_x, 'y_shift':corrected_shift_y, 'row_shift':corrected_shift_row,
                      'col_shift':corrected_shift_col, 'shift_real': shift_reliability, 'SSIM bef': img_similarity_bef,
                      'SSIM aft': img_similarity_aft})
    
    shifts.to_csv(os.path.join(directories['CoReg_Shifts'], filename), index = False)
    
    