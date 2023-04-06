#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 19:11:09 2021

This module retrieve the shoreline

@author: paolasouto
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from osgeo import osr, ogr, gdal
from spectral import *
from scipy import ndimage
from scipy import interpolate
import itertools
from skimage.transform import resize

## Packages for unmixing

import pysptools.util as util
import pysptools.eea as eea #endmembers extraction algorithms
import pysptools.abundance_maps as amap
import pysptools.classification as cls
import pysptools.material_count as cnt
import pysptools.distance as dst




from PRISMA_SDS import PRIS_tools_v2, PRIS_img_v2


def ReadSS(path, filename): # I don't know if this must be here or in tools
    
    """
    
    Load the csv storing the Spectral Signatures
    
    """
    
    SS = pd.read_csv(os.path.join(path, filename), converters={'signature': read_list})
    
    return SS

def read_list(s):
    
    return [float(x) for x in s[1:-1].replace(' ', '').split(',')]


def interface_pix(directories, plot_interface_pix, plot_mean_ss):
    
    """
    
    This function detects the interface pixels
    
    inputs: 
        
        * directories: dictionary with the structure folders path
        * plot_interface_pix: Boolean. True if the interface pix position in the Cropped PanSharpened img
        * plot_mean_ss: boolean. True if the mean and derivative reflectances per profile are desidered
    
    
    """
    
    folder_dates = [f for f in os.listdir(directories['SS_pix_csv'])]
    
    for folder in folder_dates:
        
        date = folder.split('_')[0]
        
        folder_date = os.path.join(directories['SS_pix_csv'], folder)
        
        print('1',folder_date)

        # Load the csv
        
        cs_ss = [f for f in os.listdir(folder_date)][0]
        
        print('2', cs_ss)
        
        SS = ReadSS(folder_date, cs_ss)
        
        # Obtain also the CW
        
        name_or_img =  [f for f in os.listdir(directories['scene']) if f.endswith('he5')][0]
        
        #name_or_img = 'PRS_L2D_STD_' + folder + '_0001.he5'
        print('3', folder)
        
        
        CW = PRIS_img_v2.ObtainCW(directories ,name_or_img)
    
        # Group the pixels by profile
        
        groups = SS.groupby('pr')
        
        interface_pix = pd.DataFrame({'pr':[], 'row':[], 'col':[], 'signature':[]})
        
        # OBTAIN THE INTERFACE PIXELS
        
            # Obtain the mean reflectance for each pixel in the different profiles
            
        pathS = os.path.join(directories['Interface_pix_plot'], 'MeanSS_pr', folder)
        
        if not os.path.exists(pathS):
            os.makedirs(pathS)
        
        for idd, group  in groups:
                        
            
            mean_ss = []
            
            for ii, rr in group.iterrows():
                mean_ss.append(np.mean(rr['signature']))
                
            mean_ss_old = np.array(mean_ss)
            
            ### ESTA ES LA QUE TENGO QUE INTERPOLAR
            

            
            lenss = len(mean_ss)
            #xold = np.linspace(0,lenss, num=lenss)
            xold = np.arange(0, lenss, 1)
            f = interpolate.interp1d(xold, mean_ss, kind = 'cubic')
            
            #xnew = np.linspace(0,lenss, num=lenss*4)
            xnew = np.arange(0, lenss-1,1/4)
            mean_ss_new = f(xnew)
                
            # Obtain the derivative 
            
            deriv_old = mean_ss_old[1:] - mean_ss_old[:-1]
            deriv_new = mean_ss_new[1:] - mean_ss_new[:-1]
            
            
            # Obtain the minimum value of the derivative (where the reflectance decrease faster)
            
            idd_inter_old = np.argmin(deriv_old)
            idd_inter_new = np.argmin(deriv_new)

            inter = mean_ss_new[idd_inter_new]
            
            interface_pix = interface_pix.append(group.iloc[idd_inter_old])
            
            interface_pix.reset_index(inplace = True, drop = True)
            
            
            if plot_mean_ss:
                
            
                    
                plt.figure(figsize = (12,6))
                plt.scatter(idd_inter_new, inter, s=150, c='red', alpha = 0.5, label = 'interface')
                plt.plot(mean_ss_new, 'ko', label = 'mean_reflectance')
                plt.plot(deriv_new, 'bo', label = 'derived')
                
                plt.xticks(np.arange(0, len(mean_ss_new), step=1)); # mean_diff
                plt.grid('major')
                plt.legend(loc='best')
                plt.xlabel('Pixel number')
                plt.ylabel('Mean reflectance')
                plt.title(idd, fontsize = 12)
                name = 'MeanSS_'+ idd + '_' + folder + '_new.png'
                
                plt.savefig(os.path.join(pathS, name))
                
                plt.close()
            
        
        # SAVE THE interface pixel image position
        
        name_csv = 'InterfacePix_' + folder.split('_')[0] + '.csv'
        
        interface_pix.to_csv(os.path.join(directories['Interface_pix_csv'], name_csv), index = False)
        
        
        if plot_interface_pix:
            
            # PLOT THE INTERFACE PIX ON THE PANSHARPENED RGB IMAGE
            
                # Load the image
                
            panSharp_img = [f for f in os.listdir(directories['PanSharp_Square_Cropped']) if date in f and f.endswith('tiff') ][0]
            
                
            #panSharp_img = 'PanSharpened_' + folder + '_0001.tiff'
            pansharp_img_path = os.path.join(directories['PanSharp_Square_Cropped'], panSharp_img)
            
            
            
            raster, raster_rgb = PRIS_img_v2.numpy_from_tiff(pansharp_img_path, True)
            
                # Color the interface pix
                
            for idd, rrow in interface_pix.iterrows():
                
                row = int(rrow['row'])
                col = int(rrow['col'])
                raster_rgb[row, col, :] = [1,0,0]
                
                # Genrate the image
                
            plt.figure(figsize = (12,12))
            plt.imshow(raster_rgb)
            plt.title('Interface Pixels ' + folder)
            plt.axis('off')

            
                # Store it
                
            img_name = 'Interface_pix_' + folder + '.png'
                
            pathS = os.path.join(directories['Interface_pix_plot'], 'Interface_img')
            
            if not os.path.exists(pathS):
                os.makedirs(pathS)
                
            plt.savefig(os.path.join(pathS, img_name))
                
            plt.close()
            

def interfacePix_geojson(directories):
    
    """
    
    This function transforms the 
    
    """
    
    pix_interf_csv = [f for f in os.listdir(directories['Interface_pix_csv']) if f.endswith('.csv')]
    
    panSharp_imgs = [f for f in os.listdir(directories['PanSharp_Square_Cropped']) if f.endswith('.tiff')]
    

    
    
    for file in pix_interf_csv:
        
        search = file.split('_')[1].split('.')[0]
        
        for pan_img in panSharp_imgs:
            
            if search in pan_img:
                
                # Load the image needed to retrieve the necessary information to transform from pixel coordinates
                # into geographic coordinates
                
                
                ds = gdal.Open(os.path.join(directories['PanSharp_Square_Cropped'],pan_img))

                # Load the csv file
                        
                csv_info = ReadSS(directories['Interface_pix_csv'], file)
                
                # Transfor the pixel coordinates into geographic coordinates
                
                x = []
                y = []
                
                col = []
                row = []
                
                for idd, rrow in csv_info.iterrows():
                    
                    xx, yy = PRIS_tools_v2.pixel_to_world(ds, rrow['col'], rrow['row'])
                    
                    col.append(rrow['col'])
                    row.append(rrow['row'])
                    
                    x.append(xx)
                    y.append(yy)
                    
        
                
                # Save into the geojson file
                
                geojson = PRIS_tools_v2.coor_to_geojson(x, y)
                
                # Save the geojson
                
                filename = 'SDS_derivative_' + file.split('_')[1].split('.')[0] + '.geojson'
                
                path = os.path.join(directories['SDS'], 'derivative')
                
                if not os.path.exists(path):
                
                    os.makedirs(path)
                
                PRIS_tools_v2.save_geojson(path, filename, geojson)
                
                # Save col and row in txt format
                

                filename_txt = 'SDS_derivative_' + file.split('_')[1].split('.')[0] + '.txt'
                       
                np.savetxt(os.path.join(path, filename_txt), (col,row), delimiter=' ', newline=os.linesep)
        
        
        
###############################################################################
######################### OBTAIN SENMEMBERS  ##################################
###############################################################################

def obtain_endmembers(directories, filename): # hay que añadir un chek a los scores
    
    """
    
    With this function we obtain the PanSharp beach cropped images 
    to apply the LSMA. A k-mean is applied to the image ss and the compared to 
    the 
    
    inputs:
        
        * directories: dict
        * file: str; filename storing the interface pix spectral signatures
        * plot: boolean; True to show the plot of the obtained ss ¿????
        
    outputs:
        
        * scores_BrayCrutis: pandas DataFrame with the scores obtained from comparing the endmembers
                             obtained from the k-mean applied to the image and the enmembers obtained from 
                             the images
                             
        * indx: np.array classes corrispondences. 0: SAND, 1: interface, 2: WATER
        
        * ff: np.array representing the endmembers
    
    """
    
    ##### SPECTRAL SIGNATURES INTERFACE PIX
    
    
    # Load the file
    
    ss_interface = ReadSS(directories['Interface_pix_csv'], filename)
    
    #### PROFILES SS 
    
    # Load the file
    
    # Load the Profiles SS
    
    folder_name = [f for f in os.listdir(directories['SS_pix_csv']) if filename.split('_')[1].split('.')[0] in f][0]
    #folder_name = filename.split('_')[1]+ '_' + filename.split('_')[2].split('.')[0]
    file_ss = [f for f in os.listdir(os.path.join(directories['SS_pix_csv'], folder_name))][0]
    
    # Load the ss

    SS = ReadSS(os.path.join(directories['SS_pix_csv'], folder_name), file_ss)
    
    # OBTAIN THE CW
    
    original_img =  [f for f in os.listdir(directories['scene']) if f.endswith('he5')][0]
    #original_img = 'PRS_L2D_STD_' + folder_name + '_0001.he5'
    
    CW = PRIS_img_v2.ObtainCW(directories, original_img)
    
    #### OTAIN THE SS FROM PROFILES AND INTERFACE PIXELS (it is a way to assign the correct class to the k-means ss)
    
    # Otain the Mena_SS_interface
    
    signatures = np.full([len(CW),len(ss_interface)], np.nan)

    o = 0
    
    for idd, rrow in ss_interface.iterrows():
       
       
       signatures[:,o] = rrow['signature']
       o = o+1

    Mean_SS_interface = np.nanmean(signatures, axis=1)
    
   # Obtain the mean ss for water and sand
   
    groups = SS.groupby(['pr'])
   
   # For SAND
   
    signatures_sand = np.full([len(CW), len(groups)], np.nan)

    
    o = 0
    
    for idd, group in groups:
        group.reset_index(inplace=True, drop=True)
        
        if len(group)<10: #### THIS COULD BE AVOIDED BY REMOVING THOSE PROFILES FROMM THE BEGINNING!!!!!!!
            continue
        
        signatures_sand[:,o] = group['signature'][1] # siempre el primer pixel
        
        o += 1
        
    Mean_SS_sand = np.nanmean(signatures_sand, axis=1)
    
    
    # For water
   
    signatures_water = np.full([len(CW), len(groups)], np.nan)
    
    o = 0
    for idd, group in groups:
        group.reset_index(inplace=True, drop=True)
        if len(group)<10:
            continue
        signatures_water[:,o] = group['signature'][group.index[-1]] # siempre el ultimo pixel
        o += 1
    
    Mean_SS_water = np.nanmean(signatures_water, axis=1)
    
    # Genrate the array with the ss sigantures retrived from the images
    
    ems = np.array([Mean_SS_sand, Mean_SS_interface, Mean_SS_water])
    
    
    #### THE SS FROM THE PanSharp BEACH CROPPED IMAGE USING K-MEANS 
    
    # Load the image
    
    pansharp_img = [f for f in os.listdir(directories['PanSharp_Beach_Cropped']) if f.endswith('tiff') and folder_name.split('_')[0] in f][0]
    #pansharp_img = 'BeachCrop__' + folder_name + '_0001.tiff'
    
    raster_beach, raster_beach_rgb = PRIS_img_v2.numpy_from_tiff(os.path.join(directories['PanSharp_Beach_Cropped'], pansharp_img), True)
    
    (m, ff) = kmeans(raster_beach, 4, 100);
    
    # COMPARE THE ENDMEMBERS
    
   # Compare ems with codebook
    
    dummy_array = np.full([ems.shape[0], ff.shape[0]], np.nan) # ems rows, ff col
    scores_BrayCurtis = pd.DataFrame(dummy_array)


    for i in range(ems.shape[0]):

        for ii in range(ff.shape[0]):

            scores_BrayCurtis.iloc[i,ii] = PRIS_tools_v2.Bray_Curtis(ems[i,:], ff[ii,:])
    
            # (1 - pprA.Bray_Curtis(ems[i,:], codebook[ii,:]))*100
        
    # Once the scores have been calculated  we can check how similar are the endmembers obtained from the images
    # with the control ones (the mean ones)
    

    if  np.min(scores_BrayCurtis.max(axis=0, skipna=True).values[1:]) < 80: # Este volor se debería incrementar
        
        
        #raise Exception('The endmembers are not similar to the reference ones')
        print('Warning, endmember: ', np.argmin(np.min(scores_BrayCurtis.max(axis=0, skipna=True).values[1:])), 'is not too similar')
        
    
        
    indx = scores_BrayCurtis.idxmax(axis=0, skipna=True).values
    
    # Con esta informacion tenemos que encontrar cuales son las firmas espectrales por clase

    # 1 interface
    # 0 sand
    # 2 water
    
    return scores_BrayCurtis, indx, ff, raster_beach
    
###############################################################################
############################## UNMIXING  ##################################
############################################################################### 


def unmixing_FCLS(raster_beach, ems): ## Here the mask could be added the problem is with the SAM
    
    """
    
    Apply the Full contrained Last Squares unmixing technique
    
    
    """
    
    am = amap.FCLS() #define am object using the amap 
    amaps = am.map(raster_beach,ems,normalize=False)
    
    return amaps

    


###############################################################################
############## SPATIAL ATRACTION MODEL, SAM ##################################
############################################################################### 


def ReturnNeighbourhoodPixels2(i, j, grid):

    # All possible displacements.
    deltas = [
        [1, 1],
        [1, 0],
        [1, -1],
        [0, 1],
        [0, -1],
        [-1, 1],
        [-1, 0],
        [-1, -1]
    ]
    
    if len(grid.shape) == 3:
        rows, cols, _ = grid.shape
    if len(grid.shape) == 2:
        rows, cols = grid.shape

    # All possible neighbours.
    ns = ((i + dy, j + dx) for dy, dx in deltas)
    return [
        (i, j) for i, j in ns
        # Remove neighbours outside the grid.
        if rows > i >= 0 and cols > j >= 0
    ]



def Calculate_Aij(S, R ,C, num_classes, fractionsMap):
    """
    This function allows to calculate the number of available sub-pixel for class c in pixel Pij
    This sentence in the article is not too clear:
        
        The number of sub-pixels for every class is calculated and limited to integer values only: 
            the remainder of the dub-pixels not being assigned because of the rounding operator is 
            assigned to the classes with the highest rest fractions
    
    Inputs:
        
        S = Scale factor
        R = Number of rows
        C = Nuber of columns
        num_classes : number of classes in the soft classification

    Outputs:
        
        Aij = Number of available sub-pxels per class (np.array of ints with the same size of the original image)
        Aij_2 =  Result of the calculation without the rounding
    
    """
    
    # Inizializate the arrays wich will contain the number of sub-pixels per class
    
    Aij = fractionsMap.copy()
    Aij_2 = fractionsMap.copy()
    

    # Perform the calculation
    
    for i in range(R):
        for j in range(C):
            for c in range(num_classes):
                #Aij[i,j,c] = np.floor(fractionsMap[i,j,c]**2 * (S**2)) # quitar el elevado 1
                Aij[i,j,c] = np.floor(fractionsMap[i,j,c] * (S**2))
                Aij_2[i,j,c] = fractionsMap[i,j,c]*(S**2)
                
    
     # Adjust Aij(c) with the remainder of sub-pixels (those not been assigned because of the round
     # operator, is assigned to the classes with the highest rest fractions)   

    for i in range(R):
        for j in range(C):
            
            Num_subPixel_toadd = S**2 - np.sum(Aij[i,j])
            if Num_subPixel_toadd == 0:
                continue
            else:

            # Aggiungiamo alla classe più alta
            
            
                ind_increasing_frac = np.argsort(Aij_2[i,j]) # returns the indexes of the Aij[i,j] values in increasing order
                band_to_encrease = ind_increasing_frac[-1]
                Aij[i,j, band_to_encrease] = Aij[i,j, band_to_encrease] + Num_subPixel_toadd
                
                #for hh in range(int(Num_subPixel_toadd)):
                    #band_to_encrease = ind_increasing_frac[-1-(hh)] # if you want to add in the highest after the one with more subpixels (hh+1)
                    #Aij[i,j, band_to_encrease] = Aij[i,j, band_to_encrease] + 1
                    
                    
    # There is a third interpretation of the phrase remainder...:
        # To sum the remaining sub-pixels to those fractions without any subpixel available in encreasing order
            # ind_increasing_fac = np.argsort(Aij_2[0, 5] - np.floor(Aij_2[0, 5]))[::-1]
            # band_to_increase = ind_increasing_frac[hh]
                    
                    
    if np.unique(np.sum(Aij,axis=2))[0] == S**2:
        print('All the available subpixels per pixels sum: ', np.unique(np.sum(Aij,axis=2))[0])
                    
    return Aij, Aij_2


def Calculate_pab(S, t, B, R, C, num_classes ,fractionsMap):
    """
    
    This function calculates the SPATIAL ATTRACTION based on euclidian distances 
    between subpixels and neighbouring pixels
    
    Inputs :
        
        S = Scale factor
        R = Number of rows
        C = Number of columns
        fractionMap = Result from soft classification
        
    Output:
        
        p_ab = Matrix with the spatial attarction values per pixel an per size [np.array]
        # Size = [original_image_x*S, original_image_y*S, num_classes]
    
    """
        
    Nt = [] # For store the pixels Pij neighborhoods of type t of su-pixel pab
    
    p_ab = np.full([R*S, C*S, num_classes], np.nan) # p_ab matrix

    from collections import defaultdict
    n_count = defaultdict(int)

    for i in range(R):
        for j in range(C):
            
            # neighborhood pixels (Pkm) of the pixel which contains the actual sub-pixels
            
            Neigh =  ReturnNeighbourhoodPixels2(i, j, fractionsMap) # 

            for a in range(S*i, S*(i+1)):
                for b in range(S*j, S*(j+1)):
                    
                    k = [] 
                    m = []
                    dist = []
                                        
                    
                    for n in Neigh:
                        
                        
                        dy = a + 0.5 - S * (n[0] + 0.5)
                        dx = b + 0.5 - S * (n[1] + 0.5)
                        d = np.sqrt(dx**2 + dy**2)
                        
                        if t == 1:
                            
                            if d <= 1/np.sqrt(2)*(S+1):
                                
                                k.append(n[0])
                                m.append(n[1])
                                dist.append(d)

                        if t == 2:
                            
                            if d <= 1/np.sqrt(2)*(2*S-1):
                                
                                k.append(n[0])
                                m.append(n[1])
                                dist.append(d)
                                
                        if t == 3:
                            
                            if ((n[0]/S != i) or (n[1]/S != j)) and ((np.abs(n[0]/S - i)<=1) and (np.abs(n[1]/S-j) <= 1)):
                                k.append(n[0])
                                m.append(n[1])
                                dist.append(d)
                                
                    n_count[len(k)] += 1
                    P_km_c = fractionsMap[k,m,:]
                    
                    # ll = Number of neighbouring pixels
                    # The mean attarction exercised by the neighbouring pixels (Pkm) oevr 
                    # the sub-pixel pab
                    
                    P_km_c2 = P_km_c.copy()
                    
                    for ll in range(P_km_c.shape[0]):
                        P_km_c2[ll,:] = P_km_c[ll,:]/dist[ll]
                        
                    aver = np.sum(P_km_c2, axis = 0)/P_km_c2.shape[0]
                    
                    for c in range(num_classes):
                        p_ab[a,b,c] = aver[c]
                        
                        
            if B:
                
                for c in range(num_classes):
                    
                    to_sum = []
                    
                    for a in range(S*i, S*(i+1)):
                        for b in range(S*j, S*(j+1)):
                            
                            to_sum.append(p_ab[a,b,c])
                            somma = np.nansum(to_sum)
                            
                    for a in range(S*i, S*(i+1)):
                        for b in range(S*j, S*(j+1)):
                            p_ab[a,b,c] = p_ab[a,b,c]/somma
            
    p_ab[np.isnan(np.float64(p_ab))] = 0
    
    IndMax_p_ab = np.argmax(p_ab, axis=2)
    print(n_count)
    
    return p_ab, IndMax_p_ab


def Reordering_SubPixels(p_ab, Aij, R, C, S, num_classes):
    
    p_ab_2 = np.full([R*S, C*S], np.nan) # p_ab finale
    
    for i in range(R):
        for j in range(C):
            #print(i,j)
            # TO SOLVE THE NON ATTRACTED SUBPIXELS (p_ab == 0) from the beeggining because can cause problems
            r0,c0 = np.where(np.nansum(p_ab[S*i:S*(i+1), S*j:S*(j+1)], axis=2) == 0)
            gg = np.where(np.nansum(p_ab[S*i:S*(i+1),S*j:S*(j+1)],axis=2) != 0)
            for rr0, cc0 in zip(r0,c0):
                for r1,c1 in zip(gg[0], gg[1]):
                    if r1 == rr0:
                        p_ab[S*i:S*(i+1),S*j:S*(j+1)][r0,c0,:] = p_ab[S*i:S*(i+1),S*j:S*(j+1)][r1,c1,:]
                        break
                        

                    
            while np.sum(Aij[i,j,:]) != 0:
                sub_pixels = p_ab[S*i:S*(i+1), S*j:S*(j+1),:]
                highest = np.where(sub_pixels==np.nanmax(sub_pixels)) 
                a = highest[0][0]
                b = highest[1][0]
                d = highest[2][0]
                #classes = list(np.linspace(0,num_classes-1,num_classes))    
                
                if Aij[i,j,d] > 0:
                    p_ab_2[S*i +a, S*j +b] = d
                    Aij[i,j,d] -= 1
                    # for c in range(num_classes):
                    #     p_ab[S*i +a, S*j + b, c] = 0
                    p_ab[S*i + a, S*j + b, :] = 0
        
                else:
                    #p_ab_2[S*i + a, S*j + b] = 0
                    p_ab[S*i + a, S*j + b, d] = 0
                    
    return p_ab_2


def Spatial_Attraction_Model(fractionsMap, S, t, B):
    """
    This is the MAIN function
    
    Inputs:
        
        S = Scale factor (int). It can be 2,3,4. Higher S increase the number of sup-pixels quadratically
        and make the sub-pixel mapping problem more complicated. RESULTS IN LESS FAVOURABLE ACCURACY VALUES
        
        fractionsMap = Result of the soft classification [np.array with the same size of the original image]
        
        t = Describe the type of neighbouring [int]. t=1 Touching; t=2 Quadrant; t=3 Surrounding (still to implement)
        
        Note: With S=2 there is no difference between touching and quadrante. Quadrant better results
        
        B = Normalization boolean (True/False)
        
    Output:
        
        p_ab_2 : Hard classification (1 image) containing pab
        
    
    """
    print('Iniziamo la COSA')
    # FractionMao
    
    fractionsMap = fractionsMap
    
    
    # Shapes
    
    R = fractionsMap.shape[0]
    C = fractionsMap.shape[1]
    num_classes = fractionsMap.shape[2]
    
    
    # Calculate the number of the available sub-pixels for class c in pixel Pij
    
    Aij, Aij_2 = Calculate_Aij(S, R, C, num_classes, fractionsMap)
    
    # Calculate p_ab
    
    p_ab, IndMax_p_ab = Calculate_pab(S, t, B, R, C, num_classes, fractionsMap)
    
    
    # Reassign the pixels 
    
    p_ab_2 = Reordering_SubPixels(p_ab, Aij, R, C, S, num_classes)
    
    return p_ab_2, IndMax_p_ab, Aij_2, p_ab


###############################################################################
########################## EDGES SELECTION ##################################
############################################################################### 

def edges_selection(p_ab_2, edges_sobel):
    
    """
    This function calculates the distance of the edges to the nearest water pixels
    
    Input:
        
        * p_ab_2: np.array. Class map resulting from apply the SAM
        * edges_soble: np.array. Detected borders
        
    Output:
        
        * edges_sobel_near_water: np.array. Edges nearest to water
    
    """
    
    # 1) SELECT THE WATER PIXELS
    
    Wat_ind = np.where(p_ab_2 == 1)
    coor_wat = pd.DataFrame({'row':Wat_ind[0], 'col':Wat_ind[1]})
    
    # 2) BUFFER AROUND THE EDGES 
    # Credo non sia necessario considerari gli oblicui
    
    buffer_x = []
    buffer_y = []
    
    for x_ed, y_ed in zip(edges_sobel[:,1], edges_sobel[:,0]):
        
        buffer_x.append([x_ed - i for i in range(4)])
        buffer_x.append([x_ed + i for i in range(4)])
        buffer_y.append([y_ed + i for i in range(4)])
        buffer_y.append([y_ed - i for i in range(4)])
        
        
    buffer_x = list(itertools.chain(*buffer_x))
    buffer_y = list(itertools.chain(*buffer_y))
    
    coor_bufer = pd.DataFrame({'row':buffer_y, 'col':buffer_x})
    
    # 3) BUFER PIXELS THAT ARE WATER (avoid calculate the distances to the sand pixels)
    
    water_bufer = coor_bufer.merge(coor_wat, left_on=['row','col'], right_on=['row','col'], how='inner')
    
    # 4) CALCULATE THE DISTANCES sobel-water
    

    dist_min = []
    
    for x_ed, y_ed in zip(edges_sobel[:,1], edges_sobel[:,0]):
            
        # Calculate all the distances of this edge to the pixels of the class water
        dist_min.append(np.min([np.sqrt((xx-x_ed)**2 + (yy-y_ed)**2) for xx,yy in zip(water_bufer['col'], water_bufer['row'])]))
        
        
    edges_sobel_min_dist = pd.DataFrame({'c':edges_sobel[:,1], 'r': edges_sobel[:,0], 'dist':dist_min})
    edge_sobel_near_water = edges_sobel_min_dist[edges_sobel_min_dist['dist'] < 0.5]
    
    #edges_sobel_near_water2 = edge_sobel_near_water.to_numpy()
    
    ## CLEAN ALSO POSIBLE INCONSISTENCES IN WATER
    
    edges_sobel2 = np.asarray([edge_sobel_near_water['r'], edge_sobel_near_water['c']]).T
    
    Sand_ind = np.where(p_ab_2 == 3)
    coor_sand = pd.DataFrame({'row':Sand_ind[0], 'col':Sand_ind[1]})
    
    buffer_x = []
    buffer_y = []
    
    for x_ed, y_ed in zip(edges_sobel[:,1], edges_sobel[:,0]):
        buffer_x.append([x_ed - i for i in range(15)])
        buffer_x.append([x_ed + i for i in range(15)])
        buffer_y.append([y_ed + i for i in range(15)])
        buffer_y.append([y_ed - i for i in range(15)])
        
        
    buffer_x = list(itertools.chain(*buffer_x))
    buffer_y = list(itertools.chain(*buffer_y))
    
    coor_bufer = pd.DataFrame({'row':buffer_y, 'col':buffer_x})
    
    sand_bufer = coor_bufer.merge(coor_sand, left_on=['row','col'], right_on=['row','col'], how='inner')
    
    dist_min = []

    for x_ed, y_ed in zip(edges_sobel2[:,1], edges_sobel2[:,0]):
    
        #print(x_ed,y_ed)
        
        # Calculate all the distances of this edge to the pixels of the class water
        dist_min.append(np.min([np.sqrt((xx-x_ed)**2 + (yy-y_ed)**2) for xx,yy in zip(sand_bufer['col'], sand_bufer['row'])]))
    
    edges_sobel_min_dist2 = pd.DataFrame({'c':edges_sobel2[:,1], 'r': edges_sobel2[:,0], 'dist':dist_min})

    dist_clean = edges_sobel_min_dist2['dist'].quantile(0.90) # distance that contains the 90% of the distances
    
    edge_clean = edges_sobel_min_dist2[edges_sobel_min_dist2['dist'] < dist_clean]
    edge_clean.reset_index(inplace = True, drop = True)
    
    return edge_clean


## I have to find a way to remove the outliers


###############################################################################
############################# SAVE THE SDS ##################################
############################################################################### 


def save_subpixel_SDS(directories, filename, border, S, edge_clean): #edges_sobel_near_water
    
    # 1) LOAD THE IMAGE
    
    folder_name = [f for f in os.listdir(directories['PanSharp_Square_Cropped']) if filename.split('_')[1].split('.')[0] in f][0]
    
    #folder_name = filename.split('_')[1]+ '_' + filename.split('_')[2].split('.')[0]
    
    print('El nombre del folder es:', folder_name)
    
    #pansharp_img = [f for f in os.listdir(directories['PanSharp_Square_Cropped']) if  folder_name.split('_')[0] in f][0]
    
    pansharp_img = [f for f in os.listdir(directories['PanSharp_Beach_Cropped']) if f.startswith('BeachCrop') and folder_name.split('_')[1] in f][0]
    
    print('El nombre de la imagen cropped a la beach: ', pansharp_img)
    
    #pansharp_img = 'BeachCrop__' + folder_name + '_0001.tiff'

    #raster_gdal  = gdal.Open(os.path.join(directories['PanSharp_Square_Cropped'], pansharp_img))
    raster_gdal  = gdal.Open(os.path.join(directories['PanSharp_Beach_Cropped'], pansharp_img))
    gt = raster_gdal.GetGeoTransform()
    wkt = raster_gdal.GetProjection()
    
    
    # 2) RESOLUTION OF THE RESIZED IMAGE
    
    delta_lon_sub = (5/S)
    delta_lat_sub = -(5/S)
    
    # 3) GENERATE THE NEW GeoT
    
    gt2 = [gt[0], delta_lon_sub, gt[2], gt[3], gt[4], delta_lat_sub]
    
    # 4) OBTAIN THE WORLD COORDINATES
    
    xx = []; yy = []

    for idd, rrow in edge_clean.iterrows():
        
        x, y = PRIS_tools_v2.pixel_to_world_resized(raster_gdal, gt2, rrow['c'], rrow['r'])
        
        xx.append(x)
        
        yy.append(y)
        
    cols = (edge_clean['c']/S).to_numpy()
    rows = (edge_clean['r']/S).to_numpy()
        
    # 5) GENERATE THE GEOJSON
    
    geoson = PRIS_tools_v2.coor_to_geojson(xx,yy)
    
    # 6) SAVE THE INFO
    

    
    path_store = os.path.join(directories['SDS'], 'k-means')
    path_store_txt = os.path.join(path_store, 'txt')
    
    if not os.path.exists(path_store):
        os.makedirs(path_store)
        
    if not os.path.exists(path_store_txt):
        os.makedirs(path_store_txt)
        
    print(folder_name)
        
    name_sds = folder_name.split('.')[0] +'_' + border +'_kmeans' + '.geojson'
    name_sds_txt = folder_name.split('.')[0] +'_' + border + '_kmeans' + '.txt'
    
    fileout = open(os.path.join(path_store_txt, name_sds_txt), "w")
    
    
    for pos in range(len(cols)):
        fileout.write('{:9.3f} {:9.3f}\n'.format(cols[pos], rows[pos]))
    fileout.close()
    
        
    #name_sds2 = 'SubPix_SDS_' + folder_name + '.txt'
    
    #np.savetxt(os.path.join(path_store, name_sds2), edges_sobel_near_water, fmt='%.18g', delimiter=' ', newline=os.linesep)

    
    PRIS_tools_v2.save_geojson(path_store, name_sds, geoson)
    
    

###############################################################################
############################# OBTAIN THE SDS ##################################
############################################################################### 

def obtain_SDS(directories):
    
    """
    
    This function calls the others functions to obtain the SDS
    
    """
    
    
    interface_csv = [f for f in os.listdir(directories['Interface_pix_csv']) if f[-3:] == 'csv']
    
    tot_img = len(interface_csv)
    
    u = 0
    
    for file in interface_csv:
    
        # 1) OBTAIN THE ENDMEMBERS
        
        print(file)
        
        scores_Bray, indx, ems, raster_beach = obtain_endmembers(directories, file)
        
        print(scores_Bray)
        
        print(indx)
        
        # Resize the image
        
        #subset_resized = resize(raster_beach, (raster_beach.shape[0]*S, raster_beach.shape[1]*S))
        
        # 2) UNMIXING TECHNIQUE
        
        amaps = unmixing_FCLS(raster_beach, ems)
        
        # 3 ) ASPATIAL ATTRACTION MODEL
        
        # 3.1) Define the parameters that will be used in SAM
        
        S = 4 # How many pixels are formed inside each original pixel = S^2
        t = 2 # type of neighbourhood
        B = False # Normalization boolean (no se hacer la normalizacion)
        number_classes = amaps.shape[2]
           
        # 3.2) Apply the SAM
        
        p_ab_2, IndMax_p_ab, Aij, p_ab = Spatial_Attraction_Model(amaps, S, t, B)
        
        # 4) CLASSES BETWEEN THE BORDERS MUST BE FOUND
        
        # We mantain the classes water and interface, index of the water and interface classes
        
        inter_class0 = np.where(indx == 1)[0]# esta es la clase interface [0]
        water_class0 = np.where(indx == 2)[0]
        
        inter_class0 = np.where(indx == 1)[0]# esta es la clase interface [0]
        water_class0 = np.where(indx == 2)[0]
        
        score_values = []
        
        if indx[1] != 2 or indx[2]!=1 or indx[3]!=0:
            delete = [0, 3]
        else:
        
            if len(water_class0)>1:
                for idx in water_class0:
                    score_values = scores_Bray.iloc[2, idx]
                row_max, col_max = np.where(scores_Bray == score_values)
                water_class = np.where(indx == 2)[0][col_max[0]]
            if len(water_class0)==1:
                water_class = np.where(indx == 2)[0]
        
            score_values = []    
        
            if len(inter_class0)>1:
                for idx in inter_class0:
                    score_values = scores_Bray.iloc[1, idx]
                row_max, col_max = np.where(scores_Bray == score_values)
                inter_class = np.where(indx == 1)[0][col_max[0]]
            if len(inter_class0)==1:
                inter_class = np.where(indx == 1)[0]
       
            # Index of the classes that must be deleted
            
            delete = [f for f in range(0,4) if (f != inter_class) and (f != water_class) ]
        
        P_AB = p_ab_2.copy()
        
        # Mask the class thumbnail and sand
        
        for f in delete:
            
            P_AB[P_AB == f] = np.nan
            
        # 5) OBTAIN THE BORDER WITH SOBEL
        
        sx = ndimage.sobel(P_AB, axis=0, mode='constant')
        sy = ndimage.sobel(P_AB, axis=1, mode='constant')
        sob = np.hypot(sx, sy)
        
        sob[sob == 0] = np.nan
        
        edges_sobel = np.argwhere(~np.isnan(sob))
                                  
        # 6) SELECT ONLY THE POINTS NEAR TO THE WATER CLASS 
        
        edge_clean = edges_selection(p_ab_2, edges_sobel)

        # 7) SAVE THE SDS
        
        # 7.1) Geojson Format
            
        save_subpixel_SDS(directories, file, border, S, edge_clean)
        
        print('Saved Subpixel SDS: ', (u/tot_img)*100, '%', end='\r')
        
        
        u = u +1
        
        
        