# Toolset to calibrate change attribution algorithm

import math
import numpy as np
import scipy.stats as stats
import scipy.ndimage
from scipy.ndimage.measurements import label

import biota.indices

import matplotlib.pyplot as plt
import pdb



def getChange(data_t1, data_t2, F_threshold = 10., C_threshold = 0.2, min_area = 0, output = False):
    """
    Returns pixels that meet change detection thresholds for country.
    
    Args:
        F_threshold = threshold above which a pixel is forest
        C_threshold = threshold of proportional change with which to begin change detection
    """
    
    from osgeo import gdal
    
    #TODO: Test that data_t1 and data_t2 are from the same location, with the same extent
    #TODO: Combine masks from data_t1 and data_t2
    
    AGB_t1 = data_t1.getAGB()
    AGB_t2 = data_t2.getAGB()
    
    # Pixels that move from forest to nonforest or vice versa
    F_NF = np.logical_and(AGB_t1 >= F_threshold, AGB_t2 < F_threshold)
    NF_F = np.logical_and(AGB_t1 < F_threshold, AGB_t2 >= F_threshold)
    F_F = np.logical_and(AGB_t1 >= F_threshold, AGB_t2 >= F_threshold)
    NF_NF = np.logical_and(AGB_t1 < F_threshold, AGB_t2 < F_threshold)
    
    # Change pixels
    CHANGE = np.logical_or(((AGB_t1 - AGB_t2) / AGB_t1) >= C_threshold, ((AGB_t1 - AGB_t2) / AGB_t1) < (- C_threshold))
    NOCHANGE = CHANGE == False
    
    # Trajectory
    DECREASE = AGB_t2 < AGB_t1
    INCREASE = AGB_t2 >= AGB_t1
    
    # Get a minimum pixel extent. Loss/Gain events much have a spatial extent greater than min_pixels, and not occur in nonforest.
    if min_area > 0:
        
        pixel_area = data_t1.xRes * data_t1.yRes * 0.0001
        min_pixels = int(round(min_area / pixel_area))
        
        CHANGE_INCREASE, _ = biota.indices.getContiguousAreas(CHANGE & INCREASE & (NF_F | F_F), True, min_pixels = min_pixels)
        CHANGE_DECREASE, _ = biota.indices.getContiguousAreas(CHANGE & DECREASE & (F_NF | F_F), True, min_pixels = min_pixels)
        CHANGE = np.logical_or(CHANGE_INCREASE, CHANGE_DECREASE)
        NOCHANGE = CHANGE == False
    
    metrics = {}
    
    # These are all the possible definitions. Note: they *must* add up to one.
    metrics['deforestation'] = (F_NF & CHANGE) * 1
    metrics['degradation'] = (F_F & CHANGE & DECREASE) * 1
    metrics['minorloss'] = ((F_F | F_NF) & NOCHANGE & DECREASE) * 1

    metrics['aforestation'] = (NF_F & CHANGE) * 1
    metrics['growth'] = (F_F & CHANGE & INCREASE) * 1
    metrics['minorgain'] = ((F_F | NF_F) & NOCHANGE & INCREASE) * 1
            
    metrics['notforest'] = (NF_NF * 1)
    
    if output:
        output_im = np.zeros_like(AGB_t1) + data_t1.nodata
        output_im[metrics['notforest'].data == 1] = 0
        output_im[metrics['deforestation'].data == 1] = 1
        output_im[metrics['degradation'].data == 1] = 2
        output_im[metrics['minorloss'].data == 1] = 3
        output_im[metrics['minorgain'].data == 1] = 4
        output_im[metrics['growth'].data == 1] = 5
        output_im[metrics['aforestation'].data == 1] = 6
        
        output_im[np.logical_or(data_t1.mask, data_t2.mask)] = data_t1.nodata
        
        biota.IO.outputGeoTiff(output_im, data_t1.output_pattern%('CHANGE', '%s_%s'%(str(data_t1.year), str(data_t2.year))), data_t1.geo_t, data_t1.proj, output_dir = data_t1.output_dir, dtype = gdal.GDT_Int32, nodata = data_t1.nodata)
    
    return metrics


def getChangeDownsampled(data_t1, data_t2, shrink_factor = 45, output = True):
    """
    Downsample getChange to address false positives.
    
    This is a palceholder for something more advanced.
    """
    
    print "WARNING: This function is not yet functional. Come back later!"
    
    from osgeo import gdal
    
    metrics = getChange(data_t1, data_t2, F_threshold = 15., C_threshold = 0.25, min_area = 2.)
    
    # Calculate output output size based on reduction shrink_factor
    output_size = int(round(((data_t1.xSize + data_t1.ySize) / 2.) / shrink_factor))
    
    change_downsampled = np.zeros((output_size, output_size)) + data_t1.nodata
        
    for n, ymin in enumerate(np.linspace(0, data_t1.ySize, output_size, endpoint = False)):
        ymax = ymin + (data_t1.ySize / output_size)
        for m, xmin in enumerate(np.linspace(0, data_t1.xSize, output_size, endpoint = False)):
            xmax = xmin + (data_t1.xSize / output_size)
            
            deforestation = metrics['deforestation'][int(round(ymin)):int(round(ymax)), int(round(xmin)):int(round(xmax))]
            aforestation = metrics['aforestation'][int(round(ymin)):int(round(ymax)), int(round(xmin)):int(round(xmax))]
            this_mask = np.logical_or(data_t1.mask, data_t2.mask)[int(round(ymin)):int(round(ymax)), int(round(xmin)):int(round(xmax))]
            
            change_downsampled[n,m] = aforestation.sum() - deforestation.sum()
        
    if output: biota.IO.outputGeoTiff(change_downsampled, data_t1.output_pattern%('CHANGEDOWNSAMPLED', '%s_%s'%(str(data_t1.year), str(data_t2.year))), data_t1.shrinkGeoT(shrink_factor), data_t1.proj, output_dir = data_t1.output_dir, dtype = gdal.GDT_Float32, nodata = data_t1.nodata)
