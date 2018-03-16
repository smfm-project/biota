# Toolset to calibrate change attribution algorithm

import math
import numpy as np
import scipy.stats as stats
import scipy.ndimage
from scipy.ndimage.measurements import label

import biota.indices

import matplotlib.pyplot as plt
import pdb



def getChange(data_t1, data_t2, forest_threshold = 10., intensity_threshold = 0.2, area_threshold = 0, output = False):
    """
    Returns pixels that meet change detection thresholds for country.
    
    Args:
        forest_threshold = threshold above which a pixel is forest
        intensity_threshold = threshold of proportional change with which to accept a change as real
        area_threshold
    """
    
    from osgeo import gdal
    
    AGB_t1 = data_t1.getAGB()
    AGB_t2 = data_t2.getAGB()
    
    # Combine masks, for case where they differ. TODO: Set alternative where AGB is not masked by user.
    mask = np.logical_or(AGB_t1.mask, AGB_t2.mask)
    AGB_t1.mask = mask
    AGB_t2.mask = mask
    
    # Pixels that move from forest to nonforest (F_NF) or vice versa (NF_F)
    F_NF = np.logical_and(AGB_t1 >= forest_threshold, AGB_t2 < forest_threshold)
    NF_F = np.logical_and(AGB_t1 < forest_threshold, AGB_t2 >= forest_threshold)
    
    # Pixels that remain forest (F_F) or nonforest (NF_NF)
    F_F = np.logical_and(AGB_t1 >= forest_threshold, AGB_t2 >= forest_threshold)
    NF_NF = np.logical_and(AGB_t1 < forest_threshold, AGB_t2 < forest_threshold)
    
    # Get change pixels given requirements
    CHANGE = np.logical_or(((AGB_t1 - AGB_t2) / AGB_t1) >= intensity_threshold, ((AGB_t1 - AGB_t2) / AGB_t1) < (- intensity_threshold))
    NOCHANGE = CHANGE == False
    
    # Trajectory (changes can be positive or negative)
    DECREASE = AGB_t2 < AGB_t1
    INCREASE = AGB_t2 >= AGB_t1
    
    # Get a minimum pixel extent. Loss/Gain events much have a spatial extent greater than min_pixels, and not occur in nonforest.
    if area_threshold > 0:
        
        pixel_area = data_t1.xRes * data_t1.yRes * 0.0001
        min_pixels = int(round(area_threshold / pixel_area))
        
        # Get areas of change that meet minimum area requirement
        CHANGE_INCREASE, _ = biota.indices.getContiguousAreas(CHANGE & INCREASE & (NF_F | F_F), True, min_pixels = min_pixels)
        CHANGE_DECREASE, _ = biota.indices.getContiguousAreas(CHANGE & DECREASE & (F_NF | F_F), True, min_pixels = min_pixels)
        CHANGE = np.logical_or(CHANGE_INCREASE, CHANGE_DECREASE)
        NOCHANGE = CHANGE == False
    
    metrics = {}
    
    # These are all the possible definitions. Note: they *must* add up to one.
    metrics['deforestation'] = F_NF & CHANGE
    metrics['degradation'] = F_F & CHANGE & DECREASE
    metrics['minorloss'] = (F_F | F_NF) & NOCHANGE & DECREASE
    
    metrics['aforestation'] = NF_F & CHANGE
    metrics['growth'] = F_F & CHANGE & INCREASE
    metrics['minorgain'] = (F_F | NF_F) & NOCHANGE & INCREASE
            
    metrics['nonforest'] = NF_NF
    
    if output:
        output_im = np.zeros_like(AGB_t1) + data_t1.nodata
        output_im[metrics['nonforest'].data] = 0
        output_im[metrics['deforestation'].data] = 1
        output_im[metrics['degradation'].data] = 2
        output_im[metrics['minorloss'].data] = 3
        output_im[metrics['minorgain'].data] = 4
        output_im[metrics['growth'].data] = 5
        output_im[metrics['aforestation'].data] = 6
        
        output_im[np.logical_or(data_t1.mask, data_t2.mask)] = 99
        
        biota.IO.outputGeoTiff(output_im, data_t1.output_pattern%('CHANGE', '%s_%s'%(str(data_t1.year), str(data_t2.year))), data_t1.geo_t, data_t1.proj, output_dir = data_t1.output_dir, dtype = gdal.GDT_Byte, nodata = 99)
    
    # Also return AGB in metrics (with combined masks)
    metrics['AGB_t1'] = AGB_t1
    metrics['AGB_t2'] = AGB_t2
    
    return metrics




"""

def getChangeDownsampled(data_t1, data_t2, shrink_factor = 45, output = True):
    '''
    Downsample getChange to address false positives.
    
    This is a palceholder for something more advanced.
    '''
    
    print "WARNING: This function is not yet functional. Come back later!"
        
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
"""