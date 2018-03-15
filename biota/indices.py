#!/usr/bin/env python

import numpy as np
import scipy.ndimage as ndimage

import pdb

import biota.IO

def getContiguousAreas(data, value, min_pixels = 1):
    '''
    Get pixels that come from the same contigous area.
    
    Args:
        data: A numpy array
        value: Pixel value to include in contiguous_area (e.g. True for forest)
        min_area: What minimum area should be included (number of pixels, queens move)
    
    Returns:
        A binary array of pixels that meet the conditions
    '''
    
    from scipy.ndimage.measurements import label
    
    # If masked, we use this flag to save the mask for later.   
    masked = np.ma.isMaskedArray(data)
    
    # If any pixels are masked, we give them the value of the nearest valid pixel.
    if masked:
        mask = data.mask
        ind = ndimage.distance_transform_edt(data.mask, return_distances = False, return_indices = True)
        data = data.data[tuple(ind)]
    
    # Extract area that meets condition
    binary_array = (data == value) * 1
    
    # Label contigous areas with a number
    structure = ndimage.generate_binary_structure(2,2) # This connects diagonal elements
    location_id, n_areas = label(binary_array, structure = structure)
    
    # Get count of each value in array
    label_area = np.bincount(location_id.flatten())[1:]
    
    # Find those IDs that meet minimum area requirements
    include_id = np.arange(1, n_areas + 1)[label_area >= min_pixels]
    
    # Get a binary array of location_id pixels that meet the minimum area requirement
    contiguous_area = np.in1d(location_id, include_id).reshape(data.shape).astype(np.bool)
    
    # Return an array giving values to each area
    location_id[contiguous_area == False] = 0
    
    # Re-number location_ids 1 to n, given that some unique value shave now been removed
    location_id_unique, location_id_indices = np.unique(location_id, return_inverse = True)
    location_id = np.arange(0, location_id_unique.shape[0], 1)[location_id_indices].reshape(data.shape)
       
    # Put mask back in if input was a masked array
    if masked:
        contiguous_area = np.ma.array(contiguous_area, mask = mask)
        location_id = np.ma.array(location_id, mask = mask)

    return contiguous_area, location_id


def calculateLDI(data, output = False, forest_threshold = 10., shrink_factor = 45):
    """
    Landcover Division Index (LDI) is an indicator of the degree of habitat coherence, ranging from 0 (fully contiguous) to 1 (very fragmented).
    
    LDI is defined as the probability that two randomly selected points in the landscape are situated in two different patches of the habitat (Jaeger 2000; Mcgarigal 2015).
    
    Note that fragmentation of both an undisturbed forest area and contiguous agriculture will both be low, so interpret this index carefully
    
    Args:
        shrink_factor = Number of pixels to build into a single patch.
    
    Returns:
        An array with LDI values.
    """
    
    from osgeo import gdal
    
    def _computeLDI(unique_ids):
        '''
        Calculates an index of landscape 
        '''
        
        # If masked array, separate mask
        if np.ma.isMaskedArray(unique_ids):
            mask = unique_ids.mask
            unique_ids = unique_ids.data
        else:
            mask = np.zeros_like(unique_ids, dtype=np.bool)
        
        # Calculate LDI
        selection_1 = np.zeros_like(unique_ids,dtype = np.bool).ravel()
        selection_2 = np.zeros_like(unique_ids,dtype = np.bool).ravel()
        
        # Draw 450 random pixels
        selection_1[:(selection_1.shape[0]/5)] = True
        selection_2[:(selection_2.shape[0]/5)] = True
        
        # Shuffle to get random pixels
        np.random.shuffle(selection_1)
        np.random.shuffle(selection_2)
        
        selection_1 = selection_1.reshape(unique_ids.shape)
        selection_2 = selection_2.reshape(unique_ids.shape)
        
        # This ignores IDs where one or the other of selection falls in an area of nodata
        mask_selection = np.logical_or(mask[selection_1], mask[selection_2])
        ID1 = unique_ids[selection_1][mask_selection==False]
        ID2 = unique_ids[selection_2][mask_selection==False]
        
        # Shuffle again to randomise pairs spatially
        np.random.shuffle(ID1)
        np.random.shuffle(ID2)
        
        # Number of random points from different patches
        different_patches = np.sum(ID1 != ID2)
        
        # Divide by potential matches
        if (mask_selection==False).sum() > 0:
            LDI = float(different_patches) / (mask_selection==False).sum()
            
        # Unless everything is masked, in which case return nan.
        else:
            LDI = np.nan
        
        return LDI
    
    woody_cover = data.getWoodyCover(forest_threshold = forest_threshold, min_forest_area = 0.)
        
    # Calculate output output size based on reduction shrink_factor
    output_size = int(round(((data.xSize + data.ySize) / 2.) / shrink_factor))
    
    LDI = np.zeros((output_size, output_size)) + data.nodata
            
    for n, ymin in enumerate(np.linspace(0, data.ySize, output_size, endpoint = False)):
        ymax = ymin + (data.ySize / output_size)
        for m, xmin in enumerate(np.linspace(0, data.xSize, output_size, endpoint = False)):
            xmax = xmin + (data.xSize / output_size)
            
            this_patch = woody_cover[int(round(ymin)):int(round(ymax)), int(round(xmin)):int(round(xmax))]
            this_mask = data.mask[int(round(ymin)):int(round(ymax)), int(round(xmin)):int(round(xmax))]
            
            # Get connected patches of forest and nonforest
            _, forest_id = getContiguousAreas(this_patch, True)
            _, nonforest_id = getContiguousAreas(this_patch, False)
            
            # Supply a unique ID to each habitat patch, whether forest or nonforest.
            unique_ids = forest_id.copy()
            unique_ids[nonforest_id != 0] = forest_id[nonforest_id != 0] + (nonforest_id[nonforest_id != 0] + np.max(forest_id) + 1)
            
            # Return LDI (if at least 50 % of data is present)
            if this_mask.sum() <= ((shrink_factor ** 2) * 0.5):
                LDI[n, m] = _computeLDI(unique_ids)
            
    if output: biota.IO.outputGeoTiff(LDI, data.output_pattern%('LDI', str(data.year)), data.shrinkGeoT(shrink_factor), data.proj, output_dir = data.output_dir, dtype = gdal.GDT_Float32, nodata = data.nodata)
    
    return LDI


def calculateTWC(data, shrink_factor = 45, forest_threshold = 10., output = False):
    """
    Total woody cover (TWC) describes the proportion of woody cover in a downsampled image.
    """

    from osgeo import gdal
    
    woody_cover = data.getWoodyCover(forest_threshold = forest_threshold, min_forest_area = 0.)
        
    # Calculate output output size based on reduction shrink_factor
    output_size = int(round(((data.xSize + data.ySize) / 2.) / shrink_factor))
    
    TWC = np.zeros((output_size, output_size)) + data.nodata
        
    for n, ymin in enumerate(np.linspace(0, data.ySize, output_size, endpoint = False)):
        ymax = ymin + (data.ySize / output_size)
        for m, xmin in enumerate(np.linspace(0, data.xSize, output_size, endpoint = False)):
            xmax = xmin + (data.xSize / output_size)
            
            WC = woody_cover[int(round(ymin)):int(round(ymax)), int(round(xmin)):int(round(xmax))]
            this_mask = data.mask[int(round(ymin)):int(round(ymax)), int(round(xmin)):int(round(xmax))]
            
            # Get proportion of woody cover in patch (if at least 50 % of data is present)
            if this_mask.sum() <= ((shrink_factor ** 2) * 0.5):
                TWC[n, m] = float((WC[this_mask == False]).sum()) / ((shrink_factor ** 2) - this_mask.sum())
    
    if output: biota.IO.outputGeoTiff(TWC, data.output_pattern%('TWC', str(data.year)), data.shrinkGeoT(shrink_factor), data.proj, output_dir = data.output_dir, dtype = gdal.GDT_Float32, nodata = data.nodata)
    
    return TWC


def calculateWCC(data_t1, data_t2, shrink_factor = 45, output = False):
    """
    Woody cover change (WCC) describes the loss of woody cover between two images
    """
    
    from osgeo import gdal
    
    # TODO Test that nodata values and extents are identidical for data_t1 and data_t2
    
    # Get total woody cover for time 1 and time 2
    TWC_t1 = calculateTWC(data_t1, threshold = threshold, shrink_factor = shrink_factor)
    TWC_t2 = calculateTWC(data_t2, threshold = threshold, shrink_factor = shrink_factor)
    
    # Nodata value is -1
    WCC = np.zeros_like(TWC_t1) + data_t1.nodata
    
    # Woody cover change is defined as the difference between these two
    mask = np.logical_or(TWC_t1 == data_t1.nodata, TWC_t2 == data_t2.nodata)
    WCC[mask == False] = TWC_t2[mask == False] - TWC_t1[mask == False]
    
    if output: biota.IO.outputGeoTiff(WCC, data_t1.output_pattern%('WCC', '%s_%s'%(str(data_t1.year), str(data_t2.year))), data_t1.shrinkGeoT(shrink_factor), data_t1.proj, output_dir = data_t1.output_dir, dtype = gdal.GDT_Float32, nodata = data_t1.nodata)