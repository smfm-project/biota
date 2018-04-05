#!/usr/bin/env python


# This is a set of scripts for calibration of a biomass-backscatter curve

import numpy as np
import scipy.stats

import biota
import biota.mask

import pdb



def extractGamma0(shp, plot_field, agb_field, tile_example, verbose = False, units = 'natural'):
    '''
    Extract gamma0 from ALOS tiles.
    
    Args:
        shp: A shapefile containing plot data
        plot_field: The shapefile field containing plot names
        agb_field: The shapefile field containing AGB estimates
        tile_example: An exemplar ALOS tile (from LoadTile()) which contains information on the processing chain.
        #TODO: year?
    Returns:
        A dictionary containing plot names, AGB and gamma0 values
    '''
    
    # Use example data to get processing steps
    downsample_factor = tile_example.downsample_factor
    lee_filter = tile_example.lee_filter
    year = tile_example.year
    dataloc = tile_example.dataloc
    
    # Extract relevant info from shapefile
    plot_names = biota.mask.getField(shp, plot_field) 
    agb = biota.mask.getField(shp, agb_field).astype(np.float)
    
    # Identify tiles that contain gamma0 data for the shapefile
    tiles = biota.mask.getTilesInShapefile(shp)
    
    # Check whether all the necessary tiles are all present. Decide what to do if they aren't.
    # TODO
    
    # Generate output array for backscatter
    gamma0_hv = np.empty_like(agb, dtype = np.float32)
    gamma0_hv[:] = np.nan
    
    gamma0_hh = gamma0_hv.copy()
    doy = gamma0_hv.copy()
    
    # For each tile covered by shapefile
    for lat, lon in tiles:
        
        if verbose: print 'Doing lat: %s, lon: %s'%(str(lat), str(lon))
        
        # Load tile
        tile = biota.LoadTile(dataloc, lat, lon, year, downsample_factor = downsample_factor, lee_filter = lee_filter)
        
        # Get backscatter (both polarisations) and DOY
        data_gamma0_hv = tile.getGamma0(polarisation = 'HV', units = units)
        data_gamma0_hh = tile.getGamma0(polarisation = 'HH', units = units)
        data_doy = tile.getDOY()
        
        # Mask out each plot
        plot_mask = biota.mask.rasterizeShapefile(tile, shp, location_id = True, buffer_size = 0.)
        
        # Extract values for each plot
        for n in np.unique(plot_mask[plot_mask != 0]):
            
            # Get mask for plot and tile
            this_mask = np.logical_and(plot_mask == n, tile.mask == False)
            
            # Add metrics to output array
            gamma0_hv[n-1] = np.mean(data_gamma0_hv[this_mask])
            gamma0_hh[n-1] = np.mean(data_gamma0_hh[this_mask])
            doy[n-1] = np.median(data_doy[this_mask])
    
    # Return data as a dictionary
    data_dict = {}
    data_dict['plot_name'] = plot_names
    data_dict['AGB'] = agb
    data_dict['gamma0_HH'] = gamma0_hh
    data_dict['gamma0_HV'] = gamma0_hv
    data_dict['DOY'] = doy
        
    return data_dict


def fitLinearModel(data_dict):
    '''
    Fit a linear model to relate AGB to gamma0 backscatter.
    
    Args:
        data_dict: Dictionary output from extractGamma0() function.
    Returns:
        model slope, model intercept
    '''
    
    # Select only data that have values. NaN can arise in cases of masked areas in ALOS tiles
    sel = np.logical_and(np.isfinite(data_dict['gamma0_HV']), np.isfinite(data_dict['AGB']))
    
    assert sel.sum() > 0, "No usable data in data_dict."
    
    slope, intercept, r_value, p_value, stderr = scipy.stats.linregress(data_dict['gamma0_HV'][sel], data_dict['AGB'][sel])
    
    print "r-squared:", r_value**2
    print "p value:", p_value
    
    return slope, intercept