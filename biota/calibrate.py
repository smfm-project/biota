#!/usr/bin/env python


# This is a set of scripts for calibration of a biomass-backscatter curve

import argparse
import csv
import numpy as np
import scipy.stats

import biota
import biota.mask

import pdb



def extractGamma0(dataloc, year, shp, plot_field, agb_field, buffer_size = 0, verbose = False, units = 'natural'):
    '''
    Extract gamma0 from ALOS tiles.

    Args:
        shp: A shapefile containing plot data
        year: Year to load
        plot_field: The shapefile field containing plot names
        agb_field: The shapefile field containing AGB estimates
    Returns:
        A dictionary containing plot names, AGB and gamma0 values
    '''

    # Use example data to get processing steps
    #downsample_factor = tile_example.downsample_factor
    #lee_filter = tile_example.lee_filter
    #year = tile_example.year
    #dataloc = tile_example.dataloc

    assert (year >= 2007 and year <= 2010) or year >= 2015, "Invalid year (%s) input"%str(year)

    downsample_factor = 1
    lee_filter = True
    window_size = 3

    # Extract relevant info from shapefile
    plot_names = biota.mask.getField(shp, plot_field)

    if agb_field is not None:
        agb = biota.mask.getField(shp, agb_field).astype(np.float)

    # Identify tiles that contain gamma0 data for the shapefile
    tiles = biota.mask.getTilesInShapefile(shp)

    # Check whether all the necessary tiles are all present. Decide what to do if they aren't.
    # TODO

    # Generate output array for backscatter
    gamma0_mean_hv = np.empty_like(plot_names, dtype = np.float32)
    gamma0_mean_hv[:] = np.nan

    gamma0_mean_hh = gamma0_mean_hv.copy()
    gamma0_std_hv = gamma0_mean_hv.copy()
    gamma0_std_hh = gamma0_mean_hv.copy()
    doy = gamma0_mean_hv.copy()

    # For each tile covered by shapefile
    for lat, lon in tiles:

        if verbose: print('Doing lat: %s, lon: %s'%(str(lat), str(lon)))

        # Load tile
        try:
            tile = biota.LoadTile(dataloc, lat, lon, year, downsample_factor = downsample_factor, lee_filter = lee_filter, window_size = window_size)
        except:
            continue

        # Get backscatter (both polarisations) and DOY
        data_gamma0_hv = tile.getGamma0(polarisation = 'HV', units = units)
        data_gamma0_hh = tile.getGamma0(polarisation = 'HH', units = units)
        data_doy = tile.getDOY()

        # Mask out each plot
        plot_mask = biota.mask.maskShapefile(tile, shp, location_id = True, buffer_size = buffer_size)

        # Extract values for each plot
        for n in np.unique(plot_mask[plot_mask != 0]):

            # Get mask for plot and tile
            this_mask = np.logical_and(plot_mask == n, tile.mask == False)

            if this_mask.sum() == 0: continue

            # Add metrics to output array
            gamma0_mean_hv[n-1] = np.nanmean(data_gamma0_hv[this_mask])
            gamma0_mean_hh[n-1] = np.nanmean(data_gamma0_hh[this_mask])
            gamma0_std_hv[n-1] = np.nanstd(data_gamma0_hv[this_mask])
            gamma0_std_hh[n-1] = np.nanstd(data_gamma0_hh[this_mask])
            doy[n-1] = np.nanmedian(data_doy[this_mask])

            if np.isnan(np.nanmean(data_gamma0_hv[this_mask])): pdb.set_trace()

    # Return data as a dictionary
    data_dict = {}
    data_dict['plot_name'] = plot_names.tolist()
    data_dict['gamma0_mean_HH'] = gamma0_mean_hh.tolist()
    data_dict['gamma0_mean_HV'] = gamma0_mean_hv.tolist()
    data_dict['gamma0_std_HH'] = gamma0_std_hh.tolist()
    data_dict['gamma0_std_HV'] = gamma0_std_hv.tolist()
    data_dict['DOY'] = doy.tolist()

    if agb_field is not None:
        data_dict['plot_AGB'] = agb.tolist()

    with open('gamma0_%s_by_plot.csv'%str(year), 'wb') as f:  # Just use 'w' mode in 3.x
        writer = csv.writer(f, delimiter = ',')
        writer.writerow(list(data_dict.keys()))
        for row in range(len(plot_names)):
            writer.writerow([data_dict[k][row] for k in list(data_dict.keys())])

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
    sel = np.logical_and(np.isfinite(data_dict['gamma0_HV']), np.isfinite(data_dict['plot_AGB']))

    assert sel.sum() > 0, "No usable data in data_dict."

    slope, intercept, r_value, p_value, stderr = scipy.stats.linregress(data_dict['gamma0_HV'][sel], data_dict['AGB'][sel])

    print("r-squared:", r_value**2)
    print("p value:", p_value)

    return slope, intercept


