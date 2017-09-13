#!/usr/bin/env python

import argparse
import datetime
import matplotlib.pyplot as plt
import numpy as np
import os
from osgeo import gdal

"""
This is a script to calibrate ALOS mosaic tiles (and, eventually, cross-calibrate ALOS-1/ALOS-2)
"""


def generateFilenames(lat, lon, year, data_dir):
    """
    Determines ALOS mosaic filename for a given year, latitutude and longitude.
    """

    # Test that inputs are of reasonable lats/lons/years
    assert lat < 90. or lat > -90., "Latitude must be between -90 and 90 degrees."
    assert lon < 180. or lon > -180., "Longitude must be between -180 and 180 degrees."
    assert (year >= 2007 and year <= 2010) or (year >= 2015 and year <= datetime.datetime.now().year), "Years must be in the range 2007 - 2010 and 2015 - present. Your year was %s."%str(year)
    
    # Get hemisphere
    hem_NS = 'S' if lat < 0 else 'N'
    hem_EW = 'W' if lon < 0 else 'E'
    
    # Get lat/lon for directory (lat/lon of upper-left hand corner of 5x5 tiles)
    lat_dir = hem_NS + str(abs(lat + (5 - lat) % 5)).zfill(2)
    lon_dir = hem_EW +  str(abs(lon - 5 + (5 - lon) % 5)).zfill(3)
    
    # Get lat/lon for filename
    lat_file = hem_NS + str(abs(lat)).zfill(2)
    lon_file = hem_EW + str(abs(lon)).zfill(3)

    # Directories and files have standardised pattern
    name_pattern = '%s%s_%s_%s'
    
    # Directory name patterns name patterns are different for ALOS-1/ALOS-2
    if year >= 2015:
        name_pattern += '_F02DAR'

    # Generate directory name
    directory = data_dir + '/' + name_pattern%(lat_dir, lon_dir, str(year)[-2:], 'MOS') + '/'

    # Generate file names
    HV_file = directory + name_pattern%(lat_file, lon_file, str(year)[-2:], 'sl_HV')
    mask_file = directory + name_pattern%(lat_file, lon_file, str(year)[-2:], 'mask')
    
    return HV_file, mask_file


def openTile(HV_file, mask_file):
    """
    Opens ALOS mosaic tile, returns gdal dataset and array.
    """
    
    HV_ds = gdal.Open(HV_file, 0)
    HV_data = HV_ds.ReadAsArray()

    mask_ds = gdal.Open(mask_file, 0)
    mask_data = mask_ds.ReadAsArray()
    
    data_out = np.ma.array(HV_data, mask = np.logical_or(mask_data != 255, HV_file == 0))
    
    return HV_ds, data_out


def calibrateToGamma0(DN):
    """
    Calibrates ALOS-1/ALOS-2 data from digital numbers to gamma0 backscatter.
    Follows instructions in:
    Global 25m Resolution PALSAR-2/PALSAR Mosaic and Forest/Non-Forest Map (FNF): Dataset Description
    """
    
    gamma0 = 10 * np.ma.log10(DN.astype(np.float) ** 2) - 83. # units = decibels
    gamma0 = 10 ** (gamma0 / 10.) # Convert to natural units
    
    return gamma0


def calibrateToAGB(gamma0):
    """
    Placeholder equation to calibrate backscatter (gamma0) to AGB (tC/ha).
    """
    
    AGB = 715.667 * gamma0 - 5.967
    
    return AGB


def buildMap(fig, ax, data, lat, lon, title ='', cbartitle = '', vmin = 10., vmax = 60., cmap = 'YlGn'):
    """
    Builds a standardised map for overviewFigure().
    """
    
    im = ax.imshow(data, vmin = vmin, vmax = vmax, cmap = cmap, interpolation = 'nearest')
    
    ax.set_xticks(np.arange(0,4501,450))
    ax.set_yticks(np.arange(0,4501,450))
    ax.set_xticklabels(np.arange(lon,lon + 1.01, 0.1))
    ax.set_yticklabels(np.arange(lat, lat - 1.01, - 0.1))
    ax.tick_params(labelsize = 5)
    ax.set_xlabel('Longitude', fontsize = 5)
    ax.set_ylabel('Latitude', fontsize = 5)
    #ax.grid()
    ax.set_title(title, fontsize = 8)
    
    cbar = fig.colorbar(im, ax = ax, fraction = 0.046, pad = 0.04)
    cbar.ax.tick_params(labelsize = 6)
    cbar.set_label(cbartitle, fontsize = 7)
    

def overviewFigure(AGB_t1, AGB_t2, t1, t2, lat, lon, output_dir = os.getcwd()):
    """
    Create an overview image image showing biomass and proportional biomass change for the tile being processed.
    """
    
    # Update masks to exclude areas outisde forest definition. Good for visualisation
    AGB_t1 = np.ma.array(AGB_t1, mask = np.logical_or(AGB_t1.mask, AGB_t1 < 10.))
    AGB_t2 = np.ma.array(AGB_t2, mask = np.logical_or(AGB_t2.mask, AGB_t1 < 10.))
    
    
    AGB_change = (AGB_t2 - AGB_t1) / (t2 - t1) # tC/ha/yr

    AGB_pcChange = 100 * (AGB_change / AGB_t1) # %/yr
    
    fig = plt.figure(figsize = (7, 6))
    
    # Plot a map of AGB at t1
    ax1 = fig.add_subplot(2, 2, 1)
    buildMap(fig, ax1, AGB_t1, lat, lon, title = 'AGB %s'%str(t1), cbartitle = 'tC/ha')
    
    # Plot a map of AGB at t2
    ax2 = fig.add_subplot(2, 2, 2)
    buildMap(fig, ax2, AGB_t2, lat, lon, title = 'AGB %s'%str(t2), cbartitle = 'tC/ha')    
    
    # Plot a map of absolute AGB change   
    ax3 = fig.add_subplot(2, 2, 3)
    buildMap(fig, ax3, AGB_change, lat, lon, title = 'AGB change (%s-%s)'%(str(t1),str(t2)), cbartitle = 'tC/ha/yr',
             vmin = -10., vmax = 10., cmap = 'RdBu')    
    
    # Plot a map of % AGB change
    ax4 = fig.add_subplot(2, 2, 4)
    buildMap(fig, ax4, AGB_pcChange, lat, lon, title = 'AGB change (%s-%s)'%(str(t1),str(t2)), cbartitle = '%/yr',
             vmin = -50., vmax = 50., cmap = 'RdBu')    
    
    plt.tight_layout()
    
    hem_NS = 'S' if lat < 0 else 'N'
    hem_EW = 'W' if lon < 0 else 'E'
    
    filename = '/overview_%s%s.png'%(hem_NS + str(abs(lat)), hem_EW + str(abs(lon)))
    
    plt.savefig(output_dir + filename, dpi = 150)
    plt.close()


data_dir = '/home/sbowers3/DATA/ALOS_data/ALOS_mosaic/kilwa/'
output_dir = '/home/sbowers3/DATA/ALOS_data/ALOS_mosaic/kilwa/'

t1 = 2007
t2 = 2016

lat = -11
lon = 38

HVfile_t1, maskfile_t1 = generateFilenames(lat, lon, t1, data_dir)
HVfile_t2, maskfile_t2 = generateFilenames(lat, lon, t2, data_dir)

ds_t1, data_t1 = openTile(HVfile_t1, maskfile_t1)
ds_t2, data_t2 = openTile(HVfile_t2, maskfile_t2)

gamma0_t1 = calibrateToGamma0(data_t1)
gamma0_t2 = calibrateToGamma0(data_t2)

AGB_t1 = calibrateToAGB(gamma0_t1)
AGB_t2 = calibrateToAGB(gamma0_t2)

overviewFigure(AGB_t1, AGB_t2, t1, t2, lat, lon, output_dir = output_dir)








"""

def hist_match(source, template):
    '''
    Adjust the pixel values of a grayscale image such that its histogram
    matches that of a target image

    Arguments:
    -----------
        source: np.ndarray
            Image to transform; the histogram is computed over the flattened
            array
        template: np.ndarray
            Template image; can have different dimensions to source
    Returns:
    -----------
        matched: np.ndarray
            The transformed output image
    '''

    oldshape = source.shape
    source = source.ravel()
    template = template.ravel()

    # get the set of unique pixel values and their corresponding indices and
    # counts
    s_values, bin_idx, s_counts = np.unique(source, return_inverse=True,
                                            return_counts=True)
    t_values, t_counts = np.unique(template, return_counts=True)

    # take the cumsum of the counts and normalize by the number of pixels to
    # get the empirical cumulative distribution functions for the source and
    # template images (maps pixel value --> quantile)
    s_quantiles = np.cumsum(s_counts).astype(np.float64)
    s_quantiles /= s_quantiles[-1]
    t_quantiles = np.cumsum(t_counts).astype(np.float64)
    t_quantiles /= t_quantiles[-1]

    # interpolate linearly to find the pixel values in the template image
    # that correspond most closely to the quantiles in the source image
    interp_t_values = np.interp(s_quantiles, t_quantiles, t_values)

    return interp_t_values[bin_idx].reshape(oldshape)


def normalise(source, template):
    '''
    Matches maximum value in template image.
    '''
    
    return ((source.astype(np.float) / np.max(source)) * np.max(template)).astype(source.dtype)


def medianFit(source, template):
    
    out = np.zeros_like(source)
    
    for i in np.unique(source):
        out[source==i] = np.median(template[template==i])
    
    return out


    
def linearFit(source, template):
    '''
    Matches source to template with linear regression, based on log10 of source and template
    '''
    
    import scipy.stats
    
    # Create linear regression object    
    s = np.logical_and(np.isfinite(np.log10(source)), np.isfinite(np.log10(template)))
    
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(np.log10(source[s]), np.log10(template[s]))
    
    return 10 ** ((slope * np.log10(source)) + intercept)
    
""" 

