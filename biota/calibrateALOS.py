#!/usr/bin/env python

import argparse
from osgeo import gdal
import matplotlib.pyplot as plt
import numpy as np

"""
This is a script to calibrate ALOS mosaic tiles (and, eventually, cross-calibrate ALOS-1/ALOS-2)
"""



def calibrateToGamma0(DN):
    """
    Calibrates ALOS-1/ALOS-2 data from digital numbers to gamma0 backscatter.
    Follows instructions in:
    Global 25m Resolution PALSAR-2/PALSAR Mosaic and Forest/Non-Forest Map (FNF): Dataset Description
    """
    
    gamma0 = 10 * np.log10(DN.astype(np.float) ** 2) - 83. # units = decibels
    gamma0 = 10 ** (gamma0 / 10.) # Natural units
    
    return gamma0


file_2007 = '/home/sbowers3/DATA/ALOS_data/ALOS_mosaic/kilwa/S10E035_07_MOS/S11E038_07_sl_HV'
file_2016 = '/home/sbowers3/DATA/ALOS_data/ALOS_mosaic/kilwa/S10E035_16_MOS_F02DAR/S11E038_16_sl_HV_F02DAR'


ds_2007 = gdal.Open(file_2007, 0)
ds_2016 = gdal.Open(file_2016, 0)

data_2007 = ds_2007.ReadAsArray()
data_2016 = ds_2016.ReadAsArray()



def hist_match(source, template):
    """
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
    """

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
    """
    Matches maximum value in template image.
    """
    
    return ((source.astype(np.float) / np.max(source)) * np.max(template)).astype(source.dtype)


def medianFit(source, template):
    
    out = np.zeros_like(source)
    
    for i in np.unique(source):
        out[source==i] = np.median(template[template==i])
    
    return out


    
def linearFit(source, template):
    """
    Matches source to template with linear regression, based on log10 of source and template
    """
    
    import scipy.stats
    
    # Create linear regression object    
    s = np.logical_and(np.isfinite(np.log10(source)), np.isfinite(np.log10(template)))
    
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(np.log10(source[s]), np.log10(template[s]))
    
    return 10 ** ((slope * np.log10(source)) + intercept)
    
    

#data_2016 = hist_match(data_2016, data_2007)
#data_2016 = normalise(data_2016, data_2007)
#data_2016 = medianFit(data_2016, data_2007)


data_2007 = calibrateToGamma0(data_2007)
data_2016 = calibrateToGamma0(data_2016)

#data_2016 = linearFit(data_2016, data_2007)



plt.subplot(121)
plt.imshow(data_2007, vmin = 0, vmax = 0.06)
plt.colorbar()
plt.subplot(122)
plt.imshow(data_2016, vmin = 0, vmax = 0.06)
plt.colorbar()
plt.show()


#plt.scatter(np.log10(data_2007[::10,::10]), np.log10(data_2016[::10,::10]))
#plt.show()