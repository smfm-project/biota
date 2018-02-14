#!/usr/bin/env python

import numpy as np
import signal
import scipy.ndimage as ndimage

import pdb

"""
This file contains scripts to filter ALOS data.
"""

def enhanced_lee_filter(img, window_size = 7, n_looks = 16):
    '''
    Filters a masked array with the enhanced lee filter.
    
    Based on formulation here: http://www.pcigeomatics.com/geomatica-help/concepts/orthoengine_c/Chapter_825.html
    
    Args:
        img: A masked array
    Returns:
        A masked array with a filtered verison of img
    '''
    
    assert type(window_size == int), "Window size must be an integer. You input the value: %s"%str(window_size)
    assert (window_size % 2) == 1, "Window size must be an odd number. You input the value: %s"%str(window_size)
    assert window_size > 3, "Window size must be at least 3. You input the value: %s"%str(window_size)

    # Inner function to calculate mean with a moving window and a masked array
    def _window_mean(img, window_size = 3):
        '''
        Based on https://stackoverflow.com/questions/18419871/improving-code-efficiency-standard-deviation-on-sliding-windows
        '''
        
        from scipy import signal
        
        c1 = signal.convolve2d(img, np.ones((window_size, window_size)) / (window_size ** 2), boundary = 'symm')
        
        border = window_size/2
        
        return c1[border:-border, border:-border]
        
    # Inner function to calculate standard deviation with a moving window and a masked array
    def _window_stdev(img, window_size = 3):
        '''
        Based on https://stackoverflow.com/questions/18419871/improving-code-efficiency-standard-deviation-on-sliding-windows
        and http://nickc1.github.io/python,/matlab/2016/05/17/Standard-Deviation-(Filters)-in-Matlab-and-Python.html
        '''
        
        from scipy import signal
        
        c1 = signal.convolve2d(img, np.ones((window_size, window_size)) / (window_size ** 2), boundary = 'symm')
        c2 = signal.convolve2d(img*img, np.ones((window_size, window_size)) / (window_size ** 2), boundary = 'symm')
        
        border = window_size / 2
        
        variance = c2 - c1 * c1
        variance[variance < 0] += 0.01 # Prevents divide by zero errors.
        
        return np.sqrt(variance)[border:-border, border:-border]
    
    # Damping factor, set to 1 which is adequate for most SAR images
    k = 1
        
    cu = (1./n_looks) ** 0.5
    cmax =  (1 + (2./n_looks)) ** 0.5
    
    
    # Interpolate across nodata areas. No standard Python filters understand nodata values; this is a simplification
    
    indices = ndimage.distance_transform_edt(img.mask, return_distances = False, return_indices = True)
    data = img.data[tuple(indices)]
    
    img_mean = _window_mean(data, window_size = window_size)
    img_std = _window_stdev(data, window_size = window_size)

    ci = img_std / img_mean
    ci[np.isfinite(ci) == False] = 0.
    
    W = np.zeros_like(ci)
    
    # There are three conditions in the enhanced lee filter
    W[ci <= cu] = 1.
    W[ci >= cmax] = 0.
    
    s = np.logical_and(ci > cu, ci < cmax)
    W[s] = np.exp((-k * (ci[s] - cu)) / (cmax - ci[s]))
        
    img_filtered = (img_mean * W) + (data * (1. - W))
    
    img_filtered = np.ma.array(img_filtered, mask = np.logical_or(np.isnan(img_filtered), img.mask))
    
    img_filtered.data[img_filtered.mask] = 0.
    
    return img_filtered
