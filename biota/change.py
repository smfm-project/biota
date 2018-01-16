# Toolset to calibrate change attribution algorithm

import math
import numpy as np
import biota.calibration as cal
import scipy.stats as stats
import scipy.ndimage
from scipy.ndimage.measurements import label

import pdb




def lee_filter(img, size):
    """
    # From https://stackoverflow.com/questions/39785970/speckle-lee-filter-in-python
    """
    from scipy.ndimage.filters import uniform_filter
    from scipy.ndimage.measurements import variance
    
    img_mean = uniform_filter(img, (size, size))
    img_sqr_mean = uniform_filter(img**2, (size, size))
    img_variance = img_sqr_mean - img_mean**2

    overall_variance = variance(img)

    img_weights = img_variance**2 / (img_variance**2 + overall_variance**2)
    img_output = img_mean + img_weights * (img - img_mean)
    
    return img_output



def getChange(AGB_t1, AGB_t2, F_threshold = 10., C_threshold = 0.2):
    """
    Returns pixels that meet change detection thresholds for country.
    
    Args:
        F_threshold = threshold above which a pixel is forest
        C_threshold = threshold of proportional change with which to begin change detection
    """

    # First 2007 to 2010 
    
    # Pixels that move from forest to notforest or vice versa
    F_NF = np.logical_and(AGB_t1 >= F_threshold, AGB_t2 < F_threshold)
    NF_F = np.logical_and(AGB_t1 < F_threshold, AGB_t2 >= F_threshold)
    F_F = np.logical_and(AGB_t1 >= F_threshold, AGB_t2 >= F_threshold)
    NF_NF = np.logical_and(AGB_t1 < F_threshold, AGB_t2 < F_threshold)
    
    # Change pixels
    CHANGE = np.logical_or(((AGB_t1 - AGB_t2) / AGB_t1) >= C_threshold, ((AGB_t1 - AGB_t2) / AGB_t1) < (- C_threshold))
    NOCHANGE = CHANGE==False

    # Trajectory
    DECREASE = AGB_t2 < AGB_t1
    INCREASE = AGB_t2 >= AGB_t1

    metrics = {}
    
    # These are all the possible definitions. Note: they *must* add up to one.
    metrics['deforestation'] = (F_NF & CHANGE) * 1
    metrics['degradation'] = (F_F & CHANGE & DECREASE) * 1
    metrics['minorloss'] = ((F_F | F_NF) & NOCHANGE & DECREASE) * 1

    metrics['aforestation'] = (NF_F & CHANGE) * 1
    metrics['growth'] = (F_F & CHANGE & INCREASE) * 1
    metrics['minorgain'] = ((F_F | NF_F) & NOCHANGE & INCREASE) * 1
            
    metrics['notforest'] = (NF_NF * 1)

    return metrics


def getContiguousAreas(data, value, min_pixels):
    '''
    Get pixels that come from the same contigous area.
    
    Args:
        data: A numpy array
        min_val: The minimum value to contain within the cotiguous area (inclusive)
        max_val: The maxiumum value to contain within the contigous area (inclusive)
        min_area: What minimum area should be included (number of pixels, queens move)
    
    Returns:
        A binary array of pixels that meet the conditiuons
    '''
    
    # Extract area that meets condition
    binary_array = data == value * 1
    
    # Label contigous areas with a number
    s = scipy.ndimage.generate_binary_structure(2,2) # This connects diagonal elements
    location_id, n_areas = label(binary_array, structure = s)
    
    # Get count of each value in array
    label_area = np.bincount(location_id.flatten())[1:]
    
    # Find those IDs that meet minimum area requirements
    include_id = np.arange(1, n_areas + 1)[label_area >= min_pixels]
    
    # Get a binary array of location_id pixels that meet the minimum area requirement
    contiguous_area = np.in1d(location_id, include_id).reshape(data.shape) * 1
    
    # Return an array giving values to each area
    #location_id[contiguous_area == False] = 0
    
    return contiguous_area#, location_id


def _window_mean(img, window_size = 3):
    '''
    Based on https://stackoverflow.com/questions/18419871/improving-code-efficiency-standard-deviation-on-sliding-windows
    '''
    
    from scipy import signal
    
    c1 = signal.convolve2d(img, np.ones((window_size, window_size)) / (window_size ** 2), boundary = 'symm')
    
    border = window_size/2
    
    return c1[border:-border, border:-border]
    

def _window_stdev(img, window_size = 3):
    '''
    Based on https://stackoverflow.com/questions/18419871/improving-code-efficiency-standard-deviation-on-sliding-windows
    and http://nickc1.github.io/python,/matlab/2016/05/17/Standard-Deviation-(Filters)-in-Matlab-and-Python.html
    '''
    
    from scipy import signal
    
    c1 = signal.convolve2d(img, np.ones((window_size, window_size)) / (window_size ** 2), boundary = 'symm')
    c2 = signal.convolve2d(img*img, np.ones((window_size, window_size)) / (window_size ** 2), boundary = 'symm')
    
    border = window_size / 2
    
    return np.sqrt(c2 - c1 * c1)[border:-border, border:-border]




def enhanced_lee_filter(img, window_size = 3, n_looks = 16):
    '''
    Filters a masked array with the enhanced lee filter.
    
    Args:
        img: A masked array
    Returns:
        A masked array with a filtered verison of img
    '''
    
    assert type(window_size == int), "Window size must be an integer. You input the value: %s"%str(window_size)
    assert (window_size % 2) == 1, "Window size must be an odd number. You input the value: %s"%str(window_size)
    
    k = 1. #Adequate for most SAR images
    
    # Set parameters to default
    #cu = 0.523
    #cmax = 1.73
    
    cu = (1./n_looks) ** 0.5
    cmax =  (1 + (2./n_looks)) ** 0.5

    img_mask = img.data
    img_mask[img.mask == True] = np.nan
    
    img_mean = _window_mean(img_mask, window_size = window_size)
    img_std = _window_stdev(img_mask, window_size = window_size)

    ci = img_std / img_mean
    ci[np.isfinite(ci) == False] = 0.
    
    w_t = np.zeros_like(ci)
    
    # There are three conditions in the enhanced lee filter
    w_t[ci <= cu] = 1.
    w_t[ci >= cmax] = 0.
    
    s = np.logical_and(ci > cu, ci < cmax)
    w_t[s] = np.exp((-k * (ci[s] - cu)) / (cmax - ci[s]))
        
    img_filtered = (img_mean * w_t) + (img_mask * (1. - w_t))
    
    img_filtered = np.ma.array(img_filtered, mask = np.isnan(img_filtered))
    
    img_filtered.data[img_filtered.mask] = 0.
    
    return img_filtered
    
def rebin(a, shape = (450, 450)):
    '''
    Reduce array size to 0.5 ha resolution by taking mean, to increase ENL.
    # From https://stackoverflow.com/questions/8090229/resize-with-averaging-or-rebin-a-numpy-2d-array
    '''
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)



data_dir = '/home/sbowers3/DATA/ALOS_data/ALOS_mosaic/gorongosa/'
output_dir = '/home/sbowers3/DATA/ALOS_data/ALOS_mosaic/gorongosa/'

#Nhambita = lat -18 lon 33

lat = -19
lon = 33#34#39
    
data_2007 = cal.ALOS(lat, lon, 2007, data_dir)
#data_2008 = cal.ALOS(lat, lon, 2008, data_dir)
#data_2009 = cal.ALOS(lat, lon, 2009, data_dir)
data_2010 = cal.ALOS(lat, lon, 2010, data_dir)

gamma0_2007 = data_2007.getGamma0(units='decibels')
#gamma0_2008 = data_2008.getGamma0()
#gamma0_2009 = data_2009.getGamma0()
gamma0_2010 = data_2010.getGamma0(units='decibels')


gamma0_2007 = enhanced_lee_filter(gamma0_2007, window_size = 3, n_looks = 16)
gamma0_2010 = enhanced_lee_filter(gamma0_2010, window_size = 3, n_looks = 16)
#gamma0_2007 = enhanced_lee_filter(rebin(gamma0_2007), window_size = 5, n_looks = 16*4)
#gamma0_2010 = enhanced_lee_filter(rebin(gamma0_2010), window_size = 5, n_looks = 16*4)


AGB_2007 = 715.667 * (10 ** (gamma0_2007 / 10.)) - 5.967
AGB_2010 = 715.667 * (10 ** (gamma0_2010 / 10.)) - 5.967

#AGB_2007 = data_2007.getAGB()
#AGB_2010 = data_2010.getAGB()

#AGB_2008 = lee_filter(data_2008.getAGB(), 10)
#AGB_2009 = lee_filter(data_2009.getAGB(), 10)
#AGB_2010 = lee_enhanced_filter(data_2010.getAGB())

#change_07_08 = getChange(AGB_2007, AGB_2008)
#change_08_09 = getChange(AGB_2008, AGB_2009)
#change_09_10 = getChange(AGB_2009, AGB_2010)
change_07_10 = getChange(AGB_2007, AGB_2010, F_threshold = 10., C_threshold = 0.2)

contiguous_change_area = getContiguousAreas(change_07_10['degradation'] + change_07_10['deforestation'], 1, 8)

deforestation = change_07_10['deforestation'] * contiguous_change_area
degradation = change_07_10['degradation'] * contiguous_change_area



#deforestation_07_08 = getContiguousAreas(change_07_08['deforestation'], 1, 8)
#deforestation_08_09 = getContiguousAreas(change_08_09['deforestation'], 1, 8)
#deforestation_09_10 = getContiguousAreas(change_09_10['deforestation'], 1, 8)

#deforestation = deforestation_07_10.sum() / (4500 * 4500.)
#total_deforestation = (deforestation_07_08.sum() + deforestation_08_09.sum() + deforestation_09_10.sum()) / (4500. * 4500)


# This result looks nice. TODO: calculate % change of entire contiguous change region

from osgeo import osr, gdal
srs = osr.SpatialReference()
srs.ImportFromEPSG(4326)

driver = gdal.GetDriverByName('GTiff')
ds = driver.Create(output_dir + 'forest_loss.tif', 4500, 4500, 1, gdal.GDT_UInt16)
ds.SetGeoTransform(data_2007.geo_t)
ds.SetProjection(srs.ExportToWkt())
ds.GetRasterBand(1).WriteArray(deforestation*2 + degradation)
ds = None


"""
geo_t = (33.0, 0.00022222222222222223*10, 0.0, -18.0, 0.0, -0.00022222222222222223*10)
driver = gdal.GetDriverByName('GTiff')
ds = driver.Create(output_dir + 'forest_loss.tif', 450, 450, 1, gdal.GDT_UInt16)
ds.SetGeoTransform(geo_t)
ds.SetProjection(srs.ExportToWkt())
ds.GetRasterBand(1).WriteArray(deforestation*2 + degradation)
ds = None
"""
"""

contiguous_change_area_gain = getContiguousAreas(change_07_10['aforestation'] + change_07_10['growth'], 1, 8)

aforestation = change_07_10['aforestation'] * contiguous_change_area_gain
growth = change_07_10['growth'] * contiguous_change_area_gain

# A deforestation rate from differencing?

out = np.zeros((450,450))
for x in range(450):
    for y in range(450):
        out[y,x] = aforestation[y*10:(y*10)+10,x*10:(x*10)+10].data.sum() - deforestation[y*10:(y*10)+10,x*10:(x*10)+10].data.sum()
        

print "Deforestation rate: %s"%str(np.sum(out[out<0])/(4500*4500.))
print "Aforestation rate: %s"%str(np.sum(out[out>0])/(4500*4500.))

"""
