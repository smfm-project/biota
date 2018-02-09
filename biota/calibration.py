#!/usr/bin/env python

import argparse
import datetime as dt
import itertools
import math
import numpy as np
import os
from PIL import Image, ImageDraw
from scipy import ndimage
import scipy.stats as stats

import matplotlib.pyplot as plt
import pdb



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
        
        return np.sqrt(c2 - c1 * c1)[border:-border, border:-border]
    
    k = 1. #Adequate for most SAR images
        
    cu = (1./n_looks) ** 0.5
    cmax =  (1 + (2./n_looks)) ** 0.5
    
    # Or set parameters to default?
    #cu = 0.523
    #cmax = 1.73
    
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


def outputGeoTiff(data, filename, geo_t, proj, output_dir = os.getcwd(), dtype = 6, nodata = None):
    """
    Writes a GeoTiff file to disk.
    
    Args:
        data: A numpy array.
        geo_t: A GDAL geoMatrix (ds.GetGeoTransform()).
        proj: A GDAL projection (ds.GetProjection()).
        filename: Specify an output file name.
        output_dir: Optioanlly specify an output directory. Defaults to working directory.
        dtype: gdal data type (gdal.GDT_*). Defaults to gdal.GDT_Float32.
        nodata: The nodata value for the array
    """
    
    from osgeo import osr, gdal
    
    # Get full output path
    output_path = '%s/%s.tif'%(os.path.abspath(output_dir), filename.rstrip('.tif'))
    
    # Save image with georeference info
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(output_path, data.shape[0], data.shape[1], 1, dtype, options = ['COMPRESS=LZW'])
    ds.SetGeoTransform(geo_t)
    ds.SetProjection(proj)
    
    # Set nodata
    if nodata != None:
        ds.GetRasterBand(1).SetNoDataValue(nodata)
    
    # Write data for masked and unmasked arrays
    if np.ma.isMaskedArray(data):
        ds.GetRasterBand(1).WriteArray(data.filled(nodata))
    else:
        ds.GetRasterBand(1).WriteArray(data)
    ds = None


def getContiguousAreas(data, value, min_pixels = 1):
    '''
    Get pixels that come from the same contigous area.
    
    Args:
        data: A numpy array
        value: Pixel value to include in contiguous_area (e.g. True for forest)
        min_area: What minimum area should be included (number of pixels, queens move)
    
    Returns:
        A binary array of pixels that meet the conditiuons
    '''
    
    from scipy.ndimage.measurements import label
    
    # If any pixels are masked, we give them the value of the nearest valid pixel.
    if np.ma.isMaskedArray(data):
        if data.mask.sum() > 0:
            ind = ndimage.distance_transform_edt(data.mask, return_distances = False, return_indices = True)
            data = np.ma.array(data.data[tuple(ind)], mask = data.mask)
    
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
    if np.ma.isMaskedArray(data):
        contiguous_area = np.ma.array(contiguous_area, mask = data.mask)
        location_id = np.ma.array(location_id, mask = data.mask)

    return contiguous_area, location_id


def calculateLDI(unique_ids):
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


class ALOS(object):
    """
    An ALOS mosaic tile.
    Mosaic tiles have the following properties:

    Attributes:
        lat:
        lon:
        DN: An array of uncalibrated digital numbers from the ALOS tile.
        mask:
    """
        
    def __init__(self, lat, lon, year, dataloc):
        """
        Loads data and metadata for an ALOS mosaic tile.
        """
               
        # Test that inputs are of reasonable lats/lons/years
        assert type(lat) == int, "Latitude must be an integer."
        assert lat < 90. or lat > -90., "Latitude must be between -90 and 90 degrees."
        assert type(lon) == int, "Longitude must be an integer."
        assert lon < 180. or lon > -180., "Longitude must be between -180 and 180 degrees."
        assert type(year) == int, "Year must be an integer."
        assert (year >= 2007 and year <= 2010) or (year >= 2015 and year <= dt.datetime.now().year), "Years must be in the range 2007 - 2010 and 2015 - present. Your input year was %s."%str(year)
        
        self.lat = lat
        self.lon = lon
        self.year = year
        
        # Deterine hemispheres
        self.hem_NS = 'S' if lat < 0 else 'N'
        self.hem_EW = 'W' if lon < 0 else 'E'
        
        # Determine whether ALOS-1 or ALOS-2
        self.satellite = self.__getSatellite()
        
        # Determine filenames
        self.dataloc = dataloc.rstrip('/')
        self.directory = self.__getDirectory()
        self.HH_path = self.__getHHPath()
        self.HV_path = self.__getHVPath()
        self.mask_path = self.__getMaskPath()
        self.date_path = self.__getDatePath()
        
        # Stop of the ALOS tile doesn't exist
        if not os.path.isfile(self.HV_path): 
            raise IOError('ALOS tile does not exist in the file system.')
        
        # Get GDAL geotransform and projection
        self.geo_t = self.__getGeoT()
        self.proj = self.__getProj()
        
        # Get Raster size
        self.xSize, self.ySize = self.__getSize(self.HH_path)
        
        # Determine x and y resolution in meters
        self.xRes = self.__getXRes()
        self.yRes = self.__getYRes()
        
        # Record the nodata value
        self.nodata = self.__getNodata()
        
        # Load DN, mask, and day of year
        self.mask = self.getMask()
        #self.DN = self.getDN()
        #self.DOY = self.getDOY()

    def __getSatellite(self):
        """
        Return the sensor for the ALOS mosaic tile
        """
        
        if self.year >= 2015:
            satellite = 'ALOS-2'
        else:
            satellite = 'ALOS-1'
        
        return satellite
    
    def __getDirectory(self):
        """
        Return the directory containing ALOS data for a given lat/lon. Assumes its distributed as a 5x5 tile at present.
        """
        
        # Calculate the hemisphere and lat/lon of the 5x5 tile
        lat_dir = self.lat + (5 - self.lat) % 5
        lon_dir = self.lon - (self.lon % 5)
        hem_NS_dir = 'S' if lat_dir < 0 else 'N'
        hem_EW_dir = 'W' if lon_dir < 0 else 'E'

        # Get lat/lon for directory (lat/lon of upper-left hand corner of 5x5 tiles)
        lat_dir = hem_NS_dir + str(abs(lat_dir)).zfill(2)
        lon_dir = hem_EW_dir + str(abs(lon_dir)).zfill(3)
                
        # Directories and files have standardised pattern
        name_pattern = '%s%s_%s_%s'
        
        # Directory name patterns name patterns are different for ALOS-1/ALOS-2
        if self.satellite == 'ALOS-2':
            name_pattern += '_F02DAR'
        
        # Generate directory name
        directory = self.dataloc + '/' + name_pattern%(lat_dir, lon_dir, str(self.year)[-2:], 'MOS') + '/'
        
        return directory
    
    def __getFilename(self, append_pattern):
        """
        Return the filename for ALOS data for a given lat/lon.
        """
      
        # Get lat/lon for filename
        lat_file = self.hem_NS + str(abs(self.lat)).zfill(2)
        lon_file = self.hem_EW + str(abs(self.lon)).zfill(3)
        
        name_pattern = '%s%s_%s_%s'
        
        # Directory name patterns name patterns are different for ALOS-1/ALOS-2
        if self.satellite == 'ALOS-2':
            name_pattern += '_F02DAR'

        # Generate file name
        return name_pattern%(lat_file, lon_file, str(self.year)[-2:], append_pattern)
    
    def __getHHPath(self):
        """
        Determines the filepath to HV data.
        """
        
        return self.__getDirectory() + self.__getFilename('sl_HH')
    
    def __getHVPath(self):
        """
        Determines the filepath to HV data.
        """
        
        return self.__getDirectory() + self.__getFilename('sl_HV')
    
    def __getMaskPath(self):
        """
        Determines the filepath to mask data.
        """
        
        return self.__getDirectory() + self.__getFilename('mask')
    
    def __getDatePath(self):
        """
        Determines the filepath to DOY data.
        """
        
        return self.__getDirectory() + self.__getFilename('date')
    
    def __getXRes(self):
        """
        Approximates the resolution of a pixel in meters.
        From: https://stackoverflow.com/questions/639695/how-to-convert-latitude-or-longitude-to-meters
        """
        
        m_per_degree = 111132.954 * math.cos(self.lat - 0.5)
        
        return m_per_degree / self.xSize
        
    def __getYRes(self):
        """
        Approximates the resolution of a pixel in meters.
        From: https://stackoverflow.com/questions/639695/how-to-convert-latitude-or-longitude-to-meters
        """
        
        m_per_degree = 111132.92 - 559.82 * math.cos(2 * (self.lat - 0.5)) + 1.175 * math.cos(4 * self.lat - 0.5)
        
        return m_per_degree / self.ySize
    
    def __getNodata(self):
        """
        Return the nodata value for ALOS array. Should always be 0, so this function is a placeholder should the format change.
        """
        return 0
    
    def __getGeoT(self):
        """
        Fetches a GDAL GeoTransform for tile.
        """
        
        from osgeo import gdal
         
        ds = gdal.Open(self.__getHHPath(), 0)
        geo_t = ds.GetGeoTransform()
               
        return geo_t
    
    def __getProj(self):
        """
        Fetches projection info for tile.
        """
        
        from osgeo import gdal
                
        ds = gdal.Open(self.__getHHPath(), 0)
        proj = ds.GetProjection()
        
        return proj
    
    def __getSize(self, filepath):
        """
        Determines the size of a tile.
        """
        
        from osgeo import gdal
   
        ds = gdal.Open(filepath, 0)
        x_size = ds.RasterXSize
        y_size = ds.RasterYSize
        
        return x_size, y_size
    
    def getMask(self):
        """
        Loads the mask into a numpy array.
        """
        
        from osgeo import gdal
        
        mask_ds = gdal.Open(self.mask_path, 0)
        mask = mask_ds.ReadAsArray()
        
        return mask != 255
    
    def getDOY(self):
        """
        Loads date values into a numpy array.
        """
        
        from osgeo import gdal
        
        date_ds = gdal.Open(self.date_path, 0)
        day_after_launch = date_ds.ReadAsArray()
        
        # Get list of unique dates in ALOS tile
        unique_days = np.unique(day_after_launch)
        
        if self.satellite == 'ALOS-2':
            launch_date = dt.datetime(2014,5,24)
        else:
            launch_date =  dt.datetime(2006,1,24)
        
        # Determine the Day of Year associated with each
        unique_dates = [launch_date  + dt.timedelta(days=int(d)) for d in unique_days]
        unique_doys = [(d - dt.datetime(d.year,1,1,0,0)).days + 1 for d in unique_dates]
        
        # Build a Day Of Year array
        DOY = np.zeros_like(day_after_launch)
        for day, doy in zip(unique_days, unique_doys):
            DOY[day_after_launch == day] = doy
                
        return DOY
    
    def getDN(self, polarisation = 'HV'):
        """
        Loads DN (raw) values into a numpy array.
        """
        
        from osgeo import gdal

        assert polarisation == 'HH' or polarisation == 'HV', "polarisation must be either 'HH' or 'HV'."
               
        if polarisation == 'HV':
            DN_ds = gdal.Open(self.HV_path, 0)
        else:
            DN_ds = gdal.Open(self.HH_path, 0)

        DN = DN_ds.ReadAsArray()
               
        return np.ma.array(DN, mask = self.mask)
    
    def getGamma0(self, polarisation = 'HV', units = 'natural', lee_filter = True, output = False):
        """
        Calibrates data to gamma0 (baskscatter) in decibels or natural units.
        """
        
        assert units == 'natural' or units == 'decibels', "Units must be 'natural' or 'decibels'. You input %s."%units
        
        # Calibrate DN to units of dB
        gamma0 = 10 * np.ma.log10(self.getDN(polarisation = polarisation).astype(np.float) ** 2) - 83. # units = decibels
        
        # Apply filter based on dB values
        if lee_filter:
            gamma0 = enhanced_lee_filter(gamma0)
        
        # Convert to natural units where specified
        if units == 'natural': gamma0 = 10 ** (gamma0 / 10.)
        
        # Keep masked values tidy
        gamma0.data[self.mask] = 0
        
        if output: self.__outputGeoTiff(gamma0, 'Gamma0')
        
        return gamma0
    
    def getAGB(self, lee_filter = True, output = False):
        """
        Calibrates data to aboveground biomass (AGB).
        Placeholder equation to calibrate backscatter (gamma0) to AGB (tC/ha).
        """
        
        # ALOS-1
        if self.satellite == 'ALOS-1':
            AGB = 715.667 * self.getGamma0(units = 'natural', polarisation = 'HV', lee_filter = lee_filter) - 5.967
            
        elif self.satellite == 'ALOS-2':
            AGB = 715.667 * self.getGamma0(units = 'natural', polarisation = 'HV', lee_filter = lee_filter) - 5.967
            
        else:       
            raise ValueError("Unknown satellite named '%s'. self.satellite must be 'ALOS-1' or 'ALOS-2'."%self.satellite)
        
        # Keep masked values tidy
        AGB.data[self.mask] = 0.
        
        if output: self.__outputGeoTiff(AGB, 'AGB')
        
        return AGB
    
    def getWoodyCover(self, threshold, min_area = 0., lee_filter = True, output = False):
        """
        Get woody cover, based on a threshold of AGB.
        min_area in ha
        """
        
        from osgeo import gdal
        
        woody_cover = self.getAGB(lee_filter = lee_filter) >= float(threshold)
               
        if min_area > 0:
            
            # Calculate number of pixels in min_area
            min_pixels = int(round(min_area / (self.yRes * self.xRes * 0.0001)))
            
            # Remove pixels that aren't part of a forest block of size at least min_pixels
            contiguous_area, location_id = getContiguousAreas(woody_cover, True, min_pixels = min_pixels)
            woody_cover[contiguous_area == False] = False
        
        if output: self.__outputGeoTiff(woody_cover * 1, 'WoodyCover', dtype = gdal.GDT_Byte)
        
        return woody_cover
    
    def getForestPatches(self, threshold, min_area = 0, lee_filter = True, output = False):
        """
        Get numbered forest patches, based on woody cover threshold.
        """
        
        from osgeo import gdal
        
        woody_cover = self.getWoodyCover(threshold, lee_filter = lee_filter)
        
        # Calculate number of pixels in min_area
        min_pixels = int(round(min_area / (self.yRes * self.xRes * 0.0001)))
                
        # Get areas that meet that threshold
        contiguous_area, location_id = getContiguousAreas(woody_cover, True, min_pixels = min_pixels)
                
        if output: self.__outputGeoTiff(woody_cover * 1, 'ForestPatches', dtype = gdal.GDT_Int32)
        
        return location_id
    
    def __shrinkGeoT(self, factor):
        """
        Function to modify gdal geo_transform to output at correct resolution for downsampled products.
        """
        
        geo_t = list(self.geo_t)
        
        geo_t[1] = geo_t[1] * (self.xSize / factor)
        geo_t[-1] = geo_t[-1] * (self.ySize / factor)
        
        return tuple(geo_t)
    
    def calculateLDI(self, threshold, lee_filter = True, output = False):
        """
        Landcover Division Index (LDI) is an indicator of the degree of habitat coherence, ranging from 0 (fully contiguous) to 1 (very fragmented).
        
        LDI is defined as the probability that two randomly selected points in the landscape are situated in two different patches of the habitat (Jaeger 2000; Mcgarigal 2015).
        
        Args:
        
        Returns:
            An array with LDI values.
        """
        
        from osgeo import gdal
        
        woody_cover = self.getWoodyCover(threshold, lee_filter = lee_filter)
        
        # Reduce to an image of 100 x 100 pixels
        factor = 100
        
        LDI = np.zeros((factor, factor))
        
        for ymin in np.arange(0, self.ySize, self.ySize / factor):
            ymax = ymin + (self.ySize / factor)
            for xmin in np.arange(0, self.xSize, self.xSize / factor):
                xmax = xmin + (self.xSize / factor)
                
                this_patch = woody_cover[ymin:ymax, xmin:xmax]
                
                # Get connected patches of forest and nonforest
                forest_area, forest_id = getContiguousAreas(this_patch, True)
                nonforest_area, nonforest_id = getContiguousAreas(this_patch, False)
                
                # Supply a unique ID to each habitat patch, whether forest or nonforest.
                unique_ids = forest_id.copy()
                unique_ids[nonforest_id != 0] = forest_id[nonforest_id != 0] + (nonforest_id[nonforest_id != 0] + np.max(forest_id) + 1)
                                                
                LDI[ymin / (self.ySize / factor), xmin / (self.xSize / factor)] = calculateLDI(unique_ids)
        
        if output: self.__outputGeoTiff(LDI, 'LDI', geo_t = self.__shrinkGeoT(factor), nodata = None)
        
        return LDI
    
    def calculateTWC(self, threshold, lee_filter = True, output = False):
        """
        Total woody cover (TWC) describes the propotion of woody cover in a downsampled image.
        """
        
        from osgeo import gdal
        
        woody_cover = self.getWoodyCover(threshold, lee_filter = lee_filter)
        
        # Reduce to an image of 100 x 100 pixels
        factor = 100
        
        TWC = np.zeros((factor, factor))
        
        for ymin in np.arange(0, self.ySize, self.ySize / factor):
            ymax = ymin + (self.ySize / factor)
            for xmin in np.arange(0, self.xSize, self.xSize / factor):
                xmax = xmin + (self.xSize / factor)
                
                this_patch = woody_cover[ymin:ymax, xmin:xmax]
                
                # Get proportion of woody cover in patch
                TWC[ymin / (self.ySize / factor), xmin / (self.xSize / factor)] = float(this_patch.sum()) / ((self.ySize / factor) * (self.xSize / factor))
        
        if output: self.__outputGeoTiff(TWC, 'TWC', geo_t = self.__shrinkGeoT(factor), nodata = None)
        
        return TWC
    
    def __outputGeoTiff(self, data, output_name, output_dir = os.getcwd(), dtype = 6, geo_t = '', proj = '', nodata = ''):
        """
        Output a GeoTiff file.
        """
        
        # Where nothing specified, use standard geo_t and proj definitions. This is to allow these to be overwritten for functions that reduce output size.
        if geo_t == '': geo_t = self.geo_t
        if proj == '': proj = self.proj
        if nodata == '': nodata = self.nodata
        
        # Generate a standardised filename
        filename = '%s_%s%s.tif'%(output_name, self.hem_NS + str(abs(self.lat)).zfill(2), self.hem_EW + str(abs(self.lon)).zfill(3))
        
        # Write to disk
        outputGeoTiff(data, filename, geo_t, proj, output_dir = output_dir, dtype = dtype)




if __name__ == '__main__':
    
    data_dir = '/home/sbowers3/DATA/ALOS_data/ALOS_mosaic/gorongosa/'
    output_dir = '/home/sbowers3/DATA/ALOS_data/ALOS_mosaic/gorongosa/'

    t1 = 2007
    t2 = 2010

    lat = -19#-9#-11
    lon = 33#34#39
    
    data_t1 = ALOS(lat, lon, t1, data_dir)
    LDI_t1 = data_t1.calculateLDI(15., output=True)
    """
    data_t2 = ALOS(lat, lon, t2, data_dir)
    LDI_t2 = data_t2.calculateLDI(15., lee_filter = True)
    
    plt.subplot(121)
    plt.imshow(LDI_t1, vmin=0, vmax=1, cmap='viridis', interpolation = 'nearest')
    plt.colorbar()
    plt.subplot(122)
    plt.imshow(LDI_t2, vmin=0, vmax=1, cmap='viridis', interpolation = 'nearest')
    plt.colorbar()
    plt.show()
    """
    """
    data_t2 = ALOS(lat, lon, t2, data_dir)
    
    # Build masks (optionally with buffers)
    lakes = '/home/sbowers3/DATA/GIS_data/mozambique/diva/MOZ_wat/MOZ_water_areas_dcw.shp'
    rivers = '/home/sbowers3/DATA/GIS_data/mozambique/diva/MOZ_wat/MOZ_water_lines_dcw.shp'
    mozambique = '/home/sbowers3/DATA/GIS_data/mozambique/diva/MOZ_adm/MOZ_adm0.shp'
    wdpa = '/home/sbowers3/DATA/GIS_data/mozambique/WDPA/MOZ_WDPA.shp'
    
    lake_mask = rasterizeShapefile(data_t1, lakes, buffer_size = 0.005)
    river_mask = rasterizeShapefile(data_t1, rivers, buffer_size = 0.005)
    water_mask = np.logical_or(river_mask, lake_mask)
    
    moz_mask = rasterizeShapefile(data_t1, mozambique)
    wdpa_mask = rasterizeShapefile(data_t1, wdpa)

    data_t1.mask = np.logical_or(np.logical_or(data_t1.mask, water_mask), moz_mask == False)
    data_t2.mask = np.logical_or(np.logical_or(data_t2.mask, water_mask), moz_mask == False)
        
    overviewFigure(data_t1, data_t2, output_dir = output_dir)
    """

"""



"""

"""
for lat in np.arange(-19,-14,1):
    for lon in np.arange(30,35,1):
        HVfile_t1, maskfile_t1 = generateFilenames(lat, lon, t1, data_dir)
        HVfile_t2, maskfile_t2 = generateFilenames(lat, lon, t2, data_dir)
        
        ds_t1, data_t1 = openTile(HVfile_t1, maskfile_t1)
        ds_t2, data_t2 = openTile(HVfile_t2, maskfile_t2)

        gamma0_t1 = calibrateToGamma0(data_t1)
        gamma0_t2 = calibrateToGamma0(data_t2)
        
        gamma0_change = (gamma0_t2 - gamma0_t1) / (t2 - t1)
        
        gamma0_pcChange = 100 * (gamma0_change / gamma0_t1)
        
        ##outputGeoTiff(gamma0_pcChange, ds_t1.GetGeoTransform(), output_dir + 'test_outputs/', output_name = 'gamma0_change')
"""








def _buildMap(fig, ax, data, lat, lon, title ='', cbartitle = '', vmin = 10., vmax = 60., cmap = 'YlGn'):
    """
    Builds a standardised map for overviewFigure().
    """
    
    import matplotlib.pyplot as plt
        
    im = ax.imshow(data, vmin = vmin, vmax = vmax, cmap = cmap, interpolation = 'nearest')
    
    ax.set_xticks(np.arange(0,4501,450))
    ax.set_yticks(np.arange(0,4501,450))
    ax.set_xticklabels(np.arange(lon, lon + 1.01, 0.1))
    ax.set_yticklabels(np.arange(lat, lat - 1.01, - 0.1))
    ax.tick_params(labelsize = 5)
    ax.set_xlabel('Longitude', fontsize = 5)
    ax.set_ylabel('Latitude', fontsize = 5)
    ax.set_title(title, fontsize = 8)
    
    cbar = fig.colorbar(im, ax = ax, fraction = 0.046, pad = 0.04)
    cbar.ax.tick_params(labelsize = 6)
    cbar.set_label(cbartitle, fontsize = 7)
    

def overviewFigure(data_t1, data_t2, output_dir = os.getcwd(), output_name = 'overview'):
    """overviewFigure(data_t1, data_t2, t1, t2, geo_t, output_dir = os.getcwd())
    
    Generate an overview image showing biomass and proportional biomass change for the tile being processed.
    
    Args:
        data_t1:
        data_t2: 
        output_name: Optionally specify an output string to precede output file. Defaults to 'overview'.
    """
    
    import matplotlib.pyplot as plt
    
    assert data_t1.geo_t == data_t2.geo_t, "The two ALOS tiles must be from the same location."
    
    # Get upper left longitude and latitude from GeoMatrix
    lon, lat = data_t1.lat, data_t1.lon
    
    # Update masks to exclude areas outisde forest definition. Good for visualisation
    AGB_t1 = data_t1.getAGB()
    AGB_t2 = data_t2.getAGB()
    
    AGB_t1 = np.ma.array(AGB_t1, mask = np.logical_or(AGB_t1.mask, AGB_t1 < 10.))
    AGB_t2 = np.ma.array(AGB_t2, mask = np.logical_or(AGB_t2.mask, AGB_t1 < 10.))
        
    AGB_change = (AGB_t2 - AGB_t1) / (data_t2.year - data_t1.year) # tC/ha/yr

    AGB_pcChange = 100 * (AGB_change / AGB_t1) # %/yr
    
    fig = plt.figure(figsize = (7, 6))
    
    # Plot a map of AGB at t1
    ax1 = fig.add_subplot(2, 2, 1)
    _buildMap(fig, ax1, AGB_t1, lat, lon, title = 'AGB %s'%str(t1), cbartitle = 'tC/ha')
    
    # Plot a map of AGB at t2
    ax2 = fig.add_subplot(2, 2, 2)
    _buildMap(fig, ax2, AGB_t2, lat, lon, title = 'AGB %s'%str(t2), cbartitle = 'tC/ha')    
    
    # Plot a map of absolute AGB change   
    ax3 = fig.add_subplot(2, 2, 3)
    _buildMap(fig, ax3, AGB_change, lat, lon, title = 'AGB change (%s-%s)'%(str(t1),str(t2)),
              cbartitle = 'tC/ha/yr', vmin = -10., vmax = 10., cmap = 'RdBu')    
    
    # Plot a map of % AGB change
    ax4 = fig.add_subplot(2, 2, 4)
    _buildMap(fig, ax4, AGB_pcChange, lat, lon, title = 'AGB change (%s-%s)'%(str(t1),str(t2)),
              cbartitle = '%/yr', vmin = -50., vmax = 50., cmap = 'RdBu')    
    
    plt.tight_layout()
    
    # Determine filename
    hem_NS = 'S' if lat < 0 else 'N'
    hem_EW = 'W' if lon < 0 else 'E'
    
    output_path = '%s/%s_%s%s.png'%(output_dir, output_name, data_t1.hem_NS + str(abs(lat)).zfill(2), 
                                    data_t1.hem_EW + str(abs(lon)).zfill(3))
    
    plt.savefig(output_path, dpi = 150)
    plt.close()


def dilateMask(mask, buffer_px):
    """
    Dilate a boolean (True/False) numpy array by a specified number of pixels.
        
    Args:
        mask: A boolean (True/False) numpy array, with 'True' representing locations to add a buffer.
        buffer_px: A number of pixels to add around each 'True' array element.

    Returns:
        The mask array with dilated 'True' locations.
    """
        
    mask_dilated = ndimage.morphology.binary_dilation(mask, iterations = buffer_px)
    
    return mask_dilated


def _coordinateTransformer(shp):
    """
    Generates function to transform coordinates from a source shapefile CRS to EPSG.
    
    Args:
        shp: Path to a shapefile.
    
    Returns:
        A function that transforms shapefile points to EPSG.
    """
    
    from osgeo import ogr, osr
        
    driver = ogr.GetDriverByName('ESRI Shapefile')
    ds = driver.Open(shp)
    layer = ds.GetLayer()
    spatialRef = layer.GetSpatialRef()
    
    # Create coordinate transformation
    inSpatialRef = osr.SpatialReference()
    inSpatialRef.ImportFromWkt(spatialRef.ExportToWkt())

    outSpatialRef = osr.SpatialReference()
    outSpatialRef.ImportFromEPSG(4326)

    coordTransform = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)
    
    return coordTransform


def _world2Pixel(geo_t, x, y, buffer_size = 0):
    """
    Uses a gdal geomatrix (ds.GetGeoTransform()) to calculate the pixel location of a geospatial coordinate.
    Modified from: http://geospatialpython.com/2011/02/clip-raster-using-shapefile.html.
    
    Args:
        geo_t: A gdal geoMatrix (ds.GetGeoTransform().
        x: x coordinate in map units.
        y: y coordinate in map units.
        buffer_size: Optionally specify a buffer size. This is used when a buffer has been applied to extend all edges of an image, as in rasterizeShapfile().
    
    Returns:
        A tuple with pixel/line locations for each input coordinate.
    """
    ulX = geo_t[0] - buffer_size
    ulY = geo_t[3] + buffer_size
    xDist = geo_t[1]
    yDist = geo_t[5]
    
    pixel = int((x - ulX) / xDist)
    line = int((y - ulY) / yDist)
    
    return (pixel, line)


def getField(shp, field):
    '''
    Get values from a field in a shapefile attribute table.
    
    Args:
        shp: A string pointing to a shapefile
        field: A string with the field name of the attribute of interest
    
    Retuns:
        An array containing all the values of the specified attribute
    '''
    
    import shapefile
    
    assert os.path.isfile(shp), "Shapefile %s does not exist."%shp
    
    # Read shapefile
    sf = shapefile.Reader(shp)
    
    # Get the column number of the field of interest
    for n, this_field in enumerate(sf.fields[1:]):
        
        fieldname = this_field[0]
        
        if fieldname == field:
            
            field_n = n
    
    assert 'field_n' in locals(), "Attribute %s not found in shapefile."%str(field)
    
    value_out = []
    
    # Cycle through records:
    for s in sf.records():
        value_out.append(s[field_n])
    
    return np.array(value_out)


def getBBox(shp, field, value):
    '''
    Get the bounding box of a shape in a shapefile.
    
    Args:
        shp: A string pointing to a shapefile
        field: A string with the field name of the attribute of interest
        value: The value of the field for the shape of interest
    
    Retuns:
        An list with the bounding box in the format [minlon, minlat, maxlon, maxlat]
    '''
        
    import shapefile
    
    assert os.path.isfile(shp), "Shapefile %s does not exist."%shp

    assert (np.sum(getField(shp, field) == value) > 1) == False, "The value name in a field must be unique. In the field %s there are %s records with value %s."%(str(field), str(np.sum(getField(shp, field) == value)), str(value))
    
    # Read shapefile
    sf = shapefile.Reader(shp)
    
    shapes = np.array(sf.shapes())
    
    # Get bounding box    
    bbox = shapes[getField(shp, field) == value][0].bbox
    
    return bbox
    



def rasterizeShapefile(data, shp, buffer_size = 0., field = None, value = None):
    """
    Rasterize points, lines or polygons from a shapefile to match ALOS mosaic data.
        
    Args:
        data: An ALOS object
        shp: Path to a shapefile consisting of points, lines and/or polygons. This does not have to be in the same projection as ds
        buffer_size: Optionally specify a buffer to add around features of the shapefile, in decimal degrees.

    Returns:
        A numpy array with a boolean mask delineating locations inside (True) and outside (False) the shapefile [and optional buffer].
    """
    
    import shapefile
    from osgeo import gdalnumeric
    
    assert np.logical_or(np.logical_and(field == None, value == None), np.logical_and(field != None, value != None)), "If specifying field or value, both must be defined. At present, field = %s and value = %s"%(str(field), str(value))
    
    # Determine size of buffer to place around lines/polygons
    buffer_px = int(round(buffer_size / data.geo_t[1]))
    
    # Create output image. Add a buffer around the image array equal to the maxiumum dilation size. This means that features just outside ALOS tile extent can contribute to dilated mask.
    rasterPoly = Image.new("I", (data.ySize + (buffer_px * 2), data.xSize + (buffer_px * 2)), 0)
    rasterize = ImageDraw.Draw(rasterPoly)
    
    # The shapefile may not have the same CRS as ALOS mosaic data, so this will generate a function to reproject points.
    coordTransform = _coordinateTransformer(shp)
    
    # Read shapefile
    sf = shapefile.Reader(shp) 
    
    # Get shapes
    shapes = np.array(sf.shapes())
    
    # If extracting a mask for just a single field.
    if field != None:
        
        shapes = shapes[getField(shp, field) == value]
            
    # For each shape in shapefile...
    for shape in shapes:
                
        # Get shape bounding box
        sxmin, symin, sxmax, symax = shape.bbox
        
        # Transform bounding box points
        sxmin, symin, z = coordTransform.TransformPoint(sxmin, symin)
        sxmax, symax, z = coordTransform.TransformPoint(sxmax, symax)
        
        # Go to the next record if out of bounds
        if sxmax < data.geo_t[0] - buffer_size: continue
        if sxmin > data.geo_t[0] + (data.geo_t[1] * data.ySize) + buffer_size: continue
        if symax < data.geo_t[3] + (data.geo_t[5] * data.xSize) + buffer_size: continue
        if symin > data.geo_t[3] - buffer_size: continue
        
        #Separate polygons with list indices
        n_parts = len(shape.parts) #Number of parts
        indices = shape.parts #Get indices of shapefile part starts
        indices.append(len(shape.points)) #Add index of final vertex
        
        for part in range(n_parts):
            
            start_index = shape.parts[part]
            end_index = shape.parts[part+1]
            
            points = shape.points[start_index:end_index] #Map coordinates
            pixels = [] #Pixel coordinantes
            
            # Transform coordinates to pixel values
            for p in points:
                
                # First update points from shapefile projection to ALOS mosaic projection
                lon, lat, z = coordTransform.TransformPoint(p[0], p[1])

                # Then convert map to pixel coordinates using geo transform
                pixels.append(_world2Pixel(data.geo_t, lon, lat, buffer_size = buffer_size))

            # Draw the mask for this shape...
            # if a point...
            if shape.shapeType == 0:
                rasterize.point(pixels, 1)

            # a line...
            elif shape.shapeType == 3:
                rasterize.line(pixels, 1)
  
            # or a polygon.
            elif shape.shapeType == 5:  
                rasterize.polygon(pixels, 1)
        
    #Converts a Python Imaging Library array to a gdalnumeric image.
    mask = gdalnumeric.fromstring(rasterPoly.tobytes(),dtype=np.uint32)
    mask.shape = rasterPoly.im.size[1], rasterPoly.im.size[0]
    
    # If any buffer pixels are slected, dilate the masked area by buffer_px pixels
    if buffer_px > 0:
        mask = dilateMask(mask, buffer_px)
    
    # Get rid of image buffer
    mask = mask[buffer_px:mask.shape[0]-buffer_px, buffer_px:mask.shape[1]-buffer_px]
    
    # Get rid of record numbers
    mask = mask > 0
    
    return mask




def getTilesInShapefile(shp):
    """
    Identify all the ALOS tiles that fall within a shapefile.
    
    Args:
        shp: Path to a shapefile consisting of polygons. This can be in any projection.
    
    Returns:
        The lat/lon indicators of which ALOS tiles are covered by the shapefile
    """
    
    import shapefile
    
    # The shapefile may not have the same CRS as ALOS mosaic data, so this will generate a function to reproject points.    
    coordTransform = _coordinateTransformer(shp)
    
    lats, lons = [], []
    tiles_to_include = set([])
        
    for shape in shapefile.Reader(shp).shapes():
        
        # Get the bbox for each shape in the shapefile
        lonmin, latmin, lonmax, latmax = shape.bbox
        
        # Transform points to WGS84
        lonmin, latmin, z = coordTransform.TransformPoint(lonmin, latmin)
        lonmax, latmax, z = coordTransform.TransformPoint(lonmax, latmax)
        
        # Get the tiles that cover the area of the shapefile
        latrange = range(int(math.ceil(latmin)), int(math.ceil(latmax)+1), 1)
        lonrange = range(int(math.floor(lonmin)), int(math.floor(lonmax)+1), 1)
        tiles = list(itertools.product(latrange,lonrange))
        
        # Add them to tiles_to_include if not already tere
        [tiles_to_include.add(t) for t in tiles]
    
    return sorted(list(tiles_to_include))
    






