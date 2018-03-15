#!/usr/bin/env python

import argparse
import datetime as dt
import itertools
import math
import numpy as np
import os
from scipy import ndimage
import scipy.stats as stats
import skimage.measure

import matplotlib.pyplot as plt
import pdb

import biota.filter
import biota.indices
import biota.IO
import biota.mask



class LoadTile(object):
    """
    An ALOS mosaic tile.
    Mosaic tiles have the following properties:

    Attributes:
        lat:
        lon:
        DN: An array of uncalibrated digital numbers from the ALOS tile.
        mask:
    """
        
    def __init__(self, dataloc, lat, lon, year, downsample_factor = 1, lee_filter = False, output_dir = os.getcwd()):
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
        assert downsample_factor >= 1 and type(downsample_factor) == int, "Downsampling factor must be an integer greater than 1."
        assert type(lee_filter) == bool, "Option lee_filter must be set to 'True' or 'False'."
        assert os.path.isdir(output_dir), "Specified output directory (%s) does not exist"%str(output_dir)
        
        self.lat = lat
        self.lon = lon
        self.year = year
        self.downsample_factor = downsample_factor
        self.lee_filter = lee_filter
        
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
        
        # Set up locations for file output
        self.output_pattern = self.__getOutputPattern()
        self.output_dir = output_dir
        
        # Stop of the ALOS tile doesn't exist
        if not os.path.isfile(self.HV_path): 
            raise IOError('ALOS tile does not exist in the file system.')

        # Get Raster size
        self.ySize, self.xSize = self.__getSize(self.HH_path)
                
        # Get GDAL geotransform and projection
        self.geo_t = self.__getGeoT()
        self.proj = self.__getProj()
        
        # Get the raster extent
        self.extent = self.__getExtent()
        
        # Determine x and y resolution in meters
        self.xRes = self.__getXRes()
        self.yRes = self.__getYRes()
        
        # Record the nodata value
        self.nodata = self.__getNodata()
        
        # Get Equivalent Number of Looks
        self.nLooks = self.__getnLooks()
        
        # Update image metadata where a degree of resampling is included
        if self.downsample_factor != 1:
            self.__rebin()
        
        # Load DN, mask, and day of year
        self.mask = self.getMask()
                
        
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
    
    def __getOutputPattern(self):
        """
        Generates a filename pattern for data output.
        """
        
        return '%s_%s%s_%s.tif'%('%s', self.hem_NS + str(abs(self.lat)).zfill(2), self.hem_EW + str(abs(self.lon)).zfill(3), '%s')
    
    def __getNodata(self):
        """
        Return the nodata value for ALOS array. Should always be 0, so this function is a placeholder should the format change.
        """
        
        return 999999
    
    def __getnLooks(self)    :
        """
        Define number of looks. Out of the box, this is 16 for the ALOS mosaic.
        """
        
        return 16

    def __getSize(self, filepath):
        """
        Determines the size of a tile.
        """
                
        y_size, x_size = biota.IO.loadSize(filepath)
        
        return y_size, x_size
    
    def __getGeoT(self):
        """
        Generates a GDAL GeoTransform for tile.
        """
        
        geo_t = (float(self.lon), 1./self.xSize, 0., float(self.lat), 0., -1./self.xSize)
               
        return geo_t
    
    def __getProj(self):
        """
        Fetches projection info for tile.
        """
        
        proj = biota.IO.loadProjection(self.__getHHPath())
        
        return proj
    
    def __getExtent(self):
        """
        Calculates the raster exent in the format [xmin. ymin, xmax, ymax] using the geoTransform and raster size.
        """
        
        xmin = self.geo_t[0]
        ymax = self.geo_t[3]
        xmax = self.geo_t[0] + (self.geo_t[1] * self.xSize)
        ymin = self.geo_t[3] + (self.geo_t[5] * self.ySize)
        
        return (xmin, ymin, xmax, ymax)    
    
    def __getXRes(self):
        """
        Approximates the resolution of a pixel in meters.
        From: https://en.wikipedia.org/wiki/Geographic_coordinate_system
        """
        
        lat = math.radians(self.lat - 0.5)
                
        m_per_degree = (111412.84 * math.cos(lat)) - (93.5 * math.cos(3 * lat)) + (0.118 * math.cos(5 * lat))
        
        return round(m_per_degree / self.xSize, 3)
        
    def __getYRes(self):
        """
        Approximates the resolution of a pixel in meters.
        From: https://en.wikipedia.org/wiki/Geographic_coordinate_system
        """

        lat = math.radians(self.lat - 0.5)
        
        m_per_degree = 111132.92 - (559.82 * math.cos(2 * lat)) + (1.175 * math.cos(4 * lat))
        
        return round(m_per_degree / self.ySize, 3)
    
    def shrinkGeoT(self, downsample_factor):
        """
        Function to modify gdal geo_transform to output at correct resolution for downsampled products.
        """
        
        geo_t = list(self.geo_t)
                
        geo_t[1] = 1. / int(round(self.xSize / float(downsample_factor)))
        geo_t[-1] = (1. / int(round(self.ySize / float(downsample_factor)))) * -1
        
        return tuple(geo_t)
    
    def __rebin(self):
        """
        Reduce array size to a new resolution by taking mean, to increase equivalent number of looks.
        # From https://stackoverflow.com/questions/8090229/resize-with-averaging-or-rebin-a-numpy-2d-array
        """
        
        # Update geo_t
        self.geo_t = self.shrinkGeoT(self.downsample_factor)
        
        # Take the mean of X/Y resolutions in metres.
        native_resolution = (self.xRes + self.yRes) / 2.

        # Update extents to new resolution
        new_xSize = int(round(self.xSize / float(self.downsample_factor)))
        new_ySize = int(round(self.ySize / float(self.downsample_factor)))
        
        self.xRes = (self.xSize * self.xRes) / new_xSize
        self.yRes = (self.ySize * self.yRes) / new_xSize
        self.xSize = new_xSize
        self.ySize = new_ySize
                
        # Update number of looks
        self.nLooks = self.nLooks * (self.downsample_factor ** 2)
    
    def getMask(self, masked_px_count = False):
        """
        Loads the mask into a numpy array.
        """
                
        mask = biota.IO.loadArray(self.mask_path) != 255
        
        # If resampling, this removes the mask from any pixel with > 75 % data availability.
        if self.downsample_factor != 1 and not masked_px_count:
            mask = skimage.measure.block_reduce(mask, (self.downsample_factor, self.downsample_factor), np.sum) > ((self.downsample_factor ** 2) * 0.25)
     
        # This is an option to return the sum of masked pixels. It's used to downsample the DN array.
        if self.downsample_factor != 1 and masked_px_count:
            mask = skimage.measure.block_reduce(mask, (self.downsample_factor, self.downsample_factor), np.sum)
        
        return mask
    
    def updateMask(self, shp, buffer_size = 0.):
        """
        Function to add further pixels to mask based on shapefiles and buffers.
        """
        
        # Rasterize the shapefile, optionally with a buffer
        shp_mask = biota.mask.rasterizeShapefile(self, shp, buffer_size = buffer_size)
        
        # Add the rasterized shapefile to the mask
        self.mask = np.logical_or(self.mask, shp_mask)     
        
    def resetMask(self):
        """
        Function to reset a mask to the default.
        """
        
        self.mask = self.getMask()
       
        
    def getDN(self, polarisation = 'HV'):
        """
        Loads DN (raw) values into a numpy array.
        """
        
        assert polarisation == 'HH' or polarisation == 'HV', "polarisation must be either 'HH' or 'HV'."
               
        if polarisation == 'HV':
            DN = biota.IO.loadArray(self.HV_path)
        else:
            DN = biota.IO.loadArray(self.HH_path)
                
        # Rebin DN with mean. This is acceptable as the DN values are on a linear scale.
        if self.downsample_factor != 1:
            
            # Load the sum of DNs
            DN_sum = skimage.measure.block_reduce(DN, (self.downsample_factor, self.downsample_factor), np.sum)
            
            # Amd the sum of masked pixels
            mask_sum = self.getMask(masked_px_count = True)
            
            # Divide the sum of DNs by the sum of unmasked pixels to get the mean DN value
            DN = np.zeros_like(DN_sum)
            DN[self.mask == False] = (DN_sum.astype(np.float)[self.mask == False] / ((self.downsample_factor ** 2) - mask_sum[self.mask == False])).astype(np.int)
        
        return np.ma.array(DN, mask = self.mask)
        
    def getDOY(self, output = False):
        """
        Loads date values into a numpy array.
        """
        
        day_after_launch = biota.IO.loadArray(self.date_path)
        
        # Get list of unique dates in ALOS tile
        unique_days = np.unique(day_after_launch)
        
        # Days are counted from the launch date
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
        
        # If required, downsample the DOY array      
        if self.downsample_factor != 1:
            
            # Take the latest DOY where > 1 date, as there's no numpy function for mode
            DOY = skimage.measure.block_reduce(DOY, (self.downsample_factor, self.downsample_factor), np.max)
        
        if output: self.__outputGeoTiff(DOY, 'DOY', dtype = gdal.GDT_Int16)

        return DOY
    
    def getGamma0(self, polarisation = 'HV', units = 'natural', output = False):
        """
        Calibrates data to gamma0 (baskscatter) in decibels or natural units.
        """
        
        assert units == 'natural' or units == 'decibels', "Units must be 'natural' or 'decibels'. You input %s."%units
        
        # Calibrate DN to units of dB
        gamma0 = 10 * np.ma.log10(self.getDN(polarisation = polarisation).astype(np.float) ** 2) - 83. # units = decibels
        
        # Apply filter based on dB values
        if self.lee_filter:
            gamma0 = biota.filter.enhanced_lee_filter(gamma0, n_looks = self.nLooks)
        
        # Convert to natural units where specified
        if units == 'natural': gamma0 = 10 ** (gamma0 / 10.)
        
        # Keep masked values tidy
        gamma0.data[self.mask] = self.nodata
        
        if output: self.__outputGeoTiff(gamma0, 'Gamma0')
        
        return gamma0
    
    def getAGB(self, output = False):
        """
        Calibrates data to aboveground biomass (AGB).
        Placeholder equation to calibrate backscatter (gamma0) to AGB (tC/ha).
        """
        
        # ALOS-1
        if self.satellite == 'ALOS-1':
            AGB = 715.667 * self.getGamma0(units = 'natural', polarisation = 'HV') - 5.967
        
        # ALOS-2 (to calculate)
        elif self.satellite == 'ALOS-2':
            AGB = 715.667 * self.getGamma0(units = 'natural', polarisation = 'HV') - 5.967
            
        else:       
            raise ValueError("Unknown satellite named '%s'. self.satellite must be 'ALOS-1' or 'ALOS-2'."%self.satellite)
        
        # Keep masked values tidy
        AGB.data[self.mask] = self.nodata
        
        if output: self.__outputGeoTiff(AGB, 'AGB')
        
        return AGB
    
    def getWoodyCover(self, forest_threshold = 10., min_forest_area = 0., output = False):
        """
        Get woody cover, based on a threshold of AGB.
        min_area in ha
        """
        
        woody_cover = self.getAGB() >= float(forest_threshold)
               
        if min_forest_area > 0:
            
            # Calculate number of pixels in min_area (assuming input is given in hecatres)
            min_pixels = int(round(min_forest_area / (self.yRes * self.xRes * 0.0001)))
            
            # Remove pixels that aren't part of a forest block of size at least min_pixels
            contiguous_area, _ = biota.indices.getContiguousAreas(woody_cover, True, min_pixels = min_pixels)
            woody_cover[contiguous_area == False] = False
        
        woody_cover.data[self.mask] = self.nodata
        
        if output: self.__outputGeoTiff(woody_cover * 1, 'WoodyCover', dtype = gdal.GDT_Int32)
        
        return woody_cover
    
    def getForestPatches(self, forest_threshold = 10., min_forest_area = 0., output = False):
        """
        Get numbered forest patches, based on woody cover threshold.
        """
                
        woody_cover = self.getWoodyCover(forest_threshold = forest_threshold, min_forest_area = 0.)
        
        # Calculate number of pixels in min_area
        min_pixels = int(round(min_forest_area / (self.yRes * self.xRes * 0.0001)))
                
        # Get areas that meet that threshold
        _, location_id = biota.indices.getContiguousAreas(woody_cover, True, min_pixels = min_pixels)
        
        # Tidy up masked pixels     
        location_id.data[location_id.mask] = self.nodata
        
        if output: self.__outputGeoTiff(location_id * 1, 'ForestPatches', dtype = gdal.GDT_Int32)
        
        return location_id
    
    def __outputGeoTiff(self, data, output_name, dtype = 6):
        """
        Output a GeoTiff file.
        """
                
        # Generate a standardised filename
        filename = self.output_pattern%(output_name, str(self.year))
        
        # Write to disk
        biota.IO.outputGeoTiff(data, filename, self.geo_t, self.proj, output_dir = self.output_dir, dtype = dtype, nodata = self.nodata)
