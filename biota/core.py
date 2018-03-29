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

"""
These are two classes for loading individual tiles and compliling them into a change object.
"""



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
        
    def __init__(self, data_dir, lat, lon, year, forest_threshold = 10., area_threshold = 0., downsample_factor = 1, lee_filter = False, output_dir = os.getcwd()):
        """
        Loads data and metadata for an ALOS mosaic tile.
        """
        
        from osgeo import gdal
        
        # Test that inputs are of reasonable lats/lons/years
        assert type(lat) == int, "Latitude must be an integer."
        assert lat < 90. or lat > -90., "Latitude must be between -90 and 90 degrees."
        assert type(lon) == int, "Longitude must be an integer."
        assert lon < 180. or lon > -180., "Longitude must be between -180 and 180 degrees."
        assert type(year) == int, "Year must be an integer."
        assert (year >= 2007 and year <= 2010) or (year >= 2015 and year <= dt.datetime.now().year), "Years must be in the range 2007 - 2010 and 2015 - present. Your input year was %s."%str(year)
        assert downsample_factor >= 1 and type(downsample_factor) == int, "Downsampling factor must be an integer greater than 1."
        assert type(lee_filter) == bool, "Option lee_filter must be set to 'True' or 'False'."
        assert os.path.isdir(os.path.expanduser(data_dir)), "Specified data directory (%s) does not exist"%str(data_dir)
        assert os.path.isdir(os.path.expanduser(output_dir)), "Specified output directory (%s) does not exist"%str(output_dir)
        assert type(forest_threshold) == float or type(forest_threshold) == int, "Forest threshold must be numeric."
        assert type(area_threshold) == float or type(area_threshold) == int, "Area threshold must be numeric."
        
        self.lat = lat
        self.lon = lon
        self.year = year
        self.downsample_factor = downsample_factor
        self.lee_filter = lee_filter
        self.forest_threshold = forest_threshold
        self.area_threshold = area_threshold
        
        # Deterine hemispheres
        self.hem_NS = 'S' if lat < 0 else 'N'
        self.hem_EW = 'W' if lon < 0 else 'E'
        
        # Determine whether ALOS-1 or ALOS-2
        self.satellite = self.__getSatellite()
        
        # Determine filenames
        self.data_dir = os.path.expanduser(data_dir.rstrip('/'))
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
        self.nodata_byte = self.__getNodata(dtype = gdal.GDT_Byte)
        
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
        
        assert os.path.isdir(self.data_dir), "Data location must be a directory."
        
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
        directory = self.data_dir + '/' + name_pattern%(lat_dir, lon_dir, str(self.year)[-2:], 'MOS') + '/'
        
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
        
        return '%s_%s_%s%s.tif'%('%s', str(self.year), self.hem_NS + str(abs(self.lat)).zfill(2), self.hem_EW + str(abs(self.lon)).zfill(3))
    
    def __getNodata(self, dtype = 6):
        """
        Return nodata values
        """
        
        from osgeo import gdal
        
        # Generate a nodata value
        if dtype == gdal.GDT_Byte:
            nodata = 99
        else:
            nodata = 999999
        
        return nodata
    
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
        
        # If resampling, this removes the mask from any pixel with >= 75 % data availability.
        if self.downsample_factor != 1 and not masked_px_count:
            mask = skimage.measure.block_reduce(mask, (self.downsample_factor, self.downsample_factor), np.sum) >= ((self.downsample_factor ** 2) * 0.25)
     
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
        
        from osgeo import gdal
        
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
        
        if output: self.__outputGeoTiff(DOY, 'DOY', dtype = gdal.GDT_Int32)

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
        
        # Don't rerun processing if already present in memory
        if not hasattr(self, 'AGB'):
                    
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
            
            # Save output to class
            self.AGB = AGB
        
        if output: self.__outputGeoTiff(self.AGB, 'AGB')
        
        return self.AGB

    def getWoodyCover(self, output = False):
        """
        Get woody cover, based on a threshold of AGB.
        min_area in ha
        """
        
        from osgeo import gdal
        
        # Don't rerun processing if already present in memory
        if not hasattr(self, 'WoodyCover'):
            
            WoodyCover = self.getAGB() >= float(self.forest_threshold)
                
            if self.area_threshold > 0:
                
                # Calculate number of pixels in min_area (assuming input is given in hecatres)
                min_pixels = int(round(self.area_threshold / (self.yRes * self.xRes * 0.0001)))
                
                # Remove pixels that aren't part of a forest block of size at least min_pixels
                contiguous_area, _ = biota.indices.getContiguousAreas(WoodyCover, True, min_pixels = min_pixels)
                
                WoodyCover.data[contiguous_area == False] = False
            
            WoodyCover.data[self.mask] = False
            
            # Save output to class
            self.WoodyCover = WoodyCover
        
        if output: 
            
            # Convert data to integer type for output
            WoodyCover_out = WoodyCover.astype(np.uint8)
            
            self.__outputGeoTiff(WoodyCover_out, 'WoodyCover', dtype = gdal.GDT_Byte)
        
        return self.WoodyCover
    
    def getForestPatches(self, output = False):
        """
        Get numbered forest patches, based on woody cover threshold.
        """
        
        from osgeo import gdal
        
        # Don't rerun processing if already present in memory
        if not hasattr(self, 'ForestPatches'):
            
            WoodyCover = self.getWoodyCover()
            
            # Calculate number of pixels in min_area
            min_pixels = int(round(self.area_threshold / (self.yRes * self.xRes * 0.0001)))
            
            # Get areas that meet that threshold
            _, location_id = biota.indices.getContiguousAreas(WoodyCover, True, min_pixels = min_pixels)
            
            # Tidy up masked pixels     
            location_id.data[location_id.mask] = self.nodata
            
            # Save output to class
            self.ForestPatches = location_id
        
        if output: self.__outputGeoTiff(location_id * 1, 'ForestPatches', dtype = gdal.GDT_Int32)
        
        return self.ForestPatches
    
    def __outputGeoTiff(self, data, output_name, dtype = 6):
        """
        Output a GeoTiff file.
        """
        
        from osgeo import gdal
        
        # Generate a standardised filename
        filename = self.output_pattern%output_name
        
        # Generate a nodata value appropriate for datatype
        nodata = self.__getNodata(dtype = dtype)
        
        # Write to disk
        biota.IO.outputGeoTiff(data, filename, self.geo_t, self.proj, output_dir = self.output_dir, dtype = dtype, nodata = nodata)


class LoadChange(object):
    """
    Input is two mosaic tiles from LoadTile, will output maps and change statistics.
    """
        
    def __init__(self, data_t1, data_t2, change_intensity_threshold = 0.2, change_area_threshold = 0, output_dir = os.getcwd(), output = False):
        '''
        Initialise
        '''
        
        from osgeo import gdal
        
        self.data_t1 = data_t1
        self.data_t2 = data_t2
        
        # Ensure that tiles are compatable
        self.__testTiles()
        
        # Get basic properties
        self.lat = data_t1.lat
        self.lon = data_t1.lon
        self.xRes = data_t1.xRes
        self.yRes = data_t1.yRes
        self.xSize = data_t1.xSize
        self.ySize = data_t1.ySize
        
        self.hem_NS = data_t1.hem_NS
        self.hem_EW = data_t1.hem_EW
        
        # Get GDAL geotransform and projection
        self.geo_t = data_t1.geo_t
        self.proj = data_t1.proj
                
        # Change definitions
        self.change_intensity_threshold = change_intensity_threshold
        self.change_area_threshold = change_area_threshold
        
        self.year_t1 = data_t1.year
        self.year_t2 = data_t2.year
        
        self.output_dir = output_dir
        self.output_pattern = self.__getOutputPattern()
        
        # Nodata currently hardwired to 99
        self.nodata = data_t1.nodata
        self.nodata_byte = data_t1.nodata_byte
        
        # Calculate combined mask
        self.mask = self.__combineMasks()
        
    def __testTiles(self):
        '''
        Test that input tiles are from reasonable lats/lons/years
        '''
        
        assert self.data_t1.lat == self.data_t2.lat and self.data_t1.lon == self.data_t2.lon, "Input tiles must be from the same location."
        assert self.data_t1.year <= self.data_t2.year, "Input data_t2 must be from a later year than data_t1."
        assert self.data_t1.year != self.data_t2.year, "Input data_t1 cannot be from the same year as data_t2."
        assert self.data_t1.lee_filter == self.data_t2.lee_filter, "Only one of the input tiles has been filtered. Both tiles should have the same pre-processing parameters."
        assert self.data_t1.proj == self.data_t2.proj, "Input tiles do not have the same projection."
        assert self.data_t1.xSize == self.data_t2.xSize and self.data_t1.ySize == self.data_t2.ySize, "Input tiles do not have the same resolution."
        assert self.data_t1.geo_t == self.data_t2.geo_t, "Input tiles do not have the same geo_transform."
        assert self.data_t1.forest_threshold == self.data_t2.forest_threshold, "'forest_threshold' must be identical for both input tiles."
        assert self.data_t1.area_threshold == self.data_t2.area_threshold, "'area_threshold' must be identical for both input tiles."
        assert self.data_t1.downsample_factor == self.data_t2.downsample_factor, "'downsample_facor' must be identical for both input tiles."


    def __combineMasks(self):
        '''
        Add together masks from data_t1 and data_t2
        '''
        
        return np.logical_or(self.data_t1.mask, self.data_t2.mask)
        
    def __getOutputPattern(self):
        """
        Generates a filename pattern for data output.
        """
        
        return '%s_%s_%s_%s%s.tif'%('%s', str(self.year_t1), str(self.year_t2), self.hem_NS + str(abs(self.lat)).zfill(2), self.hem_EW + str(abs(self.lon)).zfill(3))

    def __getNodata(self, dtype = 6):
        """
        Return nodata values
        """
        
        from osgeo import gdal
        
        # Generate a nodata value
        if dtype == gdal.GDT_Byte:
            nodata = 99
        else:
            nodata = 999999
        
        return nodata    
            
    def getAGBChange(self, output = False):
        '''
        '''
        
        # Only run processing if not already done
        if not hasattr(self, 'AGBChange'):
            
            AGB_change = self.data_t2.getAGB() - self.data_t1.getAGB()
            
            self.AGB_change = AGB_change
            
        if output: self.__outputGeoTiff(AGB_change, 'AGBChange')
        
        return self.AGB_change
        
        
    def getChangeType(self, output = False):
        '''
        Returns pixels that meet change detection thresholds for country.
        
        Args:
            forest_threshold = threshold above which a pixel is forest
            intensity_threshold = threshold of proportional change with which to accept a change as real
            area_threshold
        '''
        
        from osgeo import gdal
        
        # Only run processing if not already done
        if not hasattr(self, 'ChangeType'):

            # Pixels that move from forest to nonforest (F_NF) or vice versa (NF_F)
            F_NF = np.logical_and(self.data_t1.getWoodyCover(), self.data_t2.getWoodyCover() == False)
            NF_F = np.logical_and(self.data_t1.getWoodyCover() == False, self.data_t2.getWoodyCover())
            
            # Pixels that remain forest (F_F) or nonforest (NF_NF)
            F_F = np.logical_and(self.data_t1.getWoodyCover(), self.data_t2.getWoodyCover())
            NF_NF = np.logical_and(self.data_t1.getWoodyCover() == False, self.data_t2.getWoodyCover() == False)
            
            # Get pixels of change greater than intensity threshold
            CHANGE = np.logical_or((self.getAGBChange() / self.data_t1.getAGB()) >= self.change_intensity_threshold, (self.getAGBChange() / self.data_t1.getAGB()) < (- self.change_intensity_threshold))
            NOCHANGE = CHANGE == False
            
            # Trajectory (changes can be positive or negative)
            DECREASE = self.data_t2.getAGB() < self.data_t1.getAGB()
            INCREASE = self.data_t2.getAGB() >= self.data_t1.getAGB()
            
            # Get a minimum pixel extent. Loss/Gain events much have a spatial extent greater than min_pixels, and not occur in nonforest.
            if self.change_area_threshold > 0:
                
                min_pixels = int(round(self.change_area_threshold / (self.yRes * self.xRes * 0.0001)))
                
                # Get areas of change that meet minimum area requirement
                CHANGE_INCREASE, _ = biota.indices.getContiguousAreas(CHANGE & INCREASE & (NF_F | F_F), True, min_pixels = min_pixels)
                CHANGE_DECREASE, _ = biota.indices.getContiguousAreas(CHANGE & DECREASE & (F_NF | F_F), True, min_pixels = min_pixels)
                CHANGE = np.logical_or(CHANGE_INCREASE, CHANGE_DECREASE)
                NOCHANGE = CHANGE == False
            
            change_type = {}
            
            # These are all the possible definitions. Note: they *must* add up to one.
            change_type['deforestation'] = F_NF & CHANGE
            change_type['degradation'] = F_F & CHANGE & DECREASE
            change_type['minorloss'] = (F_F | F_NF) & NOCHANGE & DECREASE
            
            change_type['aforestation'] = NF_F & CHANGE
            change_type['growth'] = F_F & CHANGE & INCREASE
            change_type['minorgain'] = (F_F | NF_F) & NOCHANGE & INCREASE
                    
            change_type['nonforest'] = NF_NF
            
            self.ChangeType = change_type
            
        if output:
            
            # Image with coded change values
            output_im = np.zeros((self.ySize, self.xSize), dtype = np.int8) + 99
            
            output_im[self.ChangeType['nonforest'].data] = 0
            output_im[self.ChangeType['deforestation'].data] = 1
            output_im[self.ChangeType['degradation'].data] = 2
            output_im[self.ChangeType['minorloss'].data] = 3
            output_im[self.ChangeType['minorgain'].data] = 4
            output_im[self.ChangeType['growth'].data] = 5
            output_im[self.ChangeType['aforestation'].data] = 6
            
            # Keep things tidy
            output_im[self.mask] = self.nodata_byte
            
            self.__outputGeoTiff(output_im, 'ChangeType', dtype = gdal.GDT_Byte)
                
        return self.ChangeType
        
    
    def __sumChange(self, change_type, scale = 1):
        '''
        Function for scaling then summing and (optionally) scaling change statistics.
        '''
        
        return {k: np.sum(v * scale) for k, v in change_type.items()}
    
    def getAreaSum(self, proportion = False, output = False):
        '''
        Extract a change area in hectares. 
        '''
        
        # Classify change type
        change_type = self.getChangeType()
                        
        # Get area change in units of ha/pixel
        change_hectares = (self.mask == False) * (self.xRes * self.yRes * 0.0001)
        
        totals = self.__sumChange(change_type, scale = change_hectares)
        
        if proportion: 
            for change in totals:
                totals[change] = totals[change] / change_hectares.sum()
        
        return totals
    
    
    def getAGBSum(self, proportion = False, output = False):
        '''
        Extract a change magnitude in tonnes of carbon.
        TODO: Proportional measures need work.
        TODO: Add similar summary stats to LoadTile()?
        '''
        
        # Classify change type
        change_type = self.getChangeType()
        
        # Get AGB change in units of tC/pixel
        change_AGB = self.getAGBChange() * (self.xRes * self.yRes * 0.0001)
        
        totals = self.__sumChange(change_type, scale = change_AGB)
                
        if proportion: 
            for change_type in totals:    
                totals[change_type] = totals[change_type] / (self.data_t1.getAGB() * self.xRes * self.yRes * 0.0001).sum()
        
        return totals
    
    
    def __outputGeoTiff(self, data, output_name, dtype = 6):
        """
        Output a GeoTiff file.
        """
                
        # Generate a standardised filename
        filename = self.output_pattern%output_name
        
        nodata = self.__getNodata(dtype = dtype)
        
        # Write to disk
        biota.IO.outputGeoTiff(data, filename, self.geo_t, self.proj, output_dir = self.output_dir, dtype = dtype, nodata = nodata)
        
        