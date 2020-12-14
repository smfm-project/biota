#!/usr/bin/env python

import argparse
import csv
import datetime as dt
import itertools
import math
import numpy as np
import os
from osgeo import gdal
from scipy import ndimage
import scipy.stats as stats
import scipy.ndimage.morphology
import skimage.measure

import matplotlib.pyplot as plt
import pdb

import biota.filter
import biota.indices
import biota.IO
import biota.mask
import biota.SM

"""
These are two classes for loading individual tiles and compliling them into a change object.
"""



class LoadTile(object):
    """
    Class to load an ALOS mosaic tile, and extract properties related to properties of forest in the tile.

    Mosaic tiles require the following attributes to be loaded:
        data_dir: Directory containing data from the ALOS mosaic. Files should retain the original file structure provided by JAXA.
        lat: Latitude of the upper-left hand corner of the mosaic tile.
        lon: Longitude of the upper-left hand corner of the mosaic tile.
        year: Year of the ALOS mosaic tile
        
    Further optional attributes are:
        forest_threshold: The theshold of aboveground biomass that separates forest from non-forest in units of tonnes carbon per hectare (tC/ha). Defaults to 10 tC/ha (approximately 10% canopy cover in soutehrn Africa).
        area_threshold: Contiguous area required to meet the definiton of forest, in units of hectares. Defaults to 0 ha.
        downsample_factor:
        lee_filter: Apply a radar speckle filter to the ALOS image. Defaults to True.
        window_size: Size of lee_filter window. Must be an odd integer. Defaults to 5.
        contiguity: When applying an area_threshold, a forest area could be considered continuous when directly adjacent ("rook's move") or diagonally adjacent to another forest pixel ("queen's move"). To switch, set this parameter to either 'rook' or 'queen'. Defaults to 'queen'.
        output_dir: Directory to save output GeoTiff images. Defaults to present working directory.

    For example, to load an ALOS tile:
        tile_2015 = biota.LoadTile('/path/to/data_dir/', -15, 30, 2015)
        
    To load an ALOS tile, but applying a 20 tC/ha forest threshold and a 1 hectare forest definition:
        tile_2015 = biota.LoadTile('/path/to/data_dir/', -15, 30, 2015, forest_thresold = 20., area_threshold = 1.)
    """


    def __init__(self, data_dir, lat, lon, year, parameter_file = None, forest_threshold = 10., area_threshold = 0., downsample_factor = 1, lee_filter = True, window_size = 5, contiguity = 'queen', sm_dir = os.getcwd(), sm_interpolation = 'average', output_dir = os.getcwd()):
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
        assert type(window_size) == int, "Option window_size must be an integer."
        assert window_size % 2 == 1, "Option window_size must be an odd integer."
        assert contiguity in ['rook', 'queen'], "Contiguity constraint must be 'rook' or 'queen'."
        assert os.path.isdir(os.path.expanduser(data_dir)), "Specified data directory (%s) does not exist"%str(data_dir)
        assert os.path.isdir(os.path.expanduser(sm_dir)), "Specified soil moisture directory (%s) does not exist"%str(sm_dir)
        assert sm_interpolation in ['nearest', 'average', 'cubic'], "Soil moisture interpolation type must be one of 'nearest', 'average' or 'cubic'."
        assert os.path.isdir(os.path.expanduser(output_dir)), "Specified output directory (%s) does not exist"%str(output_dir)
        assert type(forest_threshold) == float or type(forest_threshold) == int, "Forest threshold must be numeric."
        assert type(area_threshold) == float or type(area_threshold) == int, "Area threshold must be numeric."

        self.lat = lat
        self.lon = lon
        self.year = year
        self.downsample_factor = downsample_factor
        self.lee_filter = lee_filter
        self.window_size = window_size
        self.contiguity = contiguity
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
        self.HV_path = self.__getHVPath()
        self.HH_path = self.__getHHPath()
        self.mask_path = self.__getMaskPath()
        self.date_path = self.__getDatePath()

        # Load AGB parameter file (default, or specified)
        self.AGB_parameters_path = self.__getAGBParametersPath() if parameter_file == None else os.path.expanduser(parameter_file)
        assert os.path.exists(self.AGB_parameters_path), "Specified parameter file (%s) not found."%str(self.AGB_parameters_path)
        
        # Set up soil moisture data location
        self.SM_dir = os.path.expanduser(sm_dir.rstrip('/'))
        self.SM_interpolation = sm_interpolation

        # Set up locations for file output
        self.output_pattern = self.__getOutputPattern()
        self.output_dir = os.path.expanduser(output_dir.rstrip('/'))

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

        # Load in AGB parameters
        self.loadAGBParameters(self.AGB_parameters_path)

        # Update image metadata where a degree of resampling is included
        if self.downsample_factor != 1:
            self.__rebin()

        # Load mask
        self.mask = self.__getMask()


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

        # Directories and files have standardised pattern
        name_pattern = '%s%s_%s_%s'

        # Directory name patterns name patterns are different for ALOS-1/ALOS-2
        if self.satellite == 'ALOS-2':
            name_pattern += '_F02DAR'

        # First, check if the 1x1 tile is present at data_dir
        lat_dir = self.hem_NS + str(abs(self.lat)).zfill(2)
        lon_dir = self.hem_EW + str(abs(self.lon)).zfill(3)

        # Generate directory name
        directory = self.data_dir + '/' + name_pattern%(lat_dir, lon_dir, str(self.year)[-2:], 'MOS') + '/'

        # If no 1x1 tile, check that 5x5 tile exists
        if not os.path.isdir(directory):

            # Calculate the hemisphere and lat/lon of the 5x5 tile
            lat_large = self.lat + (5 - self.lat) % 5
            lon_large = self.lon - (self.lon % 5)
            hem_NS_large = 'S' if lat_large < 0 else 'N'
            hem_EW_large = 'W' if lon_large < 0 else 'E'

            # Get lat/lon for directory (lat/lon of upper-left hand corner of 5x5 tiles)
            lat_dir = hem_NS_large + str(abs(lat_large)).zfill(2)
            lon_dir = hem_EW_large + str(abs(lon_large)).zfill(3)

            # Generate directory name
            directory = self.data_dir + '/' + name_pattern%(lat_dir, lon_dir, str(self.year)[-2:], 'MOS') + '/'

            if not os.path.isdir(directory):
                raise IOError('No tile for lat: %s, lon: %s exists in the specified data directory.'%(str(self.lat), str(self.lon)))

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

        HHpath = self.__getDirectory() + self.__getFilename('sl_HH')

        # Stop if the ALOS data don't exist
        if not os.path.isfile(HHpath):
            raise IOError('No data found for HH polarisation for lat: %s, lon: %s.'%(str(self.lat), str(self.lon)))

        return HHpath

    def __getHVPath(self):
        """
        Determines the filepath to HV data.
        """

        HVpath = self.__getDirectory() + self.__getFilename('sl_HV')

        # Stop if the ALOS data don't exist
        if not os.path.isfile(HVpath):
            raise IOError('No data found for HV polarisation for lat: %s, lon: %s.'%(str(self.lat), str(self.lon)))

        return HVpath


    def __getMaskPath(self):
        """
        Determines the filepath to mask data.
        """

        mask = self.__getDirectory() + self.__getFilename('mask')

        # Stop if the ALOS data don't exist
        if not os.path.isfile(mask):
            raise IOError('No data found for mask for lat: %s, lon: %s.')

        return mask

    def __getDatePath(self):
        """
        Determines the filepath to DOY data.
        """

        date = self.__getDirectory() + self.__getFilename('date')

        # Stop if the ALOS data don't exist
        if not os.path.isfile(date):
            raise IOError('No data found for date for lat: %s, lon: %s.')

        return date

    def __getAGBParametersPath(self):
        """
        Get default AGB model parameters
        """
        
        return '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-1]) + '/cfg/McNicol2018.csv'
        

    def __getOutputPattern(self):
        """
        Generates a filename pattern for data output.
        """

        return '%s_%s_%s%s.tif'%('%s', str(self.year), self.hem_NS + str(abs(self.lat)).zfill(2), self.hem_EW + str(abs(self.lon)).zfill(3))

    def __getNodata(self, dtype = 6):
        """
        Return nodata values
        """

        # Generate a nodata value
        if dtype == gdal.GDT_Byte:
            nodata = 255
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

    def __getMask(self, masked_px_count = False, output = False, show = False):
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

        if output: self.__outputGeoTiff(mask, 'Mask', dtype = gdal.GDT_Byte)

        if show: self.__showArray(mask, title = 'Mask', cmap = 'coolwarm')

        return mask

    def updateMask(self, filename, buffer_size = 0., classes = [], output = False, show = False):
        """
        Function to add further pixels to mask based on a numpy array, shapefile or GeoTiff file, with optional buffers.
        """

        # Load mask
        mask = biota.mask.updateMask(self, filename, buffer_size = buffer_size, classes = classes)

        # Add the new raster masks to the existing mask
        self.mask = np.logical_or(self.mask, mask)

        if output: self.__outputGeoTiff(self.mask * 1, 'Mask', dtype = gdal.GDT_Byte)

        if show: self.__showArray(self.mask * 1, title = 'Mask', cmap = 'coolwarm')

    def resetMask(self):
        """
        Function to reset a mask to the default.
        """

        self.mask = self.__getMask()

    def loadAGBParameters(self, csv_file):
        """
        Load AGB parameters for ALOS-1 and ALOS-2 from appropriately formatted .csv file.
        """
        
        assert csv_file.endswith('.csv'), "Specified parameter file (%s) must be a .csv."%csv_file
        assert os.path.exists(os.path.expanduser(csv_file)), "Specified parameter file (%s) doesn't exist."%csv_file
        
        with open(csv_file) as csvfile:
            reader = csv.reader(csvfile, delimiter = ',')
            header = next(reader)
            self.A1_gradient, self.A1_intercept = [float(i) for i in next(reader)[1:]]
            self.A2_gradient, self.A2_intercept = [float(i) for i in next(reader)[1:]]


    def getDN(self, polarisation = 'HV', output = False, show = False):
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

            # Load the sum of contributing pixels
            block_sum = skimage.measure.block_reduce(np.ones_like(DN), (self.downsample_factor, self.downsample_factor), np.sum)

            # And the sum of masked pixels
            mask_sum = self.__getMask(masked_px_count = True)

            # Divide the sum of DNs by the sum of unmasked pixels to get the mean DN value
            DN = np.zeros_like(DN_sum)

            DN[self.mask == False] = (DN_sum.astype(np.float)[self.mask == False] / (block_sum - mask_sum)[self.mask == False]).astype(np.int)

        if output: self.__outputGeoTiff(DN, 'DN', dtype = gdal.GDT_Int32)

        if show: self.__showArray(DN, title = 'DN', cbartitle = 'Digital Number', cmap = 'Spectral_r')

        return np.ma.array(DN, mask = self.mask)

    def __getDay(self):
        '''
        '''

        day_after_launch = biota.IO.loadArray(self.date_path)

        # If required, downsample the dates array
        if self.downsample_factor != 1:

            day_after_launch = skimage.measure.block_reduce(day_after_launch, (self.downsample_factor, self.downsample_factor), np.max)

        return day_after_launch

    def getDate(self, output = False, show = False):
        """
        Loads date values into a numpy array.
        """

        # Don't rerun processing if already present in memory
        if not hasattr(self, 'date'):

            day_after_launch = self.__getDay()

            # Get list of unique dates in ALOS tile
            unique_days = np.unique(day_after_launch[day_after_launch != 0])

            # Days are counted from the launch date
            if self.satellite == 'ALOS-2':
                launch_date = dt.date(2014, 5, 24)
            else:
                launch_date = dt.date(2006, 1, 24)

            # Determine the Day of Year associated with each
            unique_dates = [launch_date  + dt.timedelta(days=int(d)) for d in unique_days]

            # Build a Day Of Year array
            dates = np.zeros_like(day_after_launch, dtype='datetime64[D]')
            dates_int = np.zeros_like(day_after_launch, dtype = np.int32)
            for day, date in zip(unique_days, unique_dates):
                dates[day_after_launch == day] = np.datetime64(date,'D')
                dates_int[day_after_launch == day] = np.int(np.datetime64(date,'D').astype(dt.date).strftime('%Y%m%d'))

            # Save output to class
            self.date = dates
            self.date_int = dates_int

        if output: self.__outputGeoTiff(self.dates_int, 'Date', dtype = gdal.GDT_Int32)

        if show: print("Sorry, matplotlib doesn't support display of numpy datetime objects right now.")

        return self.date

    def getDOY(self, output = False, show = False):
        """
        Loads day of year values into a numpy array.
        """

        # Don't rerun processing if already present in memory
        if not hasattr(self, 'DOY'):

            day_after_launch = self.__getDay()
            dates = self.getDate()

            # Determine the Day of Year associated with each
            unique_dates = np.unique(dates)
            unique_doys = [(d.astype(dt.date) - dt.date(d.astype(dt.date).year,1,1)).days + 1 for d in unique_dates]
            unique_days = np.unique(day_after_launch)

            # Build a Day Of Year array
            DOY = np.zeros_like(dates, dtype = np.int16)
            for day, doy in zip(unique_days, unique_doys):
                DOY[day_after_launch == day] = doy

            # Save output to class
            self.DOY = DOY

        if output: self.__outputGeoTiff(self.DOY, 'DOY', dtype = gdal.GDT_Int32)

        if show: self.__showArray(self.DOY, title = 'DOY', cbartitle = 'Day', vmin = 0, vmax = 366, cmap = 'Spectral')

        return self.DOY

    def getYearArray(self, output = False, show = False):
        """
        Loads year of overpass into a numpy array
        """

        # Don't rerun processing if already present in memory
        if not hasattr(self, 'YearArray'):

            dates = self.getDate()

            YearArray = dates.astype('datetime64[Y]').astype(int) + 1970

            # Save output to class
            self.YearArray = YearArray

        if output: self.__outputGeoTiff(self.YearArray, 'YearArray', dtype = gdal.GDT_Int32)

        if show: self.__showArray(self.YearArray, title = 'YearArray', cbartitle = 'Year', vmin = self.year-1, vmax = self.year+1, cmap = 'Spectral')

        return self.YearArray

    def getSM(self, output = False, show = False, search_days = 7):
        """
        Loads a soil moisture map using the ESA CCI soil moisture product.
        """

        # Don't rerun processing if already present in memory
        if not hasattr(self, 'SM'):

            SM = biota.SM.getSM(self, search_days = search_days, interpolation = self.SM_interpolation)

            # Keep masked values tidy
            SM.data[self.mask] = self.nodata

            # Save output to class
            self.SM = SM

        if output: self.__outputGeoTiff(self.SM, 'SM')

        if show: self.__showArray(self.SM, title = 'Soil Moisture', cbartitle = 'm^2/m^2', vmin = 0., vmax = 0.3, cmap = 'Blues')

        return self.SM

    def getGamma0(self, polarisation = 'HV', units = 'natural', output = False, show = False):
        """
        Calibrates data to gamma0 (baskscatter) in decibels or natural units.
        """

        assert units == 'natural' or units == 'decibels', "Units must be 'natural' or 'decibels'. You input %s."%units
        assert polarisation == 'HH' or polarisation == 'HV', "Polarisation must be 'HH' or 'HV'. You input %s."%polarisation

        # Calibrate DN to units of dB
        gamma0 = 10 * np.ma.log10(self.getDN(polarisation = polarisation).astype(np.float) ** 2) - 83. # units = decibels
        # Apply filter based on dB values
        if self.lee_filter:
            gamma0 = biota.filter.enhanced_lee_filter(gamma0, n_looks = self.nLooks, window_size = self.window_size)


        # Convert to natural units where specified
        if units == 'natural': gamma0 = 10 ** (gamma0 / 10.)

        # Keep masked values tidy
        gamma0.data[self.mask] = self.nodata
        gamma0.mask = self.mask

        if output: self.__outputGeoTiff(gamma0, 'Gamma0')

        if show:
            # Different display settings depending on options
            if polarisation == 'HH' and units == 'natural': vmin, vmax = 0, 0.15
            if polarisation == 'HV' and units == 'natural': vmin, vmax = 0, 0.06
            if polarisation == 'HH' and units == 'decibels': vmin, vmax = -15, -5
            if polarisation == 'HV' and units == 'decibels': vmin, vmax = -20, -10

            self.__showArray(gamma0, title = 'Gamma0 %s'%polarisation, cbartitle = units, vmin = vmin, vmax = vmax, cmap = 'Greys_r')

        return gamma0

    def getAGB(self, output = False, show = False):
        """
        Calibrates data to aboveground biomass (AGB).
        Placeholder equation to calibrate backscatter (gamma0) to AGB (tC/ha).
        """

        # Don't rerun processing if already present in memory
        if not hasattr(self, 'AGB'):

            # ALOS-1
            if self.satellite == 'ALOS-1':
                AGB = self.A1_gradient * self.getGamma0(units = 'natural', polarisation = 'HV') + self.A1_intercept

            # ALOS-2 (to calculate)
            elif self.satellite == 'ALOS-2':
                AGB = self.A2_gradient * self.getGamma0(units = 'natural', polarisation = 'HV') + self.A2_intercept

            else:
                raise ValueError("Unknown satellite named '%s'. self.satellite must be 'ALOS-1' or 'ALOS-2'."%self.satellite)

            # Save output to class
            self.AGB = AGB

        # Keep masked values tidy
        #self.AGB.data[self.mask] = self.nodata
        self.AGB.mask = self.mask

        if output: self.__outputGeoTiff(self.AGB, 'AGB')

        if show: self.__showArray(self.AGB, title = 'AGB', cbartitle = 'tC/ha', vmin = 0, vmax = 40, cmap = 'YlGn')

        return self.AGB

    def getWoodyCover(self, output = False, show = False):
        """
        Get woody cover, based on a threshold of AGB.
        min_area in ha
        """

        # Don't rerun processing if already present in memory
        if not hasattr(self, 'WoodyCover'):

            WoodyCover = self.getAGB() >= float(self.forest_threshold)

            if self.area_threshold > 0:

                # Calculate number of pixels in min_area (assuming input is given in hecatres)
                min_pixels = int(round(self.area_threshold / (self.yRes * self.xRes * 0.0001)))

                # Remove pixels that aren't part of a forest block of size at least min_pixels
                contiguous_area, _ = biota.indices.getContiguousAreas(WoodyCover, True, min_pixels = min_pixels, contiguity = self.contiguity)

                WoodyCover.data[contiguous_area == False] = False

            # Save output to class
            self.WoodyCover = WoodyCover

        # Keep masked values tidy
        #self.WoodyCover.data[self.mask] = False
        self.WoodyCover.mask = self.mask

        if output:

            # Convert data to integer type for output
            WoodyCover_out = self.WoodyCover.astype(np.uint8)

            self.__outputGeoTiff(WoodyCover_out, 'WoodyCover', dtype = gdal.GDT_Byte)

        if show: self.__showArray(self.WoodyCover, title = 'Woody Cover', cbartitle = '', vmin = 0, vmax = 1, cmap = 'summer_r')

        return self.WoodyCover

    def getForestPatches(self, output = False, show = False):
        """
        Get numbered forest patches, based on woody cover threshold.
        """

        # Don't rerun processing if already present in memory
        if not hasattr(self, 'ForestPatches'):

            WoodyCover = self.getWoodyCover()

            # Calculate number of pixels in min_area
            min_pixels = int(round(self.area_threshold / (self.yRes * self.xRes * 0.0001)))

            # Get areas that meet that threshold
            _, ForestPatches = biota.indices.getContiguousAreas(WoodyCover, True, min_pixels = min_pixels, contiguity = self.contiguity)

            # Save output to class
            self.ForestPatches = ForestPatches

        # Tidy up masked pixels
        #ForestPatches.data[np.ma.getmaskarray(self.ForestPatches)] = self.nodata
        self.ForestPatches.mask = self.mask

        if output: self.__outputGeoTiff(self.ForestPatches * 1, 'ForestPatches', dtype = gdal.GDT_Int32)

        if show: self.__showArray(self.ForestPatches, title = 'Forest patches', cbartitle = 'Patch ID', cmap = 'Spectral')

        return self.ForestPatches

    def __outputGeoTiff(self, data, output_name, dtype = 6):
        """
        Output a GeoTiff file.
        """

        # Generate a standardised filename
        filename = self.output_pattern%output_name

        # Generate a nodata value appropriate for datatype
        nodata = self.__getNodata(dtype = dtype)

        # Write to disk
        biota.IO.outputGeoTiff(data, filename, self.geo_t, self.proj, output_dir = self.output_dir, dtype = dtype, nodata = nodata)

    def __showArray(self, data, title = '', cbartitle = '', vmin = None, vmax = None, cmap = None):
        """
        Display data from a tile.
        """

        biota.IO.showFigure(data, self.lat, self.lon, title = title, cbartitle = cbartitle, vmin = vmin, vmax = vmax, cmap = cmap)


class LoadChange(object):
    """
    Input is two mosaic tiles from LoadTile, will output maps and change statistics.
    """

    def __init__(self, tile_t1, tile_t2, change_intensity_threshold = 0.2, change_magnitude_threshold = 0., change_area_threshold = 0, deforestation_threshold = None, contiguity = 'queen', combine_areas = True, output_dir = os.getcwd()):
        '''
        Initialise
        '''

        assert contiguity in ['rook', 'queen'], "Contiguity constraint must be 'rook' or 'queen'."

        self.tile_t1 = tile_t1
        self.tile_t2 = tile_t2

        # Ensure that tiles are compatable
        self.__testTiles()

        # Get basic properties
        self.lat = tile_t1.lat
        self.lon = tile_t1.lon
        self.xRes = tile_t1.xRes
        self.yRes = tile_t1.yRes
        self.xSize = tile_t1.xSize
        self.ySize = tile_t1.ySize

        self.hem_NS = tile_t1.hem_NS
        self.hem_EW = tile_t1.hem_EW
        
        # Get GDAL geotransform and projection
        self.geo_t = tile_t1.geo_t
        self.proj = tile_t1.proj

        # Change definitions
        self.change_intensity_threshold = change_intensity_threshold
        self.change_magnitude_threshold = change_magnitude_threshold
        self.change_area_threshold = change_area_threshold
        self.deforestation_threshold = deforestation_threshold
        self.combine_areas = combine_areas
        
        # Set deforestation_treshold equal to forest_threshold if not in use
        if self.deforestation_threshold is None: self.deforestation_threshold = self.tile_t1.forest_threshold
        
        self.contiguity = contiguity

        self.year_t1 = tile_t1.year
        self.year_t2 = tile_t2.year

        self.output_dir = output_dir
        self.output_pattern = self.__getOutputPattern()

        # Nodata currently hardwired to 99
        self.nodata = tile_t1.nodata
        self.nodata_byte = tile_t1.nodata_byte

        # Calculate combined mask
        self.mask = self.__combineMasks()

    def __testTiles(self):
        '''
        Test that input tiles are from reasonable lats/lons/years
        '''

        assert self.tile_t1.lat == self.tile_t2.lat and self.tile_t1.lon == self.tile_t2.lon, "Input tiles must be from the same location."
        assert self.tile_t1.year <= self.tile_t2.year, "Input tile_t2 must be from a later year than tile_t1."
        assert self.tile_t1.year != self.tile_t2.year, "Input tile_t1 cannot be from the same year as tile_t2."
        assert self.tile_t1.lee_filter == self.tile_t2.lee_filter, "Only one of the input tiles has been filtered. Both tiles should have the same pre-processing parameters."
        assert self.tile_t1.proj == self.tile_t2.proj, "Input tiles do not have the same projection."
        assert self.tile_t1.xSize == self.tile_t2.xSize and self.tile_t1.ySize == self.tile_t2.ySize, "Input tiles do not have the same resolution."
        assert self.tile_t1.geo_t == self.tile_t2.geo_t, "Input tiles do not have the same geo_transform."
        assert self.tile_t1.forest_threshold == self.tile_t2.forest_threshold, "'forest_threshold' must be identical for both input tiles."
        assert self.tile_t1.area_threshold == self.tile_t2.area_threshold, "'area_threshold' must be identical for both input tiles."
        assert self.tile_t1.downsample_factor == self.tile_t2.downsample_factor, "'downsample_facor' must be identical for both input tiles."


    def __combineMasks(self):
        '''
        Add together masks from tile_t1 and tile_t2
        '''

        return np.logical_or(self.tile_t1.mask, self.tile_t2.mask)

    def __getOutputPattern(self):
        """
        Generates a filename pattern for data output.
        """

        return '%s_%s_%s_%s%s.tif'%('%s', str(self.year_t1), str(self.year_t2), self.hem_NS + str(abs(self.lat)).zfill(2), self.hem_EW + str(abs(self.lon)).zfill(3))

    def __getNodata(self, dtype = 6):
        """
        Return nodata values
        """

        # Generate a nodata value
        if dtype == gdal.GDT_Byte:
            nodata = 255
        else:
            nodata = 999999

        return nodata

    def updateMask(self, filename, buffer_size = 0., classes = [], output = False, show = False):
        """
        Function to add further pixels to mask based on a numpy array, shapefile or GeoTiff file, with optional buffers.
        """

        # Load mask
        mask = biota.mask.updateMask(self, filename, buffer_size = buffer_size, classes = classes)

        # Add the new raster masks to the existing mask
        self.mask = np.logical_or(self.mask, mask)

        if output: self.__outputGeoTiff(self.mask * 1, 'Mask', dtype = gdal.GDT_Byte)

        if show: self.__showArray(self.mask * 1, title = 'Mask', cmap = 'coolwarm')

    def resetMask(self):
        """
        Function to reset a mask to the default.
        """

        self.mask = self.__combineMasks()


    def getSMChange(self, output = False, show = False):
        """
        """

        # Only run processing if not already done
        if not hasattr(self, 'SMChange'):

            self.SMChange = self.tile_t2.getSM() - self.tile_t1.getSM()

        if output: self.__outputGeoTiff(self.SMChange, 'SMChange')

        if show: self.__showArray(self.SMChange, title = 'SM Change', cbartitle = 'm^2/m^2', vmin = -0.15, vmax = 0.15, cmap = 'RdBu')

        return self.SMChange

    def getGamma0Change(self, polarisation = 'HV', units = 'natural', output = False, show = False):
        '''
        '''

        # Always run processing as pol/units may change
        self.Gamma0_change = self.tile_t2.getGamma0(polarisation = polarisation, units = units) - self.tile_t1.getGamma0(polarisation = polarisation, units = units)

        # Add combined mask
        self.Gamma0_change.mask = self.mask

        if output: self.__outputGeoTiff(self.Gamma0_change, 'Gamma0Change')

        if show:
            if polarisation == 'HH' and units == 'natural': vmin, vmax = -0.05, 0.05
            if polarisation == 'HV' and units == 'natural': vmin, vmax = -0.025, 0.025
            if polarisation == 'HH' and units == 'decibels': vmin, vmax = -5, 5
            if polarisation == 'HV' and units == 'decibels': vmin, vmax = -2.5, 2.5

            self.__showArray(self.Gamma0_change, title = 'Gamma0 Change', cbartitle = units, vmin = vmin, vmax = vmax, cmap = 'RdGy')

        return self.Gamma0_change

    def getAGBChange(self, output = False, show = False):
        '''
        '''

        # Only run processing if not already done
        if not hasattr(self, 'AGBChange'):

            self.AGB_change = self.tile_t2.getAGB() - self.tile_t1.getAGB()

        self.AGB_change.mask = self.mask

        if output: self.__outputGeoTiff(self.AGB_change, 'AGBChange')

        if show: self.__showArray(self.AGB_change, title = 'AGB Change', cbartitle = 'tC/ha', vmin = -10, vmax = 10, cmap = 'RdBu')

        return self.AGB_change


    def getChangeType(self, output = False, show = False):
        '''
        Returns pixels that meet change detection thresholds for country.

        Args:
            forest_threshold = threshold above which a pixel is forest
            intensity_threshold = threshold of proportional change with which to accept a change as real
            area_threshold
        '''

        # Only run processing if not already done
        if not hasattr(self, 'ChangeType'):

            # Pixels that move from forest to nonforest (F_NF) or vice versa (NF_F)
            F_NF = np.logical_and(self.tile_t1.getWoodyCover(), self.tile_t2.getWoodyCover() == False)
            NF_F = np.logical_and(self.tile_t1.getWoodyCover() == False, self.tile_t2.getWoodyCover())

            # Pixels that remain forest (F_F) or nonforest (NF_NF)
            F_F = np.logical_and(self.tile_t1.getWoodyCover(), self.tile_t2.getWoodyCover())
            NF_NF = np.logical_and(self.tile_t1.getWoodyCover() == False, self.tile_t2.getWoodyCover() == False)

            # Get pixels of change greater than intensity threshold
            CHANGE_INTENSITY = np.logical_or((self.getAGBChange() / self.tile_t1.getAGB()) >= self.change_intensity_threshold, (self.getAGBChange() / self.tile_t1.getAGB()) < (- self.change_intensity_threshold))

            # Get pixels of change greater than magnitude threshold
            CHANGE_MAGNITUDE = np.logical_or(self.getAGBChange() >= self.change_magnitude_threshold, self.getAGBChange() < (- self.change_magnitude_threshold))

            CHANGE = np.logical_and(CHANGE_INTENSITY, CHANGE_MAGNITUDE)
            NOCHANGE = CHANGE == False

            # Trajectory (changes can be positive or negative)
            DECREASE = self.tile_t2.getAGB() < self.tile_t1.getAGB()
            INCREASE = self.tile_t2.getAGB() >= self.tile_t1.getAGB()

            # Get a minimum pixel extent. Loss/Gain events much have a spatial extent greater than min_pixels, and not occur in nonforest.
            if self.change_area_threshold > 0:

                min_pixels = int(round(self.change_area_threshold / (self.yRes * self.xRes * 0.0001)))

                # Get areas of change that meet minimum area requirement (use 'change' pixels for area measurement)
                if self.combine_areas:
                    CHANGE_INCREASE, _ = biota.indices.getContiguousAreas(CHANGE & INCREASE & (NF_F | F_F), True, min_pixels = min_pixels, contiguity = self.contiguity)
                    CHANGE_DECREASE, _ = biota.indices.getContiguousAreas(CHANGE & DECREASE & (F_NF | F_F), True, min_pixels = min_pixels, contiguity = self.contiguity)
                    CHANGE = np.logical_or(CHANGE_INCREASE, CHANGE_DECREASE)
                    NOCHANGE = CHANGE == False
                else:
                    CHANGE_DEF, _ = biota.indices.getContiguousAreas(CHANGE & DECREASE & F_NF, True, min_pixels = min_pixels, contiguity = self.contiguity)
                    CHANGE_DEG, _ = biota.indices.getContiguousAreas(CHANGE & DECREASE & F_F, True, min_pixels = min_pixels, contiguity = self.contiguity)
                    CHANGE_GRO, _ = biota.indices.getContiguousAreas(CHANGE & INCREASE & F_F, True, min_pixels = min_pixels, contiguity = self.contiguity)
                    CHANGE_AFF, _ = biota.indices.getContiguousAreas(CHANGE & INCREASE & NF_F, True, min_pixels = min_pixels, contiguity = self.contiguity)
                    CHANGE = np.logical_or(np.logical_or(CHANGE_DEF, CHANGE_DEG), np.logical_or(CHANGE_GRO, CHANGE_AFF))
                    NOCHANGE = CHANGE == False
                            
            # Deforestation can be dramatic; allow a separate treshold to specify an end-biomass for deforestation 
            DEFORESTED = self.tile_t2.getAGB() < self.deforestation_threshold
            
            change_type = {}
            
            # These are all the possible definitions. Note: they *must* add up to one.
            change_type['deforestation'] = F_NF & CHANGE & DEFORESTED
            change_type['degradation'] = (F_F & CHANGE & DECREASE) | (F_NF & CHANGE & (DEFORESTED == False))
            change_type['minorloss'] = (F_F | F_NF) & NOCHANGE & DECREASE

            change_type['afforestation'] = NF_F & CHANGE
            change_type['growth'] = F_F & CHANGE & INCREASE
            change_type['minorgain'] = (F_F | NF_F) & NOCHANGE & INCREASE

            change_type['nonforest'] = NF_NF
            
            # Also produce a code for output
            change_code = np.zeros((self.ySize, self.xSize), dtype = np.int8) + self.nodata_byte

            change_code[change_type['nonforest'].data] = 0
            change_code[change_type['deforestation'].data] = 1
            change_code[change_type['degradation'].data] = 2
            change_code[change_type['minorloss'].data] = 3
            change_code[change_type['minorgain'].data] = 4
            change_code[change_type['growth'].data] = 5
            change_code[change_type['afforestation'].data] = 6

            # Save to class
            self.ChangeType = change_type
            self.ChangeCode = change_code

        # TO FIX (include mask)
        #for ct in ['nonforest', 'deforestation', 'degradation', 'minorloss', 'minorgain', 'growth', 'afforestation']:
        #
        #self.ChangeType.mask = self.mask
        self.ChangeCode[self.mask] = 255

        if output: self.__outputGeoTiff(self.ChangeCode, 'ChangeType', dtype = gdal.GDT_Byte)

        if show:
            # Hide minor gain, minor loss and nonforest in display output
            change_code_display = np.ma.array(self.ChangeCode, mask = np.zeros_like(self.ChangeCode, dtype = np.bool))
            change_code_display.mask[np.isin(change_code_display, [0, 3, 4, 255])] = True
            self.__showArray(change_code_display, title = 'Change type', cbartitle = 'Class', vmin = 1, vmax = 6, cmap = 'Spectral')

        return self.ChangeType
    
    def getRiskMap(self, output = False, show = False, buffer_size = 50.):
        '''
        Fuction to return 'low' 'medium' and 'high' risks of deforestation, based on buffers around change locations.
        
        Authors: Samuel Bowers (University of Edinburgh) and Muri Soares (Fundo Nacional de Desenvolvimento Sustentavel)
        '''
        
        # Get change type for each pixel
        change_type = self.getChangeType()
        
        # Extract 'deforestation' and 'degradation' pixels
        deforestation = change_type['deforestation'].astype(np.int8)
        
        # Repeat change detection, but with no minimum change area threshold. Only process if deforestation_threshold exists
        if self.deforestation_threshold != self.tile_t1.forest_threshold:
            change_type_noDF = LoadChange(self.tile_t1, self.tile_t2, change_intensity_threshold = self.change_intensity_threshold, \
                    change_magnitude_threshold = self.change_magnitude_threshold, change_area_threshold = self.change_area_threshold, \
                    deforestation_threshold = None, contiguity = self.contiguity, output_dir = self.output_dir).getChangeType()
            deforestation_noDF = change_type_noDF['deforestation'].astype(np.int8)
        else:
            deforestation_noDF = change_type['deforestation'].astype(np.int8)
               
        # Set up output image
        risk_map = np.zeros_like(deforestation).astype(np.int8)
                
        # High risk of change
        risk_map[deforestation == 1] = 1
        
        # Medium risk of change (no deforestation_threshold)
        risk_map[np.logical_and(risk_map == 0, deforestation_noDF == 1)] = 2
        
        # Low risk of change
        low_dilate = scipy.ndimage.morphology.binary_dilation((np.logical_or(deforestation == 1, deforestation_noDF == 1)).astype(np.int8), iterations = int(round(buffer_size / ((self.tile_t1.xRes + self.tile_t1.yRes) / 2.), 0))) # 50 m buffer
        
        risk_map[np.logical_and(risk_map == 0, low_dilate)] = 3
                
        self.risk_map = np.ma.array(risk_map, mask = self.mask)
        
        self.risk_map.mask = self.mask
        
        if output: self.__outputGeoTiff(self.risk_map,'RiskMap')
        
        if show: self.__showArray(self.risk_map, title = 'Deforestation risk map', cbartitle = 'Low - High risk', vmin = 0, vmax = 3, cmap = 'autumn')
        
        return self.risk_map
     
    def __sumChange(self, change_type, scale = 1):
        '''
        Function for scaling then summing and (optionally) scaling change statistics.
        '''

        return {k: np.sum(v * scale) for k, v in list(change_type.items())}

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

        if output: print('TODO')

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
                totals[change_type] = totals[change_type] / (self.tile_t1.getAGB() * self.xRes * self.yRes * 0.0001).sum()

        if output: print('TODO')

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

    def __showArray(self, data, title = '', cbartitle = '', vmin = None, vmax = None, cmap = None):
        '''
        '''

        biota.IO.showFigure(data, self.lat, self.lon, title = title, cbartitle = cbartitle, vmin = vmin, vmax = vmax, cmap = cmap)
