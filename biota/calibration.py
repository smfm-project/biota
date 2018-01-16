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
import pdb


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
        self.HV_path = self.__getHVPath()
        self.mask_path = self.__getMaskPath()
        self.date_path = self.__getDatePath()

        # Stop of the ALOS tile doesn't exist
        if not os.path.isfile(self.HV_path): 
            raise IOError
        
        # Get GDAL geotransform
        self.geo_t = self.__getGeoT(self.HV_path)
        
        # Get Raster size
        self.xSize, self.ySize = self.__getSize(self.HV_path)
        
        # Load DN, mask, and day of year
        self.mask = self.getMask()
        self.DN = self.getDN()
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
        Return the directory containing ALOS data for a given lat/lon.
        """
        
        # Get lat/lon for directory (lat/lon of upper-left hand corner of 5x5 tiles)
        lat_dir = self.hem_NS + str(abs(self.lat + (5 - self.lat) % 5)).zfill(2)
        lon_dir = self.hem_EW +  str(abs(self.lon - (self.lon % 5))).zfill(3)
                
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
    
    def __getGeoT(self, filepath):
        """
        Fetches a GDAL GeoTransform for tile.
        """
        
        ds = gdal.Open(filepath, 0)
        geo_t = ds.GetGeoTransform()
               
        return geo_t
    
    def __getSize(self, filepath):
        """
        Determines the size of a tile.
        """
        
        ds = gdal.Open(filepath, 0)
        x_size = ds.RasterXSize
        y_size = ds.RasterYSize
        
        return x_size, y_size
    
    def getMask(self):
        """
        Loads the mask into a numpy array.
        """
        
        mask_ds = gdal.Open(self.mask_path, 0)
        mask = mask_ds.ReadAsArray()
        
        return mask != 255
    
    def getDOY(self):
        """
        Loads date values into a numpy array.
        """
        
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
        
    
    def getDN(self):
        """
        Loads DN (raw) values into a numpy array.
        """
        
        DN_ds = gdal.Open(self.HV_path, 0)
        DN = DN_ds.ReadAsArray()
        
        return DN
    
        
    def getGamma0(self, units = 'natural'):
        """
        Calibrates data to gamma0 (baskscatter) in decibels or natural units.
        """
        
        assert units == 'natural' or units == 'decibels', "Units must be 'natural' or 'decibels'. You input %s."%units
        
        DN = np.ma.array(self.DN, mask = np.logical_or(self.mask, self.DN == 0))
        
        gamma0 = 10 * np.ma.log10(DN.astype(np.float) ** 2) - 83. # units = decibels
        
        if units == 'decibels':
            return gamma0
        elif units == 'natural':
            return 10 ** (gamma0 / 10.) # Convert to natural units
        
    def getAGB(self):
        """
        Calibrates data to aboveground biomass (AGB).
        Placeholder equation to calibrate backscatter (gamma0) to AGB (tC/ha).
        """
        
        AGB = 715.667 * self.getGamma0() - 5.967
        
        return AGB


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


def outputGeoTiff(data, geo_t, output_dir, output_name = 'output'):
    """
    Writes a GeoTiff file to disk.
    
    Args:
        data: A numpy array containing ALOS mosaic data.
        geo_t: A GDAL geoMatrix (ds.GetGeoTransform()) for data.
        output_dir: Directory to write output file.
        output_name: Optionally specify an output string to precede output file. Defaults to 'output'.
    """
    
    from osgeo import osr, gdal
    
    lon, lat = geo_t[0], geo_t[3]
    
    hem_NS = 'S' if lat < 0 else 'N'
    hem_EW = 'W' if lon < 0 else 'E'
    
    output_path = '%s/%s_%s%s.tif'%(output_dir, output_name, hem_NS + str(abs(lat)).zfill(2), hem_EW + str(abs(lon)).zfill(3))
    
    # Get 
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(output_path, data.shape[0], data.shape[1], 1, gdal.GDT_Float32, )
    ds.SetGeoTransform(geo_t)
    ds.SetProjection(srs.ExportToWkt())
    ds.GetRasterBand(1).WriteArray(data)
    ds = None


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


def rasterizeShapefile(data, shp, buffer_size = 0., binary_mask = True):
    """
    Rasterize points, lines or polygons from a shapefile to match ALOS mosaic data.
        
    Args:
        data:
        shp: Path to a shapefile consisting of points, lines and/or polygons. This does not have to be in the same projection as ds
        buffer_size: Optionally specify a buffer to add around features of the shapefile, in decimal degrees.

    Returns:
        A numpy array with a boolean mask delineating locations inside (True) and outside (False) the shapefile [and optional buffer].
    """
    
    import shapefile
    
    # Determine size of buffer to place around lines/polygons
    buffer_px = int(round(buffer_size / data.geo_t[1]))
    
    # Create output image. Add an buffer around the image array equal to the maxiumum dilation size. This means that features just outside ALOS tile extent can contribute to dilated mask.
    rasterPoly = Image.new("I", (data.ySize + (buffer_px * 2), data.xSize + (buffer_px * 2)), 0)
    rasterize = ImageDraw.Draw(rasterPoly)
    
    # The shapefile may not have the same CRS as ALOS mosaic data, so this will generate a function to reproject points.
    coordTransform = _coordinateTransformer(shp)
    
    # Read shapefile
    sf = shapefile.Reader(shp) 
    
    # For each shape in shapefile...
    for n,shape in enumerate(sf.iterShapes()):
        
        # If shape number not requried
        if binary_mask:
            n = 1
        
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
                rasterize.point(pixels, n)

            # a line...
            elif shape.shapeType == 3:
                rasterize.line(pixels, n)
  
            # or a polygon.
            elif shape.shapeType == 5:  
                rasterize.polygon(pixels, n)
    
    #Converts a Python Imaging Library array to a gdalnumeric image.
    mask = gdalnumeric.fromstring(rasterPoly.tobytes(),dtype=np.uint32)
    mask.shape = rasterPoly.im.size[1], rasterPoly.im.size[0]
    
    # If any buffer pixels are slected, dilate the masked area by buffer_px pixels
    if buffer_px > 0:
        mask = dilateMask(mask, buffer_px)
    
    # Get rid of image buffer
    mask = mask[buffer_px:mask.shape[0]-buffer_px, buffer_px:mask.shape[1]-buffer_px]
    
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
    




if __name__ == '__main__':
    
    data_dir = '/home/sbowers3/DATA/ALOS_data/ALOS_mosaic/gorongosa/'
    output_dir = '/home/sbowers3/DATA/ALOS_data/ALOS_mosaic/gorongosa/'

    t1 = 2007
    t2 = 2010

    lat = -18#-9#-11
    lon = 33#34#39
    
    data_t1 = ALOS(lat, lon, t1, data_dir)
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







