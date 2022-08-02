
import itertools
import math
import numpy as np
import os
from PIL import Image, ImageDraw
import scipy.ndimage as ndimage

import biota.IO

import pdb


def dilateMask(mask, buffer_px, location_id = False):
    """
    Dilate a boolean (True/False) numpy array by a specified number of pixels.

    Args:
        mask: A boolean (True/False) numpy array, with 'True' representing locations to add a buffer.
        buffer_px: A number of pixels to add around each 'True' array element.

    Returns:
        The mask array with dilated 'True' locations.
    """

    if location_id == False:
        mask_dilated = ndimage.morphology.binary_dilation(mask, iterations = buffer_px)

    else:
        mask_dilated = np.zeros_like(mask)
        for i in np.unique(mask[mask > 0]):
            mask_dilated[ndimage.morphology.binary_dilation(mask == i, iterations = buffer_px)] = i

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

    pixel = (x - ulX) / xDist
    line = (y - ulY) / yDist

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

    # Extract data type from shapefile. Interprets N (int), F (float) and C (string), sets others to string.
    this_dtype = sf.fields[1:][field_n][1]

    if this_dtype == 'N':
        dtype = np.int
    elif this_dtype == 'F':
        dtype = np.float32
    elif this_dtype == 'C':
        dtype = np.str
    else:
        dtype = np.str

    value_out = []

    # Cycle through records:
    for s in sf.records():
        value_out.append(s[field_n])

    return np.array(value_out, dtype = dtype)


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

def maskArray(tile, array, classes = [], buffer_size = 0.):
    '''
    Extract a mask from a numpy array based on specified classes

    Args:
        tile: An ALOS tile (biota.LoadTile())
        array: A numpy array
        classes: A list of values to add to be masked
        buffer_size: Optionally specify a buffer to add around maked pixels, in meters.

    Returns:
        A numpy array with a boolean mask
    '''

    assert (tile.ySize, tile.xSize) == array.shape, "Numpy array shape must be identical to the tile extent."

    # If a binary mask is input, assume that True should be added to mask
    if array.dtype == np.bool and classes == []:
        classes = [True]

    # Identify pixels that are classes to be masked
    mask = np.in1d(array, classes).reshape(array.shape)

    if buffer_size > 0.:

        # Determine the size of the buffer in pixels
        buffer_px = int(buffer_size / ((tile.xRes + tile.yRes) / 2.))

        # Dilate the mask
        mask = dilateMask(mask, buffer_px)

    return mask


def maskRaster(tile, raster, classes = [], buffer_size = 0.):
    '''
    Extract a mask from a GeoTiff based on specified classes

    Args:
        tile: An ALOS tile (biota.LoadTile())
        raster: A GeoTiff or VRT file, with integer values
        classes: A list of values to add to be masked
        buffer_size: Optionally specify a buffer to add around maked pixels, in meters.

    Returns:
        A numpy array with a boolean mask
    '''

    from osgeo import gdal

    assert raster.rstrip('/').split('.')[-1] == 'tif' or raster.rstrip('/').split('.')[-1] == 'tiff' or raster.rstrip('/').split('.')[-1] == 'vrt', "raster input must be a GeoTiff or VRT file."

    raster = os.path.expanduser(raster)
    assert os.path.exists(raster), "GeoTiff file %s does not exist in the file system."%raster

    # Load raster + reporject to match tile
    reampled_image = biota.IO.loadRaster(raster, tile)

    # Identify pixels that are classes to be masked
    mask = maskArray(tile, reampled_image, classes = classes, buffer_size = buffer_size)

    return mask


def maskShapefile(tile, shp, buffer_size = 0., field = None, value = None, location_id = False):
    """
    Rasterize points, lines or polygons from a shapefile to match ALOS mosaic data.

    Args:
        tile: An ALOS tile (biota.LoadTile())
        shp: Path to a shapefile consisting of points, lines and/or polygons. This does not have to be in the same projection as ds
        buffer_size: Optionally specify a buffer to add around features of the shapefile, in meters.
        field: Optionally specify a single shapefile field to extract (you must also specify its value)
        value: Optionally specify a single shapefile field value to extract (you must also specify the field name)
        location_id: Set True to return a unique ID for each masked shape. Note: This is not zero indexed, but starts at 1.

    Returns:
        A numpy array with a boolean (or integer) mask delineating locations inside and outside the shapefile and optional buffer.
    """

    import shapefile
    from osgeo import gdalnumeric

    assert np.logical_or(np.logical_and(field == None, value == None), np.logical_and(field != None, value != None)), "If specifying field or value, both must be defined. At present, field = %s and value = %s"%(str(field), str(value))

    shp = os.path.expanduser(shp)
    assert os.path.exists(shp), "Shapefile %s does not exist in the file system."%shp

    # Determine the size of the buffer in degrees
    buffer_size_degrees = buffer_size / (((tile.xRes * tile.xSize) + (tile.yRes * tile.ySize)) / 2.)

    # Determine size of buffer to place around lines/polygons
    buffer_px = int(round(buffer_size_degrees / tile.geo_t[1]))

    # Create output image. Add a buffer around the image array equal to the maxiumum dilation size. This means that features just outside ALOS tile extent can contribute to dilated mask.
    rasterPoly = Image.new("I", (tile.xSize + (buffer_px * 2), tile.ySize + (buffer_px * 2)), 0)
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
    for n, shape in enumerate(shapes):

        # Get shape bounding box
        if shape.shapeType == 1 or shape.shapeType == 11:
            # Points don't have a bbox, calculate manually
            sxmin = np.min(np.array(shape.points)[:,0])
            sxmax = np.max(np.array(shape.points)[:,0])
            symin = np.min(np.array(shape.points)[:,1])
            symax = np.max(np.array(shape.points)[:,1])
        else:
            sxmin, symin, sxmax, symax = shape.bbox

        # Transform bounding box points
        sxmin, symin, z = coordTransform.TransformPoint(sxmin, symin)
        sxmax, symax, z = coordTransform.TransformPoint(sxmax, symax)

        # Go to the next record if out of bounds
        if sxmax < tile.geo_t[0] - buffer_size_degrees: continue
        if sxmin > tile.geo_t[0] + (tile.geo_t[1] * tile.xSize) + buffer_size_degrees: continue
        if symax < tile.geo_t[3] + (tile.geo_t[5] * tile.ySize) + buffer_size_degrees: continue
        if symin > tile.geo_t[3] - buffer_size_degrees: continue

        #Separate polygons with list indices
        n_parts = len(shape.parts) #Number of parts
        indices = shape.parts #Get indices of shapefile part starts
        indices.append(len(shape.points)) #Add index of final vertex

        # Catch to allow use of point shapefiles, which don't have parts
        if shape.shapeType == 1 or shape.shapeType == 11:
            n_parts = 1
            points = shape.points

        for part in range(n_parts):

            if shape.shapeType != 1 and shape.shapeType != 11:

                start_index = shape.parts[part]
                end_index = shape.parts[part+1]
                points = shape.points[start_index:end_index] #Map coordinates

            pixels = [] #Pixel coordinantes

            # Transform coordinates to pixel values
            for p in points:

                # First update points from shapefile projection to ALOS mosaic projection
                lon, lat, z = coordTransform.TransformPoint(p[0], p[1])

                # Then convert map to pixel coordinates using geo transform
                pixels.append(_world2Pixel(tile.geo_t, lon, lat, buffer_size = buffer_size_degrees))

            # Draw the mask for this shape...
            # if a point...
            if shape.shapeType == 0 or shape.shapeType == 1 or shape.shapeType == 11:
                rasterize.point(pixels, n+1)

            # a line...
            elif shape.shapeType == 3 or shape.shapeType == 13:
                rasterize.line(pixels, n+1)

            # or a polygon.
            elif shape.shapeType == 5 or shape.shapeType == 15:
                rasterize.polygon(pixels, n+1)

            else:
                print('Shapefile type %s not recognised!'%(str(shape.shapeType)))

    #Converts a Python Imaging Library array to a gdalnumeric image.
    mask = gdalnumeric.fromstring(rasterPoly.tobytes(),dtype=np.uint32)
    mask.shape = rasterPoly.im.size[1], rasterPoly.im.size[0]

    # If any buffer pixels are slected, dilate the masked area by buffer_px pixels
    if buffer_px > 0:
        mask = dilateMask(mask, buffer_px, location_id = location_id)

    # Get rid of image buffer
    mask = mask[buffer_px:mask.shape[0]-buffer_px, buffer_px:mask.shape[1]-buffer_px]

    if location_id == False:
        # Get rid of record numbers
        mask = mask > 0

    return mask


def getTilesInShapefile(shp, field = None, value = None):
    """
    Identify all the ALOS tiles that fall within a shapefile.

    Args:
        shp: Path to a shapefile consisting of polygons. This can be in any projection.

    Returns:
        The lat/lon indicators of which ALOS tiles are covered by the shapefile
    """

    import shapefile

    assert np.logical_or(np.logical_and(field == None, value == None), np.logical_and(field != None, value != None)), "If specifying field or value, both must be defined. At present, field = %s and value = %s"%(str(field), str(value))

    # The shapefile may not have the same CRS as ALOS mosaic data, so this will generate a function to reproject points.
    coordTransform = _coordinateTransformer(shp)

    lats, lons = [], []
    tiles_to_include = set([])

    # Get shapes
    shapes = np.array(shapefile.Reader(shp).shapes())

    if field != None:
        shapes = shapes[getField(shp, field) == value]

    for shape in shapes:
        # Get the bbox for each shape in the shapefile
        if shape.shapeType == 1 or shape.shapeType == 11:
            # Points don't have a bbox, calculate manually
            sxmin = np.min(np.array(shape.points)[:,0])
            sxmax = np.max(np.array(shape.points)[:,0])
            symin = np.min(np.array(shape.points)[:,1])
            symax = np.max(np.array(shape.points)[:,1])
        else:
            sxmin, symin, sxmax, symax = shape.bbox

        # Transform points to WGS84
        lonmin, latmin, z = coordTransform.TransformPoint(sxmin, symin)
        lonmax, latmax, z = coordTransform.TransformPoint(sxmax, symax)

        # Get the tiles that cover the area of the shapefile
        latrange = list(range(int(math.ceil(latmin)), int(math.ceil(latmax)+1), 1))
        lonrange = list(range(int(math.floor(lonmin)), int(math.floor(lonmax)+1), 1))
        tiles = list(itertools.product(latrange,lonrange))

        # Add them to tiles_to_include if not already tere
        [tiles_to_include.add(t) for t in tiles]

    return sorted(list(tiles_to_include))


def updateMask(tile, filename, buffer_size = 0., classes = []):
    """
    Function to generate a VRT, GeoTiff, shapefile, or numpy array mask to match an ALOS tile.

    Args:
        ::
        tile: An ALOS tile (biota.LoadTile())
        filename: A GeoTiff, VRT or shapefile (string)m or a numpy array
        buffer_size: Optionally specify a buffer to add around maked pixels, in meters.
        classes: A list of values to add to be masked from a raster input. May be omitted where a boolean numpy array is input.

    Returns:
        A boolean numpy array
    """

    if type(filename) == str:

        file_type = filename.split('/')[-1].split('.')[-1]

        assert file_type in ['shp', 'tif', 'tiff', 'vrt'], "Input must be a numpy array, GeoTiff, VRT, or a shapefile."

    elif type(filename) == np.ndarray or type(filename) == np.ma.core.MaskedArray:
        file_type = 'array'

    else:

        assert False, "Input must be a numpy array, GeoTiff, VRT, or a shapefile."

    if file_type == 'shp':

        # Rasterize the shapefile, optionally with a buffer
        mask = biota.mask.maskShapefile(tile, filename, buffer_size = buffer_size)

    elif file_type in ['tif', 'tiff', 'vrt']:

        assert classes != [], "If adding a GeoTiff or VRT file to the mask, you must also specify the class values to add to the mask (e.g. classes = [20, 160, 170, 190, 210])."

        # Resample and extract values from shapefile, optionally with a buffer
        mask = biota.mask.maskRaster(tile, filename, classes = classes, buffer_size = buffer_size)

    else:

        if filename.dtype != np.bool:
            assert classes != [], "If adding a non-boolean numpy array file to the mask, you must also specify the class values to add to the mask (e.g. classes = [20, 160, 170, 190, 210])."

        mask = biota.mask.maskArray(tile, filename, classes = classes, buffer_size = buffer_size)

    return mask
