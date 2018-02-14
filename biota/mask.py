
import numpy as np
from PIL import Image, ImageDraw
import scipy.ndimage as ndimage

import pdb

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
        buffer_size: Optionally specify a buffer to add around features of the shapefile, in meters.

    Returns:
        A numpy array with a boolean mask delineating locations inside (True) and outside (False) the shapefile [and optional buffer].
    """
    
    import shapefile
    from osgeo import gdalnumeric
    
    assert np.logical_or(np.logical_and(field == None, value == None), np.logical_and(field != None, value != None)), "If specifying field or value, both must be defined. At present, field = %s and value = %s"%(str(field), str(value))
    
    # Determine size of buffer to place around lines/polygons
    #buffer_px = int(round(buffer_size / data.geo_t[1]))
    
    # Determine the size of the buffer in degrees
    buffer_size_degrees = buffer_size / (((data.xRes * data.xSize) + (data.yRes * data.ySize)) / 2.)
    
    # Determine size of buffer to place around lines/polygons
    #buffer_px = int(round(buffer_size / ((data.xRes + data.yRes) / 2.)))
    buffer_px = int(round(buffer_size_degrees / data.geo_t[1]))
    
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
        if sxmax < data.geo_t[0] - buffer_size_degrees: continue
        if sxmin > data.geo_t[0] + (data.geo_t[1] * data.ySize) + buffer_size_degrees: continue
        if symax < data.geo_t[3] + (data.geo_t[5] * data.xSize) + buffer_size_degrees: continue
        if symin > data.geo_t[3] - buffer_size_degrees: continue
        
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
                pixels.append(_world2Pixel(data.geo_t, lon, lat, buffer_size = buffer_size_degrees))

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
