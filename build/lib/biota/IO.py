
import matplotlib.pyplot as plt
import numpy as np
import os


import pdb

def loadArray(filepath):
    """
    Use gdal to load a geospatial image into a numpy array
    """
    
    from osgeo import gdal
    
    ds = gdal.Open(filepath, 0)
    
    return ds.ReadAsArray()
    
    
def loadGeoTransform(filepath):
    """
    Use gdal to load a gdal geotransform (affine transform) tuple
    """
    
    from osgeo import gdal
        
    ds = gdal.Open(filepath, 0)
    
    return ds.GetGeoTransform()


def loadProjection(filepath):
    """
    Use gdal to load projection info
    """
    
    from osgeo import gdal
    
    ds = gdal.Open(filepath, 0)
    
    return ds.GetProjection()


def loadSize(filepath):
    """
    Use gdal to load raster
    """
    
    from osgeo import gdal
    
    ds = gdal.Open(filepath, 0)
    
    return ds.RasterYSize, ds.RasterXSize


def loadRaster(raster, tile, resampling = 0, dtype = 3):
    """
    Loads a raster image (e.g. GeoTiff), reprojecting to match a tile.
    
    Args:
        raster: path to raster file (string)
        tile: A tile from biota.LoadTile()
        resampling: Optionally specify a gdal resampling type. Defaults to gdal.GRA_NearestNeighbour.
        
    Returns:
        A numpy array matching tile
    """
    
    from osgeo import gdal
    
    # Open GeoTiff and get metadata
    ds_source = gdal.Open(raster)
    proj_source = loadProjection(raster)
    
    # Create output file matching ALOS tile
    gdal_driver = gdal.GetDriverByName('MEM')
    ds_dest = gdal_driver.Create('', tile.xSize, tile.ySize, 1, dtype)
    ds_dest.SetGeoTransform(tile.geo_t)
    ds_dest.SetProjection(tile.proj)
       
    # Reproject input GeoTiff to match the ALOS tile
    gdal.ReprojectImage(ds_source, ds_dest, proj_source, tile.proj, resampling)
    
    # Load resampled image into memory
    resampled = ds_dest.GetRasterBand(1).ReadAsArray()
    
    return resampled

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
    output_path = '%s/%s.tif'%(os.path.abspath(os.path.expanduser(output_dir)), filename.rstrip('.tif'))
    
    # Save image with georeference info
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(output_path, data.shape[1], data.shape[0], 1, dtype, options = ['COMPRESS=LZW'])
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
    


def buildMap(fig, ax, data, lat, lon, title ='', cbartitle = '', vmin = None, vmax = None, cmap = None, big_labels = False):
    """
    Builds a standardised map for overviewFigure() and showFigure().
    """
    
    labelsize = 5
    titlesize = 8
    cbarlabelsize = 6
    cbartitlesize = 7
    
    if big_labels:
        labelsize *= 2
        titlesize *= 2
        cbarlabelsize *= 2
        cbartitlesize *= 2
    
        
    im = ax.imshow(data, vmin = vmin, vmax = vmax, cmap = cmap, interpolation = 'nearest', aspect = 'auto')
    
    ax.set_yticks(np.arange(0, data.shape[0] + data.shape[0]/10, data.shape[0]/10))
    ax.set_xticks(np.arange(0, data.shape[1] + data.shape[1]/10, data.shape[1]/10))
    ax.set_yticklabels(["%.1f" %i for i in np.arange(lat, lat - 1.01, - 0.1)])
    ax.set_xticklabels(["%.1f" %i for i in np.arange(lon, lon + 1.01, 0.1)])
    ax.tick_params(labelsize = labelsize)
    ax.set_ylabel('Latitude', fontsize = labelsize)
    ax.set_xlabel('Longitude', fontsize = labelsize)
    ax.set_title(title, fontsize = titlesize)
    
    cbar = fig.colorbar(im, ax = ax, fraction = 0.046, pad = 0.04)
    cbar.ax.tick_params(labelsize = cbarlabelsize)
    cbar.set_label(cbartitle, fontsize = cbartitlesize)


def showFigure(data, lat, lon, title = None, cbartitle = None, vmin = None, vmax = None, cmap = None, big_labels = True):
    '''
    Show an overview image.
    
    Args:
        data: A numpy array containing the data to display
        lat: Latitude, an integer
        lon: Longitude, an integer
        title: A string to show as the figure title
        cbartitle: A string to show as the colorbar title
        vmin: Min colorbar range
        vmax: Max colorbar range
        cmap: String of a matplotlib colorbar
    '''
    
    # Set up new figure
    fig = plt.figure(figsize = (7, 6))
    ax = fig.add_subplot(1, 1, 1)
    
    # Plot map
    buildMap(fig, ax, data, lat, lon, title = title, cbartitle = cbartitle, vmin = vmin, vmax = vmax, cmap = cmap, big_labels = big_labels)
    
    # Show in window
    plt.show()
        
    # Tidy up
    plt.close()
        

def overviewFigure(tile_change, output = False, show = True):
    """
    Generate an overview image showing biomass and proportional biomass change for an ALOS tile.
    
    Args:
        tile_change: An ALOS change object from biota.LoadChange()
        output: Set to True to output a summary png image
        show: Set to True to display image on screen
    """
    
    import matplotlib.pyplot as plt
    
    # Load AGB, and update masks to exclude areas outisde forest definition. Good for visualisation
    AGB_t1 = tile_change.tile_t1.getAGB()
    AGB_t2 = tile_change.tile_t2.getAGB()
    
    # Mask out areas < 10 tC/ha
    AGB_t1 = np.ma.array(AGB_t1, mask = np.logical_or(AGB_t1.mask, AGB_t1 < 10.))
    AGB_t2 = np.ma.array(AGB_t2, mask = np.logical_or(AGB_t2.mask, AGB_t1 < 10.))
        
    AGB_change = AGB_t2 - AGB_t1
    #AGB_pcChange = 100 * (AGB_change / AGB_t1) # %
    
    # Calculate change type
    change_type = tile_change.getChangeType()
    change_code = tile_change.ChangeCode
    
    # Set minor loss and minor gain to nodata
    change_code[np.logical_or(change_code == 3, change_code == 4)] = 0
    change_code.mask[change_type.data == 0] = True
    
    fig = plt.figure(figsize = (7, 6))
    
    # Plot a map of AGB at t1
    ax1 = fig.add_subplot(2, 2, 1)
    buildMap(fig, ax1, AGB_t1, tile_change.lat, tile_change.lon, title = 'AGB %s'%str(tile_change.year_t1), cbartitle = 'tC/ha', vmin = 10., vmax = 40., cmap = 'YlGn')
    
    # Plot a map of AGB at t2
    ax2 = fig.add_subplot(2, 2, 2, sharex = ax1, sharey = ax1)
    buildMap(fig, ax2, AGB_t2, tile_change.lat, tile_change.lon, title = 'AGB %s'%str(tile_change.year_t2), cbartitle = 'tC/ha', vmin = 10., vmax = 40., cmap = 'YlGn')    
    
    # Plot a map of absolute AGB change   
    ax3 = fig.add_subplot(2, 2, 3, sharex = ax1, sharey = ax1)
    buildMap(fig, ax3, AGB_change, tile_change.lat, tile_change.lon, title = 'AGB change (%s-%s)'%(str(tile_change.tile_t1.year),str(tile_change.year_t2)),
              cbartitle = 'tC/ha', vmin = -10., vmax = 10., cmap = 'RdBu')    
    
    # Plot a map of % AGB change
    ax4 = fig.add_subplot(2, 2, 4, sharex = ax1, sharey = ax1)
    buildMap(fig, ax4, change_code, tile_change.lat, tile_change.lon, title = 'Change type (%s-%s)'%(str(tile_change.tile_t1.year),str(tile_change.year_t2)),
              vmin = 1., vmax = 6., cmap = 'Spectral')
    
    plt.tight_layout()
    
    # Output image to png
    if output:
        output_pattern = tile_change.output_pattern.replace('.tif','.png')
    
        output_path = os.path.abspath(os.path.expanduser('%s/%s'%(tile_change.output_dir, output_pattern%('OverviewFigure'))))
    
        plt.savefig(output_path, dpi = 150)
    
    # Display image on screen
    if show:
        plt.show()
        
    plt.close()


