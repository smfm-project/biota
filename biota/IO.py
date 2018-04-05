
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
    


def _buildMap(fig, ax, data, lat, lon, title ='', cbartitle = '', vmin = 10., vmax = 60., cmap = 'YlGn'):
    """
    Builds a standardised map for overviewFigure().
    """
    
    import matplotlib.pyplot as plt
        
    im = ax.imshow(data, vmin = vmin, vmax = vmax, cmap = cmap, interpolation = 'nearest')
    
    ax.set_yticks(np.arange(0,data.shape[0],data.shape[0]/10))
    ax.set_xticks(np.arange(0,data.shape[1],data.shape[1]/10))
    ax.set_yticklabels(np.arange(lat, lat - 1.01, - 0.1))
    ax.set_xticklabels(np.arange(lon, lon + 1.01, 0.1))
    ax.tick_params(labelsize = 5)
    ax.set_xlabel('Longitude', fontsize = 5)
    ax.set_ylabel('Latitude', fontsize = 5)
    ax.set_title(title, fontsize = 8)
    
    cbar = fig.colorbar(im, ax = ax, fraction = 0.046, pad = 0.04)
    cbar.ax.tick_params(labelsize = 6)
    cbar.set_label(cbartitle, fontsize = 7)
    

def overviewFigure(tile_change, output = False, show = True):
    """
    Generate an overview image showing biomass and proportional biomass change for an ALOS tile.
    
    Args:
        tile_change: An ALOS change object from biota.LoadChange()
    """
    
    import matplotlib.pyplot as plt
    
    # Load AGB, and update masks to exclude areas outisde forest definition. Good for visualisation
    AGB_t1 = tile_change.data_t1.getAGB()
    AGB_t2 = tile_change.data_t2.getAGB()
    
    # Mask out areas < 10 tC/ha
    AGB_t1 = np.ma.array(AGB_t1, mask = np.logical_or(AGB_t1.mask, AGB_t1 < 10.))
    AGB_t2 = np.ma.array(AGB_t2, mask = np.logical_or(AGB_t2.mask, AGB_t1 < 10.))
        
    AGB_change = AGB_t2 - AGB_t1

    AGB_pcChange = 100 * (AGB_change / AGB_t1) # %
    
    fig = plt.figure(figsize = (7, 6))
    
    # Plot a map of AGB at t1
    ax1 = fig.add_subplot(2, 2, 1)
    _buildMap(fig, ax1, AGB_t1, tile_change.lat, tile_change.lon, title = 'AGB %s'%str(tile_change.year_t1), cbartitle = 'tC/ha')
    
    # Plot a map of AGB at t2
    ax2 = fig.add_subplot(2, 2, 2)
    _buildMap(fig, ax2, AGB_t2, tile_change.lat, tile_change.lon, title = 'AGB %s'%str(tile_change.year_t2), cbartitle = 'tC/ha')    
    
    # Plot a map of absolute AGB change   
    ax3 = fig.add_subplot(2, 2, 3)
    _buildMap(fig, ax3, AGB_change, tile_change.lat, tile_change.lon, title = 'AGB change (%s-%s)'%(str(tile_change.data_t1.year),str(tile_change.year_t2)),
              cbartitle = 'tC/ha', vmin = -10., vmax = 10., cmap = 'RdBu')    
    
    # Plot a map of % AGB change
    ax4 = fig.add_subplot(2, 2, 4)
    _buildMap(fig, ax4, AGB_pcChange, tile_change.lat, tile_change.lon, title = 'AGB change (%s-%s)'%(str(tile_change.data_t1.year),str(tile_change.year_t2)),
              cbartitle = '%', vmin = -50., vmax = 50., cmap = 'RdBu')    
    
    plt.tight_layout()
    
    if output:
        output_pattern = tile_change.output_pattern.replace('.tif','.png')
    
        output_path = os.path.abspath(os.path.expanduser('%s/%s'%(tile_change.output_dir, output_pattern%('OverviewFigure'))))
    
        plt.savefig(output_path, dpi = 150)
    
    if show:
        plt.show()
        
    plt.close()


