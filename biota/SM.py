
import datetime as dt
import glob
import math
import numpy as np
from osgeo import gdal
import scipy.interpolate

import biota

import pdb
    


def _readTimeStack(tile, date, search_days = 7):
    '''
    Load soil moisture time stack
    '''
    
    # Load list of netCDF files
    nc_files = np.array(sorted(glob.glob(tile.SM_dir + '/*.nc')))
    
    assert len(nc_files) > 0, "No data found in soil moisture data path (%s)."%nc_files
    
    # Grab their dates
    nc_dates = np.array([dt.datetime.strptime(nc_file.split('/')[-1].split('-')[-2][:8], '%Y%m%d').date() for nc_file in nc_files])

    sel = np.logical_and(nc_dates >= date - dt.timedelta(search_days), nc_dates <= date + dt.timedelta(search_days))

    assert len(nc_files[sel]) > 0, "No data for tile date found in soil moisture data path (%s)."%nc_files

    for n, nc_file in enumerate(nc_files[sel]):
            
        # Load dataset
        ds = gdal.Open('NETCDF:%s:sm'%nc_file)
        
        # Get GeoTransform()
        geo_t = ds.GetGeoTransform()
        
        # Get netCDF extent
        ulLon, ulLat = geo_t[0], geo_t[3]
        lonDist, latDist = geo_t[1], geo_t[5]
        
        # Get UL coord for tile
        xOffset = int(math.floor((tile.lon - ulLon) / lonDist))
        yOffset = int(math.ceil((ulLat - tile.lat) / latDist)) * -1
        
        # And the extent of pixels to load
        yCount = int(abs(math.floor(1. / latDist)))
        xCount = int(math.ceil(1. / lonDist))
        
        this_sm = ds.ReadAsArray(xOffset, yOffset, xCount, yCount)
        
        if 'sm_out' not in locals():
            sm_out = np.ma.zeros((this_sm.shape[0], this_sm.shape[1], (search_days * 2) + 1), dtype = np.float32)
        
        sm_out[:,:,n] = np.ma.array(this_sm, mask = this_sm == -9999.)
            
    geo_t_out = (ulLon + (xOffset * lonDist), geo_t[1], geo_t[2], ulLat + (yOffset * latDist), geo_t[4], geo_t[5])
    
    return sm_out, geo_t_out


def _interpolateTime(sm):
    '''
    '''

    # Interpolate through time
    sm_interp = np.zeros((sm.shape[0], sm.shape[1]), dtype = np.float32)
    
    for i in range(sm.shape[0]):
        for j in range(sm.shape[1]):
            
            this_sm = sm[i,j,:]
            x = np.arange(sm.shape[-1])
            y = this_sm.data
            
            if np.sum(this_sm.mask == False) > 1:
                f = scipy.interpolate.interp1d(x[this_sm.mask == False], y[this_sm.mask == False], fill_value = 'extrapolate')
                interp = f(x)[(x[this_sm.mask == False]).shape[0]/2]
            else:
                interp = -9999.
            
            sm_interp[i, j] = interp
            
    return sm_interp

def _resampleSM(sm_interp, tile, geo_t):
    '''
    '''
        
    # Create output file matching ALOS tile
    gdal_driver = gdal.GetDriverByName('MEM')
    
    ds_source = gdal_driver.Create('', sm_interp.shape[0], sm_interp.shape[1], 1, gdal.GDT_Float32)
    ds_source.SetGeoTransform(geo_t)
    ds_source.SetProjection(tile.proj)
    ds_source.GetRasterBand(1).WriteArray(sm_interp)
    
    ds_dest = gdal_driver.Create('', tile.ySize, tile.xSize, 1, gdal.GDT_Float32)
    ds_dest.SetGeoTransform(tile.geo_t)
    ds_dest.SetProjection(tile.proj)
    
    # Reproject input GeoTiff to match the ALOS tile
    gdal.ReprojectImage(ds_source, ds_dest, tile.proj, tile.proj, gdal.GRA_NearestNeighbour)
    
    # Load resampled image into memory
    sm_resampled = ds_dest.GetRasterBand(1).ReadAsArray()
    
    return np.ma.array(sm_resampled, mask = sm_resampled == -9999.)


def _loadSM(tile, date, search_days = 7):
    '''
    Load and reproject soil moisture data
    '''
    
    # Load SM time stack (within 1 week of measurement)
    sm, geo_t = _readTimeStack(tile, date, search_days = search_days)
    
    # Interpolate over time gaps
    sm_interp = _interpolateTime(sm)
    
    # Reproject to match ALOS tile (nearest neighbor)
    sm_resampled = _resampleSM(sm_interp, tile, geo_t)
    
    return sm_resampled

            
def getSM(tile, search_days = 7):
    """
    Function to load a resampled soil moisture image from the ESA CCI soil moisture product. The function returns a masked array with estimated volumetric soil moisture content (m^2/m^2), averaged for each satellite overpass in a tile.
    
    Args:
        tile: A biota.LoadTile() object.
        
    Returns:
        A masked array of estimated soil moisture (m^2/m^2)
    """
        
    # Build an array of dates in string format, which is quicker to search
    dates = tile.getDate().astype(np.int)
    
    # Generate output array
    out = np.zeros_like(tile.mask, dtype = np.float32) + 999999.
    
    for date in np.unique(dates): #np.unique(tile.getDate()).astype(np.int):        
        
        # Nodata (0 = 1st Jan 1970)
        if date == 0:
            continue
        
        # Interpolate through time and draw out the central value, and zoom to scale of tile
        sm_resampled = _loadSM(tile, dt.datetime(1970,1,1).date() + dt.timedelta(date), search_days = search_days)
        
        #out[tile.getDate() == date] = sm_resampled[tile.getDate() == date]
        out[dates == date] = np.ma.mean(sm_resampled[dates == date])
        
        print 'TODO: min proportional coverage'
    
    # Add the mask
    out = np.ma.array(out, mask = np.logical_or(tile.mask, out == 999999.))
        
    return out
    

if __name__ == '__main__':
    """
    """
    
    sm_dir = '/home/sbowers3/guasha/sam_bowers/soil_moisture/'
    alos_dir = '/home/sbowers3/SMFM/ALOS_data/'
    output_dir = '/home/sbowers3/guasha/sam_bowers/carla_outputs/'
    
    for lat in range(-15, -5, 1):
        for lon in range(30, 40, 1):
    
            #lat, lon = -15, 30
            print lat, lon
            #lat, lon = -10, 17
            
            #lat = -15; lon = 35
            
            try:
                tile_2007 = biota.LoadTile(alos_dir, lat, lon, 2007, downsample_factor = 1)
                tile_2010 = biota.LoadTile(alos_dir, lat, lon, 2010, downsample_factor = 1)
                tile_change = biota.LoadChange(tile_2007, tile_2010)
                
                sm_2007 = getSM(tile_2007, sm_dir)
                sm_2010 = getSM(tile_2010, sm_dir)
                sm_change = sm_2010 - sm_2007
                
                sm_change_out = (sm_change * 0).astype(np.int)
                np.ma.set_fill_value(sm_change_out, 999999)
                sm_change_out[np.logical_and(sm_change.data < -0.1, sm_change.mask == False)] = 1
                sm_change_out[np.logical_and(sm_change.data > 0.1, sm_change.mask == False)] = 1
                
                #biota.IO.showFigure(sm_2007, lat, lon, title = 'Soil moisture', cbartitle = 'm^2/m^2', vmin = 0, vmax = 0.25, cmap = 'Blues')
                
            except:
                print '    No data for this tile'
                continue
            
            
            biota.IO.outputGeoTiff(sm_2007, tile_2007.output_pattern%'SM', tile_2007.geo_t, tile_2007.proj, output_dir = output_dir, dtype = 6, nodata = 999999.)
            biota.IO.outputGeoTiff(sm_2010, tile_2010.output_pattern%'SM', tile_2010.geo_t, tile_2010.proj, output_dir = output_dir, dtype = 6, nodata = 999999.)
            biota.IO.outputGeoTiff(sm_change_out, tile_change.output_pattern%'SM', tile_change.geo_t, tile_change.proj, output_dir = output_dir, dtype = 5, nodata = 999999)
            