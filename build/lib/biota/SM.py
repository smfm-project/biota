
import datetime as dt
import glob
import math
import numpy as np
from osgeo import gdal
import scipy.interpolate

import biota

import pdb


def _findFiles(tile, date, search_days = 7, data_source = "CCI"):
    """
    """

    assert data_source == "CCI" or data_source == "SMAP", "Soil moisture data_source must be CCI or SMAP."

    if data_source == "CCI":
        data_files = np.array(sorted(glob.glob(tile.SM_dir + '/*.nc')))
    elif data_source == "SMAP":
        data_files = np.array(sorted(glob.glob(tile.SM_dir + '/*.h5')))

    assert len(data_files) > 0, "No data found in soil moisture data path (%s)."%data_files

    # Grab their dates
    if data_source == "CCI":
        data_dates = np.array([dt.datetime.strptime(data_file.split('/')[-1].split('-')[-2][:8], '%Y%m%d').date() for data_file in data_files])

    elif data_source == "SMAP":
        data_dates = np.array([dt.datetime.strptime(data_file.split('/')[-1].split('_')[-3], '%Y%m%d').date() for data_file in data_files])


    sel = np.logical_and(data_dates >= date - dt.timedelta(search_days), data_dates <= date + dt.timedelta(search_days))

    assert len(data_files[sel]) > 0, "No data for tile date found in soil moisture data path (%s)."%tile.SM_dir

    return data_files[sel]

def _readSM(sm_file, data_source = "CCI"):
    """
    """

    assert data_source == "CCI" or data_source == "SMAP", "Soil moisture data_source must be CCI or SMAP."

    if data_source == "CCI":

        # Load dataset
        ds = gdal.Open('NETCDF:%s:sm'%sm_file)

        # Get GeoTransform()
        geo_t = ds.GetGeoTransform()

    else:

        # Load dataset
        ds = gdal.Open('HDF5:"%s"://Soil_Moisture_Retrieval_Data_AM/soil_moisture'%sm_file)

        # Build geo_t equivalent
        latitudes = gdal.Open('HDF5:"%s"://Soil_Moisture_Retrieval_Data_AM/latitude'%sm_file).ReadAsArray()
        longitudes = gdal.Open('HDF5:"%s"://Soil_Moisture_Retrieval_Data_AM/longitude'%sm_file).ReadAsArray()

        # Get GeoTransform()
        geo_t = (-180, 360./3856, 0.0, 84.65642 + ((360./3856)/2.), 0.0, -(360./3856))

    return ds, geo_t


def _readTimeStack(tile, date, search_days = 7, data_source = "CCI"):
    '''
    Load soil moisture time stack
    '''

    # Load list of netCDF files
    sm_files = _findFiles(tile, date, search_days = search_days, data_source = data_source)

    for n, sm_file in enumerate(sm_files):

        ds, geo_t = _readSM(sm_file, data_source = data_source)

        # Get file extent
        ulLon, ulLat = geo_t[0], geo_t[3]
        lonDist, latDist = geo_t[1], geo_t[5]

        # Get UL coord for tile
        xOffset = int(math.floor((tile.lon - ulLon) / lonDist))
        yOffset = int(math.floor((ulLat - tile.lat) / latDist * -1))

        # And the extent of pixels to load. Add one in case of un-aligned pixels
        yCount = int(abs(math.floor(1. / latDist))+1)
        xCount = int(math.ceil(1. / lonDist)+1)

        this_sm = ds.ReadAsArray(xOffset, yOffset, xCount, yCount)

        if 'sm_out' not in locals():
            sm_out = np.ma.zeros((this_sm.shape[0], this_sm.shape[1], (search_days * 2) + 1), dtype = np.float32)

        # Both datasets have -9999. as a nodata value
        sm_out[:,:,n] = np.ma.array(this_sm, mask = this_sm == -9999.)

        # Close dataset
        ds = None

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
                interp = f(x)[int(round((x[this_sm.mask == False]).shape[0] / 2))]
            else:
                interp = -9999.

            sm_interp[i, j] = interp

    return sm_interp

def _resampleSM(sm_interp, tile, geo_t, interpolation = 'avearge'):
    '''
    '''

    # Create output file matching ALOS tile
    gdal_driver = gdal.GetDriverByName('MEM')

    ds_source = gdal_driver.Create('', sm_interp.shape[0], sm_interp.shape[1], 1, gdal.GDT_Float32)
    ds_source.SetGeoTransform(geo_t)
    ds_source.SetProjection(tile.proj)
    ds_source.GetRasterBand(1).SetNoDataValue(-9999)
    ds_source.GetRasterBand(1).WriteArray(sm_interp)

    ds_dest = gdal_driver.Create('', tile.ySize, tile.xSize, 1, gdal.GDT_Float32)
    ds_dest.SetGeoTransform(tile.geo_t)
    ds_dest.SetProjection(tile.proj)

    # Select interpolation type
    if interpolation == 'average' or interpolation == 'nearest':
        interp_type = gdal.GRA_NearestNeighbour
    elif interpolation == 'cubic':
        interp_type = gdal.GRA_CubicSpline

    # Reproject input GeoTiff to match the ALOS tile
    gdal.ReprojectImage(ds_source, ds_dest, tile.proj, tile.proj, interp_type)

    # Load resampled image into memory
    sm_resampled = ds_dest.GetRasterBand(1).ReadAsArray()

    return np.ma.array(sm_resampled, mask = sm_resampled == -9999.)


def _loadSM(tile, date, search_days = 7, interpolation = 'average', data_source = "CCI"):
    '''
    Load and reproject soil moisture data
    '''

    # Load SM time stack (within 1 week of measurement)
    sm, geo_t = _readTimeStack(tile, date, search_days = search_days, data_source = data_source)

    # Interpolate over time gaps
    sm_interp = _interpolateTime(sm)

    # Reproject to match ALOS tile (nearest neighbor)
    sm_resampled = _resampleSM(sm_interp, tile, geo_t, interpolation = interpolation)

    return sm_resampled


def getSM(tile, search_days = 7, interpolation = 'average'):
    """
    Function to load a resampled soil moisture image from the ESA CCI soil moisture product. The function returns a masked array with estimated volumetric soil moisture content (m^2/m^2), averaged for each satellite overpass in a tile.

    Args:
        tile: A biota.LoadTile() object.

    Returns:
        A masked array of estimated soil moisture (m^2/m^2)
    """

    assert interpolation in ['average', 'nearest', 'cubic'], "Soil moisture interpolation type must be 'average', 'nearest', or 'cubic'."

    # Build an array of dates in string format, which is quicker to search
    dates = tile.getDate().astype(np.int)

    # Generate output array
    out = np.zeros_like(tile.mask, dtype = np.float32) + 999999.

    for date in np.unique(dates): #np.unique(tile.getDate()).astype(np.int):

        # Nodata (0 = 1st Jan 1970)
        if date == 0:
            continue

        # Determine source data type
        if len(glob.glob(tile.SM_dir + "/ESACCI-SOILMOISTURE-*.nc")) > 0:
            data_source = "CCI"
        elif len(glob.glob(tile.SM_dir + "/SMAP_L3*.h5")) > 0:
            data_source = "SMAP"
        else:
            assert False, "No data from CCI or SMAP soil moisture products found in sm_dir (%s)"%tile.SM_dir

        # Interpolate through time and draw out the central value, and zoom to scale of tile
        sm_resampled = _loadSM(tile, dt.datetime(1970,1,1).date() + dt.timedelta(date), search_days = search_days, interpolation = interpolation, data_source = data_source)

        # Only output soil moisture where > 25 % data exists
        if float((sm_resampled.mask[dates == date]).sum()) / (sm_resampled.mask[dates == date]).shape[0] <= 0.25:

            if interpolation == 'average':
                # Calculate the mean average of values in overpass
                out[dates == date] = np.ma.mean(sm_resampled[dates == date])
            elif interpolation == 'nearest' or interpolation == 'cubic':
                out[dates == date] = sm_resampled[dates == date]

    # Add the mask
    out = np.ma.array(out, mask = np.logical_or(tile.mask, out == 999999.))

    return out
