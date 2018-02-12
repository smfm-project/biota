# Toolset to calibrate change attribution algorithm

import math
import numpy as np
import biota.calibration as cal
import scipy.stats as stats
import scipy.ndimage
from scipy.ndimage.measurements import label

import matplotlib.pyplot as plt
import pdb




def getChange(AGB_t1, AGB_t2, F_threshold = 10., C_threshold = 0.2, min_pixels = 1):
    """
    Returns pixels that meet change detection thresholds for country.
    
    Args:
        F_threshold = threshold above which a pixel is forest
        C_threshold = threshold of proportional change with which to begin change detection
    """
       
    # Pixels that move from forest to nonforest or vice versa
    F_NF = np.logical_and(AGB_t1 >= F_threshold, AGB_t2 < F_threshold)
    NF_F = np.logical_and(AGB_t1 < F_threshold, AGB_t2 >= F_threshold)
    F_F = np.logical_and(AGB_t1 >= F_threshold, AGB_t2 >= F_threshold)
    NF_NF = np.logical_and(AGB_t1 < F_threshold, AGB_t2 < F_threshold)
    
    # Change pixels
    CHANGE = np.logical_or(((AGB_t1 - AGB_t2) / AGB_t1) >= C_threshold, ((AGB_t1 - AGB_t2) / AGB_t1) < (- C_threshold))
    NOCHANGE = CHANGE == False
    
    # Trajectory
    DECREASE = AGB_t2 < AGB_t1
    INCREASE = AGB_t2 >= AGB_t1
    
    # Get a minimum pixel extent
    if min_pixels != 1:
        CHANGE_INCREASE, _ = cal.getContiguousAreas(CHANGE & INCREASE, True, min_pixels = min_pixels)
        CHANGE_DECREASE, _ = cal.getContiguousAreas(CHANGE & DECREASE, True, min_pixels = min_pixels)
        CHANGE = np.logical_or(CHANGE_INCREASE, CHANGE_DECREASE)
    
    metrics = {}
    
    # These are all the possible definitions. Note: they *must* add up to one.
    metrics['deforestation'] = (F_NF & CHANGE) * 1
    metrics['degradation'] = (F_F & CHANGE & DECREASE) * 1
    metrics['minorloss'] = ((F_F | F_NF) & NOCHANGE & DECREASE) * 1

    metrics['aforestation'] = (NF_F & CHANGE) * 1
    metrics['growth'] = (F_F & CHANGE & INCREASE) * 1
    metrics['minorgain'] = ((F_F | NF_F) & NOCHANGE & INCREASE) * 1
            
    metrics['notforest'] = (NF_NF * 1)

    return metrics



class CHANGE(object):
    """
    A tree cover change object, based on ALOS tiles.
    Change objects have the following properties:

    Attributes:
        lat:
        lon:
        DN:
        mask:
    """
        
    def __init__(self, dataloc, lat, lon, year1, year2):
        """
        Loads data and metadata for an ALOS mosaic tile.
        """
        
        # Test that inputs are of reasonable lats/lons/years
        assert type(lat) == int, "Latitude must be an integer."
        assert lat < 90. or lat > -90., "Latitude must be between -90 and 90 degrees."
        assert type(lon) == int, "Longitude must be an integer."
        assert lon < 180. or lon > -180., "Longitude must be between -180 and 180 degrees."
        assert type(year1) == int, "year1 must be an integer."
        assert (year1 >= 2007 and year2 <= 2010) or (year1 >= 2015 and year1 <= dt.datetime.now().year1), "year1 must be in the range 2007 - 2010 and 2015 - present. Your input year1 was %s."%str(year1)
        assert type(year2) == int, "year2 must be an integer."
        assert (year2 >= 2007 and year2 <= 2010) or (year2 >= 2015 and year2 <= dt.datetime.now().year2), "year2 must be in the range 2007 - 2010 and 2015 - present. Your input year2 was %s."%str(year2)
        
        # Set up the two ALOS tiles to compare
        self.data_t1 = cal.ALOS(dataloc, lat, lon, year1, factor = 3)
        self.data_t2 = cal.ALOS(dataloc, lat, lon, year2, factor = 3)
        
    def getChange(self, F_threshold = 10., C_threshold = 0.2, min_pixels = 1):
        '''
        '''
        
        AGB_t1 = self.data_t1.getAGB(lee_filter = True)
        AGB_t2 = self.data_t2.getAGB(lee_filter = True)
        
        metrics = getChange(AGB_t1, AGB_t2, F_threshold = F_threshold, C_threshold = C_threshold, min_pixels = min_pixels)
        
        return metrics
        




if __name__ == '__main__':
    '''
    '''

    data_dir = '/home/sbowers3/DATA/ALOS_data/ALOS_mosaic/gorongosa/'
    output_dir = '/home/sbowers3/DATA/ALOS_data/ALOS_mosaic/gorongosa/'

    #Nhambita = lat -18 lon 33

    lat = -18#-19
    lon = 33#34#39
    
    data_2007 = cal.ALOS(data_dir, lat, lon, 2007)
    #data_2008 = cal.ALOS(lat, lon, 2008, data_dir)
    #data_2009 = cal.ALOS(lat, lon, 2009, data_dir)
    data_2010 = cal.ALOS(data_dir, lat, lon, 2010 )
    
    change_07_10 = CHANGE(data_dir, lat, lon, 2007, 2010)
    metrics = change_07_10.getChange(F_threshold = 15, C_threshold = 0.2, min_pixels = 6)
    
    plt.imshow(metrics['deforestation'] * 2 + metrics['degradation'])
    plt.colorbar()
    plt.show()

    #change_07_10 = getChange(AGB_2007, AGB_2010, F_threshold = 10., C_threshold = 0.2)
    #contiguous_change_area = getContiguousAreas(change_07_10['degradation'] + change_07_10['deforestation'], 1, 8)

    #deforestation = change_07_10['deforestation'] * contiguous_change_area
    #degradation = change_07_10['degradation'] * contiguous_change_area

    

    #deforestation_07_08 = getContiguousAreas(change_07_08['deforestation'], 1, 8)
    #deforestation_08_09 = getContiguousAreas(change_08_09['deforestation'], 1, 8)
    #deforestation_09_10 = getContiguousAreas(change_09_10['deforestation'], 1, 8)

    #deforestation = deforestation_07_10.sum() / (4500 * 4500.)
    #total_deforestation = (deforestation_07_08.sum() + deforestation_08_09.sum() + deforestation_09_10.sum()) / (4500. * 4500)


    # This result looks nice. TODO: calculate % change of entire contiguous change region
    """
    lakes = '/home/sbowers3/DATA/GIS_data/mozambique/diva/MOZ_wat/MOZ_water_areas_dcw.shp'
    rivers = '/home/sbowers3/DATA/GIS_data/mozambique/diva/MOZ_wat/MOZ_water_lines_dcw.shp'
    lake_mask = cal.rasterizeShapefile(data_2007, lakes, buffer_size = 0.005)
    river_mask = cal.rasterizeShapefile(data_2007, rivers, buffer_size = 0.005)
    water_mask = np.logical_or(river_mask, lake_mask)
    
    output_im = np.ma.array(deforestation*2+degradation,mask=water_mask)
    output_im.data[output_im.mask] = 0
    
    from osgeo import osr, gdal
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)

    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(output_dir + 'forest_loss.tif', 4500, 4500, 1, gdal.GDT_UInt16)
    ds.SetGeoTransform(data_2007.geo_t)
    ds.SetProjection(srs.ExportToWkt())
    ds.GetRasterBand(1).WriteArray(output_im)
    ds.GetRasterBand(1).SetNoDataValue(0)
    ds = None
    """
    

    """
    geo_t = (33.0, 0.00022222222222222223*10, 0.0, -18.0, 0.0, -0.00022222222222222223*10)
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(output_dir + 'forest_loss.tif', 450, 450, 1, gdal.GDT_UInt16)
    ds.SetGeoTransform(geo_t)
    ds.SetProjection(srs.ExportToWkt())
    ds.GetRasterBand(1).WriteArray(deforestation*2 + degradation)
    ds = None
    """
    """

    contiguous_change_area_gain = getContiguousAreas(change_07_10['aforestation'] + change_07_10['growth'], 1, 8)

    aforestation = change_07_10['aforestation'] * contiguous_change_area_gain
    growth = change_07_10['growth'] * contiguous_change_area_gain

    # A deforestation rate from differencing?

    out = np.zeros((450,450))
    for x in range(450):
        for y in range(450):
            out[y,x] = aforestation[y*10:(y*10)+10,x*10:(x*10)+10].data.sum() - deforestation[y*10:(y*10)+10,x*10:(x*10)+10].data.sum()
            

    print "Deforestation rate: %s"%str(np.sum(out[out<0])/(4500*4500.))
    print "Aforestation rate: %s"%str(np.sum(out[out>0])/(4500*4500.))

    """
