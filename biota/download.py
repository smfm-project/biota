#!/usr/bin/env python

import argparse
import datetime
import ftplib
import math
import os
import subprocess
import tarfile

import pdb

"""
This is a simple script to assist in the downloading of data from the ALOS mosaic product.
"""


def generateURL(lat, lon, year, large_tile = False):
    """
    Generates a URL to fetch ALOS mosaic data from JAXA FTP server.
    """
    
    # Test that inputs are reasonable lats/lons
    if large_tile:
        assert lat % 5 == 0, "Latitudes must be a multiple of 5 degrees if downloading a 5x5 degree tile."
        assert lon % 5 == 0, "Longitudes must be a multiple of 5 degrees if downloading a 5x5 degree tile."
    
    # Test that years are reasonable
    assert (year >= 2007 and year <= 2010) or (year >= 2015 and year <= datetime.datetime.now().year), "Years must be in the range 2007 - 2010 and 2015 - present. Your year was %s."%str(year)
    
    # Test that lat/lon is reasonable
    assert (lat >= -90 and lat <= 90) and (lon >= -180 and lon <= 180), "Latitudes must be between -90 and 90 degrees, and lontitudes between -180 and 180 degrees."
    
    # Get hemisphere
    hem_NS = 'S' if lat < 0 else 'N'
    hem_EW = 'W' if lon < 0 else 'E'
    
    tile_name = '%s%s%s%s'%(hem_NS, str(abs(lat)).zfill(2), hem_EW, str(abs(lon)).zfill(3))

    # Get URL for a 5x5 tile
    if large_tile:
        
        # Filename patterns are different for ALOS-1/ALOS-2
        if year <= 2010:
            url = 'ftp://ftp.eorc.jaxa.jp/pub/ALOS/ext1/PALSAR_MSC/25m_MSC/%s/%s_%s_MOS.tar.gz'
        else:
            url = 'ftp://ftp.eorc.jaxa.jp/pub/ALOS-2/ext1/PALSAR-2_MSC/25m_MSC/%s/%s_%s_MOS_F02DAR.tar.gz'
                
        url = url%(str(year), tile_name, str(year)[-2:])

    # Get URL for a 1x1 tile
    else:
        
        # Filename patterns are different for ALOS-1/ALOS-2
        if year <= 2010:
            url = 'ftp://ftp.eorc.jaxa.jp/pub/ALOS/ext1/PALSAR_MSC/25m_MSC/%s/%s/%s_%s_MOS.tar.gz'
        else:
            url = 'ftp://ftp.eorc.jaxa.jp/pub/ALOS-2/ext1/PALSAR-2_MSC/25m_MSC/%s/%s/%s_%s_MOS_F02DAR.tar.gz'
        
        # Special consideration for directory name of 1x1 tile       
        lat_dir = int(5 * math.ceil(float(lat)/5))
        lon_dir = int(5 * math.floor(float(lon)/5))
        hem_NS_dir = 'S' if lat_dir < 0 else 'N'
        hem_EW_dir = 'W' if lon_dir < 0 else 'E'
        
        directory_name = '%s%s%s%s'%(hem_NS_dir, str(abs(lat_dir)).zfill(2), hem_EW_dir, str(abs(lon_dir)).zfill(3))
        
        url = url%(str(year), directory_name, tile_name, str(year)[-2:])
        
    return url


def download(url, output_dir = os.getcwd()):
    """
    Download data from JAXA FTP server
    """
    
    # Check that output directory exists
    assert os.path.isdir(output_dir), "The output directory (%s) does not exist. Create it, then try again."%output_dir
    
    # Check that output file doesn't already exist
    this_file = '%s/%s'%(output_dir,url.split('/')[-1])
    
    if os.path.exists(this_file):
        print "WARNING: File %s already exists. Skipping."%this_file
    elif os.path.exists(this_file[:-7]):
        print "WARNING: File %s already exists. Skipping."%this_file[:-7]

    else:
        # Download
        exit_status = subprocess.call(['wget', '-nc', url, '-P', output_dir])
        
        if exit_status == 8:
            raise ValueError('Download failed: a tile for that location was not found on the FTP server.')
        if exit_status != 0:
            raise ValueError('Download failed with wget error code %s.'%str(exit_status))
    
    # Determine absolute path of downloaded file
    filepath = '%s/%s'%(output_dir.rstrip('/'), url.split('/')[-1])
    
    return filepath


def decompress(targz_file, remove = False):
    '''
    Unzips .tar.gz ALOS mosaic files downloaded from JAXA, and removes original where requested.
    '''
    
    assert targz_file.endswith("tar.gz"), "File name must end with .tar.gz to be decompressed."
    
    # Check that output file doesn't already exist
    if os.path.exists(targz_file.split('/')[-1].split('.')[0]):
        print 'WARNING: File %s already exists at output location. Not extracting.'%targz_file
    
    else:
        print 'Extracting %s'%targz_file
            
        tar = tarfile.open(targz_file, "r:gz")
        tar.extractall(path = targz_file[:-7])
        tar.close()
        
        if remove: removeTarGz(targz_file)


def removeTarGz(targz_file):
    """
    Deletes ALOS-1/ALOS-2 .tar.gx files from disk.
    Input is a compress ALOS-1/ALOS-2 file from JAXA.
    """
    
    assert targz_file.endswith('_MOS.tar.gz') or targz_file.endswith('_MOS_F02DAR.tar.gz'), "removeTarGz function should only be used to delete ALOS-1/ALOS2 .tar.gz files"
    
    os.remove(targz_file)


def checkYears(years):
    """
    Reduces input years to those available for the ALOS mosaic (or may be available in future).
    """
    
    years_cleansed = []
    
    for y in years:
        
        # Remove years before ALOS-1
        if y < 2007:   
            raise ValueError("Can't download data for year %s; no data from ALOS are available before 2007."%str(y))
            
        # Remove years between ALOS-1 and ALOS-2
        elif y > 2010 and y < 2015:
            raise ValueError("Can't download data for year %s; no data from ALOS are available 2011 to 2014 (inclusive)."%str(y))
        
        # Remove years from the future, which can't possibly exist yet.
        elif y > datetime.datetime.now().year:
            raise ValueError("Can't download data for %s; this year is in the future."%str(y))
        


def main(lat, lon, year, large_tile = False, output_dir = os.getcwd(), remove = False):
    '''
    Run through data download and preparation chain
    '''
    
    # Generate download URL
    url = generateURL(lat, lon, year, large_tile = large_tile)
    
    # Get file, sending it to output_dir
    filepath = download(url, output_dir = output_dir)
    
    # Decompress downloaded file, and delete original if remove == True
    decompress(filepath, remove = remove)

    

if __name__ == '__main__':

    # Set up command line parser
    parser = argparse.ArgumentParser(description = 'Download ALOS-1/2 data from JAXA, specifying a particular year and latitude/longitude.')

    # Required arguments
    parser.add_argument('-lat', '--latitude', type = int, help = "Latitude of tile upper-left corner.")
    parser.add_argument('-lon', '--longitude', type = int, help = "Longitude of tile upper-left corner.")
    
    # Optional arguments
    parser.add_argument('-l', '--large', action = 'store_true', default = False, help = "Download large tiles. ALOS mosaic tiles are available in 1x1 or 5x5 degree tiles. If downloading large volumes of data, it's usually better to use the latter. If this option is chosen, you must select a lat and lon that's a multiple of 5 degrees.")
    parser.add_argument('-y', '--years', type = int, nargs = '+', default = [2007, 2008, 2009, 2010, 2015, 2016], help = "Year of datat o download. Defaults to downloading all data.")
    parser.add_argument('-o', '--output_dir', type = str, default = os.getcwd(), help = "Optionally specify an output directory. Defaults to the present working directory.")
    parser.add_argument('-r', '--remove', action='store_true', default = False, help = "Optionally remove downloaded .tar.gz files after decompression.")

    # Get arguments from command line
    args = parser.parse_args()
    
    # Cleanse input years
    checkYears(args.years)
    
    # Run through entire processing sequence
    for year in args.years:
        
        main(args.latitude, args.longitude, year, large_tile = args.large, output_dir = args.output_dir, remove = args.remove)
