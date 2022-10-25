#!/usr/bin/env python

import argparse
import datetime
import ftplib
import math
import os
import subprocess
import sys
import tarfile
import tqdm
from zipfile import ZipFile

import pdb

"""
This is a script to assist in the downloading of data from the ALOS mosaic product.
"""



def checkYears(years):
    """
    Checks input years, and throws an error for those available for the ALOS mosaic (or may be available in future).

    Args:
        years: A list of years
    """

    for year in years:

        # Check years before ALOS-1
        if year < 2007:
            raise ValueError("Can't download data for year %s; no data from ALOS are available before 2007."%str(year))

        # Check years between ALOS-1 and ALOS-2
        elif year > 2010 and year < 2015:
            raise ValueError("Can't download data for year %s; no data from ALOS are available 2011 to 2014 (inclusive)."%str(year))

        # Check years from the future, which can't possibly exist yet.
        elif year > datetime.datetime.now().year:
            raise ValueError("Can't download data for %s; this year is in the future."%str(year))



def generateURL(lat, lon, year, large_tile = False):
    """
    Generates a URL to fetch ALOS mosaic data from JAXA FTP server.

    Args:
        lat: Latitude
        lon: Longitude
        year: Year
        large_tile: Set True to return the URL for a large (5x5 tile). Defaults to False (1x1)
    Returns:
        A URL
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
            url = 'ftp://ftp.eorc.jaxa.jp/pub/ALOS-2/ext1/PALSAR-2_MSC/25m_MSC/%s/%s_%s_MOS_F02DAR.zip'

        url = url%(str(year), tile_name, str(year)[-2:])

    # Get URL for a 1x1 tile
    else:

        # Filename patterns are different for ALOS-1/ALOS-2
        if year <= 2010:
            url = 'ftp://ftp.eorc.jaxa.jp/pub/ALOS/ext1/PALSAR_MSC/25m_MSC/%s/%s/%s_%s_MOS.tar.gz'
        else:
            url = 'ftp://ftp.eorc.jaxa.jp/pub/ALOS-2/ext1/PALSAR-2_MSC/25m_MSC/%s/%s/%s_%s_MOS_F02DAR.zip'

        # Special consideration for directory name of 1x1 tile
        lat_dir = int(5 * math.ceil(float(lat)/5))
        lon_dir = int(5 * math.floor(float(lon)/5))
        hem_NS_dir = 'S' if lat_dir < 0 else 'N'
        hem_EW_dir = 'W' if lon_dir < 0 else 'E'

        directory_name = '%s%s%s%s'%(hem_NS_dir, str(abs(lat_dir)).zfill(2), hem_EW_dir, str(abs(lon_dir)).zfill(3))

        url = url%(str(year), directory_name, tile_name, str(year)[-2:])

    return url


def download(lat, lon, year, large_tile = False, output_dir = os.getcwd(), verbose = False):
    """
    Download data from JAXA FTP server.

    Args:
        lat: Latitude
        lon: Longitude
        year: Year
        large_tile: Set True to return the URL for a large (5x5 tile). Defaults to False (1x1)
        output_dir: Output data to a specified directory. Defaults to current working directory.
        verbose: Set True to print progress.
    Returns:
        A location of the output file
    """

    if verbose: print('Downloading %stile for year: %s, lat: %s, lon: %s'%('large ' if large_tile else '', str(year), str(lat), str(lon)))

    # Generate download URL
    url = generateURL(lat, lon, year, large_tile = large_tile)
    
    # Check that output directory exists
    output_dir = os.path.abspath(os.path.expanduser(output_dir))
    assert os.path.isdir(output_dir), "The output directory (%s) does not exist. Create it, then try again."%output_dir

    output_file = '%s/%s'%(output_dir,url.split('/')[-1])

    # Check that output file doesn't already exist
    if os.path.exists(output_file) or os.path.exists(output_file[:-7]):
        raise ValueError("File %s already exists. Skipping."%output_file)

    # Login to FTP
    try:
        ftp = ftplib.FTP('ftp.eorc.jaxa.jp')
        login = ftp.login()
        cmd = ftp.cwd('/'.join(url.split('/')[3:-1]))
    except ftplib.all_errors as e:
        errorcode_string = str(e).split(None, 1)[0]
        raise OSError('Failed to connect to JAXA FTP server with error code %s.'%errorcode_string)

    # Test if file exists on remote server
    if url.split('/')[-1] not in ftp.nlst():
        ftp.quit()
        raise ValueError('Download failed: a tile for that location was not found on the FTP server.')

    try:
        with open(output_file, 'wb') as f:
            ftp.sendcmd("TYPE i")
            total = ftp.size(url.split('/')[-1])

            with tqdm.tqdm(total = total, unit = 'B', unit_scale = True, disable = not verbose) as pbar:
                def cb(data):
                    pbar.update(len(data))
                    f.write(data)
                ftp.retrbinary('RETR %s'%url.split('/')[-1], cb)

    except (Exception, KeyboardInterrupt) as e:
        #Tidy up in case of interrupted file transfer
        os.remove(output_file)
        ftp.close()
        raise

    # Logout politely
    exit = ftp.quit()

    return output_file
    

def decompress_tar(targz_file):
    '''
    Unzips a .tar.gz file.
    '''

    with tarfile.open(targz_file, "r:gz") as tar:
        def is_within_directory(directory, target):
            
            abs_directory = os.path.abspath(directory)
            abs_target = os.path.abspath(target)
        
            prefix = os.path.commonprefix([abs_directory, abs_target])
            
            return prefix == abs_directory
        
        def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
        
            for member in tar.getmembers():
                member_path = os.path.join(path, member.name)
                if not is_within_directory(path, member_path):
                    raise Exception("Attempted Path Traversal in Tar File")
        
            tar.extractall(path, members, numeric_owner=numeric_owner) 
            
        
        safe_extract(tar, path=targz_file[:-"7"])

def decompress_zip(zip_file):
    '''
    Unzips a .zip file.
    '''

    with ZipFile(zip_file, 'r') as zip:
        zip.extractall(path = zip_file.split('.')[0])

def decompress(file_path, remove=False):
    """wraps around decompress_tar and decompress_zip to decompress yearly ALOS and ALOS-2 mosaics, and removes the compressed files if remove==True"""
    
    assert os.path.exists(file_path), "File %s not found for decompression"%file_path

    # Check that output file doesn't already exist
    if os.path.exists(file_path.split('.')[0]):
        raise ValueError('File %s already exists at output location. Not extracting.'%file_path)

    # select right decompress function
    if file_path.endswith('.zip'):
        print('Extracting %s'%file_path)
        decompress_zip(file_path)

    elif file_path.endswith('.tar.gz'): 
        print('Extracting %s'%file_path)
        decompress_tar(file_path)

    else: 
        raise ValueError("File name must end with .zip or .tar.gz to be decompressed.")

    if remove:
        assert file_path.endswith('_MOS_F02DAR.tar.gz') or file_path.endswith('_MOS.tar.gz') or file_path.endswith('_MOS_F02DAR.zip'), "remove function should only be used to delete ALOS-1/ALOS2 .tar.gz or .zip files"
        os.remove(file_path)
