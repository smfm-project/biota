#!/usr/bin/env python

import argparse
import datetime
import ftplib
import os
import subprocess
import tarfile

"""
This is a simple script to assist in the downloading of data from the ALOS mosaic product.
"""


def generateURL(lat, lon, year):
    """
    Generates a URL to fetch ALOS mosaic data from JAXA FTP server.
    """
    
    # Test that inputs are reasonable lats/lons/years
    assert lat % 5 == 0, "Latitudes must be a multiple of 5 degrees."
    assert lon % 5 == 0, "Longitudes must be a multiple of 5 degrees."
    assert (year >= 2007 and year <= 2010) or (year >= 2015 and year <= datetime.datetime.now().year), "Years must be in the range 2007 - 2010 and 2015 - present. Your year was %s."%str(year)
    
    # Get hemisphere
    hem_NS = 'S' if lat < 0 else 'N'
    hem_EW = 'W' if lon < 0 else 'E'
    
    # Filename patterns are different for ALOS-1/ALOS-2
    if year <= 2010:
        url = 'ftp://ftp.eorc.jaxa.jp/pub/ALOS/ext1/PALSAR_MSC/25m_MSC/%s/%s%s%s%s_%s_MOS.tar.gz'
    else:
        url = 'ftp://ftp.eorc.jaxa.jp/pub/ALOS-2/ext1/PALSAR-2_MSC/25m_MSC/%s/%s%s%s%s_%s_MOS_F02DAR.tar.gz'
    
    url = url%(str(year), hem_NS, str(abs(lat)).zfill(2), hem_EW, str(abs(lon)).zfill(3), str(year)[-2:])
    
    return url


def download(url, output_dir = os.getcwd()):
    """
    Download data from JAXA FTP server
    """
    
    subprocess.call(['wget', url, '-P', output_dir])
    
    return url.split('/')[-1]


def decompress(filename, dataloc = os.getcwd(), remove = False):
    '''
    Unzips .tar.gz ALOS mosaic files downloaded from JAXA, and removes original where requested.
    '''
    
    # Get a list of zip files matching the Level 1C file pattern
    targz_file = dataloc + filename
    
    assert targz_file.endswith("tar.gz"), "File name must end with .tar.gz to be decompressed."
    
    tar = tarfile.open(targz_file, "r:gz")
    tar.extractall()
    tar.close()
    
    if remove: removeTarGz(targz_file)


def removeTarGz(targz_file):
    """
    Deletes ALOS-1/ALOS-2 .tar.gx files from disk.
    Input is a compress ALOS-1/ALOS-2 file from JAXA.
    """
    
    assert targz_file.endswith('_MOS.tar.gz') or targz_file.endswith('_MOS_F02DAR.tar.gz'), "removeTarGz function should only be used to delete ALOS-1/ALOS2 .tar.gz files"
    
    os.remove(targx_file)
    

def main(lat, lon, year, output_dir = os.getcwd(), remove = False):
    '''
    Run through data download and preparation chain
    '''
    
    # Generate download URL
    url = generateURL(lat, lon, year)
    
    # Get file, sending it to output_dir
    filename = download(url, output_dir = output_dir)
    
    # Decompress downloaded file, and delete original if remove == True
    decompress(url.split('/')[-1], remove = remove)

    
def getYears(years):
    """
    Reduces input years to those available for the ALOS mosaic (or may be available in future).
    """
    
    years_cleansed = []
    
    for y in years:
        
        # Remove years before ALOS-1
        if y < 2007:   
            print 'WARNING: Not downloading data for %s; no data from ALOS are available before 2007'%str(y)
            
        # Remove years between ALOS-1 and ALOS-2
        elif y > 2010 and y < 2015:
            print 'WARNING: Not downloading data for %s; no data from ALOS are available 2011 to 2014 inclusive.'%str(y)
        
        # Remove years from the future, which can't possibly exist yet.
        elif y > datetime.datetime.now().year:
            print 'WARNING: Not downloading data for %s; this year is in the future.'%str(y)
        
        # We will attempt to download all remaining years
        else:
            years_cleansed.append(y)
        
    return sorted(years_cleansed)


if __name__ == '__main__':

    # Set up command line parser
    parser = argparse.ArgumentParser(description = 'Download ALOS-1/2 data from JAXA, specifying a particular year and latitude/longitude.')

    # Required arguments
    parser.add_argument('-lat', '--latitude', type = int, help = "Latitude of tile lower-left corner. Must be a multiple of 5 degrees.")
    parser.add_argument('-lon', '--longitude', type = int, help = "Longitude of tile lower-left corner. Must be a multiple of 5 degrees.")
    
    # Optional arguments
    parser.add_argument('-y', '--years', type = int, nargs = '+', default = [2007,2010,2016], help = "Year of data to download. Defaults to downloading all data.")
    parser.add_argument('-o', '--output_dir', type = str, default = os.getcwd(), help = "Optionally specify an output directory. Defaults to the present working directory.")
    parser.add_argument('-r', '--remove', action='store_true', default = False, help = "Optionally remove downloaded .zip files after decompression.")

    # Get arguments from command line
    args = parser.parse_args()
    
    # Cleanse input years
    years = getYears(args.years)
    
    # Run through entire processing sequence
    for year in years:
        main(args.latitude, args.longitude, year, output_dir = args.output_dir, remove = args.remove)