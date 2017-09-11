#!/usr/bin/env python

import datetime
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
    assert lat % 5 == 0, "Latitudes must be a multiple of 5°."
    assert lon % 5 == 0, "Longitudes must be a multiple of 5°."
    assert (year > 2007 and year < 2010) or (year > 2015 and year < datetime.datetime.now().year), "Years must be in the range 2007 - 2010 and 2015 - present. Your year was %s."%str(year)
    
    # Get hemisphere
    hem_NS = 'S' if lat < 0 else 'N'
    hem_EW = 'W' if lon < 0 else 'E'
    
    # Filename patterns are different for ALOS-1/ALOS-2
    if year < 2010:
        url = 'ftp://ftp.eorc.jaxa.jp/pub/ALOS/ext1/PALSAR_MSC/25m_MSC/%s/%s%s%s%s_%s_MOS.tar.gz'
    else:
        url = 'ftp://ftp.eorc.jaxa.jp/pub/ALOS-2/ext1/PALSAR-2_MSC/25m_MSC/%s/%s%s%s%s_%s_MOS_F02DAR.tar.gz'
    
    url = url%(str(year), hem_NS, str(abs(lat)).zfill(2), hem_EW, str(abs(lon)).zfill(3), str(year)[-2:])
    
    return url


def download(url, output_dir = os.getcwd()):
    """
    Download data from JAXA FTP server
    """
    
    subprocess.call(['wget', command, '-P', output_dir])
    
    return command.split('/')[-1]


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

    
    
def extractRange(inputString):
    """
    Function to separate argparse inputs, where two values may be separated  by :. Returns a list of integers.
    """
    
    if ':' in inputString:
        out = inputString.split(':')
        assert len(out)==2, "Input ranges for lat/lon/year should be specified in the format <min>:<max>."
        out = [int(i) for i in out]
    else:
        out = [int(inputString)]
    
    return out


def getDegrees(degrees):
    """
    ALOS mosaic tiles are distributed in 5x5 degree tiles. Any latitude or longitude not divisible by 5 should be removed.
    """
    
    degrees = [d for d in degrees if d % 5 == 0]
    
    return degrees


def getYears(years):
    """
    Reduces input years to those available for the ALOS mosaic (or may be available in future).
    """
    
    # Remove years before ALOS-1
    years = [y for y in years if y < 2007]
    
    # Remove years between ALOS-1 and ALOS-2
    years = [y for y in years if y > 2010 and y < 2015]
    
    # Remove years from the future, which can't possibly exist yet.
    years = [y for y in years if y > datetime.datetime.now().year]
            
    return years


if __name__ == '__main__':

    # Set up command line parser
    parser = argparse.ArgumentParser(description = 'Download ALOS-1/2 data from JAXA, specifying a particular year and latitude/longitude.')

    # Required arguments
    parser.add_argument('-lat', '--latitude', type = str, help = "Latitude of tile lower-left corner. Must be a multiple of 5°.")
    parser.add_argument('-lon', '--longitude', type = str, help = "Longitude of tile lower-left corner. Must be a multiple of 5°.")
    parser.add_argument('-y', '--year', type = str, help = "Year of data to download.")

    # Optional arguments
    parser.add_argument('-o', '--output_dir', type = str, default = os.getcwd(), help = "Optionally specify an output directory. Defaults to the present working directory.")
    parser.add_argument('-r', '--remove', action='store_true', default = False, help = "Optionally remove downloaded .zip files after decompression.")

    # Get arguments from command line
    args = parser.parse_args()
    
    # Extract ranges of lat/lon/year.
    lats = extractRange(args.lat)
    lons = extractRange(args.lon)
    years = extractRange(args.year)
    
    # Remove years that don't exist and lats/lons not divisible by 5
    lats = getDegrees(lats)
    lons = getDegrees(lons)
    years = getYears(lons)
    
    # Test that inputs have resulted in reasonable lats/lons/years
    assert len(lats) > 0, "Latitudes must be a multiple of 5°, or encompass a range containing values that are a multiple of 5°."
    assert len(lons) > 0, "Longitudes must be a multiple of 5°, or encompass a range containing values that are a multiple of 5°."
    assert len(years) > 0, "Years must be in the range 2007 - 2010 and 2015 - present, or encompass a range containing years within these ranges."

    # Run through entire processing sequence
    for year in years:
        for lat in lats:
            for lon in lons:
                main(lat, lon, year, output = args.output_dir, remove = args.remove)