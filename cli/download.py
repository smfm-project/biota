import argparse
import os
import sys

import biota.download

import pdb

"""
This is the script that runs the command line for download.py
"""



def main(lat, lon, years, 
         large_tile = False,
		 output_dir = os.getcwd(),
		 remove = False):
    '''
    Run through data download and preparation chain
    '''

    # Allow single year input or list
    if type(years) != list: years = [years]

    # Cleanse input years
    biota.download.checkYears(years)

    for year in years:

        # Download file, provided it exists, else continue. Exit with KeyboardInterrupt.
        try:
            filepath = biota.download.download(lat, lon, year, large_tile = large_tile, output_dir = output_dir, verbose = True)
        except KeyboardInterrupt:
            sys.exit(0)
        except Exception as e:
            print(e)
            continue

        # Decompress downloaded file, and delete original if remove == True
        biota.download.decompress(filepath, remove = remove)



if __name__ == '__main__':

    # Set up command line parser
    parser = argparse.ArgumentParser(description = 'Download ALOS-1/2 data from JAXA, specifying a particular year and latitude/longitude.')

    parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')

    # Required arguments
    required.add_argument('-lat', '--latitude', metavar = 'DEG', type = int, help = "Latitude of tile upper-left corner.")
    required.add_argument('-lon', '--longitude', metavar = 'DEG', type = int, help = "Longitude of tile upper-left corner.")

    # Optional arguments
    optional.add_argument('-l', '--large', action = 'store_true', default = False, help = "Download large tiles. ALOS mosaic tiles are available in 1x1 or 5x5 degree tiles. If downloading large volumes of data, it's usually better to use the latter. If this option is chosen, you must select a lat and lon that's a multiple of 5 degrees.")
    optional.add_argument('-y', '--years', metavar = 'Y', type = int, nargs = '+', default = [2007, 2008, 2009, 2010, 2015, 2016], help = "Year of data to download. Defaults to downloading all data.")
    optional.add_argument('-o', '--output_dir',metavar = 'DIR',  type = str, default = os.getcwd(), help = "Optionally specify an output directory. Defaults to the present working directory.")
    optional.add_argument('-r', '--remove', action='store_true', default = False, help = "Optionally remove downloaded .tar.gz files after decompression.")

    # Get arguments from command line
    args = parser.parse_args()

    # Run through entire processing sequence
    try:
        main(args.latitude, args.longitude, args.years, large_tile = args.large, output_dir = args.output_dir, remove = args.remove)
    except KeyboardInterrupt:
        sys.exit(0)
