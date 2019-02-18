import argparse
import datetime
import ftplib
import math
import os
import subprocess
import sys
import tarfile
import tqdm

import pdb

import biota3

"""
This is the script that runs the command line to extract and plot properties for one year.
"""




def main(dir, lat, lon, years,
        output = 'all',
        speckle = False,
        downsample = 1,
        polarisation = 'HV',
        forest_threshold = 10,
        area_threshold = 0,
        output_dir = os.getcwd()):
    '''
    Comment this meaningfully
    '''

    # Allow single year input or list
    if type(years) != list: years = [years]


    for year in years:

        # Process  tile, provided it exists, else continue. Exit with KeyboardInterrupt.
        try:
            # Load the tile
            tile = biota3.LoadTile(dir, lat, lon, year, lee_filter = speckle, downsample_factor = downsample, forest_threshold = forest_threshold, area_threshold = area_threshold, output_dir = output_dir)

            # Here come the choices
            if output == 'Gamma0' or output == 'all':
                print ("Calculating Gamma0")
                G0_tile = tile.getGamma0(polarisation = polarisation, units = 'decibels', output = True)
            if output == 'AGB' or output == 'all':
                print ("Calculating Above-Ground Biomass")
                AGB_tile =tile.getAGB(output = True)
            if output == 'WoodyCover' or output == 'all':
                print ("Calculating Woody Cover")
                WC_tile =tile.getWoodyCover(output = True)




        except KeyboardInterrupt:
            sys.exit(0)
        except Exception as e:
            print(e)
            continue




if __name__ == '__main__':

    # Set up command line parser
    parser = argparse.ArgumentParser(description = "Process downloaded ALOS-1/2 data to output biomass and forest cover.")

    parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')
    specific = parser.add_argument_group('Output-specific arguments')

    # Required arguments
    required.add_argument('-dir', '--data_directory', metavar = 'DIR', type = str, help = "absolute path to data directory")
    required.add_argument('-lat', '--latitude', metavar = 'DEG', type = int, help = "Latitude of tile upper-left corner.")
    required.add_argument('-lon', '--longitude', metavar = 'DEG', type = int, help = "Longitude of tile upper-left corner.")
    required.add_argument('-y', '--years', metavar = 'Y', type = int, nargs = '+', help = "Years of data to process.")

    # Optional arguments
    optional.add_argument('-o', '--output', choices = ['Gamma0', 'AGB', 'WoodyCover', 'all'], default = 'all', help = "Choose which kind of output you want. Defaults to all possible outputs.")
    optional.add_argument('-lf', '--speckle', action = 'store_false', default = True, help = "Apply speckle filtering. Defaults to True.")
    optional.add_argument('-ds', '--downsample', metavar = 'FACTOR', action = 'store', type = int, default = 1, help = "Apply downsampling. Defaults to 1.")
    optional.add_argument('-od', '--output_dir',metavar = 'DIR',  type = str, default = os.getcwd(), help = "Optionally specify an output directory. Defaults to the present working directory.")

    # Arguments specific to a type of output
    specific.add_argument('-pz', '--polarisation', choices = ['HV', 'HH', 'VH', 'VV'], default = 'HV', help = "If you have selected Gamma0 as an output, choose the polarisation. Defaults to HV.")
    specific.add_argument('-ft', '--forest_threshold', metavar = 'THRESHOLD', action = 'store', type = float, default = 10, help = "If you have selected WoodyCover as an output, choose the miminum forest biomass threshold. Defaults to 10tC/ha.")
    specific.add_argument('-at', '--area_threshold', metavar = 'THRESHOLD', action = 'store', type = float, default = 0, help = "If you have selected WoodyCover as an output, choose the minimum forest area threshold. Defaults to 0ha.")



    # Get arguments from command line
    args = parser.parse_args()

    # Run through entire processing sequence
    try:
        main(args.data_directory, args.latitude, args.longitude, args.years,
        output = args.output,
        speckle = args.speckle,
        downsample = args.downsample,
        polarisation = args.polarisation,
        forest_threshold = args.forest_threshold,
        area_threshold = args.area_threshold,
        output_dir = args.output_dir)
    except KeyboardInterrupt:
        sys.exit(0)
