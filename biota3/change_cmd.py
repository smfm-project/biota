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
This is the script that runs the command line to extract and plot property variations between two years.
"""




def main(dir, lat, lon, year1, year2,
        output = 'all',
        speckle = False,
        downsample = 1,
        forest_threshold = 10,
        area_threshold = 0,
        change_area_threshold = 0,
        change_magnitude_threshold = 0,
        change_intensity_threshold = 0,
        output_dir = os.getcwd()):
    '''
    Comment this meaningfully
    '''

    # Allow single year input or list
    #if type(year1) != list: year1 = [year1]


    # Cleanse input years - Maybe not useful here
    #dld.checkYears(args.years)

    #for year in years:

    # Process  tile, provided it exists, else continue. Exit with KeyboardInterrupt.
    try:
        # Load the tiles
        year1 = int(year1[0])
        year2 = int(year2[0])

        print ('loading tile:', year1)
        tile1 = biota3.LoadTile(dir, lat, lon, year1, lee_filter = speckle, downsample_factor = downsample, forest_threshold = forest_threshold, area_threshold = area_threshold, output_dir = output_dir)
        print ('loading tile:', year2)
        tile2 = biota3.LoadTile(dir, lat, lon, year2, lee_filter = speckle, downsample_factor = downsample, forest_threshold = forest_threshold, area_threshold = area_threshold, output_dir = output_dir)

        tile_change = biota3.LoadChange(tile1, tile2, change_area_threshold = change_area_threshold, change_magnitude_threshold = change_magnitude_threshold, change_intensity_threshold = change_intensity_threshold)

        # Here come the choices
        if output == 'AGB' or output == 'all':
            print ("Calculating Biomass Change")
            AGB_tile =tile_change.getAGBChange(output = True)
        if output == 'ChangeType' or output == 'all':
            print ("Calculating Change Type")
            WC_tile =tile_change.getChangeType(output = True)




    except KeyboardInterrupt:
        sys.exit(0)
    except Exception as e:
        print(e)
        #continue




if __name__ == '__main__':

    # Set up command line parser
    parser = argparse.ArgumentParser(description = "Process downloaded ALOS-1/2 data to output biomass and forest cover change between 2 years.")

    parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')
    specific = parser.add_argument_group('Output-specific arguments')

    # Required arguments
    required.add_argument('-dir', '--data_directory', metavar = 'DIR', type = str, help = "absolute path to data directory")
    required.add_argument('-lat', '--latitude', metavar = 'DEG', type = int, help = "Latitude of tile upper-left corner.")
    required.add_argument('-lon', '--longitude', metavar = 'DEG', type = int, help = "Longitude of tile upper-left corner.")
    required.add_argument('-y1', '--year1', metavar = 'Y', type = int, nargs = '+', help = "First year of data to process.")
    required.add_argument('-y2', '--year2', metavar = 'Y', type = int, nargs = '+', help = "Second year of data to process.")

    # Optional arguments
    optional.add_argument('-o', '--output', choices = ['AGB', 'ChangeType', 'all'], default = 'all', help = "Choose which kind of output you want. Defaults to all possible outputs.")
    optional.add_argument('-lf', '--speckle', action = 'store_false', default = True, help = "Apply speckle filtering. Defaults to True.")
    optional.add_argument('-ds', '--downsample', metavar = 'FACTOR', action = 'store', type = int, default = 1, help = "Apply downsampling. Defaults to 1.")
    optional.add_argument('-od', '--output_dir',metavar = 'DIR',  type = str, default = os.getcwd(), help = "Optionally specify an output directory. Defaults to the present working directory.")

    # Arguments specific to a type of output
    specific.add_argument('-ft', '--forest_threshold', metavar = 'THRESHOLD', action = 'store', type = float, default = 10, help = "If you have selected WoodyCover as an output, choose the miminum forest biomass threshold. Defaults to 10tC/ha.")
    specific.add_argument('-at', '--area_threshold', metavar = 'THRESHOLD', action = 'store', type = float, default = 0, help = "If you have selected WoodyCover as an output, choose the minimum forest area threshold. Defaults to 0ha.")
    specific.add_argument('-cat', '--change_area_threshold', metavar = 'CAT', action = 'store', type = float, default = 0, help = "If you have selected ChangeType as an output, choose the minimum change in forest area threshold. Defaults to 0ha.")
    specific.add_argument('-cmt', '--change_magnitude_threshold', metavar = 'CMT', action = 'store', type = float, default = 0, help = "If you have selected ChangeType as an output, choose the minimum change in biomass threshold. Defaults to 0tC/ha.")
    specific.add_argument('-cit', '--change_intensity_threshold', metavar = 'CIT', action = 'store', type = float, default = 0, help = "If you have selected ChangeType as an output, choose the minimum relative change in forest biomass threshold. Defaults to 0.")

    # Get arguments from command line
    args = parser.parse_args()

    # Run through entire processing sequence
    try:
        main(args.data_directory, args.latitude, args.longitude, args.year1, args.year2,
        output = args.output,
        speckle = args.speckle,
        downsample = args.downsample,
        forest_threshold = args.forest_threshold,
        area_threshold = args.area_threshold,
        change_area_threshold = args.change_area_threshold,
        change_magnitude_threshold = args.change_magnitude_threshold,
        change_intensity_threshold = args.change_intensity_threshold,
        output_dir = args.output_dir)
    except KeyboardInterrupt:
        sys.exit(0)
