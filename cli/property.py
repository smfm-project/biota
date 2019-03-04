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

import biota
import biota.download

"""
This is the script that runs the command line to extract and plot properties for one year.
"""




def main(dir, lat, lon, years,
        output = 'all',
        lee_filter = False,
        downsample_factor = 1,
        polarisation = 'HV',
		units = 'natural',
        forest_threshold = 10,
        area_threshold = 0,
        output_dir = os.getcwd(),
		verbose = False):
    '''
    Load an ALOS tile with biota, and generate specified output
    '''

    # Allow single year input or list
    if type(years) != list: years = [years]
            
    for year in years:

        # Process  tile, provided it exists, else continue. Exit with KeyboardInterrupt.
        try:
            # Load the tile
            tile = biota.LoadTile(dir, lat, lon, year, lee_filter = lee_filter, downsample_factor = downsample_factor, forest_threshold = forest_threshold, area_threshold = area_threshold, output_dir = output_dir)

            # Here come the choices
            if output == 'Gamma0' or output == 'all':
                if verbose: print ("Calculating Gamma0...")
                gamma0 = tile.getGamma0(polarisation = polarisation, units = units, output = True)
            if output == 'AGB' or output == 'all':
                if verbose: print ("Calculating Above-Ground Biomass...")
                AGB = tile.getAGB(output = True)
            if output == 'WoodyCover' or output == 'all':
                if verbose: print ("Calculating Woody Cover...")
                WoodyCover = tile.getWoodyCover(output = True)
            
            if verbose: print ("Done!")

        except KeyboardInterrupt:
            sys.exit(0)
			
        except Exception as e:
            print(e)
            continue




if __name__ == '__main__':

    # Set up command line parser
    parser = argparse.ArgumentParser(description = "Process ALOS-1/2 mosaic data to prpoduce estimates of forest cover and biomass.")

    parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')
    specific = parser.add_argument_group('Output-specific arguments')

    # Required arguments
    required.add_argument('-dir', '--data_directory', metavar = 'DIR', type = str, help = "Path to directory containing ALOS mosaic data.")
    required.add_argument('-lat', '--latitude', metavar = 'DEG', type = int, help = "Latitude of tile to process (upper-left corner).")
    required.add_argument('-lon', '--longitude', metavar = 'DEG', type = int, help = "Longitude of tile to process (upper-left corner).")
    required.add_argument('-y', '--years', metavar = 'Y', type = int, nargs = '+', help = "Years of data to process.")

    # Optional arguments
    optional.add_argument('-o', '--output', choices = ['Gamma0', 'AGB', 'WoodyCover', 'all'], default = 'all', help = "Choose which kind of output you want. Defaults to all possible outputs.")
    optional.add_argument('-nf', '--nofilter', action = 'store_true', default = True, help = "Use this flag if you don't want to apply a speckle filter.")
    optional.add_argument('-ds', '--downsample_factor', metavar = 'N', action = 'store', type = int, default = 1, help = "Apply downsampling to inputs by specifying an integer factor to downsample by. Defaults to no downsampling.")
    optional.add_argument('-od', '--output_dir', metavar = 'DIR', type = str, default = os.getcwd(), help = "Optionally specify an output directory. Defaults to the present working directory.")
    optional.add_argument('-v', '--verbose', action = 'store_true', default = False, help = "Print progress to terminal. Defaults to False.")

    # Arguments specific to a type of output
    specific.add_argument('-p', '--polarisation', metavar = 'POL', choices = ['HV', 'HH'], default = 'HV', help = "If you have selected Gamma0 as an output, choose the polarisation. Defaults to HV.")
    specific.add_argument('-u', '--units', metavar = 'units', choices = ['natural', 'decibels'], default = 'natural', help = "If you have selected Gamma0 as an output, choose the outputs units. Defaults to 'natural' units.")    
    specific.add_argument('-ft', '--forest_threshold', metavar = 'tC/ha', action = 'store', type = float, default = 10, help = "If you have selected WoodyCover as an output, choose the miminum forest biomass threshold. Defaults to 10 tC/ha.")
    specific.add_argument('-at', '--area_threshold', metavar = 'ha', action = 'store', type = float, default = 0, help = "If you have selected WoodyCover as an output, choose the minimum forest area threshold. Defaults to 0 ha.")


    # Get arguments from command line
    args = parser.parse_args()

    # Run through entire processing sequence
    try:
        main(args.data_directory, args.latitude, args.longitude, args.years,
        output = args.output,
        lee_filter = args.nofilter,
        downsample_factor = args.downsample_factor,
        polarisation = args.polarisation,
		units = args.units,
        forest_threshold = args.forest_threshold,
        area_threshold = args.area_threshold,
        output_dir = args.output_dir,
		verbose = args.verbose)
    except KeyboardInterrupt:
        sys.exit(0)
    except Exception as e:
        print(e)
        sys.exit(0)
