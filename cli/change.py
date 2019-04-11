import argparse
import os
import sys

import pdb

import biota

"""
This is the script that runs the command line to extract and plot property variations between two years.
"""


def main(dir, lat, lon, year1, year2,
        output = 'all',
        lee_filter = False,
        downsample_factor = 1,
        forest_threshold = 10,
        area_threshold = 0,
        change_area_threshold = 0,
        change_magnitude_threshold = 0,
        change_intensity_threshold = 0,
        deforestation_threshold = None,
        output_dir = os.getcwd(),
		verbose = False):
    '''
    Comment this meaningfully
    '''
    
    # Process  tile, provided it exists, else continue. Exit with KeyboardInterrupt.
    try:
	    
        # Load the tiles
        if verbose: print ('Loading tile:', year1)
        tile1 = biota.LoadTile(dir, lat, lon, year1, lee_filter = lee_filter, downsample_factor = downsample_factor, forest_threshold = forest_threshold, area_threshold = area_threshold, output_dir = output_dir)
		
        if verbose: print ('Loading tile:', year2)
        tile2 = biota.LoadTile(dir, lat, lon, year2, lee_filter = lee_filter, downsample_factor = downsample_factor, forest_threshold = forest_threshold, area_threshold = area_threshold, output_dir = output_dir)
        
		# Load a change object
        tile_change = biota.LoadChange(tile1, tile2, change_area_threshold = change_area_threshold, change_magnitude_threshold = change_magnitude_threshold, change_intensity_threshold = change_intensity_threshold, deforestation_threshold = deforestation_threshold, output_dir = output_dir)

        # Here come the choices
        if output == 'AGBChange' or output == 'all':
            if verbose: print ("Calculating Biomass Change...")
            AGBChange = tile_change.getAGBChange(output = True)
        if output == 'ChangeType' or output == 'all':
            if verbose: print ("Calculating Change Type...")
            ChangeType = tile_change.getChangeType(output = True)
        if output == 'RiskMap' or output == 'all':
            if verbose: print ("Calculating Deforestation Risk Map...")
            RiskMap = tile_change.getRiskMap(output = True)
        if verbose: print ("Done!")
    
    except KeyboardInterrupt:
        sys.exit(0)
    
    except Exception as e:
        print (e)
        sys.exit(0)




if __name__ == '__main__':

    # Set up command line parser
    parser = argparse.ArgumentParser(description = "Process ALOS-1/2 moasic data to output biomass and woody cover change between 2 years.")

    parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')
    specific = parser.add_argument_group('Output-specific arguments')

    # Required arguments
    required.add_argument('-dir', '--data_directory', metavar = 'DIR', type = str, help = "Path to directory containing ALOS mosaic data.")
    required.add_argument('-lat', '--latitude', metavar = 'DEG', type = int, help = "Latitude of tile to process (upper-left corner).")
    required.add_argument('-lon', '--longitude', metavar = 'DEG', type = int, help = "Longitude of tile to process (upper-left corner).")
    required.add_argument('-y1', '--year1', metavar = 'YR', type = int, help = "First year of data to process.")
    required.add_argument('-y2', '--year2', metavar = 'YR', type = int, help = "Second year of data to process.")
    
    # Optional arguments
    optional.add_argument('-o', '--output', choices = ['AGBChange', 'ChangeType', 'RiskMap', 'all'], default = 'all', help = "Choose which kind of output you want. Defaults to all possible outputs.")
    optional.add_argument('-nf', '--nofilter', action = 'store_false', default = True, help = "Use this flag if you don't want to apply a speckle filter.")
    optional.add_argument('-ds', '--downsample_factor', metavar = 'N', action = 'store', type = int, default = 1, help = "Apply downsampling to inputs by specifying an integer factor to downsample by. Defaults to no downsampling.")
    optional.add_argument('-od', '--output_dir', metavar = 'DIR', type = str, default = os.getcwd(), help = "Optionally specify an output directory. Defaults to the present working directory.")
    optional.add_argument('-v', '--verbose', action = 'store_true', default = False, help = "Print progress to terminal. Defaults to False.")
    
    # Arguments specific to a type of output
    specific.add_argument('-ct', '--change_area_threshold', metavar = 'ha', action = 'store', type = float, default = 0, help = "If you have selected ChangeType as an output, choose a threshold for a minimum change in forest area required to be flagged as a change. Defaults to 0 ha.")
    specific.add_argument('-mt', '--change_magnitude_threshold', metavar = 'tC/ha', action = 'store', type = float, default = 0, help = "If you have selected ChangeType as an output, choose the minimum absolute change in biomass to be flagged as a change. Defaults to 0 tC/ha.")
    specific.add_argument('-it', '--change_intensity_threshold', metavar = 'PC', action = 'store', type = float, default = 0, help = "If you have selected ChangeType as an output, choose the minimum relative change in biomass to be flagged as a change. Defaults to 0 percent.")
    specific.add_argument('-dt', '--deforestation_threshold', metavar = 'tC/ha', action = 'store', type = float, default = None, help = "If you have selected ChangeType as an output, choose the maxmimum residual biomass accepted in a deforested area. Remaining 'deforestation' will be allocated to 'degradation'. Defaults to the forest_threshold.")
	
    specific.add_argument('-ft', '--forest_threshold', metavar = 'tC/ha', action = 'store', type = float, default = 10, help = "If you have selected ChangeType as an output, choose the miminum forest biomass threshold in each input image. Defaults to 10 tC/ha.")
    specific.add_argument('-at', '--area_threshold', metavar = 'ha', action = 'store', type = float, default = 0, help = "If you have selected ChangeType as an output, choose the minimum forest area threshold in each input image. Defaults to 0 ha.")

    # Get arguments from command line
    args = parser.parse_args()
	
	# Warn about units
    if args.change_intensity_threshold < 1 and args.change_intensity_threshold > 0: "WARNING: Change intensity threshold should be specified in units of %. Did you specify it as a proportion by accident?"
    change_intensity_threshold = args.change_intensity_threshold / 100.

    # Run through entire processing sequence
    try:
        main(args.data_directory, args.latitude, args.longitude, args.year1, args.year2,
        output = args.output,
        lee_filter = args.nofilter,
        downsample_factor = args.downsample_factor,
        forest_threshold = args.forest_threshold,
        area_threshold = args.area_threshold,
        change_area_threshold = args.change_area_threshold,
        change_magnitude_threshold = args.change_magnitude_threshold,
        change_intensity_threshold = change_intensity_threshold,
        deforestation_threshold = args.deforestation_threshold,
        output_dir = args.output_dir,
        verbose = args.verbose)
    
    except KeyboardInterrupt:
        sys.exit(0)
