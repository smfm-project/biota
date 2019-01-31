import argparse
import biota
import os


def main(data_dir, latitude, longitude, year, lee_filter = False, output_dir = os.getcwd()):
    '''
    '''
    
    # Print progress
    print 'Doing latitude: %s, longitude: %s'%(str(latitude), str(longitude))

    # Load the ALOS tile with specified options
    tile = biota.LoadTile(data_dir, latitude, longitude, year, lee_filter = lee_filter, output_dir = output_dir)
    
    # Calculate AGB and output to GeoTiff
    AGB = tile.getAGB(output = True)
    
    return AGB


if __name__ == '__main__':
    '''
    Command line utility to produce maps of AGB using biota
    '''
    
    # Set up command line parser
    parser = argparse.ArgumentParser(description = 'Build a map of aboveground biomass using ALOS-mosaic data.')

    parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')

    # Required arguments
    required.add_argument('-d', '--data_dir', type = str, metavar = 'PATH', required = True, help = "Directory containing ALOS mosaic data.")
    required.add_argument('-lat', '--latitude', type = int, metavar = 'LON', required = True, help = "Latitude to generate.")
    required.add_argument('-lon', '--longitude', type = int, metavar = 'LAT', required = True, help = "Longitude to generate.")
    required.add_argument('-y', '--year', type = int, metavar = 'Y', required = True, help = "Year to generate.")
    
    # Optional arguments
    optional.add_argument('-f', '--filter', action = 'store_true', default = False, help = "Apply lee filter to input data.")
    optional.add_argument('-o', '--output_dir', type = str, metavar = 'PATH', default = os.getcwd(), help = "Specify an output directory. Defaults to the present working directory.")
    
    # Get arguments from command line 
    args = parser.parse_args()
    
    # Run through entire processing sequence
    AGB = main(args.data_dir, args.latitude, args.longitude, args.year, lee_filter = args.filter, output_dir = args.output_dir)
