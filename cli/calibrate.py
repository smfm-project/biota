import argparse
import os
import sys

import biota
import biota.calibrate

import pdb

"""
This is the script that runs the command line to extract data for validation purposes.
"""

def main(dataloc, years, shp, plot_field, agb_field,
         verbose = False):
    '''
    Comment this meaningfully
    '''
    
    try:
        
        # Allow processing of one or multiple year
        if type(years) != list: years = [years]

        years = [int(year) for year in years]

        for year in years:
            if verbose: print('Doing year: %s'%str(year))
            data_dict = biota.calibrate.extractGamma0(dataloc, year, shp, plot_field, agb_field, 
                                                      buffer_size = 25., verbose = verbose)

        #if agb_field is not None:
        #    slope, intercept = biota.calibrate.fitLinearModel(data_dict)

    except KeyboardInterrupt:
        sys.exit(0)
    
    except Exception as e:
        print (e)
        sys.exit(0)


if __name__ == '__main__':

    # Set up command line parser
    parser = argparse.ArgumentParser(description = 'Calibrate a biomass-backscatter relationship with ALOS mosaic data and a shapefile.')

    parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')

    # Required arguments
    required.add_argument('-y', '--years', metavar = 'YEAR', type = str, nargs = '*', help = 'Year of ALOS mosaic to load.')
    required.add_argument('-p', '--plot_field', metavar = 'NAME', type = str, help = 'Shapefile field containing a unique plot ID.')
    required.add_argument('-d', '--dataloc', metavar = 'DIR', type = str, help = 'Directory containing ALOS mosaic tiles')
    required.add_argument('shapefile', metavar = 'SHP', type = str, nargs = 1, help = "A shapefile containing plot and AGB data.")

    # Optional arguments
    optional.add_argument('-a', '--agb_field', metavar = 'NAME', type = str, default = None, help = 'Shapefile field containing an estimate of AGB.')

    # Get arguments from command line
    args = parser.parse_args()

    try:
        main(args.dataloc, args.years, args.shapefile[0], args.plot_field, args.agb_field)
    
    except KeyboardInterrupt:
        sys.exit(0)





