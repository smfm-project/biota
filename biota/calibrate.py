#!/usr/bin/env python


# This is a set of scripts for calibration of a biomass-backscatter curve

import numpy as np

import biota
import biota.mask
import pdb

data_dir = '/home/sbowers3/guasha/sam_bowers/for_Geoff/alos_data/'
output_dir = '/home/sbowers3/guasha/sam_bowers/for_Geoff/'
output_filemname = 'data_out.csv'
shp = '/home/sbowers3/guasha/sam_bowers/for_Geoff/fwmxshpfile/scolel_te_plots_w_buffer_final.shp'

years = [2007, 2008, 2009, 2010, 2015, 2016]


### The code

plot_out = []
year_out = []
data_out = []

# Extract plots from shapefile
plot_names = biota.mask.getField(shp, 'Name') 
agb = biota.mask.getField(shp, 'BUFF_DIST') 

# Identify tiles that contain gamma0 data for the shapefile
tiles = biota.mask.getTilesInShapefile(shp)

# Check whether the necessary tiles are all present
# TODO

gamma0 = np.zeros_like(agb, dtype = np.float32)

for lat, lon in tiles:#set(tuple(row) for row in tiles):
    
    data = biota.LoadTile(data_dir, lat, lon, 2007)
    data_gamma0 = data.getGamma0()
    
    plot_mask = biota.mask.rasterizeShapefile(data, shp, location_id = True, buffer_size = 0.)
    
    for n in np.unique(plot_mask[plot_mask != 0]):
        print n
        
        gamma0[n-1] = np.mean(np.logical_and(plot_mask == n, data.mask == False))
    
# And return