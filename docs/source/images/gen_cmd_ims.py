import biota

# Gamma0 backscatter

tile_t1 = biota.LoadTile('~/SMFM/ALOS_data/',-8, 38, 2007, lee_filter = True)

tile_t1.getGamma0(show=True)

tile_t1.getGamma0(polarisation='HH',units='decibels', show=True)

tile_t1 = biota.LoadTile('~/SMFM/ALOS_data',-8,38,2007, lee_filter=False)

tile_t1.getGamma0(show=True)

# AGB
tile_t1 = biota.LoadTile('~/SMFM/ALOS_data/',-8, 38, 2007, lee_filter = True)

tile_t1.getAGB(show=True)

# Woody cover
tile_t1.getWoodyCover(show = True)

tile_t1 = biota.LoadTile('~/SMFM/ALOS_data/',-8, 38, 2007, lee_filter = True, forest_threshold = 20)

tile_t1.getWoodyCover(show=True)

tile_t1 = biota.LoadTile('~/SMFM/ALOS_data/',-8, 38, 2007, lee_filter = True, forest_threshold = 20, area_threshold = 1)

tile_t1.getWoodyCover(show=True)

# Change type
tile_t1 = biota.LoadTile('~/SMFM/ALOS_data/',-8, 38, 2007, lee_filter = True)
tile_t2 = biota.LoadTile('~/SMFM/ALOS_data/',-8, 38, 2010, lee_filter = True)

tile_change = biota.LoadChange(tile_t1, tile_t2)

tile_change.getAGBChange(show=True)

tile_change.getChangeType(show=True)

tile_change = biota.LoadChange(tile_t1, tile_t2, change_area_threshold = 1)

tile_change.getChangeType(show=True)

tile_change = biota.LoadChange(tile_t1, tile_t2, change_area_threshold = 1, change_magnitude_threshold = 5)

tile_change.getChangeType(show=True)

tile_change = biota.LoadChange(tile_t1, tile_t2, change_area_threshold = 1, change_magnitude_threshold = 5, change_intensity_threshold = 0.25)

tile_change.getChangeType(show=True)
