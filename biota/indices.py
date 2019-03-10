#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndimage
from scipy.ndimage.measurements import label

import pdb

import biota.IO

def getContiguousAreas(data, value, min_pixels = 1, contiguity = 'queen'):
    '''
    Get pixels that come from the same contigous area.

    Args:
        data: A numpy array
        value: Pixel value to include in contiguous_area (e.g. True for forest)
        min_area: What minimum area should be included (number of pixels)
        contuguity: Set to rook (4-way) or queen (8-way) connectivity constraint. Defaults to 'queen'.

    Returns:
        A binary array of pixels that meet the conditions
    '''

    assert contiguity in ['rook', 'queen'], "Contiguity must be either 'rook' or 'queen'. Input recieved was <%s>."%str(contiguity)

    # Extract area that meets condition
    binary_array = (data == value) * 1

    # If masked, we use this flag to save the mask for later.
    masked = np.ma.isMaskedArray(binary_array)

    # Set masked areas to non-contiguous value
    if masked:
        mask = np.ma.getmaskarray(binary_array)
        binary_array = binary_array.filled(0)

    # Label contigous areas with a number
    if contiguity == 'rook':
        structure = ndimage.generate_binary_structure(2,1) # 4-way connectivity
    elif contiguity == 'queen':
        structure = ndimage.generate_binary_structure(2,2) # 8-way connectivity

    location_id, n_areas = label(binary_array, structure = structure)

    # Get count of each value in array
    label_area = np.bincount(location_id.flatten())[1:]

    # Find those IDs that meet minimum area requirements
    include_id = np.arange(1, n_areas + 1)[label_area >= min_pixels]

    # Get a binary array of location_id pixels that meet the minimum area requirement
    contiguous_area = np.in1d(location_id, include_id).reshape(data.shape).astype(np.bool)

    # Return an array giving values to each area
    location_id[contiguous_area == False] = 0

    # Re-number location_ids 1 to n, given that some unique value shave now been removed
    location_id_unique, location_id_indices = np.unique(location_id, return_inverse = True)
    location_id = np.arange(0, location_id_unique.shape[0], 1)[location_id_indices].reshape(data.shape)

    # Put mask back in if input was a masked array
    if masked:
        contiguous_area = np.ma.array(contiguous_area, mask = mask)
        location_id = np.ma.array(location_id, mask = mask)

    return contiguous_area, location_id


def _calculatePatchSize(tile, patch_size):
    """
    Function to automatically calculate a patch_size

    Args:
        tile: Either an ALOS tile (biota.LoadTile()) or an ALOS change object (biota.LoadChange())
        patch_size = Number of pixels to build into a single patch. Set to 'auto' for an approx 100 x 100 output image

    Returns:
        The patch size
    """

    assert patch_size == 'auto' or type(patch_size) == float or type(patch_size) == int, "Patch size must be numeric or set to 'auto'."

    # Calculate patch_size (which aims for a 100 x 100 output image if set to auto.)
    if patch_size == 'auto': patch_size = int(round(((tile.xSize + tile.ySize) / 2.) / 100.))

    return patch_size


def _buildOutputArray(tile, patch_size, dtype = np.int):
    """
    Function to generate an appropriately sized output array

    Args:
        tile: Either an ALOS tile (biota.LoadTile()) or an ALOS change object (biota.LoadChange())
        patch_size = Number of pixels to build into a single patch. Set to 'auto' for an approx 100 x 100 output image

    Returns:
        An empty array
    """

    # Get nodata value
    if dtype == np.int:
        nodata = tile.nodata_byte
    else:
        nodata = tile.nodata

    # Calculate output output size based on reduction patch_size
    output_size = int(round(((tile.xSize + tile.ySize) / 2.) / patch_size))

    # Empty array for output
    output_array = np.ma.array(np.zeros((output_size, output_size), dtype = dtype), mask = np.zeros((output_size, output_size), dtype = np.bool)) + nodata

    return output_array


def _getBlocks(tile, array):
    """
    Calculate blocks based on an ALOS tile and a downsampled array.

    Args:
        tile: Either an ALOS tile (biota.LoadTile()) or an ALOS change object (biota.LoadChange())
        array: The downsampled numpy array (from _buildOutputArray()).

    Returns:
        A list of tuples with the y, x coordinates of the output image, and ymin, ymax, xmin, xmax of input image.
    """

    ns, ms, ymins, ymaxs, xmins, xmaxs = [], [], [], [], [], []

    for n, ymin in enumerate(np.linspace(0, tile.ySize, array.shape[0], endpoint = False)):
        ymax = ymin + (tile.ySize / array.shape[0])
        for m, xmin in enumerate(np.linspace(0, tile.xSize, array.shape[1], endpoint = False)):
            xmax = xmin + (tile.xSize / array.shape[1])

            ns.append(n); ms.append(m)
            ymins.append(int(round(ymin))); ymaxs.append(int(round(ymax)))
            xmins.append(int(round(xmin))); xmaxs.append(int(round(xmax)))

    return list(zip(ns, ms, ymins, ymaxs, xmins, xmaxs))



def calculateLDI(tile, patch_size = 'auto', output = False, show = False):
    """
    Landcover Division Index (LDI) is an indicator of the degree of habitat coherence, ranging from 0 (fully contiguous) to 1 (very fragmented).

    LDI is defined as the probability that two randomly selected points in the landscape are situated in two different patches of the habitat (Jaeger 2000; Mcgarigal 2015).

    Note that fragmentation of both an undisturbed forest area and contiguous agriculture will both be low, so interpret this index carefully

    Args:
        tile: An ALOS tile from niota.LoadTile()
        patch_size: Number of pixels to build into a single patch. Set to 'auto' for an approx 100 x 100 output image

    Returns:
        An numpy array.
    """

    from osgeo import gdal

    def _computeLDI(unique_ids):
        '''
        Calculates an index of landscape heterogeneity
        '''

        # If masked array, separate data from mask
        if np.ma.isMaskedArray(unique_ids):
            mask = np.ma.getmask(unique_ids)
            unique_ids = np.ma.getdata(unique_ids)
        else:
            mask = np.zeros_like(unique_ids, dtype=np.bool)

        # Calculate LDI
        selection_1 = np.zeros_like(unique_ids,dtype = np.bool).ravel()
        selection_2 = np.zeros_like(unique_ids,dtype = np.bool).ravel()

        # Draw 1/25th of pixels at random
        selection_1[:(selection_1.shape[0]/5)] = True
        selection_2[:(selection_2.shape[0]/5)] = True

        # Shuffle to get random pixels
        np.random.shuffle(selection_1)
        np.random.shuffle(selection_2)

        selection_1 = selection_1.reshape(unique_ids.shape)
        selection_2 = selection_2.reshape(unique_ids.shape)

        # This ignores IDs where one or the other of selection falls in an area of nodata
        mask_selection = np.logical_or(mask[selection_1], mask[selection_2])
        ID1 = unique_ids[selection_1][mask_selection==False]
        ID2 = unique_ids[selection_2][mask_selection==False]

        # Shuffle again to randomise pairs spatially
        np.random.shuffle(ID1)
        np.random.shuffle(ID2)

        # Number of random points from different patches
        different_patches = np.sum(ID1 != ID2)

        # Divide by potential matches
        if (mask_selection==False).sum() > 0:
            LDI = int(round((float(different_patches) / (mask_selection==False).sum()) * 100))

        # Unless everything is masked, in which case return nan.
        else:
            LDI = np.nan

        return LDI

    # Create an output array
    patch_size = _calculatePatchSize(tile, patch_size)
    LDI = _buildOutputArray(tile, patch_size)

    # Load woody cover
    woody_cover = tile.getWoodyCover()

    # For each patch...
    for n, m, ymin, ymax, xmin, xmax in _getBlocks(tile, LDI):

        # Extract the data
        this_patch = woody_cover[ymin:ymax, xmin:xmax]

        # Get connected patches of forest and nonforest
        _, forest_id = getContiguousAreas(this_patch, True, contiguity = tile.contiguity)
        _, nonforest_id = getContiguousAreas(this_patch, False, contiguity = tile.contiguity)

        # Supply a unique ID to each habitat patch, whether forest or nonforest.
        unique_ids = forest_id.copy()
        unique_ids[nonforest_id != 0] = forest_id[nonforest_id != 0] + (nonforest_id[nonforest_id != 0] + np.max(forest_id) + 1)

        #  If at least 50 % of data is present...
        if this_patch.mask.sum() <= ((patch_size ** 2) * 0.5):

            ## Compute LDI
            LDI.data[n, m] = _computeLDI(unique_ids)

        else:

            LDI.mask[n, m] = True

    # Output GeoTiff
    if output: biota.IO.outputGeoTiff(LDI, tile.output_pattern%'LDI', tile.shrinkGeoT(patch_size), tile.proj, output_dir = tile.output_dir, dtype = gdal.GDT_Int32, nodata = tile.nodata)

    # Display
    if show: biota.IO.showFigure(LDI, tile.lat, tile.lon, title = 'LDI', cbartitle = '%', vmin = 0, vmax = 100, cmap = 'Oranges')

    return LDI


def calculateLDIChange(tile_change, patch_size = 'auto', output = False, show = False):
    """
    Landcover Division Index (LDI) is an indicator of the degree of habitat coherence, ranging from 0 (fully contiguous) to 1 (very fragmented).

    This function calculates the change in LDI between two tiles.

        Args:
        tile_change: An ALOS change object from biota.LoadChange()
        patch_size: Number of pixels to build into a single patch. Set to 'auto' for an approx 100 x 100 output image

    Returns:
        An numpy array.

    """

    from osgeo import gdal

    patch_size = _calculatePatchSize(tile_change, patch_size)

    LDI_t1 = calculateLDI(tile_change.tile_t1, patch_size = patch_size)
    LDI_t2 = calculateLDI(tile_change.tile_t2, patch_size = patch_size)

    LDI_change = LDI_t2 - LDI_t1

    # Output GeoTiff
    if output: biota.IO.outputGeoTiff(LDI_change, tile_change.output_pattern%'LDIChange', tile_change.shrinkGeoT(patch_size), tile_change.proj, output_dir = tile_change.output_dir, dtype = gdal.GDT_Int32, nodata = tile_change.nodata)

    # Display
    if show: biota.IO.showFigure(LDI_change, tile_change.lat, tile_change.lon, title = 'LDI Change', cbartitle = '%', vmin = -25, vmax = 25, cmap = 'RdBu_r')

    return LDI_change



def calculateTWC(tile, patch_size = 'auto', output = False, show = False):
    """
    Total woody cover (TWC) describes the proportion of woody cover in a downsampled image.

    Args:
        tile: An ALOS tile from biota.LoadTile()
        patch_size: Number of pixels to build into a single patch. Set to 'auto' for an approx 100 x 100 output image

    Returns:
        A numpy array,
    """

    from osgeo import gdal

    # Create an output array
    patch_size = _calculatePatchSize(tile, patch_size)
    TWC = _buildOutputArray(tile, patch_size)

    # Load woody cover
    woody_cover = tile.getWoodyCover()

    # For each patch...
    for n, m, ymin, ymax, xmin, xmax in _getBlocks(tile, TWC):

        # Extract the data
        WC = woody_cover[ymin:ymax, xmin:xmax]

        #  If at least 50 % of data is present...
        if WC.mask.sum() <= ((patch_size ** 2) * 0.5):

            # Calculate proportion of woody cover in patch
            TWC.data[n, m] = int(round((float(WC.sum()) / ((patch_size ** 2) - WC.mask.sum())) * 100))

        else:

            TWC.mask[n, m] = True

    # Output GeoTiff
    if output: biota.IO.outputGeoTiff(TWC, tile.output_pattern%'TWC', tile.shrinkGeoT(patch_size), tile.proj, output_dir = tile.output_dir, dtype = gdal.GDT_Int32, nodata = tile.nodata)

    # Display
    if show: biota.IO.showFigure(TWC, tile.lat, tile.lon, title = 'TWC', cbartitle = '%', vmin = 0, vmax = 100, cmap = 'YlGn')

    return TWC


def calculateWCC(tile_change, patch_size = 'auto', output = False, show = False):
    """
    Woody cover change (WCC) describes the loss of woody cover between two images.

    Args:
        tile_change: An ALOS change object from biota.LoadChange()
        patch_size: Number of pixels to build into a single patch. Set to 'auto' for an approx 100 x 100 output image

    Returns:
        A numpy array.
    """

    from osgeo import gdal

    # Create an output array
    patch_size = _calculatePatchSize(tile_change, patch_size)
    WCC = _buildOutputArray(tile_change, patch_size)

    # Get total woody cover for time 1 and time 2
    TWC_t1 = calculateTWC(tile_change.tile_t1, patch_size = patch_size)
    TWC_t2 = calculateTWC(tile_change.tile_t2, patch_size = patch_size)

    # Woody cover change is defined as the difference between the two
    WCC = TWC_t2 - TWC_t1

    # Output GeoTiff
    if output: biota.IO.outputGeoTiff(WCC, tile_change.output_pattern%'WCC', tile_change.shrinkGeoT(patch_size), tile_change.proj, output_dir = tile_change.output_dir, dtype = gdal.GDT_Int32, nodata = tile_change.nodata_byte)

    # Display
    if show: biota.IO.showFigure(WCC, tile_change.lat, tile_change.lon, title = 'WCC', cbartitle = '%', vmin = -30, vmax = 30, cmap = 'RdBu')


def calculateProportionalChange(tile_change, change_type, patch_size = 'auto', output = False, show = False):
    """
    Calculate the proportion of a downsampled regions subject to a class of change.

    Args:
        tile_change: An ALOS change object from biota.LoadChange()
        change_type: String representing a change from tile_change.getChangeType(). e.e. 'deforestation' or 'degradation'.
        patch_size: Number of pixels to build into a single patch. Set to 'auto' for an approx 100 x 100 output image

    Returns:
        a numpy array
    """

    from osgeo import gdal

     # Create an output array
    patch_size = _calculatePatchSize(tile_change, patch_size)
    change_downsampled = _buildOutputArray(tile_change, patch_size)

    change = tile_change.getChangeType()[change_type]

    # For each patch...
    for n, m, ymin, ymax, xmin, xmax in _getBlocks(tile_change, change_downsampled):

        # Extract the data
        C = change[ymin:ymax, xmin:xmax]

        #  If at least 50 % of data is present...
        if C.mask.sum() <= ((patch_size ** 2) * 0.5):

            # Calculate proportion of patch subject to change_type
            change_downsampled.data[n, m] = int(round((float(C.sum()) / ((patch_size ** 2) - C.mask.sum())) * 100))

        else:
            change_downsampled.mask[n, m] = True

    # Output GeoTiff
    if output: biota.IO.outputGeoTiff(change_downsampled, tile_change.output_pattern%'%sDownsampled'%change_type.title(), tile_change.shrinkGeoT(patch_size), tile_change.proj, output_dir = tile_change.output_dir, dtype = gdal.GDT_Int32, nodata = tile_change.nodata)

    # Display
    if show: biota.IO.showFigure(change_downsampled, tile_change.lat, tile_change.lon, title = '%s Downsampled'%change_type.title(), cbartitle = '%', vmin = 0, vmax = 50, cmap = 'Spectral_r')


    return change_downsampled
