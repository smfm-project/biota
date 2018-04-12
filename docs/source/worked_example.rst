


Here we'll run through an example based on an area of Southeastern Tanzania called Kilwa District (Fig). This is a location where the University of Edinburgh coordinates a network of forest plots where we have good data on aboveground biomass.


Downloading data
----------------

To follow.

Open Python and import biota
----------------------------

Open a terminal window (right-click the Desktop, and select 'Open Terminal'), and run the command ``python`` or ``ipython``. We recommend use of ``ipython`` if available, which has a range of features that make it more user-friendly than standard Python. If successful, you should see something that looks like the following:

.. code-block::
    
    Python 2.7.12 |Anaconda custom (64-bit)| (default, Jul  2 2016, 17:42:40) 
    Type "copyright", "credits" or "license" for more information.

    IPython 5.1.0 -- An enhanced Interactive Python.
    ?         -> Introduction and overview of IPython's features.
    %quickref -> Quick reference.
    help      -> Python's own help system.
    object?   -> Details about 'object', use 'object??' for extra details.

    In [1]: 

To load the biota module, type:

.. code-block:: python
    
    >> import biota

Loading a tile
--------------

Data from the ALOS mosaic is provided as a series of 1 x 1 degree tiles. To load a tile in memory, we need to tell ``biota`` what directory the ALOS mosaic data are stored in, and what latitude and longitude we want to load. To save us from writing them out repeatedly, we can store these as variables:

.. code-block:: python
    
    >>> data_dir = '~/DATA/'
    >>> latitude = -9
    >>> longitude = 38

The biota function to load an ALOS tile can be called with the function ``biota.loadTile()``, which takes inputs of (i) the data directory, (ii) the latitude, (iii) the longitude, and (iv) the year (in this order). Here we'll load in data for 2007 using the three variables we previously defined:

.. code-block:: python
    
    >>> tile_2007 = biota.LoadTile(data_dir, latitude, longitude, 2007)

The new object called ``tile_2007`` has a range of attributes. These can be accessed as follows:

.. code-block:: python

    >>> tile_2007.year
    2007
    >>> tile_2007.lat
    -9
    >>> tile_2007.lon
    38
    >>> tile_2007.directory
    '~/DATA/S05E035_07_MOS/'
    >>> tile_2007.satellite
    'ALOS-1'
    >>> tile_2007.xSize, data_2007.ySize # Raster size, in pixels
    (4500, 4500)
    >>> tile_2007.xRes, data_2007.yRes # Pixel resolution in meters
    (24.401, 24.579)

*Advanced:* The tile also contains projection information for interaction with ``GDAL``:

.. code-block:: python
    
    >>> tile_2007.extent # Extent in the format minlon, minlat, maxlon, maxlat
    (38.0, -10.0, 39.0, -9.0)
    >>> tile_2007.geo_t # A geo_transform object
    (38.0, 0.00022222222222222223, 0.0, -9.0, 0.0, -0.00022222222222222223)
    >>> tile_2007.proj # Projection wkt
    'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]'

There are a range of other options that can be specified when opening an ALOS tile, but we'll return to these.
    
Extracting backscatter information
----------------------------------

The ``biota`` module is programmed to calibrate ALOS mosaic data to interpretable units of backscatter. This is performed with the ``getGamma0()`` function. The data are returned as a masked ``numpy`` array:

.. code-block:: python
    
    >>> gamma0_2007 = tile_2007.getGamma0()
    >>> gamma0_2007
    masked_array(data =
    [[0.0669537278370757 0.04214984634805357 0.05141784577914017 ...,
    0.029133617952838833 0.024789602664736045 0.040281545637899534]
    [0.031600461516752214 0.04214984634805357 0.05141784577914017 ...,
    0.03435099209051573 0.028222499657083098 0.03354230142969638]
    [0.031600461516752214 0.04050920492690238 0.06216969020533775 ...,
    0.037654602824076254 0.04403078198836734 0.025848435873858728]
    ..., 
    [0.0900164548052426 0.0662958895217059 0.07768386584418481 ...,
    0.049509525268380976 0.0346139149132766 0.021227103665645366]
    [0.08548700525257016 0.0888309264753313 0.11198792676214335 ...,
    0.08441404357533155 0.06655132961906884 0.05196509713141002]
    [0.07134665398730806 0.05708835833035639 0.07595717558689226 ...,
    0.021496125937039534 0.027866832136739485 0.0629132766445086]],
                mask =
    [[False False False ..., False False False]
    [False False False ..., False False False]
    [False False False ..., False False False]
    ..., 
    [False False False ..., False False False]
    [False False False ..., False False False]
    [False False False ..., False False False]],
        fill_value = 1e+20)

By default the image loaded is 'HV' polarised in 'natural' units. It's also possible to specify options for the polarisation ('HV' or 'HH') and the units ('natural' or 'decibels') as follows:

.. code-block:: python
    
    >>> gamma0_HH_2007 = tile_2007.getGamma0(polarisation = 'HH', units = 'decibels')
    >>> gamma0_HV_2007 = tile_2007.getGamma0(polarisation = 'HV', units = 'decibels')

If we want to visualise this data, we can run:

.. code-block:: python
    
    >>> gamma0_2007 = tile_2007.getGamma0(polarisation = 'HV', units = 'decibels', show = True)

Which produces the following image:

.. figure:: images/gamma0.png

If we want to save this data to a geoTiff, we can run:

.. code-block:: python
    
    >>> gamma0_2007 = tile_2007.getGamma0(polarisation = 'HV', units = 'decibels', output = True)

which will write a GeoTiff file to the current working directory. This file can be processed and visualised in standard GIS and remote sensing software (e.g. QGIS, GDAL).

Building a map of AGB
---------------------

In a similar way to loading gamma0 backscatter, we can generate images of AGB. Note: by default ``biota`` uses an equation calibrated for ALOS-1 in miombo woodlands of Southern Africa. It's advisable to have a locally calibrated biomass-backscatter equation to improve on predictions.

.. code-block:: python

        agb_2007 = tile_2007.getAGB(show = True) # To display AGB
        agb_2007 = tile_2007.getAGB(output = True) # To output AGB map to a GeoTiff

.. figure:: images/agb.png

Building a forest cover map
---------------------------

With a few modifications to the above scipt, we can generate a map of forest cover baed on thresholds of minimum AGB and minimum area. Here we assume a forest definition consistent with a minimum AGB of 10 tC/ha over a minimum area of 1 hectare.

.. code-block:: python

    # Import the biota module
    import biota

    # Define a variable with the location of ALOS tiles
    data_dir = '~/SMFM/ALOS_data/'

    # Loop through the latitudes/longitudes that cover Kilwa District
    for longitude in range(35, 40):
        for latitude in range(-10, -5):
                    
            # Load the tile for the year 2010. Apply a speckle filter (lee_filter = True).
            # Apply thresholds for forest cover (10 tC/ha) with a minimum area (1 hectare)
            tile = biota.LoadTile(data_dir, longitude, latitude, 2010, lee_filter = True, forest_threshold = 10., area_threshold = 1)
            
            # Calculate forest cover and output to a GeoTiff
            tile.getWoodyCover(output = True)
