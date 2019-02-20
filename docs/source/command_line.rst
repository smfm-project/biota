Command line tools
==================

``biota`` was initially developed as a tool for use in Python. However, command line functionalities are available to download data and produce rasters of biomass and forest cover, both from single images and over a given period of time.

Downloading is available through the command ``biota``.
Downloading and raster production are available through the command ``biota3``.


Getting ALOS mosaic data
------------------------

Data from the ALOS mosaic product can be accessed from JAXA through a graphical interface `online`_ after `signing up`_. The data is delivered in either 1x1 degree tiles or 5x5 degree collections of tiles.

.. _online: http://www.eorc.jaxa.jp/ALOS/en/palsar_fnf/data/index.htm
.. _signing up: http://www.eorc.jaxa.jp/ALOS/en/palsar_fnf/registration.htm

To obtain national-scale data using this interface is possible, but for large-scale applications ``biota`` contains a command line interface to automate the download and decompression process. ``download.py`` takes latitude, longitude and years as inputs, and returns a downloaded and optionally decompressed product.

The ``biota`` tool includes a script to download ALOS mosaic tiles directly from the JAXA FTP server. The user specifies either a 1x1 or 5x5 degree tile by its upper-left corner latitude and longitude and a year or set of years. The downloader handles access of both the older ALOS-1 filename structure and the ALOS-2 filename structure. The script also optionally decompress images after download.

Help for ``download.py`` can be viewed by typing ``biota download --help`` or ``biota3 download --help``:

.. code-block:: console

    usage: download.py [-h] [-lat DEG] [-lon DEG] [-l] [-y Y [Y ...]] [-o DIR]
                    [-r]

    Download ALOS-1/2 data from JAXA, specifying a particular year and
    latitude/longitude.

    Required arguments:
    -lat DEG, --latitude DEG
                            Latitude of tile upper-left corner.
    -lon DEG, --longitude DEG
                            Longitude of tile upper-left corner.

    Optional arguments:
    -l, --large           Download large tiles. ALOS mosaic tiles are available
                            in 1x1 or 5x5 degree tiles. If downloading large
                            volumes of data, it's usually better to use the
                            latter. If this option is chosen, you must select a
                            lat and lon that's a multiple of 5 degrees.
    -y Y [Y ...], --years Y [Y ...]
                            Year of data to download. Defaults to downloading all
                            data.
    -o DIR, --output_dir DIR
                            Optionally specify an output directory. Defaults to
                            the present working directory.
    -r, --remove          Optionally remove downloaded .tar.gz files after
                            decompression.

For example, to download data for the 1x1 tile at 38 degrees longitude and -8 degrees latitude for the years 2007 and 2010, input:

.. code-block:: console

    biota download -lon 38 -lat -8 -y 2007 2010

To download all tiles for the 5x5 degree area covering 35 - 40 degrees longitude and -5 - -10 degrees latitude for the years 2007 and 2010, input:

.. code-block:: console

    biota download -lon 35 -lat -5 --large -y 2007 2010

Producing a snapshot of vegetation properties
---------------------------------------------

The ``biota3`` tool features a command line option to produce a raster of vegetation properties for a given year. The user specifies the directory where the data are stored, then specifies the designated location and year like for the ``download``. The user may choose to produce rasters for any or all of the following properties: Gamma0, Above-Ground Biomass, Woody Cover. Filtering, downsampling and options specific to each property are available.

Help for this functionality can be viewed by typing ``biota3 snapshot --help`` or ``biota3 snapshot -h``:

.. code-block:: console

    usage: snapshot_cmd.py [-h] [-dir DIR] [-lat DEG] [-lon DEG] [-y Y [Y ...]]
                           [-o {Gamma0,AGB,WoodyCover,all}] [-lf] [-ds FACTOR]
                           [-od DIR] [-pz {HV,HH,VH,VV}] [-ft THRESHOLD]
                           [-at THRESHOLD]

    Process downloaded ALOS-1/2 data to output biomass and forest cover.

    Required arguments:
      -dir DIR, --data_directory DIR
                            absolute path to data directory
      -lat DEG, --latitude DEG
                            Latitude of tile upper-left corner.
      -lon DEG, --longitude DEG
                            Longitude of tile upper-left corner.
      -y Y [Y ...], --years Y [Y ...]
                            Years of data to process.

    Optional arguments:
      -o {Gamma0,AGB,WoodyCover,all}, --output {Gamma0,AGB,WoodyCover,all}
                            Choose which kind of output you want. Defaults to all
                            possible outputs.
      -lf, --speckle        Apply speckle filtering. Defaults to True.
      -ds FACTOR, --downsample FACTOR
                            Apply downsampling. Defaults to 1.
      -od DIR, --output_dir DIR
                            Optionally specify an output directory. Defaults to
                            the present working directory.

    Output-specific arguments:
      -pz {HV,HH,VH,VV}, --polarisation {HV,HH,VH,VV}
                            If you have selected Gamma0 as an output, choose the
                            polarisation. Defaults to HV.
      -ft THRESHOLD, --forest_threshold THRESHOLD
                            If you have selected WoodyCover as an output, choose
                            the miminum forest biomass threshold. Defaults to
                            10tC/ha.
      -at THRESHOLD, --area_threshold THRESHOLD
                            If you have selected WoodyCover as an output, choose
                            the minimum forest area threshold. Defaults to 0ha.


For example, to produce a speckle-filtered map of biomass for the downloaded 1x1 tile at 38 degrees longitude and -8 degrees latitude for the year 2007, run:

.. code-block:: console

    biota3 snapshot -dir /path/to/data/ -lon 38 -lat -8 -y 2007 -o AGB


Producing vegetation change rasters
---------------------------------------------

The ``biota3`` tool features a command line option to produce a raster of vegetation change between two given years. The user specifies the directory where the data are stored, then specifies the designated location and years. The user may choose to produce rasters for any or all of the following properties: Above-Ground Biomass, Woody Cover. Filtering, downsampling and options specific to each property are available.

Help for this functionality can be viewed by typing ``biota3 change --help`` or ``biota3 change -h``:

.. code-block:: console

    usage: change_cmd.py [-h] [-dir DIR] [-lat DEG] [-lon DEG] [-y1 Y [Y ...]]
                         [-y2 Y [Y ...]] [-o {AGB,ChangeType,all}] [-lf]
                         [-ds FACTOR] [-od DIR] [-ft THRESHOLD] [-at THRESHOLD]
                         [-cat CAT] [-cmt CMT] [-cit CIT]

    Process downloaded ALOS-1/2 data to output biomass and forest cover change
    between 2 years.

    Required arguments:
      -dir DIR, --data_directory DIR
                            absolute path to data directory
      -lat DEG, --latitude DEG
                            Latitude of tile upper-left corner.
      -lon DEG, --longitude DEG
                            Longitude of tile upper-left corner.
      -y1 Y [Y ...], --year1 Y [Y ...]
                            First year of data to process.
      -y2 Y [Y ...], --year2 Y [Y ...]
                            Second year of data to process.

    Optional arguments:
      -o {AGB,ChangeType,all}, --output {AGB,ChangeType,all}
                            Choose which kind of output you want. Defaults to all
                            possible outputs.
      -lf, --speckle        Apply speckle filtering. Defaults to True.
      -ds FACTOR, --downsample FACTOR
                            Apply downsampling. Defaults to 1.
      -od DIR, --output_dir DIR
                            Optionally specify an output directory. Defaults to
                            the present working directory.

    Output-specific arguments:
      -ft THRESHOLD, --forest_threshold THRESHOLD
                            If you have selected WoodyCover as an output, choose
                            the miminum forest biomass threshold. Defaults to
                            10tC/ha.
      -at THRESHOLD, --area_threshold THRESHOLD
                            If you have selected WoodyCover as an output, choose
                            the minimum forest area threshold. Defaults to 0ha.
      -cat CAT, --change_area_threshold CAT
                            If you have selected ChangeType as an output, choose
                            the minimum change in forest area threshold. Defaults
                            to 0ha.
      -cmt CMT, --change_magnitude_threshold CMT
                            If you have selected ChangeType as an output, choose
                            the minimum change in biomass threshold. Defaults to
                            0tC/ha.
      -cit CIT, --change_intensity_threshold CIT
                            If you have selected ChangeType as an output, choose
                            the minimum relative change in forest biomass
                            threshold. Defaults to 0.

For example, to produce a speckle-filtered map of biomass change for the downloaded 1x1 tile at 38 degrees longitude and -8 degrees latitude between 2007 and 2010, run:

.. code-block:: console

    biota3 change -dir /path/to/data/ -lon 38 -lat -8 -y1 2007 -y2 2010 -o AGB
