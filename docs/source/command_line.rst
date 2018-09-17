Command line tools
==================

``biota`` is mostly a tool for use in Python. However, there is some command line functionality for downloading, which over time we may add to.

Getting ALOS mosaic data
------------------------

Data from the ALOS mosaic product can be accessed from JAXA through a graphical interface `online`_ after `signing up`_. The data is delivered in either 1x1 degree tiles or 5x5 degree collections of tiles.

.. _online: http://www.eorc.jaxa.jp/ALOS/en/palsar_fnf/data/index.htm
.. _signing up: http://www.eorc.jaxa.jp/ALOS/en/palsar_fnf/registration.htm

To obtain national-scale data using this interface is possible, but for large-scale applications ``biota`` contains a command line interface to automate the download and decompression process. ``download.py`` takes latitude, longitude and years as inputs, and returns a downloaded and optionally decompressed product.

The ``biota`` tool includes a script to download ALOS mosaic tiles directly from the JAXA FTP server. The user specifies either a 1x1 or 5x5 degree tile by its upper-left corner latitude and longitude and a year or set of years. The downloader handles access of both the older ALOS-1 filename structure and the ALOS-2 filename structure. The script also optionally decompress images after download.

Help for ``download.py`` can be viewed by typing ``biota download --help``:

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

Other tools
-----------

Other command-line tools may follow, but for the moment any further processing using ``biota`` needs to be conducted with Python.