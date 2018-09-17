Worked example: part 1
======================

In this worked example we'll run through an example based on an area of Southeastern Tanzania called Kilwa District (Fig). This is a location where the University of Edinburgh coordinates a network of forest plots where there is good data on aboveground biomass.

In this section we'll use ``biota`` to download ALOS mosaic data for the Kilwa region.

Preparation
-----------

First ensure that you’ve followed the setup instructions successfully.

Open a terminal, and use ``cd`` to navigate to the location you’d like to store data.

.. code-block:: console

    cd /home/user
    mkdir worked_example_biota
    cd worked_example_biota

Use ``mkdir`` to make a directory to contain the ALOS mosaic data:

.. code-block:: console
    
    mkdir DATA

Here we’ll demonstrate the process for downloading the large tile (5x5 degrees) covering Kilwa District.

To begin, navigate to the DATA folder.

.. code-block:: console
    
    cd DATA

Downloading data
----------------

The first step is to download Sentinel-2 level 1C data from JAXA.

For this we use the ``download.py`` tool. To do this we need to specify the latitude/longitude, and years to download. Kilwa District is contained within the large tile at longitude 35 to 40 degrees and latitude -5 - -10 degrees. For the purposes of this demonstration, we’ll download ALOS-1 data from 2007 and 2010.

These options can be encoded as follows:

.. code-block:: console
    
    biota download -lon 35 -lat -5 --large -y 2007 2010 -r

As we didn’t specify the option ``-o`` (``--output``), data will output to the current working directory. We also included ``-r`` (``--remove``) flag, meaning that intermediate .zip files downloaded from the internet will be deleted after download.

Wait for the two files (2007 and 2010) to finish downloading before proceeding to the next step. By the time the processing is complete, your ``DATA/`` directory should contain the following files (show files in the currenty working directory with the command ``ls``).

.. code-block:: console

    FILE1
    FILE2



