Setup instructions
==================

Requirements
------------

The ``biota`` module is written for Python in Linux. It should be possible to run on most Desktop PCs or a Linux server.

Installing Anaconda Python
--------------------------

We recommend running the ``biota`` module in Anaconda Python.

To install Anaconda Python, open a terminal window, change directory to the location you'd like to install Anaconda Python, and run the following commands:

.. code-block:: console
    
    wget https://repo.continuum.io/archive/Anaconda2-4.2.0-Linux-x86_64.sh
    bash Anaconda2-4.2.0-Linux-x86_64.sh

    
Once complete, you'll need to add this version of Python to your .bashrc file as follows:

.. code-block:: console
    
    # Substitute root for the path to your system's installation and .bashrc file.
    echo 'export PATH="~/anaconda2/bin:$PATH"' >> ~/.bashrc
    exec -l $SHELL


If this has functioned, on executing ``python`` in a terminal window, you should ssee the following:

.. code-block:: console

    Python 2.7.12 |Anaconda custom (64-bit)| (default, Jul  2 2016, 17:42:40) 
    [GCC 4.4.7 20120313 (Red Hat 4.4.7-1)] on linux2
    Type "help", "copyright", "credits" or "license" for more information.
    Anaconda is brought to you by Continuum Analytics.
    Please check out: http://continuum.io/thanks and https://anaconda.org
    >>> 

``biota`` requires two further modules to function. Install them as in the terminal window as follows, accepting all prompts:

.. code-block:: console
    
    conda install -c anaconda gdal
    conda install -c conda-forge pyshp

Installing biota
----------------

To install ``biota``, enter the following in a terminal window:

.. code-block:: console
    
    pip install git+https://bitbucket.org/sambowers/biota

If successful, you should now be able to import ``biota`` in Python:

.. code-block:: python
    
    import biota

To avoid having to reference the full path of the Python scripts in biota when using command line tools, add the following line to your .bashrc file: 

.. code-block:: console
    
    echo "alias biota='_biota() { python ~/biota/biota/"$1".py $(shift; echo "$@") ;}; _biota'" >> ~/.bashrc
