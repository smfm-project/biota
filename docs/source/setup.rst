Setup instructions
==================

Requirements
------------

The ``biota`` module is written for Python in Linux. It should be possible to run on most Desktop PCs or a Linux server.

.. NOTE::
    ``biota`` now requires Python 3. It may still work with Python 2.7, but this will no longer be supported.

Installing Anaconda Python
--------------------------

We recommend running the ``biota`` module in Anaconda Python.

To install Anaconda Python, open a terminal window, change directory to the location you'd like to install Anaconda Python, and run the following commands:

.. code-block:: console

    wget https://repo.anaconda.com/archive/Anaconda3-2018.12-Linux-x86_64.sh
    bash Anaconda3-2018.12-Linux-x86_64.sh

Once complete, you'll need to add this version of Python to your .bashrc file as follows (replacing ``~/`` with your installation directory):

.. code-block:: console

    # Substitute root for the path to your system's installation and .bashrc file.
    echo 'export PATH="~/anaconda3/bin:$PATH"' >> ~/.bashrc
    exec -l $SHELL

If this has functioned, on executing ``python`` in a terminal window, you should see the following:

.. code-block:: console

    Python 2.7.12 |Anaconda custom (64-bit)| (default, Jul  2 2016, 17:42:40)
    [GCC 4.4.7 20120313 (Red Hat 4.4.7-1)] on linux2
    Type "help", "copyright", "credits" or "license" for more information.
    Anaconda is brought to you by Continuum Analytics.
    Please check out: http://continuum.io/thanks and https://anaconda.org
    >>>

To ensure you are working with the appropriate version of Python as well as the correct modules, we recommend that you create an Anaconda virtual environment set up for running ``biota``. Our recommended procedure of creating the virtual environment and installation is as follows (accepting all prompts):

.. code-block:: console

    conda create -n biota python=3.7 tqdm scikit-image pyshp gdal

Activate the ``biota`` environment whenever opening a new terminal window by running this command:

.. code-block:: console

    conda activate biota

..NOTE:: Remember to activate the ``biota`` environment whenever you want to use ``biota``.

Installing biota
----------------

To install ``biota``, you will need to use the version control software ``git`` (if you don't have ``git``, follow the instructions `here <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`_ ). You can collect the ``biota``  source code with the command:

.. code-block:: console

    git clone https://bitbucket.org/sambowers/biota.git

To install ``biota``, run the following command:

.. code-block:: console

    python setup.py install

If successful, you should now be able to import ``biota`` in Python:

.. code-block:: python

    import biota

Using biota from the command line
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For most applications, the command line interface will be the most straightforward way of using ``biota``.

To avoid having to reference the full path of the Python scripts in biota when using command line tools, add a line to your .bashrc file as follows:

.. code-block:: console

    echo "alias biota='_biota() { python ~/full/path/to/biota/cli/"$1".py $(shift; echo "$@") ;}; _biota'" >> ~/.bashrc

This creates a function that enables you to call ``biota`` just by typing ``biota`` in your terminal window. To run this function, restart your terminal or run ``bash`` (you will only need to do this once). You will then need to activate the ``biota`` environment once again.

You are now ready to start using biota!

What if my install fails?
~~~~~~~~~~~~~~~~~~~~~~~~~

We've not yet anticipated all installation issues with ``biota``. If you encounter issues, please don't hesitate to get in touch with sam.bowers@ed.ac.uk.
