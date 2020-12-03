Setup instructions
==================

Preamble
--------

`biota` is written in Python and was developped on a Linux platform. Provided you can use Python and can install modules on your computer, it should be possible to run on most OS or a Linux server.

On this page, we explain how to set up your Linux or Windows machine to use `biota`.

Installing Anaconda
-------------------

We recommend running ``biota`` in Anaconda.

If you are using a Linux machine, open a terminal window, change directory to the location you downloaded Anaconda Python, and run the following commands:

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


For Windows users, go to the `Anaconda website <https://www.anaconda.com/distribution/>`_ ) and download the installer for your version of Windows (64-bit or 32-bit). Once the download is finished, run the installer. This may take a while, but when it is done you will be able to open the `Anaconda Prompt`.


Setting up your Anaconda environment
-----------------------------------

To ensure you are working with the appropriate version of Python as well as the correct modules, we recommend that you create an Anaconda virtual environment set up for running ``biota``. This is done by running the following commands in your terminal or the `Anaconda prompt` (recommended procedure):

.. code-block:: console

    conda create -n biota -c conda-forge python=3.7 tqdm scikit-image pillow pyshp gdal

.. NOTE::
  the GDAL package is notoriously temperamental. If this step fails, try again and add ` openssl=1.0` at the end of the line


Activate the ``biota`` environment whenever opening a new terminal window by running this command:

.. code-block:: console

    conda activate biota

.. NOTE::
  Remember to activate the ``biota`` environment whenever you want to use ``biota``.


If you are SURE you won't use anything else than `biota`, you can do without virtual environments. In this case, just type:

.. code-block:: console

    conda install -c conda-forge python=3.7 tqdm scikit-image pillow pyshp gdal


If you want to use the `biota` graphical interface, you need an extra package called `PyQt5`. To install it, type:

.. code-block:: console

    pip install pyqt5

Installing biota
----------------

Navigate to the folder where you want to install `biota`. To do this in both Linuw and Windows, type the following:

.. code-block:: console

    cd /full/path/to/your/favorite/folder/

.. NOTE::

  If you don't want to type the full path (and really, you souldn't), copy it from your file explorer into the terminal after the `cd ` command.

To install ``biota``, you will need to use the version control software ``git`` (if you don't have ``git``, follow the instructions `here <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`_ ). You can collect the ``biota``  source code with the command:

.. code-block:: console

    git clone https://github.com/smfm-project/biota.git

To install ``biota``, run the following command:

.. code-block:: console

    python setup.py install

Congratulations, you are now ready to use `biota`.
