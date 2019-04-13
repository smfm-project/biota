
The various forms of biota
===========================

`biota` comes in 3 forms: a Python module, a command line function, and a Graphical User Interface (GUI). We recommend the use of the Python module for users with good knowledge of Python who wish to produce and adapt their own results. The command line form of `biota` is particularly adapted for batch jobs, and has been used to produce country-scale maps of biomass loss and probably deforestation with little scripting. Finally, the `biota` GUI provides an accessible tool for users who are less comfortable with scripting and terminals.




Using biota in Python
---------------------

If your installation was successful, you should be able to import ``biota`` in Python:

.. code-block:: python

    import biota

You may now use biota in your Python scripts! See the worked examples to explore the different methods and attributes of `biota` in Python. Using `biota` in python will give you the most flexibility.



Using biota from the command line
---------------------------------

For most applications, the command line interface will be the most straightforward way of using ``biota``. It's also the best way to run batch jobs.

To avoid having to reference the full path of the Python scripts in biota when using command line tools, add the following line to your .bashrc file:

.. code-block:: console

    alias biota='_biota() { python ~/full/path/to/biota/cli/"$1".py $(shift; echo "$@") ;}; _biota'

This creates a function that enables you to call ``biota`` just by typing ``biota`` in your terminal window. To run this function, restart your terminal or run ``bash`` (you will only need to do this once). You will then need to activate the ``biota`` environment once again.

Take a look at the worked examples to explore the possibilities when you using `biota` from the command line.


Using the biota Graphical User Interface (GUI)
---------------------------------------------

The ``biota`` Graphical User Interface (GUI) allows you to use `biota` with minimal usage of the command line. In your terminal or the `Anaconda prompt` (if you're using Windows), navigate to the folder where you installed `biota` by typing:

.. code-block:: console

    cd /full/path/to/biota/

This folder should have `/biota/` at the end of its path. In this folder you will find a folder called `Biota_app/`. Navigate to it by typing:

.. code-block:: console

    cd Biota_app/

To open the GUI, type

.. code-block:: console

    python Biota.py

The GUI window should then open. Now you can use `biota` several times with minimal typing.
