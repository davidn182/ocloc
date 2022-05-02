Installation
^^^^^^^^^^^^


Requirements
============
Before installing OCloC you will need to have SAC installed on your computer. Sac can be requested at: http://ds.iris.edu/ds/nodes/dmc/forms/sac/

Installation steps
==================


First install `Anaconda <https://www.anaconda.com/products/individual#Downloads>`_ following the instructions on their site

After installing Anaconda please add the ``conda-forge`` channel as:

.. code-block:: console

   $ conda config --add channels conda-forge

Now let's create a conda environment where we will install some dependecies.

.. code-block:: console

   $ conda create -n ocloc obspy pandas seaborn
   $ conda activate ocloc  # activate the ocloc environment
   (ocloc) $

Then clone the repository in your favorite location. Now you can run the jupyter notebooks in your own laptop.

In the future we will be adding the package to Anaconda.
