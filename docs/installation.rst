Installation
^^^^^^^^^^^^

This section describes the steps to follow to install OCloC in your computer.

Requirements
============

Before installing OCloC you will need to have SAC installed on your computer. 
Sac can be requested at: http://ds.iris.edu/ds/nodes/dmc/forms/sac/

Ensure you can run Python from the command line
-----------------------------------------------

Before you go any further, make sure you have Python and that the expected
version is available from your command line. You can check this by running:

.. code-block:: console

   $ python3 --version

Installation steps
==================
1. Open your terminal and navigate to the directory where you want to clone 
ocloc and then clone the repository

.. code-block:: console

   $ git clone https://github.com/davidn182/ocloc.git
   $ cd ocloc

2. Enter the following command to create a new virtual environment named ``ocloc_env``

.. code-block:: console

   $ python3 -m venv ocloc_env
   $ source ocloc_env/bin/activate

3. Install the required packages by running the following command:

.. code-block:: console

   (ocloc_env) $ pip install -r requirements-dev.txt

4. To use this environment in the Jupyter notebook, you have to register it first:

.. code-block:: console

   (ocloc_env) $ python -m ipykernel install --user --name ocloc_env


Whenever you want to use ocloc, you will need to activate the virtual environment
by running the following command:

.. code-block:: console

   $ source ocloc_env/bin/activate

and to deactivate it, run the following command:

.. code-block:: console

   (ocloc_env) $ deactivate

4. Now you can run the jupyter notebooks located at ``/ocloc/tutorials/`` in your own laptop.
