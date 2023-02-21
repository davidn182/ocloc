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


1. Open your terminal and navigate to the directory where you want to create
 your virtual environment.

2. Enter the following command to create a new virtual environment named "myenv" 

.. code-block:: console

   $ python3 -m venv ocloc_env

3. Activate the virtual environment by running the following command:

.. code-block:: console

   $ source ocloc_env/bin/activate

Now let's clone the project. Notice that the project is called ocloc 
and the virtual environment is called ocloc_env. If you want to use ocloc for
the name of the virtual environment, make sure to install in a different
directory. For example, you can install the virtual environment in your home directory
and clone the project in your working directory.

.. code-block:: console

   (ocloc_env) $ git clone https://github.com/davidn182/ocloc.git
   (ocloc_env) $ cd ocloc

4. Install the required packages by running the following command:

.. code-block:: console
   (ocloc_env) $ pip install -r requirements-dev.txt

5. Once you have finished working in your virtual environment, you can 
deactivate it by running the following command: 

.. code-block:: console

   (ocloc_env) $ deactivate

Now you can run the jupyter notebooks located at ``/ocloc/tutorials/`` in your own laptop.