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


Developer Notes
===============

.. admonition:: Developer Note

   If you are contributing to development, follow these setup instructions to ensure a consistent environment and workflow.

**Step 1: Create a Virtual Environment**

Create a new Python virtual environment. This isolates your project dependencies from your system Python environment.

.. code-block:: console

   $ python3 -m venv ocloc_env

Make sure that no Anaconda environment is activated:

.. code-block:: console

   $ conda deactivate

**Step 2: Activate the Virtual Environment**

Activate the newly created virtual environment.

.. code-block:: console

   $ source ocloc_env/bin/activate  # On Windows, use `ocloc_env\\Scripts\\activate`

**Step 3: Clone the Repository**

Navigate to the directory where you want to clone the repository and run the following commands:

.. code-block:: console

   $ git clone https://github.com/davidn182/ocloc.git
   $ cd ocloc

**Step 4: Add Dependencies (Optional)**

If you need to add a new dependency to the project, make sure you're in the `ocloc` directory and then use `poetry`:

.. code-block:: console

   $ poetry add <your_dependency>  # Replace <your_dependency> with the package you want to add, e.g., jupyter

**Step 5: Update Your Environment**

To reflect the new changes in your virtual environment, activate it (if not already activated) and install the package in "editable" mode:

.. code-block:: console

   $ cd ocloc  # If not already in this directory
   $ pip install -e .

This will link the changes you make in the source code to your virtual environment, allowing for real-time testing and debugging.

