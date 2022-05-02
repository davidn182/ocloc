.. ocloc documentation master file, created by
   sphinx-quickstart on Wed Oct 13 09:47:23 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to OCloC's documentation!
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

OCloC is an open-source project dedicated to detect and correct timing errors when using passive seismic records. This Python Package is particularly relevant for projects that use Ocean-Bottom Seismographs because their lack of GPS connection can, and most certainly will, induce a timing error. Needless to say, interpretations of the subsurface might be wrong if the seismic stations contain timing errors. The aim of this project is to prevent this from hapening.

The main goal of the OCloC project is to rapidly detect and correct clock errors of seismic networks. Please note that for doing so, this package requires an understanding of the basics of seismic interferometry. If you are not familiar with this method please take a look at the Introduction section.

Our code is stored in github, the development branches, or the latest stable release can be found here.

.. _my-reference-label:

Contents
--------------------------

.. toctree::
   :maxdepth: 2
   
   about.rst
   installation.rst
   tutorial/tutorial.rst
   API.rst
   rst_commands_david.rst


Indices and tables
^^^^^^^^^^^^^^^^^^

* :ref:`genindex`
* :ref:`modindex` Graphs etc.

Citation
^^^^^^^^
If you use this package in your work, please refer to the following works:

* Naranjo, D., Parisi, L., Jousset, P., Weemstra, C., and Jónsson, S.: Determining OBS clock drift using ambient seismic noise , EGU General Assembly 2021, online, 19–30 Apr 2021, EGU21-13999, `<https://doi.org/10.5194/egusphere-egu21-13999>`_, 2021.
* Cornelis Weemstra, Janneke I de Laat, Arie Verdel, Pieter Smets, Systematic recovery of instrumental timing and phase errors using interferometric surface-waves retrieved from large-N seismic arrays, Geophysical Journal International, Volume 224, Issue 2, February 2021, Pages 1028–1055, `<https://doi.org/10.1093/gji/ggaa504>`_


External links
^^^^^^^^^^^^^^

* `GitHub <https://github.com/davidn182>`_
