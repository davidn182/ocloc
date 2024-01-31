.. ocloc documentation master file, created by
   sphinx-quickstart on Wed Oct 13 09:47:23 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to OCloC's documentation!
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

OCloC (OBS Clock Correction) is an open-source project dedicated to detect and correct timing errors when using passive seismic records. This Python Package is particularly relevant for projects that use Ocean-Bottom Seismometers because their lack of GPS connection can, and most certainly will, induce a timing error. Needless to say, interpretations of the subsurface might be wrong if the seismic stations contain timing errors. The aim of this project is to prevent this from hapening.

The main goal of the OCloC project is to rapidly detect and correct clock errors of seismic networks. Please note that for doing so, this package requires an understanding of the basics of seismic interferometry. If you are not familiar with this method please take a look at the Introduction section.

Our code is stored in github, the development branches, or the latest stable release can be found here.

.. _my-reference-label:

Contents
--------------------------

.. toctree::
   :maxdepth: 2
   
   installation.rst
   tutorial/tutorial.rst

Citation
^^^^^^^^
If you use this package in your work, please refer to the following works:

* Naranjo, D., Parisi, L., Jónsson, S., Jousset, P., Werthmüller, D., & Weemstra, C. (2024). Ocean Bottom Seismometer Clock Correction using Ambient Seismic Noise. Seismica, 3(1). `<https://doi.org/10.26443/seismica.v3i1.367>`_

Funding
^^^^^^^^

This work received funding from the Competitive Research Grant ZAFRAN: URF/1/4076-01-01 from KAUST granted to Sigurjon Jonsson.

This project has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement No 956965.

The seismic network and data used for the tutorials was funded by the European Community's Seventh Framework Programme under grant agreement No. 608553 (Project IMAGE). Instruments for the IMAGE project were provided by the GIPP (Geophysical Instrument pool Potsdam) and the DEPAS (German instrument pool for amphibian seismology), and can be retrieved from http://www.image-fp7.fr/Pages/default.aspx.

External links
^^^^^^^^^^^^^^

* `OCloC's GitHub <https://github.com/davidn182/ocloc>`_
