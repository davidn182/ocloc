Tutorial
###############

:Release: |version|
:Date: |today|

----

If you are here it means that you have made a huge effort to collect passive seismic records and you want to be sure that the different stations are correctly synchronized in time. Surprisingly enough, in spite the tremendous amount of work that is done when acquiring passive seismic records, most people jump into analyzing the data before checking wether or not their waveforms have the correct time. As we know, seismic data is nothing more than a series of samples measured over time. Therefore, we better make sure to have accurate timing of our waveforms after spending significant efforts in seismic campaigns.

It is important to mention here that we do not provide the routines for computing the interferometric responses that are necesary for using ocloc. The reason is tha there are many ways of calculating interferometric responses and we would like to give some freedom to the users on how to pre-process their own data

In order to detect and correct timing errors using ambient noise seismic inteferometry, we developed an algorithm that follows the processing sheme shown in figure . The functional workflow comprises five main processing steps: Interferometric response retrieval, data filtering, time-symmetry shift measurement, time-symmetry verification, construction of the linear system of equations, and results refinement and verification.

Given that Github has a limited capacity to upload data, please request the dataset to David Naranjo (d.f.naranjohernandez@tudelft.nl).

.. toctree::
   :maxdepth: 2

   loading_data.nblink
   
