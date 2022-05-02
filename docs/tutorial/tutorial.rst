Tutorial
###############

:Release: |version|
:Date: |today|

----

If you are here it means that you have made a huge effort to collect passive seismic records and you want to be sure that the different stations are correctly synchronized in time. Surprisingly enough, in spite the tremendous amount of work that is done when acquiring passive seismic records, most people jump into analyzing the data before checking wether or not their waveforms have the correct time. As we know, seismic data is nothing more than a series of samples measured over time. Therefore, we better make sure to have accurate timing of our waveforms after spending significant efforts in seismic campaigns.

This introduction section is a simple guide that will cover the basics of seismic interferometry and how to use it for detecting clock errors. If you are already familiar with the principles behind seismic interferometry please skip the introduction section.

It is important to mention here that we do not provide the routines for computing the interferometric responses that are necesary for using ocloc. The reason is tha there are many ways of calculating interferometric responses and we would like to give some freedom to the users on how to pre-process their own data

.. toctree::
   :maxdepth: 2

   0. Preface
   introduction.md
   getting_started.md
   data_requirements.md
   
