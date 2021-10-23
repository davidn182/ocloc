Tutorial
###############

The main modelling routine is :func:`empymod.model.bipole`, which can compute
the electromagnetic frequency- or time-domain field due to arbitrary finite
electric or magnetic bipole sources, measured by arbitrary finite electric or
magnetic bipole receivers. The model is defined by horizontal resistivity and
anisotropy, horizontal and vertical electric permittivities and horizontal and
vertical magnetic permeabilities. By default, the electromagnetic response is
normalized to source and receiver of 1 m length, and source strength of 1 A.


Basic example
-------------

A simple frequency-domain example, where we keep most of the parameters left at
the default value:
