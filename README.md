# OCloC
[![License](https://img.shields.io/badge/License-Apache%202.0-yellowgreen.svg)](https://opensource.org/licenses/Apache-2.0)  

`OCloC` is an open-source Python-based package to correct clock errors of ocean bottom seismometers.

## Citation:

Naranjo, D., Parisi, L., Jónsson, S., Jousset, P., Werthmüller, D., & Weemstra, C. (2024). Ocean Bottom Seismometer Clock Correction using Ambient Seismic Noise. *Seismica*, 3(1).[https://doi.org/10.26443/seismica.v3i1.367](https://doi.org/10.26443/seismica.v3i1.367)

> [!note]
> This project is under active development.

## Installation
OCloC works in OS X & Linux systems. 

Installation instructions are provided at https://ocloc.readthedocs.io/en/latest/index.html



### Side note for Windows systems:
As some subroutines were written in Fortran it can also be run on Windows but it will be necessary to compile the codes in the directory ./ocloc/recover_timing_errors-master/ by using the makefile. Note that you will need to have SAC installed in your computer.



## Contributing

1. Fork it on your own github (<https://github.com/davidn182/ocloc/fork>)
2. Create your feature branch (`git checkout -b feature/new_feature`)
3. Commit your changes (`git commit -am 'Add some new feature'`)
4. Push to the branch (`git push origin feature/new_feature`)
5. Create a new Pull Request

### Funding
This project has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement No 956965.

This work received funding from the Competitive Research Grant ZAFRAN: URF/1/4076-01-01 from KAUST granted to Sigurjon Jonsson.

The seismic network and data used for the tutorials was funded by the European Community's Seventh Framework Programme under grant agreement No. 608553 (Project IMAGE). Instruments for the IMAGE project were provided by the GIPP (Geophysical Instrument pool Potsdam) and the DEPAS (German instrument pool for amphibian seismology), and can be retrieved from http://www.image-fp7.fr/Pages/default.aspx.
