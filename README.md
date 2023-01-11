# MPSTransferFunction Example


This folder contains example code for applying a MPS transfer function to measured data. 
The method includes the representation of MPS signal responses in time- and frequency domain, as well as the MPI point-spread function and the hysteresis curve.

The method is described in the associated publication

F. Thieben, T. Knopp, M. Boberg, F. Foerger, M. Graeser, and M. Möddel. (2022) On the receive path calibration of magnetic particle imaging
systems. IEEE Transactions on Instrumentation and Measurement, 2022. doi: [10.1109/TIM.2022.3219461](https://ieeexplore.ieee.org/document/9939022).

## Installation

In order to use this code one first has to download [Julia](https://julialang.org/) (version 1.8 or later), clone this repository and navigate to the folder in the command line. The example scripts automatically activate the environment and install all necessary packages.

## Execution
After installation the example code can be executed by running `julia` and entering
```julia
include("example.jl")
```
for the initial method.


## Open MPI Data

The measurement data associated to this project will be downloaded and stored automatically, when the code is executed for the first time.
It is published under a [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/legalcode) license and can be found here:


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7457743.svg)](https://doi.org/10.5281/zenodo.7457743)
