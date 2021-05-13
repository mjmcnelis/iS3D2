**iS3D2 (c) Derek Everett, Mike McNelis, Sameed Pervaiz, Matthew Golden and Lipei Du**

Created on 2/13/2018 by Derek Everett\
Last edited on 5/13/2021 by Mike McNelis


## Summary
A Monte Carlo simulation that samples hadrons from the particlization stage of heavy-ion collisions. 

This repository is my updated version of the [iS3D](https://github.com/derekeverett/iS3D) code, which was developed from the particle sampler [iSS](https://github.com/chunshen1987/iSS).

The code reads in a freezeout surface from the preceding fluid dynamic stage and samples particles from the Cooper-Frye formula using one of five df corrections:

    1) Grad 14-moment approximation
    2) RTA Chapman-Enskog expansion
    3) PTM equilibrium distribution
    4) PTB equilibrium distribution
    5) PTMA anisotropic distribution

The shear stress, bulk pressure and baryon chemical potential / diffusion can be turned on or off during runtime.

One can also integrate the Cooper-Frye formula to obtain the continuous momentum spectra or spacetime distributions. 


## References

If you use this code, please cite the following papers:

    M. McNelis, D. Everett and U. Heinz, Comput. Phys. Commun. 228 (2021) 107604
    C. Shen, Z. Qiu, H. Song, J. Bernhard, S. Bass and U. Heinz, Comput. Phys. Commun. 199 (2016) 61-85


## Running the code

The default makefile in `src/cpp/GNUmakefile` uses the g++ compiler, which is what we use for particle sampling.

To compile and run iS3D, do either

    sh cleanMakeCPU.sh
    ./iS3D.e
    
or

    sh particlization.sh

The results from the simulation are stored in `results`.

To parallelize the continuous Cooper-Frye formula on OpenMP, edit the makefile to use the icpc compiler with the -qopenmp flag. Then compile and run iS3D on `s` threads by doing

    sh cleanMakeCPU.sh
    sh runCPU.sh s

Note that the particle sampler does not use OpenMP acceleration. 




## Freezeout surface

The input freezeout surface file in `input/surface.dat` must have one of the three formats to be read in correct

The 


Turning on the will require







This code can read in a freeze out surface from 3+1D viscous hydro or anisotropic
viscous hydro and calculate 3D smooth particle spectra or a sampled particle list.
The structure is based on iSpectra, the Cooper Frye code in the iEBE heavy ion
event generator (Chun Shen, Zhi Qiu).

to build iS3D, one can do

$ mkdir build && cd build
$ cmake ..
$ make
$ make install

To run iS3D, do

$ ./iS3D

or

$ sh runCPU.sh num_threads

where num_threads is the number of cpu threads.

The freezeout surface is read from input/surface.dat, or from memory depending on how the wrapper is called.
By default input/surface.dat contains a toy freezeout surface with one cell.
See parameters.dat for a list of compatible formats.

The results will be written in the results/ directory, so this directory must exist at runtime.
