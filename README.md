**iS3D2 (c) Mike McNelis, Derek Everett, Sameed Pervaiz, Matthew Golden and Lipei Du**

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

    sh particlization.sh

The results from the simulation are stored in `results`.

To parallelize the continuous Cooper-Frye formula on OpenMP, edit the makefile to use the icpc compiler with the -qopenmp flag. Then compile and run iS3D on `X` threads by doing

    sh cleanMakeCPU.sh
    sh runCPU.sh X

Note: the particle sampler does not use OpenMP acceleration. 


## Freezeout surface

The freezeout surface file in `input/surface.dat` must have one of the following formats to be read in correctly:

    1 = VAH or CPU VH
    5 = CPU VH (w/ thermal vorticity)
    6 = MUSIC (public)
    7 = HIC-EventGen (or VISHNU)

The code expects to read in the file format set by the parameter `mode`.

Note: I axed several file formats used by the previous version.

During runtime, the code prints the columns required by the freezeout surface reader. If your columns do not match, the code will probably crash. 

If you turn on `include_baryon`, the reader assumes there are additional columns related to net baryon charge.

