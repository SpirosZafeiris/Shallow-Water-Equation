# Shallow-Water-Equation
A simple finite-difference program written in Modern Fortran that 
solves the linearized shallow-water equations with explicit time marching
in variable bathymetry.

The domain is a rectangular structured grid domain where a portion near the
edges is sacrificed for the PML region where a specific sink function is 
applied.

The output is in tecplot .dat format and the initial conditions are user-defined.


Basic Usage :

- git clone git@github.com:SpirosZafeiris/Shallow-Water-Equation.git
- cd Shallow-Water-Equation
- mkdir build
- cd build
- Init cmake with :
    cmake ../
- Build with :
    make 
- Run with :
    ./shallow_water_equation.out


To use a fast I/O with cgns/hdf5 to create animative fields do the following:
- cd Shallow-Water-Equation
- ./install-hdf5.sh
- ./install-cgns.sh
- in main.f90 uncomment call to write_cgns in main time loop
 

