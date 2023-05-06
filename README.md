# Shallow-Water-Equation
A simple finite-difference program written in Modern Fortran that 
solves the linearized shallow-water equations with explicit time marching
in variable bathymetry.

The domain is a rectangular structured grid domain where a portion near the
edges is sacrificed for the PML region where a specific sink function is 
applied.

The output is in tecplot .dat format and the initial conditions are user-defined.


## Build

* Requires intel compilers 

1.  Clone the repository
    
    ```bash
    git clone git@github.com:SpirosZafeiris/Shallow-Water-Equation.git
    ```
    ```bash
    cd Shallow-Water-Equation
    ```
    ```bash
    mkdir build
    ```
    ```bash
    cd build/
    ```
2.  Initiate cmake:
    ```bash
    cmake ../
    ```
3. Build with Make:
    ```bash
    make
    ```

## Run

1. Run with:
    ```bash
    ./shallow_water_equation.out
    ```

## Advanced

1. To use a fast I/O with cgns/hdf5 to create animative fields do the following:
    ```bash
    cd Shallow-Water-Equation
    ./install-hdf5.sh
    ./install-cgns.sh
    ```
2. In main.f90 uncomment call to write_cgns in main time loop
 
