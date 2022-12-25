# Shallow-Water-Equation
A simple finite-difference program written in Modern Fortran that 
solves the linearized shallow-water equations with explicit time marching
in variable bathymetry.

The domain is a rectangular structured grid domain where a portion near the
edges is sacrificed for the PML region where a specific sink function is 
applied.

The output is in tecplot .dat format and the initial conditions are user-defined.

