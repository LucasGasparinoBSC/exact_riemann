# exact_riemann

## Description

Exact Riemann solver based on the code of E. Toro. Given 2 discontinuous data states, will find the pressure and velocity in the star region, then create a sampling of the 3 main variables at a given time.

## Usage

Initial data has to be modified inside the coode for now. Other data is modified in mod_constants. Program will output a .out file containing the samples,
as well as print the star values at the end of the iterative process.

## Compilation

Built using CMake and any fortran compiler. Create a buildd directory, then cmake [options] ..

