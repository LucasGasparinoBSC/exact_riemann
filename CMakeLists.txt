# Add cmake mminnimum required
cmake_minimum_required(VERSION 3.15)

# Add project name, version and language
project(exact_riemann VERSION 1.0 LANGUAGES Fortran)

# Add files for compilation
file(GLOB sources src/*.f90)

# Add executable name
add_executable(exact_riemann.x ${sources})