![logo](doc/flecsph_logo_bg.png)

# Kelvin-Helmholtz instability test

This test implements the Kelvin-Helmholtz instability test in 2 and 3D

## Generating initial data
Use the generator as follows:

    % mpirun -np X ./kh_generator kh_sqn100.par

This will produce an h5part file named `h5part_kh.h5part` to be input into
the evolution code.

## Running the evoluiton app

    % mpirun -np X ./KH

As long as FleCSI does not provide a way to read the file name, it is hardcoded
in the main_driver.cc file to be `h5part_kh.h5part`.

The current version generates a file `output_kh.h5part` with the result. 

## Visualize resutls

Paraview can be used after loading the h5part module.
