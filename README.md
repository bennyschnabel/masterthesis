# masterthesis
Matlab/ GNU Octave implementation of the MIL tensor

## MIL 2D

## MIL 3D

1. Run the data_conversion.m file to convert the image file data to a voxel file.
2. Run the main.m file to generate the Ellipsoid with the input parameters: file name, number of different angles to scan the voxel file ...

## MIL 3D Fortran

1. Add .vtk-file to folder
2. Calculate MIL Tensor with the following commands
   - gfortran -c kinds.f90
   - gfortran -c stringmod.f90
   - gfortran -c aux_routines.f90
   - gfortran -c calculate_mil.f90
   - gfortran -c main.f90
   - gfortran aux_routines.o kinds.o stringmod.o calculate_mil.o main.o -o main
3. Run main.m file
