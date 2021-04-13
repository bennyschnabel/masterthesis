# masterthesis
Matlab/ GNU Octave implementation of the MIL tensor

## MIL 2D

## MIL 3D

1. Run the data_conversion.m file to convert the image file data to a voxel file.
2. Run the main.m file to generate the Ellipsoid with the input parameters: file name, number of different angles to scan the voxel file ...

## MIL 3D Fortran

### Calculation of the MIL tensor

1. Add the *.vtk* file to the folder
2. Open and edit the *main.f90* file
   - `fileName`: Name of the *.vtk* file 
   - `noOrientations`: Number of randomly generated orientation for every subdomain
   - `dpi`: DPI of image stack (*.vtk* file)
   - `domainSize`:  Domiain size per subvolume [mm]
3. Save the file *main.f90*
4. Calculate the MIL tensor with the following commands
   - gfortran -c kinds.f90
   - gfortran -c stringmod.f90
   - gfortran -c aux_routines.f90
   - gfortran -c calculate_mil.f90
   - gfortran -c main.f90
   - gfortran aux_routines.o kinds.o stringmod.o calculate_mil.o main.o -o main
5. The result of the calculation is saved in a file with the following structure: *xxx_xxx_M.dat*

### Calculation of the fabric tensor and the stiffness matrix

