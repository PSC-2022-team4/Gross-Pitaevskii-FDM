# Gross-Pitaevskii-FDM
- Overview 
  - Equation
  - Physics 
- Installation 
  - Dependencies 
  - Build  
- Execution 
  - Generate configuration file 
  - Execute single condition 
  - Execute multiple conditions with MPI
- Example 

## Overview
### Equation
The time-dependent Gross-Pitaevskii equation is 
$$𝑖ℏ \frac{𝜕Ψ(𝒓, 𝑡)}{𝜕𝑡}=(−\frac{ ℏ^2}{2𝑚} ∇^2+𝑉(𝒓)+𝑔 |Ψ(𝒓,𝑡)|^2 )Ψ(𝒓, 𝑡)$$
g parameter is correlated to interactions of particles.  
In this project, the equation is solved numerically with various initial condition, potential values, and g values in a parallel manner.   


### Physics 
Bose-Einstein condensate(BEC) is a state of matter of a dilute gas of bosons cooled to temperatures very close to absolute zero.   
Under such conditions, bosons condensate to the same ground state.  
The Gross-Pitaevskii equation is an approximation model of BEC.   
In Hartree-Fock approximation, many body equation is turned to one body equation. Also, in the partial-wave analysis, scattering process between each bosons are approximated by the s-wave scattering.  
** TODO sample image

## Installation 
### Dependencies 
- C++ 17 
- [Cmake](https://cmake.org/) version >= 3.5
- [MPICH2](https://github.com/pmodels/mpich) 
- [Nvidia CUDA toolkit](https://github.com/NVIDIA/cuda-samples) version == 10.1
- [gtest](https://github.com/google/googletest)

### Build 
The build configuration is created with CMake.   
To create build configuration in this project, create a build directory in the project directory.  
Then, in the build directory, use cmake command to create build configuraiton as below. 
```
mkdir build 
cd build 
cmake ..
```

To build, execute:
```
make
```

Check that if the two files are generated.  
```
GrossPitaevskiiFDM_run  
GrossPitaevskiiFDM_tst   
```

## Execution 
### Generate configuration file 
In the inputs directory, there are three example configuration files.
```
benchmarking.inp
example_single.inp
example_sweep.inp
```
If you want to run the codes in a single condition(initial condition, potential, parameter g), use example_single.inp as a template.   
There are following parameters to set. 
- Main Configuration 
- Domain Configuration 
- Initial Condition Configuration 
- Equation Configuration 
- Solver Configuration 
 
For the details, refer to the code decription. #TODO  
Especially, if you want to use GPU to accelerate the execution, set parallel as true and set cuda_device = {device number} in Solver Configuration.   
For example, 
```
  parallel = true
   ...
  cuda_device = 0 
```

### Execute single condition 
To execute the program, in the build directory run ./GrossPitaevskiiFDM_run with the generated configurated file. 
```
{PATH_TO_PROGRAM}/GrossPitaevskiiFDM_run {PATH_TO_INPUT_CONFIGURATION}.inp
```
For example, 
```
./build/GrossPitaevskiiFDM_run ./inputs/example_single.inp
```

### Execute multiple conditions with MPI

To execute the program, in the build directory run ./GrossPitaevskiiFDM_run with the generated configurated file. 
If you want to execute the program with various condition parallely, use MPI.   
You have to run the program with mpiexec and set the {PROCESSOR_NUMBER} which is the number of processors that you want to use.  

```
mpiexec -np {PROCESSOR_NUMBER} {PATH_TO_PROGRAM}/GrossPitaevskiiFDM_run {PATH_TO_INPUT_CONFIGURATION}.inp
```
For example, 
```
mpiexec -np 4 ./build/GrossPitaevskiiFDM_run ./inputs/example_sweep.inp
```
