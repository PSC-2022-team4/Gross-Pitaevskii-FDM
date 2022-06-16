# Gross-Pitaevskii-FDM
- Overview 
  - Equation
  - Physics 
- Installation 
  - Dependencies 
  - Build  
- Execution 
  - Generate input file 
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
