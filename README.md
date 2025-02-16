# MLIP

MLIP is a software for Machine Learning Interatomic Potentials.
It has been developed at Skoltech (Moscow) by
Alexander Shapeev, Evgeny Podryabinkin, Konstantin Gubaev, and Ivan Novikov

## Licence
See [LICENSE](LICENSE)

## Prerequisties
* g++, gcc, gfortran, mpicxx
* Alterlatively: the same set for Intel compilers (tested with the 2017 version)
* make

## Compile
For full instructions see [INSTALL.md](INSTALL.md) or [INSTALL-cmake.md](INSTALL-cmake.md).

You might also be interested in LAMMPS-MLIP interface distributed here:
[https://gitlab.com/ashapeev/interface-lammps-mlip-2](https://gitlab.com/ashapeev/interface-lammps-mlip-2)


## Getting Started

Have a look at `doc/manual/` first.

Check the usage examples at `test/examples/`

Go to the repository [https://gitlab.com/ashapeev/mlip-2-tutorials](https://gitlab.com/ashapeev/mlip-2-tutorials) with tutorials

## Have questions?

Note that we'll not be able to answer all of your questions.
As a rule, we are supporting only the documented functionality of MLIP.
If you think you found a bug or an inconsistency in the documentation or usage examples,
please create a Gitlab.com issue.  

-----------------------------------------------------------------------------------
# MTP_AER
This is a modified [MLIP-2](https://gitlab.com/ashapeev/mlip-2/-/tree/master?ref_type=heads) package supporting atomic energy regularazition (AER).
AER is defined as:
![image](https://github.com/user-attachments/assets/a2d3a5d2-f6be-471a-a55c-dfa4c6632a19)  
AER builds up the relationship between local decompostion and global atomic-environmental information,thereby effectively constraining the atomic energy distribution. By considering AER, we reformulate the model optimization of MTP as: 
![image](https://github.com/user-attachments/assets/36616404-0fb6-4b24-b880-44ce747651e9)
The 𝜆_1 and 𝜆_2 are non-negative weights. In program, ther are determined by following command:
```bash
mlp train untrained.mtp trainset.cfg --std-weight 𝜆_1 --stdd-weight 𝜆_2 ...
```
The global term of AER is determined by chemical potential of simple substance.    
![image](https://github.com/user-attachments/assets/4beb1170-7af2-4c03-8571-1d8448d42408)
μ_t is calculated by DFT and read by program from .mtp file, like this:
![image](https://github.com/user-attachments/assets/ab84bdcc-019f-44c2-b142-12cfe26d1484)
If item "mu_ " is missed in .mtp file, the default values are set to 0.
