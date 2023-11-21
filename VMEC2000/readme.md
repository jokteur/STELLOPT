# How to build and install VMEC

There are two ways to compile VMEC:
- Follow the compilation guide on https://princetonuniversity.github.io/STELLOPT/STELLOPT%20Compilation ; You will probably need to modify one of the files in the `SHARE` folder of this git. (not tested on this version of VMEC)
- Compile with CMake [RECOMMENDED]

## Compiling with CMake
Make sure that you have a recent installation of CMake on your computer. 

### Libraries 
VMEC requires the following libraries:
- BLAS / LaPACK
- netcdf (C and Fortran libs)
- hdf5 (> 1.10), with HL, MPI and Fortran enabled
- mpi (will be optional in the future)

You can find different guides on how to install the libraries on your system:
- [Apple](https://princetonuniversity.github.io/STELLOPT/STELLOPT%20Compilation%20OSX)
- [Cineca/Marconi](https://princetonuniversity.github.io/STELLOPT/STELLOPT%20Compilation%20CINECA)
- [Redhat/Centos](https://princetonuniversity.github.io/STELLOPT/STELLOPT%20Compilation%20CentOS)
- [Ubuntu/Debian](https://princetonuniversity.github.io/STELLOPT/STELLOPT%20Compilation%20Ubuntu)

An alternative way is to use the [spack](https://spack.io/) package manager. It allows to easily compile and install specifics libraries for your system. Once installed on your system, you can:
```bash
spack install netcdf-fortran
spack install lapackpp
```

When compiling and using the library, you can simply
```bash
spack load netcdf-fortran
spack load lapackpp
```
or create your own environnement.

### Compiling
Create a build folder in the root of the git project:
```bash
mkdir build
cd build
```

Then configure the CMake project:
```bash
cmake ..
```
You can optionally activate the ANIMEC (anisotropy) or FLOW version of VMEC by adding `-DANIMEC` or `-DFLOW` respectively. It is only possible to chose one or the other, not both.
Example:
```bash
cmake .. -DFLOW
```

Once CMake is set-up, you can compile the project:
```bash
make -j
```

The binaries will be in the `build/VMEC2000` folder.