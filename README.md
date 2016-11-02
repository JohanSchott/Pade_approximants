# Pade_approximants

## Aim
The purpose is to perform an analytical continuation of a (Green's) function in the complex plane.  

## What is contained?
- Both a Python and a Fortran program.
- Parameter file `pade.par` (with standard settings) for the Fortran program, in the `tests` folder.
- Test models:
    - `betheU0`
    - `betheU4`
    - `two-poles_wc1_dw0.5`
    - `Sm7`

## Python script
- The Python notebook script is short and easy to use.
- The analytical continuation can be done by either Beach's matrix formulation, Thiele's recursive algorithm or a nonlinear Least Square minimization.
- The matrix formulation by Beach is here calculated using double precision, which usually is too little.
- (Green's) function data somewhere in the complex plane (e.g. at Matsubara frequency points) is used as input. 
- (Green's) function data somewhere else in the complex plane is the output.

## Fortran program
- The Fortran program using modified LAPACK routines to work with quadruple numerical precision.
- The Fortran prints more information during execution than the Python script. 
  It offers a more low level control of the inversion routines in the Beach algorithm.

### How to use 
- A parameter input file `pade.par` has to exist in the current/simulation directory.
- A input file `pade.in`has to exist in the current/simulation directory.
- Execute the binary `Pade_approximants` 

### Compile
- Move to the `fortran` folder.
- Copy the example Makefile: `Makefile_example` to `Makefile` and adjust it to fit to the current machine. 
- Then run: `make` and the binary `Pade_approximants` should be created.
- Access to modified LAPACK libraries using quadruple precision is required. It is stored in the `quad` folder.
- The program requires the machine has LAPACK installed somewhere. 

### Future improvements
- The parameters below are on the 2-do list to implement, but not of great importance.
0          # Impose spectral symmetry. 0: no, 1: even, 2: odd
.false.    # Shift real part on Matsubara axis before continuation starts to remove Re[f(z_inf)]
.false.    # Which inversion routines to use: .false.: lapack, .true.: mpack c++.

- Parallelize with MPI. Maybe not needed, fast enough...

### Notes
- Comparing ZGELSD with ZGELS (using double precision), the spectra shows more features using ZGELS. 
  ZGELSD gives too smooth spectra, with many features washed out or absent.
  ZGELS assumes full rank but usually we work will rank deficient cases, where it's common to also minimize the norm of the solution vector. 
  ZGELSD uses SVD and estimates the rank which can become very low, hence the smooth behavior. 
- Support for `MPACK`'s arbitrary precsion is removed to simplify compilation. For high precision routines, instead use e.g. the Mathematica software.
