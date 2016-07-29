# Pade_approximants

## Aim
The purpose is to perform analytical continuations of Green's functions.  

## What is contained?
- Both a Python and a Fortran program.

## Usability
- The Python script is short and easy to use but only support double precision calculations. 
- The Fortran program uses LAPACK routines which are modified to quadruple precision or if you want, MPACK's arbitrary precision routines.

## How to use 

### Python script
- First an input file has to be read, Green's function data somewhere in the complex plane (e.g. at Matsubara frequency points). 
- Output in the Green's function somewhere in the complex plane.
- The code is short and written so that a user easily can modify it for her/his purpose. 

### Fortran program
- A input file `pade.inp` has to exist in the current/simulation directory.
- Execute the binary `Pade_approximants` 

## Comments towards users
- The Fortran code has more features and does slightly different things than the Python script.

## Future improvements
- Make the compilation with MPACK easier.
- Consider off-diagonal Green's functions and self-energies. Especially for negative Matsubara points.

