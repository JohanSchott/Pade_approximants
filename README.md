# Pade_approximants

## Aim
The purpose is to perform an analytical continuation of a (Green's) function in the complex plane.  

## What is contained?
- Both a Python and a Fortran program.
- Parameter file `pade.par` (with standard settings) for the Fortran program, in the `tests` folder.
- Test models:
    - `Sm7`
    - `betheU0`
    - `haldane_model`
    - `two-poles_wc1_dw0.5`

## Python notebook
- The Python notebook is short and easy to use.
- The analytical continuation can be done by either Beach's matrix formulation, Thiele's recursive algorithm or a nonlinear Least Square minimization.
- The matrix formulation by Beach is here calculated using double precision, which usually is too little.
- (Green's) function data somewhere in the complex plane (e.g. at Matsubara frequency points) is used as input. 
- (Green's) function data somewhere else in the complex plane is the output.

## Fortran program
- The analytical continuation uses Beach's matrix formulation.
- A more low-level control of the inversion routines in Beach's matrix formulation is possible. 
  For example, both double and quadruple precision inversion routines exist.
- The Fortran program prints more information during execution than the Python script. 

### How to use 
- A parameter input file `pade.par` has to exist in the simulation directory.
- A input file `pade.in` has to exist in the simulation directory. 
  Four columns are expected: `Re[zin], Im[zin], Re[f(zin)] and Im[f(zin)]`
  where `zin` are the input points and `f(zin)` the corresponding function values.
- Execute the binary `pade_approximants` 
- The following files will be generated:
    - `pade_info`, gives information about the performed continuations
    - `pade_fout_all`, all the Pade approximants, evaluated on the `zout` points
    - `pade_fout`, the Pade approximant average 
    - `pade_A`, first column and -1/pi times the third column of `pade_fout`

### Parameters in `pade.par` 
Below follows a description about each line in the input parameter file `pade.par`.

1)  Works as a header (or can be left empty).

2)  Settings for the output real-axis mesh.

3)  Sets the distance above the real-axis for the output mesh.

4)  Lowest index in the file `pade.in` to use for the continuations. 
    Three parameters: nminstart, nminfinish, nminstep.
    For example, input: 0 6 2 means the first, third, fifth and seventh point 
    in `pade.in` are all used as starting points for the continuations. 
    Mirror symmetry $f(z)^* = f(z^*)$ can be enforced by using mirror points, which are added with nmin < 0.

5)  Specify how many input points to use in the continuations.
    Three parameters: Mmin, Mmax, Mstep.

6)  Specify how many Pade approximant coefficients to use in the continuations.
    Three parameters: Nmin, Nmax, Nstep.

7)  Whether to only consider continuations with equally many input points as coefficients. Fortran boolean.

8)  Select analytical continuation method.
    0: Beach's matrix formulation

9)  If to print the poles of the Pade approximants. Fortran boolean.

10) The numerical precision in the inversion routine. 64: double precision, 128: quadruple precision. 

11) Select algorithm for solving matrix equation. 
    - 0: Let LAPACK solve equation $A x = b$ in LS sense
    - 1: Let LAPACK solve normal equation: $A^{\dagger} A x = A^{\dagger} b$
    - 2: explicit SVD

12) Averaging criteria parameters for the Pade approximants on the real axis.

### Compile
- Change directory to the `fortran` folder.
- Copy the example Makefile: `Makefile_example` to `Makefile` and adjust it to fit to the current machine. 
- The program requires standard LAPACK to be installed on the computer.
- Access to modified LAPACK library using quadruple precision is required. 
  Change directory to the `quad/zgels/zgels_quad` folder,
  copy the example Makefile: `Makefile_example` to `Makefile` and adjust it to fit to the current machine. 
  Then run `make` to generate the library `libzgelsquad.a`.
- Go back to the `fortran` folder and run `make`. 
  The binary `pade_approximants` should be created.

### Future possible improvements
- The parameters below are on the 2-do list to implement, but not of great importance.
    - 0          # Impose spectral symmetry. 0: no, 1: even, 2: odd
    - .false.    # Shift real part of input data by `Re[f(z_inf)]` before the analytical continuation.
- Any suggestions?

### Notes
- Comparing ZGELSD with ZGELS (using double precision), the spectra shows more features using ZGELS. 
  ZGELSD sometimes gives too smooth spectra, with many features washed out or absent.
  ZGELS assumes full rank and gives more spectral features. Usually Beach's matrix is rank deficient. 
  For rank deficient matrices it is common to also minimize the norm of the solution vector. 
  ZGELSD uses SVD. Using rcond=-1 and rcond=10^(-40) gives similarly smooth spectra for the Sm7 test model, but the estimated effective rank became very different. 
  Hence, the different output spectra seems originate from the different algorithms in ZGELS and ZGELSD. 
  In the folder `tests/Sm7/solving_Beach_system` spectra for varous setups are shown to illustrate the differences.
- Speed: For the `Sm7` test model a simululation took 44 seconds with quadruple precision in the inversion routines and 15 seconds with double precision in the inversion routines, using the Fortran program. 
The same setup with the python notebook, but using only double precision numerics, took 1.5 seconds. 
The slow performance in the Fortran case is because of a few conversions inbetween quadruple and double variables and some calculations still done in quadruple precision.
But this enables comparisons of analytical continuation spectra where only the inversion precision is varied. 
Small note, the simulation time using double precision in Fortran was independent of if used ZGELS or ZGELSD. 
- Support for `MPACK`'s arbitrary precsion is removed to simplify compilation. For high precision routines, instead use e.g. the Mathematica software.
