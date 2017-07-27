! Program performing analytic continuations using an average of different Pade approximants.
! Written by Johan Schott and Elin Lundin.
! Given input data zin, f(zin) and complex points zout, the main output is f(zout).

! The file pade.in has four columns: Re[zin], Im[zin], Re[f(zin)] and Im[f(zin)]
! The file pade.par contains the input parameters.

! Computational time considerations:
! The number of independent continuations is denoted by the variable ns.
! The computational task for doing the continuations scales O(ns).
! Constructing deviation matrix between different continuations scales O(ns^2).
! Both performing the continuations and the deviation matrix can be parallelized.
! But it is not done since the serial program is not very slow. 
!------------------------------------------------------------------------------------------------------

program main
use loadsavetxt
use padem
implicit none
complex(kind=16),allocatable :: zin(:),fin(:)
complex(kind=16),allocatable :: zout(:)
complex(kind=8),allocatable  :: fout(:)

! Read zin and f(zin) from file pade.in
call readinputdata(zin,fin)

! Construct zout by reading from parameter file pade.par 
call readzoutdata(zout)

! Perform the analytical continuation and return fout. 
call pade(zin,fin,zout,fout)

! Save the output in file pade_fout and the spectral function 
call savefout(zout,fout)

! Deallocate variables
if(allocated(zin))  deallocate(zin)
if(allocated(fin))  deallocate(fin)
if(allocated(zout)) deallocate(zout)
if(allocated(fout)) deallocate(fout)

write(*,'(a)') "Averaging Pade approximant program successfully finished."

end program 
