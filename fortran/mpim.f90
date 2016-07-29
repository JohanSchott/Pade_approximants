module mpim
implicit none

contains

subroutine distjobs(N,a,b,nt)
implicit none
include "mpif.h"
integer,intent(in) :: N
integer,intent(out) :: a,b,nt
!dummy
integer :: nr,r
!MPI variables
integer myrank,ranks,ierr
integer status(MPI_STATUS_SIZE) 
call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)
call MPI_Comm_size(MPI_COMM_WORLD,ranks,ierr)

nr=floor(1d0*N/ranks)
r=N-nr*ranks
if(myrank<ranks-r) then
   a=nr*myrank+1
   b=nr*(myrank+1)
   nt=nr
else
   a=N-(ranks-1-myrank)*(nr+1)-nr
   b=N-(ranks-1-myrank)*(nr+1)
   nt=nr+1
endif
!write(*,*) "rank:",myrank," has a=",a,", b=",b," and nt=",nt
end subroutine

end module 
