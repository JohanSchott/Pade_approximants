module diffm
!Module config contains c and Nv, which are the only implicit input variables to this module

use config
use mpim
!use MPIroutines
implicit none

contains

subroutine getdiffmatrix(ind,emesh,dm)
!Compares all continuations in ind with each other, returning matrix dm
!Input variables ind and emesh are defined only for master so need to copy to slaves
implicit none

include "mpif.h"
integer,allocatable :: ind(:) !list with indices to continuations to consider
complex(kind=16) :: emesh(:)  !points just above real axis
real(kind=16) :: dm(:,:) !diff matrix, matrix representation of vector dv
!internal
integer :: ns !nbr of continuations
integer :: good !nbr of continuations to consider
integer :: N !nbr of tasks
integer :: nt !nbr of tasks for a rank. It's individual
!integer :: nemesh !nbr of emesh points
real(kind=16),allocatable :: dr(:),dv(:) !diff vector. Total diff vector
!dummy
complex(kind=16),allocatable :: f1(:),f2(:)
real(kind=16),allocatable :: x(:),y(:)
integer,allocatable :: iv(:),jv(:)
integer :: nr,r,a,b
integer :: i,j,k

!MPI variables
integer myrank,ranks,ierr
integer tag
integer status(MPI_STATUS_SIZE) 

call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)
call MPI_Comm_size(MPI_COMM_WORLD,ranks,ierr)

tag=myrank
ns=size(Nv,1) !since Nv is only allocated for master rank, only master gets correct ns value here
!copying data to slaves
call MPI_BCAST(ns,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
good=size(ind,1) !since ind is only allocated for master rank, only master gets correct good value here
call MPI_BCAST(good,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if(myrank /= 0 ) then
   allocate(ind(good))
endif
call MPI_BCAST(ind,good,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!Copy Nv and c to all slaves
if(myrank/=0) then
   allocate(Nv(ns))
   allocate(c(ns))
endif
call MPI_BCAST(Nv,ns,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if(myrank==0) then
   do i=1,ranks-1
      do j=1,ns
         allocate(x(Nv(j)),y(Nv(j)))
         x=real(c(j)%x(:))
         y=aimag(c(j)%x(:))
         call MPI_SEND(x,Nv(j),MPI_REAL16,i,j,MPI_COMM_WORLD,ierr) 
         call MPI_SEND(y,Nv(j),MPI_REAL16,i,ns+j,MPI_COMM_WORLD,ierr) 
         deallocate(x,y)
      enddo
   enddo
else
   !loop over all continuations, allocating and receiving Pade coefficients
   do j=1,ns
      allocate(c(j)%x(Nv(j)),x(Nv(j)),y(Nv(j)) )
      call MPI_RECV(x,Nv(j),MPI_REAL16,0,j,MPI_COMM_WORLD,status,ierr) !recive Nv from master
      call MPI_RECV(y,Nv(j),MPI_REAL16,0,ns+j,MPI_COMM_WORLD,status,ierr) !recive Nv from master
      c(j)%x=cmplx(x,y,kind=16)
      deallocate(x,y)
   enddo
endif
!distribute N tasks over all processes
N=good*(good-1)/2 !nbr of tasks (nbr of elements the lower part of triangle of dimension (good,good))
call distjobs(N,a,b,nt)
!Calculate corresponding matrix elements for each task
allocate(iv(nt),jv(nt))
k=0
r=0
!do j=1,ns
do j=1,good
   do i=j+1,good
      k=k+1
      if(a<=k .and. k<=b) then
         r=r+1
         iv(r)=i
         jv(r)=j
         !dr(r)= (i,j)
      endif
   enddo
enddo
!Calculate tasks for every process
allocate(dr(nt))
allocate(f1(size(emesh,1)),f2(size(emesh,1)))
do k=1,nt
   if(k==1) then
      call eval_pade(c(ind(iv(k)))%x,emesh,f1)
      call eval_pade(c(ind(jv(k)))%x,emesh,f2)
   elseif(iv(k)==iv(k-1)) then
      call eval_pade(c(ind(jv(k)))%x,emesh,f2)
   elseif(jv(k)==jv(k-1)) then
      call eval_pade(c(ind(iv(k)))%x,emesh,f1)
   else
      call eval_pade(c(ind(iv(k)))%x,emesh,f1)
      call eval_pade(c(ind(jv(k)))%x,emesh,f2)
   endif
   dr(k)=padediff(f1,f2)
enddo
if(allocated(f1)) deallocate(f1)
if(allocated(f2)) deallocate(f2)
!All ranks except master sends their version of dr to the master rank
if(myrank==0) then
   allocate(dv(N))
   dv(a:b)=dr  !rank 0 puts its information to big vector
   do i=1,ranks-1 !loop over all slave processes
      call MPI_RECV(j,1,MPI_INTEGER,i,i,MPI_COMM_WORLD,status,ierr)
      allocate(x(j))
      call MPI_RECV(x,j,MPI_REAL16,i,ranks+i,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(k,1,MPI_INTEGER,i,2*ranks+i,MPI_COMM_WORLD,status,ierr)
      dv(k:(k+j-1))=x
      deallocate(x)
   enddo
   !write(*,*) "dv="
   !do i=1,N
   !   write(*,*) ,dv(i)
   !enddo

   !allocate(dm(good,good))
   !write(*,*) "rank:",myrank
   dm=0q0
   k=0
   do j=1,good-1
      do i=j+1,good
         k=k+1
         dm(i,j)=dv(k)
      enddo
   enddo
   deallocate(dv)
else
   call MPI_SEND(nt,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,ierr) 
   call MPI_SEND(dr,nt,MPI_REAL16,0,ranks+tag,MPI_COMM_WORLD,ierr) 
   call MPI_SEND(a,1,MPI_INTEGER,0,2*ranks+tag,MPI_COMM_WORLD,ierr) 
endif
deallocate(dr)
if(myrank/=0) then
   deallocate(Nv,c)
endif

end subroutine

end module

