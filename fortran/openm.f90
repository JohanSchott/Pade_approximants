module openm
!This module contains interface openf.
!openf open file and save data to variables
!If first row starts with symbol: # , it is treated as a comment. 
use matlab

implicit none

interface openf
module procedure open1,open2,open3,open4
end interface

contains

integer function readstart(filen)
implicit none
character(len=800),intent(in) :: filen
character(len=800)  :: s
open(44,file=trim(filen))
read(44,'(a)') s
close(44)
if(s(1:1)=="#") then
   readstart=2
else
   readstart=1
endif
end function

subroutine open1(filen,x)
implicit none
character(len=800),intent(in) :: filen
real(kind=16),allocatable,intent(out) :: x(:)
integer :: k,i,j,ok
k=readstart(filen)
i=filelen(filen)
if(allocated(x)) deallocate(x)
allocate(x(i-k+1))
open(31,file=trim(filen),iostat=ok)
do j=1,k-1
   read(31,*)
enddo
do j=1,i-k+1
   read(31,*) x(j)
enddo
close(31)
end subroutine

subroutine open2(filen,x,y)
implicit none
character(len=800),intent(in) :: filen
real(kind=16),allocatable,intent(out) :: x(:),y(:)
integer :: k,i,j,ok
k=readstart(filen)
i=filelen(filen)
if(allocated(x)) deallocate(x)
if(allocated(y)) deallocate(y)
allocate(x(i-k+1),y(i-k+1))
open(31,file=trim(filen),iostat=ok)
do j=1,k-1
   read(31,*)
enddo
do j=1,i-k+1
   read(31,*) x(j),y(j)
enddo
close(31)
end subroutine

subroutine open3(filen,x,y,z)
implicit none
character(len=800),intent(in) :: filen
real(kind=16),allocatable,intent(out) :: x(:),y(:),z(:)
integer :: k,i,j,ok
k=readstart(filen)
i=filelen(filen)
if(allocated(x)) deallocate(x)
if(allocated(y)) deallocate(y)
if(allocated(z)) deallocate(z)
allocate(x(i-k+1),y(i-k+1),z(i-k+1))
open(31,file=trim(filen),iostat=ok)
do j=1,k-1
   read(31,*)
enddo
do j=1,i-k+1
   read(31,*) x(j),y(j),z(j)
enddo
close(31)
end subroutine

subroutine open4(filen,x,y,z,a)
implicit none
character(len=800),intent(in) :: filen
real(kind=16),allocatable,intent(out) :: x(:),y(:),z(:),a(:)
integer :: k,i,j,ok
k=readstart(filen)
i=filelen(filen)
if(allocated(x)) deallocate(x)
if(allocated(y)) deallocate(y)
if(allocated(z)) deallocate(z)
if(allocated(a)) deallocate(a)
allocate(x(i-k+1),y(i-k+1),z(i-k+1),a(i-k+1))
open(31,file=trim(filen),iostat=ok)
do j=1,k-1
   read(31,*)
enddo
do j=1,i-k+1
   read(31,*) x(j),y(j),z(j),a(j)
enddo
close(31)
end subroutine

end module 
