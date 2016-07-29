module savem

implicit none

interface save2f
module procedure save1r2f,save1c2f,save2r2f,save2c2f,save3r2f,save4r2f,save5r2f
end interface

contains

subroutine save1r2f(filen,x)
implicit none
character(len=800),intent(in) :: filen
real(kind=16),intent(in) :: x(:)
integer :: i,ok
open(31,file=trim(filen),iostat=ok)
do i=1,size(x,1)
   write(31,'(E13.4)') x(i)
enddo
close(31)
end subroutine

subroutine save1c2f(filen,x)
implicit none
character(len=800),intent(in) :: filen
complex(kind=16),intent(in) :: x(:)
integer :: i,ok
open(31,file=trim(filen),iostat=ok)
do i=1,size(x,1)
   write(31,'(2E13.4)') real(x(i)),aimag(x(i))
enddo
close(31)
end subroutine

subroutine save2r2f(filen,x,y)
implicit none
character(len=800),intent(in) :: filen
real(kind=16),intent(in) :: x(:),y(:)
integer :: i,ok
open(31,file=trim(filen),iostat=ok)
do i=1,size(x,1)
   write(31,'(2E13.4)') x(i),y(i)
enddo
close(31)
end subroutine

subroutine save2c2f(filen,x,y)
implicit none
character(len=800),intent(in) :: filen
complex(kind=16),intent(in) :: x(:),y(:)
integer :: i,ok
open(31,file=trim(filen),iostat=ok)
do i=1,size(x,1)
   write(31,'(4E13.4)') real(x(i)),aimag(x(i)),real(y(i)),aimag(y(i))
enddo
close(31)
end subroutine

subroutine save3r2f(filen,x,y,z)
implicit none
character(len=800),intent(in) :: filen
real(kind=16),intent(in) :: x(:),y(:),z(:)
integer :: i,ok
open(31,file=trim(filen),iostat=ok)
do i=1,size(x,1)
   write(31,'(3E13.4)') x(i),y(i),z(i)
enddo
close(31)
end subroutine

subroutine save4r2f(filen,x,y,z,g)
implicit none
character(len=800),intent(in) :: filen
real(kind=16),intent(in) :: x(:),y(:),z(:),g(:)
integer :: i,ok
open(31,file=trim(filen),iostat=ok)
do i=1,size(x,1)
   write(31,'(4E13.4)') x(i),y(i),z(i),g(i)
enddo
close(31)
end subroutine

subroutine save5r2f(filen,x,y,z,g,h)
implicit none
character(len=800),intent(in) :: filen
real(kind=16),intent(in) :: x(:),y(:),z(:),g(:),h(:)
integer :: i,ok
open(31,file=trim(filen),iostat=ok)
do i=1,size(x,1)
   write(31,'(5E13.4)') x(i),y(i),z(i),g(i),h(i)
enddo
close(31)
end subroutine

end module 
