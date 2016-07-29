module asym
use leastsquare
use openm
implicit none

contains

function getrealasym(fn,n,q)
!fit first,x, vs second,y, column in file fn with: y(x)=-a/x^2+b
!Let's say n_tot is the number of points in file fn.
!Fit to the last min(n,q*n_tot) points
!Return a and b as [a,b] = getrealasym(:)
implicit none
character(len=800),intent(in) :: fn
integer,intent(in) :: n
real(kind=16),intent(in) :: q

real(kind=16) :: getrealasym(2)

real(kind=16),allocatable :: x(:),y(:),k(:,:),f(:)
integer :: i,m,ntot

ntot=filelen(fn)
m=min(real(n,kind=16),ntot*q)
write(*,*) "Number of Matsubara points used for finding asymptotic:",m
call openf(fn,x,y)
allocate(f(m))
f=x(ntot-m+1:ntot)
deallocate(x)
allocate(x(m))
x=f
f=y(ntot-m+1:ntot)
deallocate(y)
allocate(y(m))
y=f
deallocate(f)
allocate(k(m,2))
do i=1,m
   k(i,1)=-1q0/(x(i)**2)
   k(i,2)=1q0
enddo
call ls(k,y,0,128,f) !min |k*[a b]-y|^2 wrt [a b]
getrealasym=f
deallocate(k,x,y,f)
end function

real(kind=16) function getsweightm(fn,n,q) !  find b in Im[G(i*w_n)] \approx -s/w_n
!fit first,x, vs third,y, column in file fn with: y(x)=-a/x
!Let's say n_tot is the number of points in file fn.
!Fit to the last min(n,q*n_tot) points
implicit none
character(len=800),intent(in) :: fn
integer,intent(in) :: n
real(kind=16),intent(in) :: q

real(kind=16),allocatable :: x(:),y(:),k(:,:),f(:)
integer :: i,m,ntot

ntot=filelen(fn)
m=min(real(n,kind=16),ntot*q)
write(*,*) "Number of Matsubara points used for finding asymptotic:",m
call openf(fn,x,f,y)
deallocate(f)
allocate(f(m))
f=x(ntot-m+1:ntot)
deallocate(x)
allocate(x(m))
x=f
f=y(ntot-m+1:ntot)
deallocate(y)
allocate(y(m))
y=f
deallocate(f)
allocate(k(m,1))
do i=1,m
   k(i,1)=-1q0/x(i)
enddo
allocate(f(1))
call ls(k,y,0,128,f) !min |k*a-y|^2 wrt f
getsweightm=f(1)
deallocate(f,k,x,y)
end function

end module 
