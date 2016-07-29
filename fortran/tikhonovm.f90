module tikhonovm
use leastsquare
implicit none

interface tikhonov
   module procedure dtikhonov,dtikhonovDefaultModel,ztikhonov
end interface tikhonov

contains

subroutine dtikhonov(A,bo,mtr,pr,nonegx,x)
!Finds Tikhonov SVD solution (using real variables):
!  Filter away singular value components, using parameter alpha.
!  x=sum_i=1^r s_i/(s_i^2+alpha^2) (u_i^dagger*b) v_i
!  x is actually just a LS solution to slightly modified problem:
!  |A      |     |b|
!  |       |*x = | |=k*x=g, namely solves normal equation: k^T*k*x=k^T*g. So instead of LS solution for A and b, find LS solution for k and g.
!  |alpha*I|     |0|
!  Another way of formulating it: min_x |A*x-b|^2+|alpha*x|^2, for a fix alpha
!  Optimal parameter alpha value can be obtained by minimizing ||res||_2+||x||_2 
! If option nonegx=.true. is LS found under constrain of x=>0
implicit none
real(kind=16),intent(in) :: A(:,:)
real(kind=16),intent(in) :: bo(:)
integer,intent(in) :: mtr,pr
logical,intent(in) :: nonegx
real(kind=16),allocatable :: x(:)
!internal 
integer :: N,M,i,j
real(kind=16),allocatable :: b(:),k(:,:),g(:),xt(:,:)
integer :: na
real(kind=16),allocatable :: alpha(:),res(:),norm(:)
character(len=800) :: str

N=size(A,1)
M=size(A,2)
allocate(b(size(bo,1)))
b=bo !copy
allocate(k(N+M,M),g(N+M))
na=50
allocate(alpha(na),xt(M,na),res(na),norm(na))
alpha= [ (10**(-(3q0+(i-1q0)*(9q0-3)/(na-1q0))),i=1,na) ]
write(*,'(a,E10.3,a,E10.3)') "Try alphas between:",alpha(1)," and",alpha(na)
do i=1,na !loop over different alpha
   write(*,'(a,E10.3)') "alpha=",alpha(i)
   g=0
   g(1:N)=b
   k=0
   k(1:N,1:M)=A
   do j=1,M
      k(N+j,j)=alpha(i)
   enddo
   if(nonegx) then
      call nnls(k,g,mtr,pr,x) ! Pick your favourite LS solver routine, under constrain of x>0
   else
      !call ls(k,g,0,128,x) ! This LS solver since it seams to be the best.
      call ls(k,g,mtr,pr,x) ! Pick your favourite LS solver routine.
   endif
   xt(:,i)=x
   res(i)=sqrt(sum((matmul(A,x)-b)**2))   ! |A*x-b|
   norm(i)=sqrt(sum(x**2))  ! |x|
enddo
open(54,file="output_tikhonov.dat")
write(*,'(a)') "Tikhonov output:"
write(*,'(a)') "alpha      |A*x-b|         |x|"
write(54,'(a)') "alpha      |A*x-b|         |x|"
do i=1,na
   write(*,'(3E15.5)') alpha(i),res(i),norm(i)
   write(54,'(3E15.5)') alpha(i),res(i),norm(i)
enddo
write(*,*)
close(54)
str=int2str(na)
str="("//trim(str)//"E13.4)"
open(54,file="output_tikhonov_a.dat")
do i=1,M
   write(54,trim(str)) xt(i,:)
enddo
i=minloc(log(res*norm),1)
write(*,'(a,E15.5,a,I4)') "Pick alpha=",alpha(i)," index=",i
write(54,'(a,E15.5,a,I4)') "Pick alpha=",alpha(i)," index=",i
close(54)
x=xt(:,i)
deallocate(b,alpha,xt,k,g,res,norm)
end subroutine

subroutine dtikhonovDefaultModel(A,bo,xm,mtr,pr,nonegx,x)
!Finds Tikhonov solution (using real variables) to fit also to default model xm:
!  x is just a LS solution to problem:
!  |A      |     |b|
!  |       |*x = | |=k*x=g, namely solves normal equation: k^T*k*x=k^T*g. So instead of LS solution for A and b, find LS solution for k and g.
!  |alpha*I|     |alpha*xm|
!  Another way of formulating it: min_x |A*x-b|^2+|alpha*(x-xm)|^2, for a fix alpha
!  Optimal parameter alpha value can be obtained by minimizing |A*x-b|^2+|x|^2 
! If option nonegx=.true. is LS found under constrain of x=>0
implicit none
real(kind=16),intent(in) :: A(:,:),bo(:),xm(:)
integer,intent(in) :: mtr,pr
logical,intent(in) :: nonegx
real(kind=16),allocatable :: x(:)
!internal 
integer :: N,M,i,j
real(kind=16),allocatable :: b(:),k(:,:),g(:),xt(:,:)
integer :: na
real(kind=16),allocatable :: alpha(:),res(:),norm(:)
character(len=800) :: str
N=size(A,1)
M=size(A,2)
allocate(b(size(bo,1)))
b=bo !copy
allocate(k(N+M,M),g(N+M))
na=50
allocate(alpha(na),xt(M,na),res(na),norm(na))
alpha= [ (10**(-(3q0+(i-1q0)*(9q0-3)/(na-1q0))),i=1,na) ]
write(*,'(a,E10.3,a,E10.3)') "Try alphas between:",alpha(1)," and",alpha(na)
do i=1,na !loop over different alpha
   write(*,'(a,E10.3)') "alpha=",alpha(i)
   !Extended matrix
   k=0q0
   k(1:N,1:M)=A
   do j=1,M
      k(N+j,j)=alpha(i)
   enddo
   !Extended right hand side
   g(1:N)=b
   g(N+1:N+M)=alpha(i)*xm
   if(nonegx) then
      call nnls(k,g,mtr,pr,x) ! Pick your favourite LS solver routine, under constrain of x=>0
   else
      !call ls(k,g,0,128,x) ! This LS solver since it seams to be the best.
      call ls(k,g,mtr,pr,x) ! Pick your favourite LS solver routine.
   endif
   xt(:,i)=x
   res(i)=sqrt(sum((matmul(A,x)-b)**2))   ! |A*x-b|
   norm(i)=sqrt(sum((x-xm)**2))  ! |x-xm|
enddo
open(54,file="output_tikhonov.dat")
write(*,'(a)') "Tikhonov output:"
write(*,'(a)') "alpha      |A*x-b|        |x-xm|              |x|"
write(54,'(a)') "alpha      |A*x-b|        |x-xm|              |x|"
do i=1,na
   write(*,'(4E13.5)') alpha(i),res(i),norm(i),sqrt(sum(xt(:,i)**2))
   write(54,'(4E13.5)') alpha(i),res(i),norm(i),sqrt(sum(xt(:,i)**2))
enddo
write(*,*)
close(54)
str=int2str(na)
str="("//trim(str)//"E13.4)"
open(54,file="output_tikhonov_a.dat")
do i=1,M
   write(54,trim(str)) xt(i,:)
enddo
i=minloc(log(res*norm),1)
write(*,'(a,E15.5,a,I4)') "Pick alpha=",alpha(i)," index=",i
write(54,'(a,E15.5,a,I4)') "Pick alpha=",alpha(i)," index=",i
close(54)
x=xt(:,i)
deallocate(b,alpha,xt,k,g,res,norm)
end subroutine

subroutine ztikhonov(A,bo,mtr,pr,x)
!Finds Tikhonov SVD solution (using complex variables):
!  Filter away singular value components, using parameter alpha.
!  x=sum_i=1^r s_i/(s_i^2+alpha^2) (u_i^dagger*b) v_i
!  x is actually just a LS solution to slightly modified problem:
!  |A      |     |b|
!  |       |*x = | |=k*x=g, namely solves normal equation: k^dagger*k*x=k^dagger*g. So instead of LS solution for A and b, find LS solution for k and g.
!  |alpha*I|     |0|
!  Optimal parameter alpha value can be obtained by minimizing ||res||_2+||x||_2 
implicit none
complex(kind=16),intent(in) :: A(:,:)
complex(kind=16),intent(in) :: bo(:)
complex(kind=16),allocatable :: x(:)
integer,intent(in) :: mtr,pr
!internal 
integer :: N,M,i,j
complex(kind=16),allocatable :: b(:),k(:,:),g(:),xt(:,:)
integer :: na
real(kind=16),allocatable :: alpha(:),res(:),norm(:)
N=size(A,1)
M=size(A,2)
allocate(b(size(bo,1)))
b=bo !cop,y
allocate(k(N+M,M),g(N+M))
na=50
allocate(alpha(na),xt(M,na),res(na),norm(na))
do i=1,na
   g=0
   g(1:N)=b
   k=0
   k(1:N,1:M)=A
   alpha(i)=10q0**(-(1+(i-1)*(30q0-2)/(na-1)))
   do j=1,M
      k(N+j,j)=alpha(i)
   enddo
   !call ls(k,g,0,128,x) ! This LS solver since it seams to be the best.
   call ls(k,g,mtr,pr,x) ! Pick your favourite LS solver routine.
   xt(:,i)=x
   res(i)=sum(real(matmul(A,x)-b)**2+aimag(matmul(A,x)-b)**2)   ! ||A*x-b||^2
   norm(i)=sum(real(x)**2+aimag(x)**2)  ! ||x||^2
enddo
write(*,'(a)') "Tikhonov output:"
write(*,'(a)') "alpha    ||A*x-b||^2    ||x||^2"
do i=1,na
   write(*,'(3E13.5)') alpha(i),res(i),norm(i)
enddo
i=minloc(res+norm,1)
!i=minloc(abs(alpha-10**(-30q0)),1)
write(*,*) "choose alpha=",alpha(i)
x=xt(:,i)
deallocate(b,alpha,xt,k,g,res,norm)
end subroutine

end module 
