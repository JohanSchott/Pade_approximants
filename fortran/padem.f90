module padem
use leastsquare
!use mpim
implicit none

contains

subroutine readcof(proj,contin,cof)
implicit none
character(len=800),intent(in) :: proj
integer,intent(in) :: contin
complex(kind=16),intent(inout) :: cof(:)

real(kind=16) :: re,im
character(len=800) :: str
integer :: i,n
str=int2str(contin)
open(97,file=trim(proj)//"_cof_"//trim(str))
n=size(cof,1)
do i=1,n
   read(97,*) re,im
   cof(i)=cmplx(re,im,kind=16)
enddo
call system("rm "//trim(proj)//"_cof_"//trim(str))
end subroutine

subroutine getPadeCof(z,f,cof,err)
!Fit Pade approximant to f(z) and return coefficients in Pade approximant
   implicit none
   complex(kind=16),intent(in) :: z(:),f(:)
   complex(kind=16),intent(out) :: cof(:)
   real(kind=16),intent(out) :: err(:)
   !internal
   integer :: M,N,r
   integer :: i,j
   complex(kind=16),allocatable :: A(:,:),b(:)
   complex(kind=16),allocatable :: sol(:),errv(:),fpade(:)
   M=size(z,1)
   N=size(cof)
   allocate(A(M,N),b(M))
   if(mod(N,2)==0 ) then
      r=N/2
      !construct A and b
      do i=1,M
         do j=1,r
            A(i,j) = z(i)**(j-1)
            A(i,r+j) = -z(i)**(j-1)*f(i)
         enddo
         b(i) = f(i)*z(i)**r
      enddo
   else
      r=(N-1)/2
      !construct A and b
      do i=1,M
         do j=1,r
            A(i,j) = z(i)**(j-1)
            A(i,r+1+j) = -z(i)**(j-1)*f(i)
         enddo
         A(i,r+1)=z(i)**r
         b(i) = f(i)*z(i)**r
      enddo
   endif
   call ls(A,b,0,128,sol)
   allocate(errv(M))
   errv=matmul(A,sol)-b
   err(1)=sqrt(sum(real(errv)**2+aimag(errv)**2)) ! |k*x-b|
   allocate(fpade(M))
   call eval_pade(sol,z,fpade)
   err(2)=sqrt(sum(abs(f-fpade)**2)) ! |f-fpade|
   deallocate(A,b)
   deallocate(errv,fpade)
   cof=sol !copy solution
   deallocate(sol)

end subroutine

subroutine eval_pade(p,z,sig)
!----------------------------------------------------------------
!        p(nom=2r) - coifficients of Pade approximation
!        z  - complex points
!        
!        f(z)=(p(1)+p(2)z+...+p(r)z**(r-1))/
!                        (p(r+1)+...+p(2r)z**(r-1)+z**r)
!
!        output: sig    self energy at points z
!------------------------------------------------------------------
IMPLICIT NONE 
complex(kind=16),intent(in)   :: p(:),z(:) 
complex(kind=16),intent(out)  :: sig(:)

complex(kind=16),allocatable :: a(:),b(:)
INTEGER :: j
integer :: ne,nom,r
ne=size(z)  !nbr of complex points
nom=size(p) !nbr of coefficients
allocate(a(ne),b(ne))
if(mod(nom,2)==0) then
   r = nom/2   !2r Matsubara points creates r poles
   a(:)=p(1)
   b(:)=p(r+1)
   do j=2,r
      a(:)=a(:)+p(j)*z(:)**(j-1)
      b(:)=b(:)+p(r+j)*z(:)**(j-1)
   enddo
   b(:)=b(:)+z(:)**r
   sig(:)=a(:)/b(:) !Self energy as fraction of two polynomials
else
   r = (nom-1)/2  !r nbr of poles 
   a(:)=p(1)
   b(:)=p(r+2)
   do j=2,r
      a(:)=a(:)+p(j)*z(:)**(j-1)
      b(:)=b(:)+p(r+1+j)*z(:)**(j-1)
   enddo
   a(:)=a(:)+p(r+1)*z(:)**r
   b(:)=b(:)+z(:)**r
   sig(:)=a(:)/b(:) !Self energy as fraction of two polynomials
endif
deallocate(a,b)
!do i=1,ne
!  a = p(1)
!  b = p(r+1)
!  do j = 2,r
!     a = a + p(j)*z(i)**(j-1)
!     b = b + p(r+j)*z(i)**(j-1)
!  enddo
!  b = b + z(i)**r
!  sig(i)=a/b !Self energy as fraction of two polynomials
!enddo

END SUBROUTINE  

function padediff(sigma1,sigma2)
!   Subroutine to calculate the difference between two 
!   pade approximants.
!       Input: sigma1, sigma2 - complex vectors
!       Output: diff - the sum of the lengths between
!                      each point of the input vectors
implicit none
complex(kind=16),intent(in)   :: sigma1(:),sigma2(:)
real(kind=16)   :: padediff
padediff=sum(abs(sigma1-sigma2)) !/size(sigma1)
end function

subroutine getDisSerial(good,A,dis)
implicit none
integer,intent(in) :: good(:)
real(kind=16) :: A(:,:),dis(:)

!comp variables
integer :: ns,ne
integer :: g,N,nt,s,f
integer :: i,j,k,r
integer,allocatable :: ind(:),iv(:),jv(:)
real(kind=16),allocatable :: dr(:),drtmp(:),dv(:),dm(:,:),di(:)
   
ne=size(A,1)
ns=size(A,2)
g=sum(good)
dis=0
if(g==1) then
   i=minloc(abs(good-1),1)
   dis(i)=1
else
   !Copy good data to smaller matrix
   allocate(di(g),ind(g))
   j=0
   ind=0
   do i=1,ns
       if(good(i)==1) then
           j=j+1
           ind(j)=i !stores indices of the good functions 
       endif
   enddo
   N=g*(g-1)/2 !nbr of tasks (nbr of elements the lower part of triangle of dimension (g,g))
   allocate(dv(N))
   !call distjobs(N,s,f,nt) . Instead use three lines below
   s=1
   f=N
   nt=N
   !Calculate corresponding matrix elements for each task
   allocate(iv(nt),jv(nt))
   k=0
   r=0
   do j=1,g
      do i=j+1,g
         k=k+1
         if(s<=k .and. k<=f) then
            r=r+1
            iv(r)=i
            jv(r)=j
         endif
      enddo
   enddo
   !Calculate tasks for every process
   allocate(dr(nt))
   do k=1,nt
      dr(k)=sum(abs( A(:,ind(iv(k)))-A(:,ind(jv(k))) ))
   enddo
   deallocate(iv,jv)
   !All ranks except master sends their version of dr to the master rank
   !if(myrank==0) then
      dv(s:f)=dr  !rank 0 puts its information to big vector
      !do i=1,ranks-1 !loop over all slave processes
      !   call MPI_RECV(j,1,MPI_INTEGER,i,i,MPI_COMM_WORLD,status,ierr)
      !   allocate(drtmp(j))
      !   call MPI_RECV(drtmp,j,MPI_REAL16,i,ranks+i,MPI_COMM_WORLD,status,ierr)
      !   call MPI_RECV(k,1,MPI_INTEGER,i,2*ranks+i,MPI_COMM_WORLD,status,ierr)
      !   dv(k:(k+j-1))=drtmp
      !   deallocate(drtmp)
      !enddo
   !else
   !   call MPI_SEND(nt,1,MPI_INTEGER,0,myrank,MPI_COMM_WORLD,ierr) 
   !   call MPI_SEND(dr,nt,MPI_REAL16,0,ranks+myrank,MPI_COMM_WORLD,ierr) 
   !   call MPI_SEND(s,1,MPI_INTEGER,0,2*ranks+myrank,MPI_COMM_WORLD,ierr) 
   !endif
   !deallocate(dr)
   !call MPI_BCAST(dv,N,MPI_REAL16,0,MPI_COMM_WORLD,ierr) 
   allocate(dm(g,g))
   dm=0q0
   k=0
   do j=1,g-1
      do i=j+1,g
         k=k+1
         dm(i,j)=dv(k)
      enddo
   enddo
   deallocate(dv)
   !Use the fact dm(i,j) is symmetric
   do i=1,g
       do j=i+1,g
           dm(i,j) = dm(j,i) 
       enddo
   enddo
   !Calculate di
   di=0q0
   do i=1,g
       do j=1,g
           di(i) = di(i)+dm(i,j) 
       enddo
   enddo
   deallocate(dm)
   di=di/(g-1)
   !copy back to all continuations, including unphysical ones   
   do i=1,g
      dis(ind(i))=di(i)
   enddo
   deallocate(di)
endif
end subroutine

subroutine getDisWeight(acc1,acc2,dis,good,weight) !calculates weight from dis-vector
!Calculate weight to be used to average the Pade continuations
implicit none
real(kind=16) :: acc1,acc2,dis(:)
integer :: good(:)
real(kind=16) :: weight(:)

integer :: i,j,ns,g
real(kind=16),allocatable :: di(:),disort(:),w(:)
real(kind=16) :: da,d1,d2

write(*,*) "inside getDisWeight"

g=sum(good)
allocate(di(g),disort(g),w(g))
ns=size(good)
j=0
do i=1,ns
   if(good(i)==1) then 
      j=j+1
      di(j)=dis(i) 
   endif
enddo
da=sum(di)/g
disort=sort(di)
d1=da*acc1
d2=disort(ceiling(acc2*g))
deallocate(disort)

!METHOD 1: w_i = (d_i<da*acc)
!j=0
!w=0
!accd=acc
!do while(j<1)
!    do i=1,good
!        if(di(i)<=da*accd) w(i) = 1
!    enddo
!    j=nint(sum(w))
!    accd=accd+acc
!enddo

!METHOD 2: Take the acc2 procent "best" solutions
! if acc2=0.51 and ns=2 will include both functions
!w=0
!do i=1,good
!    if(di(i)<=disort(ceiling(acc2*good))) w(i) = 1
!enddo
!j=nint(sum(w))

! METHOD 3: combine method 1 and method 2
j=0
i=0
w=0q0
do while(j<1)
    i=i+1
    if(i>1) d1=2*d1
    do i=1,g
        if( di(i)<=d1 .and. di(i)<=d2 ) w(i) = 1q0
    enddo
    j=nint(sum(w))
enddo


! METHOD 4: w_i = exp(-d_i/std)

!check how many continuations that contribute to the average.
j=0
do i=1,g
    if(w(i)>0) j=j+1
enddo
write(33,'(I4,a,I4,a)') j," out of ",g," good continuations contribute to default average."
!Normalize weights
w=w/sum(w)

!Go back to full problem
weight=0q0
j=0
do i=1,ns
   if(good(i)==1) then
      j=j+1
      weight(i)=w(j)
   endif
enddo
deallocate(di,w)
end subroutine

end module 
