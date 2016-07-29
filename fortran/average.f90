module average
use config
!use diffm
implicit none


interface faverage
   module procedure faverageN,faverageImPr0
end interface 

interface padeAverage
   module procedure padeDistsAverage,padeDistAverage
end interface 

contains

subroutine faverageImPr0(weights,asym,z,G,std)
!**************************************************
! Calculate average of the ns continuations
! Using the distribution from variable weights.
! This routine also calculates the standard deviation. 
!**************************************************
!input/output variables
implicit none
real(kind=16),intent(in)        :: weights(:)
logical,intent(in)              :: asym
complex(kind=16),intent(in)     :: z(:)
complex(kind=16),intent(out)    :: G(:) !average output
real(kind=16),intent(inout)     :: std(:,:)  !output
!computational variables
integer                         :: j,r
complex(kind=16),allocatable    :: cof(:),sigma(:),sigmaS(:)
real(kind=16)                   :: s
integer                         :: ns,ne

ns=size(weights,1)
ne=size(z)
allocate(sigma(ne),sigmaS(ne))
sigmaS=0q0
std=0q0
s=sum(weights)
if(s .le. 0q0) stop "weighted sum has to be larger than 0"
do j=1,ns   
   if(asym) then  ! Sigma calculated with coefficient Im[P_r]=0
      allocate(cof(c(j)%M))
      cof=c(j)%x
      if(mod(size(cof),2)==0) then
         r=size(cof)/2
         cof(r)=real(cof(r))
      else
         r=(size(cof)-1)/2
         cof(r)=real(cof(r))
      endif
      call eval_pade(cof,z,sigma)
      deallocate(cof)
   else
      call eval_pade(c(j)%x,z,sigma)
   endif
   sigmaS = sigmaS + sigma*weights(j)
   std(:,1) = std(:,1)+weights(j)*real(sigma)**2
   std(:,2) = std(:,2)+weights(j)*aimag(sigma)**2
enddo
sigmaS = sigmaS/s
std(:,1)=std(:,1)/s-(real(sigmaS))**2  
std(:,2)=std(:,2)/s-(aimag(sigmaS))**2 
where(std<0q0) 
    std=0q0
endwhere
std=sqrt(std)
G=SigmaS
deallocate(sigma,sigmaS)
end subroutine

subroutine faverageN(weights,z,G,std)         
!input/output variables
implicit none
real(kind=16),intent(in)        :: weights(:)
complex(kind=16),intent(in)     :: z(:)
complex(kind=16),intent(out)    :: G(:) !average output
real(kind=16),intent(inout)     :: std(:,:)  !output
call faverageImPr0(weights,.false.,z,G,std)
end subroutine

subroutine averagePadeCoeff(emesh,good,r)
!Calculate average of Pade coefficients among continuations having same number of coefficients
!The average coefficients are calculated and stored:
!  ac(i)%x(:) = 1/#cont * sum_{j with N fix and has physical spectra} c(j)%x(:)
!The corresponding spectra is calculated and stored:
!  call eval_pade(ac(i)%x(:),emesh,sigma)
!  ac(i)%A(:) = -1/pi*aimag(sigma)
!Return an average of these spectrums.
!This average is at first done simple by:
!  r=0q0
!  j=0
!  do i=1,size(ac)
!     if(minval(ac(i)%A(:))>-0.00001) then
!        r=r+ac(i)%A(:) 
!        j=j+1
!     endif
!  enddo
!  r=r/j
!Input variables: 
!  emesh - points in complex plane where want to calculate spectrums
!  c     - All the continuations performed 
!  good  - which of continuations in c are physical
!Output variables:
!  ac    - All the averaged Pade coefficients averaged continuations 
!  r     - The average of the spectrums in ac
implicit none
complex(kind=16),intent(in) :: emesh(:)
integer,intent(in) :: good(:)
real(kind=16),allocatable :: r(:)
!dummy variables
integer :: ns,i,j,k
integer,allocatable :: n(:),ncopy(:)
complex(kind=16),allocatable :: sigma(:)

!calculate how many different N-approximants there are
ns=size(c)
do i=1,ns
   if(i==1) then
      allocate(n(1))
      n(1)=c(1)%M
   elseif(minval(abs(n-c(i)%M))>0) then
      allocate(ncopy(size(n)))
      ncopy=n
      deallocate(n)
      allocate(n(size(ncopy)+1))
      n(1:size(ncopy))=ncopy
      n(size(n))=c(i)%M
      deallocate(ncopy)
   endif
enddo

!Allocate ac 
allocate(ac(size(n)))
do i=1,size(n)
   ac(i)%N=n(i) ! number of pade coefficients
   allocate(ac(i)%A(size(emesh)),ac(i)%x(ac(i)%N))
enddo
deallocate(n)

!Perform average over Pade coefficients
allocate(sigma(size(emesh)))
do i=1,size(ac)
   ac(i)%x=0q0
   k=0
   do j=1,ns
      if(c(j)%M==ac(i)%N .and. good(j)==1) then
         ac(i)%x=ac(i)%x+c(j)%x   
         k=k+1
      endif
   enddo
   ac(i)%x=(ac(i)%x)/k
   call eval_pade(ac(i)%x(:),emesh,sigma)
   ac(i)%A(:) = -1/pi*aimag(sigma)
enddo
deallocate(sigma)

!Perform average over spectrums
allocate(r(size(emesh)))
r=0q0
j=0
do i=1,size(ac)
   if(minval(ac(i)%A(:))>-0.00001) then
      r=r+ac(i)%A(:) 
      j=j+1
   endif
enddo
if(j==0) then
   write(*,*) "Warning: No Pade coefficient averaged continuations are physical"
   write(*,*) "Take an average of all instead" 
   do i=1,size(ac)
      r=r+ac(i)%A(:)
   enddo
   r=r/size(ac)
   !stop 
else
   r=r/j
   write(*,*) j," out of ",size(ac)," Pade coefficient averaged continuations are physical"
endif
end subroutine 

subroutine padeDistsAverage(weights,Gre,Gm,Gstd)
!Calculate averages and standard deviations.
!w(config type,average type), Gre(w,config type)
!Gm(w,average type),Gstd(w,average type,real and imaginary part)
implicit none
real(kind=16) :: weights(:,:)
complex(kind=16) :: Gre(:,:)
complex(kind=16),allocatable :: Gm(:,:)
real(kind=16),allocatable :: Gstd(:,:,:)
!
integer :: nw,na,ns,i,j 
nw=size(Gre,1) !number of frequencies
ns=size(Gre,2) !number of configurations
na=size(weights,2) !number of different averaging distributions
if(allocated(Gm)) deallocate(Gm)
if(allocated(Gstd)) deallocate(Gstd)
allocate(Gm(nw,na),Gstd(nw,na,2))
Gm=0q0
Gstd=0q0
do i=1,na 
   do j=1,ns
      Gm(:,i)=Gm(:,i)+Gre(:,j)*weights(j,i)
      Gstd(:,i,1)=Gstd(:,i,1)+real(Gre(:,j))**2*weights(j,i)
      Gstd(:,i,2)=Gstd(:,i,2)+aimag(Gre(:,j))**2*weights(j,i)
   enddo
   Gstd(:,i,1)=Gstd(:,i,1)-real(Gm(:,i))**2
   Gstd(:,i,2)=Gstd(:,i,2)-aimag(Gm(:,i))**2
enddo
where(Gstd<0q0 ) Gstd=0q0
Gstd=sqrt(Gstd)
end subroutine
subroutine padeDistAverage(weights,Gre,Gm,Gstd)
!Calculate averages and standard deviations.
!w(config type,average type), Gre(w,config type)
!Gm(w,average type),Gstd(w,average type,real and imaginary part)
implicit none
real(kind=16) :: weights(:)
complex(kind=16) :: Gre(:,:)
complex(kind=16),allocatable :: Gm(:,:)
real(kind=16),allocatable :: Gstd(:,:,:)
!
real(kind=16),allocatable :: w(:,:)
allocate(w(size(weights),1))
w(:,1)=weights
call padeDistsAverage(w,Gre,Gm,Gstd)
end subroutine

subroutine spectrumErr(w,Ae,A,err)
!Calculate and returns spectrum errors, spectrumerr(errors type,average type)
implicit none
real(kind=16) :: w(:),Ae(:),A(:,:)
real(kind=16),allocatable :: err(:,:)
!
integer :: i
if(allocated(err)) deallocate(err)
allocate(err(4,size(A,2)))
do i=1,size(A,2) !loop over different averages
      err(1,i)=integ(w,abs(Ae-A(:,i)))
      err(2,i)=err(1,i)/integ(w,abs(Ae))
      err(3,i)=integ(w,(Ae-A(:,i))**2)
      err(4,i)=err(3,i)/integ(w,Ae**2)
enddo
end subroutine

end module 
