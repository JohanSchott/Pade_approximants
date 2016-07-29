module config
use padem
!use savem
implicit none

type conf
   integer :: N,first,M !nbr of input Matsubara points, index of first matsubara point and number of Pade coefficients. Sometimes meaning of N and M are interchanged
   complex(kind=16),allocatable :: x(:) !Pade coefficients
   complex(kind=16),allocatable :: G(:),Gx(:) !Pade evaluated on Matsubara and real axis
   real(kind=16),allocatable :: A(:) !spectrum
end type

type(conf),allocatable       :: c(:),ac(:)
integer,allocatable :: Nv(:) !Number of Pade coefficients for each continuations

contains

subroutine confsetup(npade,Nmax,Nmin,Nstep,Mmax,Mmin,Mstep,firstmin,firstmax,firststep,ns)
implicit none
logical,intent(in) :: npade
integer,intent(in) :: Nmax,Nmin,Nstep,Mmax,Mmin,Mstep,firstmin,firstmax,firststep
integer,intent(inout) :: ns

integer :: first,i,j,N,M
type(conf),allocatable       :: ctmp(:)

! check if normal Pade or LS-Pade and create c-vector
if(npade) then
    i=0
    do N=Nmax,Nmin,-Nstep
        i=i+1
    enddo
    ns=i !nbr of different continuations/configurations (N,N)
    if(allocated(c)) deallocate(c)
    allocate(c(ns)) 
    i=0
    do N=Nmax,Nmin,-Nstep
        i=i+1
        c(i)%N=N
        c(i)%M=N
    enddo
else
    i=0
    do N=Nmax,Nmin,-Nstep
        do M=min(N,Mmax),Mmin,-Mstep
            i=i+1
        enddo
    enddo
    ns=i 
    if(allocated(c)) deallocate(c)
    allocate(c(ns)) 
    i=0
    do N=Nmax,Nmin,-Nstep
        do M=min(N,Mmax),Mmin,-Mstep
            i=i+1
            c(i)%N=N
            c(i)%M=M
        enddo
    enddo
endif

allocate(ctmp(ns))
ctmp=c
deallocate(c)
j=0
do first=firstmin,firstmax,firststep
   j=j+1
enddo
allocate(c(ns*j),Nv(ns*j))
i=0
do first=firstmax,firstmin,-firststep
   write(*,*) "first=",first
   c(ns*i+1:ns*(i+1))%N=ctmp(:)%N
   c(ns*i+1:ns*(i+1))%M=ctmp(:)%M
   c(ns*i+1:ns*(i+1))%first=first
   Nv(ns*i+1:ns*(i+1))=ctmp(:)%M
   i=i+1
enddo
ns=ns*j
deallocate(ctmp)
do i=1,ns
   allocate(c(i)%x(c(i)%M))  
enddo
end subroutine

subroutine getlambda(lambda)
implicit none
real(kind=16),intent(out) :: lambda(:)
integer :: j,r
do j=1,size(c,1)
   if(mod(c(j)%M,2)==0) then
      r=c(j)%M/2
   else
      r=(c(j)%M-1)/2
   endif
   lambda(j) = -log(abs(aimag(c(j)%x(r))))/log(10q0)
enddo
end subroutine

subroutine getSpectrums(ftype,Gre,emesh,good) ! c is implicitly an input variable
implicit none
integer,intent(in) :: ftype
complex(kind=16) :: Gre(:,:),emesh(:)
integer :: good(:)
!
integer :: i,j,k,u
real(kind=16) :: dx

if(ftype==3) then
   good=1
   do i=1,size(c)
      call eval_pade(c(i)%x,emesh,Gre(:,i))
   enddo
else
   good=0
   j=0
   do while(sum(good)<1)
      j=j+1
      if(j>1) then
         emesh=emesh+cmplx(0q0,aimag(emesh(1)),kind=16)
         write(*,*) "Warning:"
         write(*,*) "No physical continuations found."
         write(*,*) "Increase distance to real axis."
      endif
      do i=1,size(c)
         call eval_pade(c(i)%x,emesh,Gre(:,i))
         if(ftype==0 .or. ftype==1 .or. ftype==3) then
            if(maxval(aimag(Gre(:,i)))<=0) then
              good(i)=1 
            else
              good(i)=0
            endif
         elseif(ftype==2 .or. ftype==4 ) then
            dx=0.0001
            u=0
            do k=1,size(Gre(:,i))
               if(real(emesh(k))<dx .and. aimag(Gre(k,i))<0q0) then
                  u=u+1
               elseif(real(emesh(k))>dx .and. aimag(Gre(k,i))>0q0) then
                  u=u+1
               endif
            enddo
            if(u==0) then
               good(i)=1
            else 
               good(i)=0
            endif
            !if(maxval(real(emesh)*aimag(Gre(:,i)))<=0) then
            !  good(i)=1 
            !else
            !  good(i)=0
            !endif
         else
            stop "This ftype is not supported"
         endif
      enddo
   enddo
   if(j>1) then
      write(33,'(a,E13.3,a)') "Use eim=",aimag(emesh(1))," to get physical solutions."
   endif
   write(33,'(I4,a,I4,a)') sum(good)," physical continuations found (out of",size(good)," possible)."
endif

end subroutine

end module 
