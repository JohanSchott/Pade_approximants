PROGRAM main
use rw
use meshes
use input
use padeSetup
use poles   
use asym
use average
use openm
use savem
IMPLICIT NONE

!**********************************************************************************
!   Main program to perform analytic continuation by Pade approximant
!   Written by Johan Schott and Elin Lundin.
!   Version number and implementation history is documented in file rw.f90 and output file _info.dat
!   Most time is spent on constructing difference matrix between different solutions.
!   MPI parallelization is done only there.
!
!***********************************************************************************

!include "mpif.h"
type Gmean 
    integer :: n
    complex(kind=16),allocatable :: z(:),f(:)
    real(kind=16),allocatable :: err(:),aerr(:),std(:,:)
end type Gmean

!main variables
integer                      :: ns !nbr of configurations
complex(kind=16),allocatable :: pole(:),res(:) !poles and residues.
integer,allocatable          :: good(:) !if configurations are physical or not
complex(kind=16),allocatable :: Gre(:,:) !Pade approximants G(w+i*eim)
real(kind=16),allocatable    :: A(:,:)  !spectrums -1/pi*Im[G(w+i*eim)]

complex(kind=16),allocatable :: Gm(:,:) !Several different average Pade approximants G(w+i*eim)
real(kind=16),allocatable    :: Gstd(:,:,:) !Deviation in averages, both real and imaginary part 
real(kind=16),allocatable    :: Gerr(:,:) !real axis error

real(kind=16),allocatable    :: weights(:,:),dis(:)  !configuration weight and internal deviation between continuation on real axis
real(kind=16),allocatable    :: lambda(:) !pade coefficient, the higher the better is continuation
real(kind=16),allocatable    :: dev(:,:) !deviation. dev(j,1)=sqrt(sum |A*cof-b|^2). dev(j,2)=sqrt(sum |G_input-G_pade|^2)
real(kind=16)                :: sigma0  !Real constant Sigma_0 in Sigma(z) \approx Simga_0 + s/z for large |z|
real(kind=16),allocatable    :: finp(:,:) !Contains whole input file, real part shifted by sigma0

!dummy variables
character(len=800)           :: str,string
integer                      :: N,M,first,i,j,k
real(kind=8)                 :: tmpd
real(kind=16),allocatable    :: x(:),tmp(:),tmp2(:)
complex(kind=16),allocatable :: z(:),zr(:),zm(:),Ga(:,:),Gz(:) 
real                         :: start,finish !counting computational time
type(Gmean)                  :: G  


call readsettings()
call setfilenames()
call saveusedinput()
allocate(x(2),tmp(1))
tmp=0.2q0
x=getrealasym(f_mats,30,tmp(1)) !Reads input file and fit: -a/wn^2+b to last points
if(shifted) then
    sigma0=x(2)
else
    sigma0=0q0
endif
!If exact solution exists, shifting is done differently so care must be taken.. 
open(unit=33,file=trim(f_info),action="write",position="append")
write(33,*) "Re[G(i*w_n)]->",sigma0," for n->\inf"  
if(shifted) then
    write(33,'(a)') "Shift input function with this value, do continuations and shift back in the end."
else
    write(33,'(a)') "Do NOT Shift input function with this value. Directly do continuations."
endif
deallocate(x,tmp)
call openf(f_mats,x,tmp,tmp2)
!sigma0=tmp(size(tmp))
allocate(finp(size(x),3))
finp(:,1)=x
finp(:,2)=tmp-sigma0
finp(:,3)=tmp2
deallocate(x,tmp,tmp2)
open(32,file=f_mats_s,action="write")
do i=1,size(finp,1)
    write(32,'(3E80.70)') finp(i,1),finp(i,2),finp(i,3)
enddo
close(32)
!call save2f(f_mats_s,finp(:,1),finp(:,2),finp(:,3))
deallocate(finp)
call getrmesh(ftype,exact,eim,ne,f_ex,estart,estop) !constructs mesh eim above real axis
call checksettings() !checks if values in input file are valid before proceeding.
call confsetup(npade,Nmax,Nmin,Nstep,Mmax,Mmin,Mstep,firstmin,firstmax,firststep,ns) !builds c and returns ns, the number of configurations
close(33)
allocate(Gre(ne,ns),A(ne,ns))

call cpu_time(start)
if(exact) allocate(dexact(ns,4))
allocate(dev(ns,2)) 
open(unit=33,file=trim(f_info),action="write",position="append")
if(verbosity>5) then
    open(34,file=trim(trim(proj)//"_cof"),action="write")
    open(35,file=trim(f_pole),action="write")
    write(34,'(a)') "# Re[cof]       Im[cof]      NbrMats   NbrCof "
    write(35,'(a)') "# Re[pole]      Im[pole]     Re[residue]    Im[residue]     NbrMats   NbrCof"
endif
if(exact .and. real2 ) then
    open(36,file=trim(proj)//"_real2real",action="write")
    open(37,file=trim(proj)//"_real2mats",action="write")
    open(38,file=trim(proj)//"_real2pole",action="write")
    write(33,'(a)') "#  |A*c-b|   |G(w)-Gpade(w)|"
    write(36,'(a)') "# w     Re[Gpade]    Im[Gpade]    NbrMats   NbrCof"
    write(37,'(a)') "# w     Re[Gpade]    Im[Gpade]    NbrMats   NbrCof"
    write(38,'(a)') "# Re[pole]      Im[pole]     Re[residue]    Im[residue]     NbrMats   NbrCof"
endif
write(*,*) "Loop over configurations (N,M)..."
write(*,*)
do j=1,ns
    str=int2str(j)
    N=c(j)%N
    M=c(j)%M
    first=c(j)%first
    allocate(ind(N),wn(N),Gn(N))
    if(pickm==0) then
        call getIndices(f_mats_s,first,ind,wn,Gn)
    elseif(pickm==1) then
        call getIndices(f_mats_s,first,q,ind,wn,Gn)
    elseif(pickm==2) then
        call getIndices(f_mats_s,first,q,wsplit,wmax,p,ind,wn,Gn)     
    endif
    write(*,'(a,I5,a,I5,a,I3,a,I3,a,I3,a)') "continuation nbr ",j," of ",ns," (first=",first,", N=",N," and M=",M,") "
    call pade(wn,Gn,mtr,pr,rep,solver,c(j)%x,dev(j,:)) !fortran routine

    !continuation from real axis to both real and Matsubara axis
    if(exact .and. real2 .and. ftype/=4) then
        allocate(z(M),x(2))
        call getPadeCof(emesh,fexact-sigma0,z,x) !output: z as Pade coefficients, x as solver errors
        write(33,'(2E13.4)') x(1),x(2)
        allocate(zr(size(emesh)))
        call eval_pade(z,emesh,zr) !evaluate function on the input points  
        do i=1,size(emesh)
            write(36,'(3E13.4,2I4)') real(emesh(i)),real(zr(i))+sigma0,aimag(zr(i)),N,M
        enddo
        write(36,*)
        call openf(f_mats_s,x)
        allocate(zm(size(x)))
        call eval_pade(z,cmplx(0q0,x,kind=16),zm) !evaluate function on Matsubara points
        do i=1,size(x)
            write(37,*) x(i),real(zm(i))+sigma0,aimag(zm(i)),N,M
        enddo
        write(37,*)
        !calculate poles and residues
        if(mod(M,2)==0) then
            allocate(pole(M/2),res(M/2))
        else
            allocate(pole((M-1)/2),res((M-1)/2))
        endif
        call getpoles(z,pole,res)
        do k=1,size(pole)
            write(38,'(4E15.5,2I4)') real(pole(k)),aimag(pole(k)),real(res(k)),aimag(res(k)),N,M
        enddo
        deallocate(pole,res)
        write(38,*)
        deallocate(z,x,zr,zm)
    endif

    if(verbosity>5) then
        !save coefficients to one file
        do i=1,M
            write(34,'(2E70.60,2I4)') real(c(j)%x(i)),aimag(c(j)%x(i)),N,M
        enddo
        write(34,*) 
        !write poles to one file
        if(mod(M,2)==0) then
            allocate(pole(M/2),res(M/2))
        else
            allocate(pole((M-1)/2),res((M-1)/2))
        endif
        call getpoles(c(j)%x,pole,res)
        do k=1,size(pole)
            write(35,'(4E15.5,2I4)') real(pole(k)),aimag(pole(k)),real(res(k)),aimag(res(k)),N,M
        enddo
        deallocate(pole,res)
        write(35,*)
    endif

    deallocate(ind,wn,Gn)
enddo

call cpu_time(finish)

write(33,*)
write(33,'(I10,a)') int(finish-start)," seconds to do all continuations."
write(*,'(I10,a)') int(finish-start)," seconds to do all continuations."
close(34)
close(35)
if(exact .and. real2 ) then
    close(36)
    close(37)
    close(38)
endif

!******************************************
!Print information about exact solution real axis mesh
if(exact) then
  write(33,*)
  write(33,'(a)') "ne  estart  estop"
  write(33,'(I5,2F9.3)') ne,estart,estop
  write(33,*)
  open(42,file=trim(proj)//"_err",action="write")
  write(42,'(a)') "# 1:Int abs(A-Aex)dw  2:Int abs(A-Aex)dw/Int abs(Aex)dw  &
  3:Int (A-Aex)^2dw  4:Int (A-Aex)^2/Int (Aex)^2dw  -ln(E1  E2  E3  E4)/log(10)"
endif

!****************
!Calculate Pade approximants above real axis.
!Also check which are physical continuations, 
!i.e. have Im[G(w+i*eim)]<=0 if ftype=0,1,3 and
!have w*Im[G(w+i*eim)]<=0 if ftype=2,4
!If necessary to get at least one physical continuation, increase eim
allocate(good(ns))
call getSpectrums(ftype,Gre,emesh,good) ! c is implicitly an input variable

A=-1/pi*aimag(Gre)
!**************************************
!calculate dis
write(33,*)
write(33,'(a)') "calculating dis vector..."
write(*,'(a)') "calculating dis vector..."
close(33)
call cpu_time(start)

!Measure the difference inbetween spectrums
!Returning dis(:), storing devation with the other physical spectrums. 
!If continuation number "i" is unphysical, dis(i)=0
allocate(dis(ns))
!call getDis(good,A,dis)
call getDisSerial(good,A,dis)

open(unit=33,file=trim(proj)//"_info",action="write",position="append")
call cpu_time(finish)
write(33,'(I10,a)') int(finish-start)," seconds to calculate dis vector."
write(*,'(I10,a)') int(finish-start)," seconds to calculate dis vector."
write(33,*)
allocate(weights(ns,3))
call cpu_time(start)
call getDisWeight(acc,acc2,dis,good,weights(:,1)) !calculates weight from dis-vector
call cpu_time(finish)
write(*,*) "after disweight routine"
write(33,'(I10,a)') int(finish-start)," seconds to calculate dis weights."
weights(:,2)=(1q0*good)/sum(good) !physical continuations have equal weight
weights(:,3)=good/dev(:,2) !weight inverse proportional to sqrt(sum_n(|Gpade(i*w_n)-G(i*w_n)|^2)). Consider only physical continuations.
weights(:,3)=weights(:,3)/sum(weights(:,3)) !normalize to 1

!Three averges
call padeAverage(weights,Gre,Gm,Gstd) !calculate averages and standard deviations. Gm(w,average type),Gstd(w,average type,real and imaginary part)
if(ftype==4) then !symmetrize G: G(w+i*delta) = (G(w+i*delta)+G(-w+i*delta)^*)/2
  allocate(Ga(size(Gm,1),size(Gm,2))) 
  do i=1,size(Gm,2)
     Ga(:,i)=(Gm(:,i)+conjg(Gm(ne:1:-1,i)))/2
  enddo
  Gm=Ga
  deallocate(Ga)
endif
Gm=Gm+sigma0 !static shift. 
!Print averages and standard deviations
open(unit=50,file=trim(f_re),action="write")
open(unit=51,file=trim(trim(f_re)//"_p"),action="write")
open(unit=52,file=trim(trim(f_re)//"_devMats"),action="write")
str="#      w           <Re[G]>_w      <Im[G]>_w     -1/pi*<Im[G]>_w &
  std[Re[G]]_w     std[Im[G]]_w"
string=trim(str)//"  Im[Gexact]"
if(exact) then 
  write(50,'(a)') trim(string)
  write(51,'(a)') trim(string)
  write(52,'(a)') trim(string)
  do i=1,ne
         write(50,'(8E15.5)') real(emesh(i)),real(Gm(i,1)),aimag(Gm(i,1)),-1/pi*aimag(Gm(i,1)),&
         Gstd(i,1,1),Gstd(i,1,2),aimag(fexact(i))
         write(51,'(8E15.5)') real(emesh(i)),real(Gm(i,2)),aimag(Gm(i,2)),-1/pi*aimag(Gm(i,2)),&
         Gstd(i,2,1),Gstd(i,2,2),aimag(fexact(i))
         write(52,'(8E15.5)') real(emesh(i)),real(Gm(i,3)),aimag(Gm(i,3)),-1/pi*aimag(Gm(i,3)),&
         Gstd(i,3,1),Gstd(i,3,2),aimag(fexact(i))
  enddo
  !Print errors
  call spectrumErr(real(emesh),-1/pi*aimag(fexact),-1/pi*aimag(Gm),Gerr) !spectrum errors. Gerr(errors type,average type)
  write(42,'(8E15.7,a)') Gerr(:,1),-log(Gerr(:,1))/log(10d0)," # acc + acc2 weights"
  write(42,'(8E15.7,a)') Gerr(:,2),-log(Gerr(:,2))/log(10d0)," # physical weights"
  write(42,'(8E15.7,a)') Gerr(:,3),-log(Gerr(:,3))/log(10d0)," # weights proportional to sqrt(sum_n |Gpade(i*w_n)-G(i*w_n)|^2)"
else
  write(50,'(a)') trim(str)
  write(51,'(a)') trim(str)
  write(52,'(a)') trim(str)
  do i=1,ne
         write(50,'(8E15.5)') real(emesh(i)),real(Gm(i,1)),aimag(Gm(i,1)),-1/pi*aimag(Gm(i,1)),&
         Gstd(i,1,1),Gstd(i,1,2)
         write(51,'(8E15.5)') real(emesh(i)),real(Gm(i,2)),aimag(Gm(i,2)),-1/pi*aimag(Gm(i,2)),&
         Gstd(i,2,1),Gstd(i,2,2)
         write(52,'(8E15.5)') real(emesh(i)),real(Gm(i,3)),aimag(Gm(i,3)),-1/pi*aimag(Gm(i,3)),&
         Gstd(i,3,1),Gstd(i,3,2)
  enddo
endif
close(50)
close(51)
close(52)

!Diagonal physical continuation average
allocate(x(ns))
x=0
j=0
do i=1,ns
  if(c(i)%N==c(i)%M .and. good(i)==1) then
     x(i)=good(i)
     j=j+1
  endif
enddo
if(j>0) x=x/j
if(j>0) then
  call padeAverage(x,Gre,Gm,Gstd) !calculate averages and standard deviations. Gm(w,average type),Gstd(w,average type,real and imaginary part)
  if(ftype==4) then !symmetrize G: G(w+i*delta) = (G(w+i*delta)+G(-w+i*delta)^*)/2
     allocate(Ga(size(Gm,1),size(Gm,2))) 
     do i=1,size(Gm,2)
        Ga(:,i)=(Gm(:,i)+conjg(Gm(ne:1:-1,i)))/2
     enddo
     Gm=Ga
     deallocate(Ga)
  endif
  Gm=Gm+sigma0 !static shift. 
  !Print averages and standard deviations
  open(unit=53,file=trim(trim(f_re)//"_pdiag"),action="write")
  if(exact) then 
     write(53,'(a)') trim(string)
     do i=1,ne
            write(53,'(8E15.5)') real(emesh(i)),real(Gm(i,1)),aimag(Gm(i,1)),-1/pi*aimag(Gm(i,1)),&
            Gstd(i,1,1),Gstd(i,1,2),aimag(fexact(i))
     enddo
     call spectrumErr(real(emesh),-1/pi*aimag(fexact),-1/pi*aimag(Gm),Gerr) !spectrum errors. Gerr(errors type,average type)
     write(42,'(8E15.7,a)') Gerr(:,1),-log(Gerr(:,1))/log(10d0)," # physical diagonal weights" !print errors
  else
     write(53,'(a)') trim(str)
     do i=1,ne
            write(53,'(8E15.5)') real(emesh(i)),real(Gm(i,1)),aimag(Gm(i,1)),-1/pi*aimag(Gm(i,1)),&
            Gstd(i,1,1),Gstd(i,1,2)
     enddo
  endif
  close(53)
endif
deallocate(x)
!**********************************
if(verbosity>5) then
  !Calculate spectrum by averaging Pade coefficients from different
  !continuations, having the same number of coefficients. 
  !Use, in the averages, only continuations being physical.
  !The obtained spectrums are then averaged (only physical ones), to attain final spectrum.
  call averagePadeCoeff(emesh,good,x)
  if(ftype==4) then !symmetrize G: G(w+i*delta) = (G(w+i*delta)+G(-w+i*delta)^*)/2
     allocate(tmp(size(x,1)))
     tmp=(x-x(ne:1:-1))/2
     x=tmp
     deallocate(tmp)
  endif
  !write spectrum to file
  str=trim(proj)//"_re_coeffa"
  call save2f(str,real(emesh),x)
  !calculate real axis errors
  if(exact) then
     allocate(G%aerr(4))
     x=abs(-1/pi*aimag(fexact)-x)
     G%aerr(1)=integ(real(emesh),x)
     G%aerr(2)=G%aerr(1)/integ(real(emesh),1/pi*abs(aimag(fexact)))
     G%aerr(3)=integ(real(emesh),x**2)
     G%aerr(4)=G%aerr(3)/integ(real(emesh),(1/pi*abs(aimag(fexact)))**2)
     write(42,'(8E15.7,a)') G%aerr,-log(G%aerr)/log(10d0)," # average pade coefficients"
     deallocate(G%aerr)
  endif
  deallocate(x)
endif

!***************************************
!Extract Lambda from pade coeffients
allocate(lambda(ns))
call getlambda(lambda)
where(lambda<0) 
  lambda=0q0
end where
lambda=good*lambda  ! so only physical continuations have non zero lambda

!Lambda average
if(sum(lambda)>0 .and. verbosity>=2) then
  call padeAverage(lambda/sum(lambda),Gre,Gm,Gstd) !calculate averages and standard deviations. Gm(w,average type),Gstd(w,average type,real and imaginary part)
  if(ftype==4) then !symmetrize G: G(w+i*delta) = (G(w+i*delta)+G(-w+i*delta)^*)/2
     allocate(Ga(size(Gm,1),size(Gm,2))) 
     do i=1,size(Gm,2)
        Ga(:,i)=(Gm(:,i)+conjg(Gm(ne:1:-1,i)))/2
     enddo
     Gm=Ga
     deallocate(Ga)
  endif
  Gm=Gm+sigma0 !static shift. 
  !Print averages and standard deviations
  open(unit=53,file=trim(trim(f_re)//"_alambda"),action="write")
  if(exact) then 
     write(53,'(a)') trim(string)
     do i=1,ne
            write(53,'(8E15.5)') real(emesh(i)),real(Gm(i,1)),aimag(Gm(i,1)),-1/pi*aimag(Gm(i,1)),&
            Gstd(i,1,1),Gstd(i,1,2),aimag(fexact(i))
     enddo
     call spectrumErr(real(emesh),-1/pi*aimag(fexact),-1/pi*aimag(Gm),Gerr) !spectrum errors. Gerr(errors type,average type)
     write(42,'(8E15.7,a)') Gerr(:,1),-log(Gerr(:,1))/log(10d0)," # physical lambda weights" !print errors
  else
     write(53,'(a)') trim(str)
     do i=1,ne
            write(53,'(8E15.5)') real(emesh(i)),real(Gm(i,1)),aimag(Gm(i,1)),-1/pi*aimag(Gm(i,1)),&
            Gstd(i,1,1),Gstd(i,1,2)
     enddo
  endif
  close(53)
elseif(sum(lambda)<=0) then
  write(33,'(a)') "sum(lambda)=0"
endif


!*****************************************************************************************************************************
if(verbosity>9) then
  !Save all continuations on Matsubara axis
  !Shift back real part to origonal function by using sigma0 as last argument in getSigma or adding sigma0 just before printing
  call cpu_time(start)
  write(33,'(a)') "Save all Matsubara axis continuations to one file..."
  open(unit=21,file=trim(trim(f_im)//"_all"),action="write")
  write(21,'(a)') "# ind        wn            Re[Gpade]         Im[Gpade]       &
  Re[Ginput]        Im[Ginput]     N   M   first"
  do j=1,ns
      N=c(j)%N
      first=c(j)%first
      allocate(ind(N),wn(N),Gn(N))
      if(pickm==0) then
        call getIndices(f_mats_s,first,ind,wn,Gn)
      elseif(pickm==1) then
        call getIndices(f_mats_s,first,q,ind,wn,Gn)
      elseif(pickm==2) then
        call getIndices(f_mats_s,first,q,wsplit,wmax,p,ind,wn,Gn)     
      endif
      allocate(z(N))
      call eval_pade(c(j)%x,wn,z)
      z=z+sigma0
      do i=1,N
          write(21,'(I4,5E18.10,3I4)') ind(i),aimag(wn(i)),real(z(i)),aimag(z(i)),&
          real(Gn(i))+sigma0,aimag(Gn(i)),c(j)%N,c(j)%M,c(j)%first
      enddo
      write(21,*) " "
      deallocate(z)
      deallocate(ind,wn,Gn)
  enddo
  close(21)
  call cpu_time(finish)
  write(33,'(I10,a)') int(finish-start)," seconds to calculate and print Matsubara axis values _im_all."
endif


!Save all continuations on the real axis.
!*****************************************
if( ( verbosity>5 .and. exact )  .or. verbosity>9 ) then
  call cpu_time(start)
  if(verbosity>9) then
     write(33,'(a)') "Save all real axis continuations to one file..."
     open(unit=21,file=trim(proj)//"_re_all",action="write")
  endif
  if(exact .and. verbosity>9) then
     write(21,'(a)') "#         w            Re[G]            Im[G]     N   M &
  first good  weights_def    weights_phys    weights_Mats    Im[Gexact]" 
  elseif(verbosity>9) then
     write(21,'(a)') "#         w            Re[G]            Im[G]     N   M &
  first good  weights_def    weights_phys    weights_Mats"
  endif
  if(exact) then
     if(allocated(x)) deallocate(x)
     allocate(x(ne))
  endif
  do j=1,ns
     if(exact) then
        x=abs(A(:,j)+1/pi*aimag(fexact))
        dexact(j,1)=integ(real(emesh),x)
        dexact(j,2)=dexact(j,1)/integ(real(emesh),1/pi*abs(aimag(fexact)))
        dexact(j,3)=integ(real(emesh),x**2)
        dexact(j,4)=dexact(j,3)/integ(real(emesh),(1/pi*abs(aimag(fexact)))**2)
     endif
     if(verbosity>9) then
        if(exact) then
           do i=1,ne
              write(21,'(3E16.5,4I4,4E16.5)') real(emesh(i)),real(Gre(i,j))+sigma0,aimag(Gre(i,j)),&
              c(j)%N,c(j)%M,c(j)%first,good(j),weights(j,:),aimag(fexact(i))
           enddo
        else
           do i=1,ne
              write(21,'(3E16.5,4I4,3E16.5)') real(emesh(i)),real(Gre(i,j))+sigma0,aimag(Gre(i,j)),&
              c(j)%N,c(j)%M,c(j)%first,good(j),weights(j,:)
           enddo
        endif
        write(21,*)
     endif
  enddo
  if(allocated(x)) deallocate(x)
  if(verbosity>9) then
      close(21)
     call cpu_time(finish)
     write(33,'(I10,a)') int(finish-start)," seconds to calculate and print real axis spectrums in _re_all."
  endif
endif



!*******************************************************************************************************

!Matsubara average
call cpu_time(start)
write(33,*)
write(33,'(a)') "calculate matsubara average..."
call getmmesh(f_mats_s,wmax,firstmin,ind,G%z,z) !Obtain w_n=G%z and G(i*w_n)=z from input file for all Matsubara points between w_firstmin to w_n=wmax
G%n=size(G%z,1) !Number of Matsubara points considered 
allocate(G%f(G%n)) 
allocate(G%std(G%n,2))
str="# ind     w_n           <Re[G(i*w_n)]>_w <Im[G(i*w_n)]>_w &
std[Re[G(i*w_n)]]_w std[Im[G(i*w_n)]]_w Re[Ginput(i*w_n)] Im[Ginput(i*w_n)]" 
!Default weights
call faverage(weights(:,1),G%z,G%f,G%std) !calculate Pade G%f=G(i*w_n) and Pade standard deviation G%std=std[G(i*w_n]          
write(33,'(a,E15.6,a)') "|Gn-<Gn>_w|=sqrt(sum_n |Gn-<Gn>_w|^2 )=",&
  sqrt(sum(abs(G%f-z)**2)),". Using default weights."
open(unit=24,file=trim(f_im),action="write")
write(24,'(a)') trim(str)
G%f=G%f+sigma0 !shift G(i*w_n). G(i*w_n)=G(i*w_n)+sigma0
do i=1,G%n
   write(24,'(I4,7E18.9)') ind(i),aimag(G%z(i)),real(G%f(i)),aimag(G%f(i)),G%std(i,:),real(z(i))+sigma0,aimag(z(i))
enddo
close(24)
!Physical continautions have equal weights."
call faverage(weights(:,2),G%z,G%f,G%std) !calculate Pade G%f=G(i*w_n) and Pade standard deviation G%std=std[G(i*w_n]          
write(33,'(a,E15.6,a)') "|Gn-<Gn>_w|=sqrt(sum_n |Gn-<Gn>_w|^2 )=",&
  sqrt(sum(abs(G%f-z)**2)),". Physical continautions have equal weights."
open(unit=24,file=trim(f_im)//"_p",action="write")
write(24,'(a)') trim(str)
G%f=G%f+sigma0 !shift G(i*w_n). G(i*w_n)=G(i*w_n)+sigma0
do i=1,G%n
   write(24,'(I4,7E18.9)') ind(i),aimag(G%z(i)),real(G%f(i)),aimag(G%f(i)),G%std(i,:),real(z(i))+sigma0,aimag(z(i))
enddo
close(24)
!Weights inverse proportional to sqrt(sum(|G(i*w_n)-Gpade(i*w_n)|^2)). Consider only physical continautions
call faverage(weights(:,3),G%z,G%f,G%std) !calculate Pade G%f=G(i*w_n) and Pade standard deviation G%std=std[G(i*w_n]          
write(33,'(a,E15.6,a)') "|Gn-<Gn>_w|=sqrt(sum_n |Gn-<Gn>_w|^2 )=",&
  sqrt(sum(abs(G%f-z)**2)),&
  ". Weights inverse proportional to sqrt(sum(|G(i*w_n)-Gpade(i*w_n)|^2)). &
  Consider only physical continautions"
open(unit=24,file=trim(f_im)//"_devMats",action="write")
write(24,'(a)') trim(str)
G%f=G%f+sigma0 !shift G(i*w_n). G(i*w_n)=G(i*w_n)+sigma0
do i=1,G%n
   write(24,'(I4,7E18.9)') ind(i),aimag(G%z(i)),real(G%f(i)),aimag(G%f(i)),G%std(i,:),real(z(i))+sigma0,aimag(z(i))
enddo
close(24)
deallocate(G%z,G%f,G%std)
deallocate(ind,z)
call cpu_time(finish)
write(33,'(I10,a)') int(finish-start)," seconds to calculate and print averaged Matsubara values."
write(33,*)

!******************************************************************************
!Some other real axis averages
!Print default average on sparse real axis grid
str="#      w           <Re[G]>_w      <Im[G]>_w   -1/pi*<Im[G]>_w &
  std[Re[G]]_w   std[Im[G]]_w"
string=trim(str)//" Im[Gexact]"
if(verbosity>8) then
  G%n=2000
  allocate(G%z(G%n),G%f(G%n),G%std(G%n,2))
  do i=1,G%n
     G%z(i)=cmplx(-max(abs(estart),abs(estop))+(i-1)*2*max(abs(estart),abs(estop))/(G%n-1),aimag(emesh(1)))
  enddo
  call faverage(weights(:,1),G%z,G%f,G%std)  ! using default weights 
  if(ftype==4) then !symmetrize G: G(w+i*delta) = (G(w+i*delta)+G(-w+i*delta)^*)/2
     allocate(Gz(size(G%f))) 
     Gz=(G%f(:)+conjg(G%f(ne:1:-1)))/2
     G%f=Gz
     deallocate(Gz)
  endif
  G%f=G%f+sigma0
  open(unit=24,file=trim(f_re)//"_spare",action="write")
  write(24,'(a)') trim(str)
  do i=1,G%n
      write(24,'(6E16.5)') real(G%z(i)),real(G%f(i)),aimag(G%f(i)),-1/pi*aimag(G%f(i)),&
      G%std(i,1),G%std(i,2)
  enddo
  close(24)
  deallocate(G%z,G%f,G%std)
endif
call cpu_time(start)

allocate(x(ne))
G%n=ne
if(exact) allocate(G%err(G%n),G%aerr(4))
allocate(G%z(G%n),G%f(G%n),G%std(G%n,2))
G%z=emesh
if(verbosity>8) then
  call cpu_time(start)
  call faverage(weights(:,1),.true.,G%z,G%f,G%std) ! Default weights but put Im[P_r]=0 
  if(ftype==4) then !symmetrize G: G(w+i*delta) = (G(w+i*delta)+G(-w+i*delta)^*)/2
     allocate(Gz(size(G%f))) 
     Gz=(G%f(:)+conjg(G%f(ne:1:-1)))/2
     G%f=Gz
     deallocate(Gz)
  endif
  G%f=G%f+sigma0
  open(unit=24,file=trim(f_re)//"_ImPr0",action="write")
  if(exact) then
     G%err=1/pi*abs(aimag(G%f-fexact))
     G%aerr(1)=integ(real(G%z),G%err)
     G%aerr(2)=G%aerr(1)/integ(real(G%z),1/pi*abs(aimag(fexact)))
     G%aerr(3)=integ(real(G%z),G%err**2)
     G%aerr(4)=G%aerr(3)/integ(real(G%z),(1/pi*abs(aimag(fexact)))**2 )

     write(42,'(8E15.7,a)') G%aerr,-log(G%aerr)/log(10d0)," # acc + acc2 weights,but Im[P_r]=0"

     write(24,'(a)') trim(string)
     do i=1,G%n
         write(24,'(7E15.5)') real(G%z(i)),real(G%f(i)),aimag(G%f(i)),-1/pi*aimag(G%f(i)),&
         G%std(i,1),G%std(i,2),aimag(fexact(i))
     enddo
  else
     write(24,'(a)') trim(str)
     do i=1,G%n
        write(24,'(6E15.5)') real(G%z(i)),real(G%f(i)),aimag(G%f(i)),-1/pi*aimag(G%f(i)),&
        G%std(i,1),G%std(i,2)
     enddo
  endif
  close(24)
  call cpu_time(finish)
  write(33,'(I10,a)') int(finish-start)," seconds to calculate and print real axis average spectrum _re_ImPr0."
endif
deallocate(x)

if(verbosity>8 .and. sum(lambda)>0) then
  !Save continuation with the highest lambda 
  j=maxloc(lambda,1)
  write(33,'(a,E15.6,a,I4)') "max lambda= ",lambda(j),", from continuation nr:",j
  allocate(x(ns))
  x=0q0
  x(j)=1q0
  call faverage(x,G%z,G%f,G%std)         
  deallocate(x)
  if(ftype==4) then !symmetrize G: G(w+i*delta) = (G(w+i*delta)+G(-w+i*delta)^*)/2
     allocate(Gz(size(G%f))) 
     Gz=(G%f(:)+conjg(G%f(ne:1:-1)))/2
     G%f=Gz
     deallocate(Gz)
  endif
  G%f=G%f+sigma0
  open(unit=24,file=trim(trim(f_re)//"_higestlambda"),action="write")
  if(exact) then
     G%err=1/pi*abs(aimag(G%f-fexact))
     G%aerr(1)=integ(real(G%z),G%err)
     G%aerr(2)=G%aerr(1)/integ(real(G%z),1/pi*abs(aimag(fexact)))
     G%aerr(3)=integ(real(G%z),G%err**2)
     G%aerr(4)=G%aerr(3)/integ(real(G%z),(1/pi*abs(aimag(fexact)))**2 )
     write(42,'(8E15.7,a)') G%aerr,-log(G%aerr)/log(10d0)," # best lambda weights"

     write(24,'(a)') trim(string)
     do i=1,G%n
         write(24,'(7E15.5)') real(G%z(i)),real(G%f(i)),aimag(G%f(i)),-1/pi*aimag(G%f(i)),&
         G%std(i,1),G%std(i,2),aimag(fexact(i))
     enddo
  else
     write(24,'(a)') trim(str)
     do i=1,G%n
        write(24,'(6E15.5)') real(G%z(i)),real(G%f(i)),aimag(G%f(i)),-1/pi*aimag(G%f(i)),&
        G%std(i,1),G%std(i,2)
     enddo
  endif
  close(24)
endif
deallocate(G%z,G%f,G%std)
if(exact) then
  deallocate(G%err,G%aerr)
  close(42)
endif


!******************************************************
! Save file with columns: N M dis weights good dev lam dexact_1 dexact_2 dexact_3 dexact_4 F_1 F_2 F_3 F_4
call cpu_time(start)
open(unit=24,file=trim(f_weight),action="write")
if(exact .and. verbosity>5) then
  write(24,'(a)') "#  N    M  first  good    dis           weight_def      weights_phys    weights_Mats    &
  dev_matrix      dev_mats       lambda              err1           err2              err3              err4      &
             F1             F2             F3             F4"
  do j=1,ns
     write(24,'(4I5,15E16.6)') c(j)%N,c(j)%M,c(j)%first,good(j),dis(j),weights(j,:),dev(j,:),lambda(j),&
     dexact(j,:),-log(dexact(j,:))/log(10d0)
  enddo
else
  write(24,'(a)') "#  N    M  first  good    dis           weight_def      weights_phys    weights_Mats    &
  dev_matrix      dev_mats       lambda    " 
  do j=1,ns
    write(24,'(4I5,7E16.6)') c(j)%N,c(j)%M,c(j)%first,good(j),dis(j),weights(j,:),dev(j,:),lambda(j)
  enddo
endif
close(24)
call cpu_time(finish)
write(33,'(I10,a)') int(finish-start)," seconds to calculate and print _w file."

write(*,*) "after 1_w" 

!***************************************************
!save weighted averages of some observables, using default weights
!average lambda
if(sum(lambda)>0) then
  write(33,'(a,E12.3)') "<lambda>_w=",sum(lambda*weights(:,1))
endif
write(33,'(a,E12.3)') "<dis>_w=",sum(dis*weights(:,1)) !average dis, the deviation between solutions.
write(33,'(a,E12.3)') "<N>_w=",sum(c(:)%N*weights(:,1))  
write(33,'(a,E12.3)') "<M>_w=",sum(c(:)%M*weights(:,1)) 
write(33,'(a,E12.3)') "<M/N>_w=",sum(c(:)%M*weights(:,1)/c(:)%N) 
write(33,'(a,E12.3)') "<first>_w=",sum(c(:)%first*weights(:,1)) 
write(33,*)

!endif

!********************************************
!   deallocate
!********************************************
if(allocated(fexact))  deallocate(fexact)
if(allocated(dexact))  deallocate(dexact)
if(allocated(weights)) deallocate(weights)
if(allocated(good))    deallocate(good)
if(allocated(dis))     deallocate(dis)
if(allocated(emesh))   deallocate(emesh)
if(allocated(dev))     deallocate(dev)
if(allocated(lambda))  deallocate(lambda)
if(allocated(Nv))      deallocate(Nv)
if(allocated(c))       deallocate(c)
if(allocated(ac))      deallocate(ac)

!if(myrank==0) then
   write(33,*) 
   write(33,'(a)') "Pade program finished!"
   close(33)
!endif

!call MPI_Finalize(ierr)


end program

