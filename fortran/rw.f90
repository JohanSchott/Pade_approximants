module rw
use matlab
implicit none

character(len=800)  ::  f_in  !input variable required when calling binary
character(len=800)  ::  f_weight,f_done,f_info,f_re,f_std,f_im,f_pole !file names
character(len=800)  ::  f_mats_s !Name of shifted Matsubara file.
character(len=800)  ::  f_mats_c !Needed of shifted Matsubara file for c++

!varibles read from input file
character(len=800)  :: proj,projc
character(len=800)  :: f_mats
integer             :: ftype !0 for fermions, 1 for self-energy, 2 for bosonic spectrum, 3 for weird fermion and 4 for odd bosonic spectrum 
integer             :: verbosity !Determines how much is printed to files (1 little, 10 much)
logical             :: exact !if exact solution exists. If .true. mesh for exact solution is used.
logical             :: real2 !If do analytic continuation from real axis to real and Matsubara axis
character(len=800)  :: f_ex  ! Needed if has exact solution to compare with.
logical             :: mpack
integer             :: rep
integer             :: solver !1 for Least Square solution, 2 for Tikhonov solution 
integer             :: mtr ! Either 0,1,2 or 3. For npade=.false. and mpack=.true, need mtr=1.
logical             :: npade,oddc  !If N==M and if mod(M,2)==1
logical             :: shifted
integer             :: Nmax,Nmin,Nstep,Mmin,Mmax,Mstep !pade settings
integer             :: ne !nbr of points on real axis
real(kind=16)       :: eim,estart,estop  !real axis variables
integer             :: pickm
real(kind=16)       :: wmax,wsplit,q  !matsubara axis variables
integer             :: p !polynomial order of matsubara tail
integer             :: firstmin,firstmax,firststep !first index to use in input file 
real(kind=16)       :: acc !procentage of average deviation needed to be smaller than to contribute in the average
real(kind=16)       :: acc2 !procentage of number of continuations to include in the average
integer*8           :: pr  !precision

contains

subroutine readsettings()
implicit none
!read(*,'(a)') f_in !Input file to read all parameter values from.
f_in="pade.inp"
!Read from in-file
open(22,file=trim(f_in))
read(22,*) 
read(22,*) 
read(22,'(a)') proj
read(22,*) 
read(22,'(a)') f_mats
read(22,*)
read(22,*) ftype
read(22,*) 
read(22,*) verbosity
read(22,*) 
read(22,*) exact
read(22,*) 
read(22,*) real2
read(22,*) 
read(22,'(a)') f_ex
read(22,*) 
!if( .not. exact) then
    read(22,*) 
    read(22,*) estart
    read(22,*) 
    read(22,*) estop
    read(22,*) 
    read(22,*) ne
!endif
read(22,*) 
read(22,*) eim
read(22,*) 
read(22,*) 
read(22,*) mpack 
read(22,*)
read(22,*) rep
read(22,*)
read(22,*) solver
read(22,*) 
read(22,*) mtr
read(22,*) 
read(22,*) shifted
read(22,*) 
read(22,*) npade
read(22,*)
read(22,*) oddc
read(22,*) 
read(22,*) Nmax
read(22,*) 
read(22,*) Nmin
read(22,*) 
read(22,*) Nstep
read(22,*) 
read(22,*) Mmin
read(22,*) 
read(22,*) Mmax
read(22,*) 
read(22,*) Mstep
read(22,*) 
read(22,*) 
read(22,*) firstmin,firstmax,firststep
read(22,*) 
read(22,*) pickm
read(22,*) 
read(22,*) wsplit
read(22,*) 
read(22,*) wmax
read(22,*)
read(22,*) p
read(22,*)
read(22,*) q 
read(22,*)
read(22,*)
read(22,*) acc
read(22,*)
read(22,*) acc2
read(22,*)
read(22,*) pr
close(22)

write(*,'(A)') "Done reading "//trim(f_in)
end subroutine

subroutine saveusedinput()
implicit none
!Write parameters read and used from input file  
open(23,file=trim(f_info))
write(23,'(a)') "Program version: 4.00"
write(23,'(a)') "[29/06/2015] Store also spectrums, making program much faster"
write(23,'(a)') "[22/06/2015] Now calculate also average of Pade coefficients, not only spectrums"
write(23,'(a)') "[05/06/2015] MPI version, but had to sacrifice beyond quad precison possibility"
write(23,'(a)') "[02/06/2015] Switch parameter, shifted, for shifting real part or not"
write(23,'(a)') "[02/06/2015] Odd nbr of Pade coefficients possible with oddc input parameter "
write(23,'(a)') "[01/06/2015] Better estimation of sigma0 in: f(z) \approx sigma0-a/wn^2"
write(23,'(a)') "[10/05/2015] Cleaned up output-files"
write(23,'(a)') "[09/05/2015] Bugg fixed. Had unwanted quad to double precision conversion in rutine for finding poles"
write(23,'(a)') "[09/05/2015] Better definitions of deviations (on Matsubara and on matrix problem)"
write(23,'(a)') "[26/03/2015] Average over several first Matsubara index" 
write(23,'(a)') "[26/03/2015] In new file: *_wa save sum(abs(Gn-Gn_pade_average))"
write(23,'(a)') "[24/03/2015] Lapack's dgelsd routine implemented in quad precision"
write(23,'(a)') "[03/03/2015] Deviation on matsubara axis for mpack implemented"
write(23,'(a)') "[19/02/2015] Availability to use real value formulation to find Pade coefficients."
write(23,'(a)') "[17/02/2015] Store Pade coefficients in memory instead of disk for speed up."
write(23,'(a)') "[17/02/2015] Use fact: G(-i*w)=G(i*w)^* to add negative Matsubara points to fitting. Use by firstmin<=0"
write(23,'(a)') "[16/02/2015] sum_i |G(i)-G(i)_pade|^2 seams to be quality parameter."
write(23,'(a)') "[15/02/2015] Tikhonov dlb and quad precision solution available."
write(23,'(a)') "[15/02/2015] Untrucated and truncated SVD LS double precision solution available."
write(23,'(a)') "[05/02/2015] Uniform average of all physical continuations on the diagonal are now stored in 1_re_diag." 
write(23,'(a)') "[03/02/2015] Calculate |A*x-b|^2/|b|^2 (in high precision). Useful for evaluating accuarcy of matrix solvers."
write(23,'(a)') "[03/02/2015] Normal SVD Lapack routines zgelss and zgelsd with double precsion are now available."
write(23,'(a)') "[28/01/2015] Two errors estimates calculated when comparing to exact solution."  
write(23,'(a)') "[27/01/2015] Normal Lapack routine zgels with double precision is now available."
write(23,'(a)') "[26/01/2015] Easy comparison with exact solution is now implemented, one more input variable."
write(23,'(a)') "[21/01/2015] Memory leakage in c++ routine solved. Had wrong deallocate syntax."
write(23,'(a)') "[18/01/2015] Add average with weights inversly proportional to deviation to matsubara points."
write(23,'(a)') "[15/01/2015] Continuation now done on function with Re[G(i*w_n)]=0 for last Matsubara point."
write(23,'(a)') "[10/01/2015] Possibility to also pick first matsubara points and tail with uniform mesh added."
write(23,'(a)') "[09/01/2015] Program also saves average of continuations with Im[P_r]=0, with weights unchanged."
write(23,'(a)') "[08/01/2015] Possiblity to use transpose multiplied matrix for both mpack and lapack."
write(23,'(a)') "[08/01/2015] Possiblity to use mpack or lapack."
write(23,*)
write(23,*) '#Parameter values read and used from input file:'
write(23,*) '#inputfile: ',trim(f_mats)
write(23,*) '#ftype: ',ftype
write(23,*) '#verbosity=',verbosity
write(23,*) '#exact=',exact
write(23,*) '#real2=',real2
write(23,*) '#f_ex=',trim(f_ex)
write(23,*) '#eim=', eim
write(23,*) '#mpack=',mpack
if(.not. mpack) then
   write(23,*) '#solver=',solver
   write(23,*) '#rep=',rep
endif
write(23,*) '#mtr=',mtr
write(23,*) '#shifted=',shifted
write(23,*) '#npade=',npade
write(23,*) '#oddc=',oddc
if( .not. exact ) then
    write(23,*) '#estart=',estart
    write(23,*) '#estop=',estop
    write(23,*) '#ne=', ne
endif
write(23,*) '#Nmax=', Nmax
write(23,*) '#Nmin=', Nmin
write(23,*) '#Nstep=', Nstep
if( .not. npade ) then
    write(23,*) '#Mmin=', Mmin
    write(23,*) '#Mmax=', Mmax
    write(23,*) '#Mstep=', Mstep
endif
write(23,'(a,I3,a,I3,a,I3)') '#firstmin=', firstmin," firstmax=",firstmax," firststep=",firststep
write(23,*) '#pickm=',pickm
if(pickm==1) then
   write(23,*) '#q=',q
elseif(pickm==2) then
   write(23,*) '#wsplit=', wsplit
   write(23,*) '#wmax=', wmax
   write(23,*) '#p=',p
   write(23,*) '#q=',q
endif
write(23,*) '#acc=',acc
write(23,*) '#acc2=',acc2
write(23,*) '#pr=', pr
write(23,*) " "
write(23,*) "Output information:"
close(23)
end subroutine

subroutine setfilenames()
use iso_c_binding, only: C_CHAR, C_NULL_CHAR
implicit none
!files to use
f_weight=trim(proj)//"_w"
f_done=trim(proj)//"_done"
f_info=trim(proj)//"_info"
f_re=trim(proj)//"_re"
f_std=trim(proj)//"_std"
f_im=trim(proj)//"_im"
f_pole=trim(proj)//"_pole"
f_mats_s=trim(proj)//"_mats_shift"
f_mats_c = C_CHAR_""//trim(f_mats_s)//CHAR(0) 
projc=C_CHAR_""//trim(proj)//CHAR(0)
end subroutine

subroutine checksettings()
implicit none
if(mpack .and. rep .ne. 1) stop "real value formulation for mpack is not implemented yet!"
if( .not. (solver==1 .or. solver==2)) stop "LS or Tikhonov solver"
if( .not. (mtr==0 .or. mtr==1 .or. mtr==2 .or. mtr==3) ) then
   stop "mtr needs to be either 0,1,2 or 3"
endif

if(firstmin<-Nmax) stop "firstmin can't be smaller than -Nmax."

if(.not. npade ) then
    if(Mmin>Nmin .or. Mmax>Nmax) stop "Need: Mmin<Nmin and Mmax>Nmax."
    if(Mmin>Mmax .or. Nmin>Nmax) stop "Need: Mmin<Mmax and Nmin<Nmax."
    if(oddc) then
       if(mod(Mmin,2)==0) stop "Mmin has to be odd"
    else
       if(mod(Mmin,2)/=0) stop "Mmin has to be even"
       !if(mod(Mmax,2)/=0) stop "Mmax has to be even"
    endif
    if(mod(Mstep,2)/=0) stop "Mstep has to be even"
else
    Mmin=Nmin
    Mmax=Nmax
    Mstep=Nstep
    if(oddc) then
       if(mod(Mmin,2)==0) stop "Nmin has to be odd"
    else
       if(mod(Mmin,2)/=0) stop "Nmin has to be even"
    endif
    if(mod(Mstep,2)/=0) stop "Nstep has to be even"
endif
call checkmatsmesh()
call checkrealmesh()
end subroutine

subroutine checkmatsmesh()
implicit none
integer :: ntotal
integer :: ok,i
real(kind=16),allocatable ::  x(:)
ntotal=filelen(f_mats) !number of rows in f_mats
!write(*,*) "ntotal=",ntotal
!write(*,*) "firstmin=",firstmin
!write(*,*) "Nmax-1=",Nmax-1
if(ntotal<firstmin+Nmax-1) stop "Input file contains fewer points than nbr of points wanting in pade."
allocate(x(ntotal))
open(11,file=trim(f_mats))
do i=1,ntotal
    read(11,*,iostat=ok) x(i)
enddo
close(11)
if(pickm==0) then
   wmax=x(firstmin-1+Nmax)
   if(x(ntotal)<wmax) then
      stop "Too high Nmax"
   endif
elseif(pickm==1) then
   wmax=2*x(firstmin-1+floor(Nmax*q))-x(firstmin) 
   if(x(ntotal)<wmax) then
      stop "Too high Nmax"
   endif
   if( q<0.51) stop "q>0.51 is needed"
elseif(pickm==2) then
   if(x(firstmin-1+Nmax)>wmax) stop "Too low wmax value."
   if(x(firstmin-1+floor(Nmax*q))>wsplit) stop "Too low wsplit or too high q"
   wmax = min(x(ntotal),wmax)  !highest considered frequency, if not only picking the first points in input file
   open(31,file=trim(f_info),position="append")
   write(31,*) "Matsubara points picked uniformly below:", wsplit
   write(31,*) "Highest considered freq:",wmax
   close(31)
endif
deallocate(x)
write(*,*) "Successfully checked Matsubara file:",trim(f_mats)
end subroutine

subroutine checkrealmesh()
implicit none
if(eim<=0) then
   stop "Need to evaluate above real axis." 
end if
if(estop<=estart) then
   stop "Real mesh ill defined."
end if
if(ne<=0) then
   stop "Need positive number of real mesh points."
end if
end subroutine

end module 
