module meshes
use input
implicit none

!configuration dependent matsubara mesh
integer,allocatable         ::  ind(:) !indices to matsubara points for particular configuration
complex(kind=16),allocatable :: wn(:),Gn(:) !matsubara freq
complex(kind=16),allocatable :: emesh(:) !w+i*delta
complex(kind=16),allocatable :: fexact(:) ! exact solution
real(kind=16),allocatable    :: dexact(:,:) ! deviation of real axis from exact solution

contains

!subroutine getrmesh(exact,eim,ne,f_ex,fexact,emesh,estart,estop)
subroutine getrmesh(ftype,exact,eim,ne,f_ex,estart,estop)
!Construct mesh eim above real axis.
!If exact=false, use ne points uniformly between estart and estop.
!If exact=true and ftype/=4 copy from exact function.
!If exact=true and ftype=4, copy from exact function. If only contains and mirror also to negative frequencies. 
!solution is stored in fexact. Also ne,estart and estop is then determined by exact solution mesh
implicit none
integer,intent(in) :: ftype
logical,intent(in) :: exact
real(kind=16),intent(in) :: eim
integer,intent(inout) :: ne
character(len=800),intent(in) :: f_ex
!complex(kind=16),intent(out),allocatable :: fexact(:)
!complex(kind=16),intent(inout),allocatable :: emesh(:)
real(kind=16),intent(inout) :: estart,estop

integer :: i,ok
real(kind=16) :: tmp,tmp2,tmp3

if(exact) then
    !use mesh from exact solution, so can easily compare
    open(54,file=trim(f_ex))
    i=0
    do
        read(54,*,iostat=ok) tmp
        if(ok>0) stop "read error"
        if(ok<0) exit !no more lines to read
        i=i+1
    end do
    rewind(54)
    if(ftype/=4) then
       ne = i !nbr of freq in file
       allocate(fexact(ne),emesh(ne))
       do i=1,ne
           read(54,*) tmp,tmp2,tmp3 
           emesh(i)=cmplx(tmp,eim,kind=16)
           fexact(i)=cmplx(tmp2,tmp3,kind=16)
       enddo
    elseif(ftype==4) then
      !In this case we assume the exact solution to have odd spectrum.
      !If the file with exact function has only have positive frequencies,
      !double the number of frequencies by mirroring the frequencies.
      !If the file with exact function has both positive and negative
      !frequencies just read the file 
      read(54,*) tmp
      rewind(54)
      if(tmp<0) then
         ne=i
         allocate(fexact(ne),emesh(ne))
         do i=1,ne
             read(54,*) tmp,tmp2,tmp3 
             emesh(i)=cmplx(tmp,eim,kind=16)
             fexact(i)=cmplx(tmp2,tmp3,kind=16)
         enddo
      else
         ne=2*i !nbr of freq in mirrored mesh
         allocate(fexact(ne),emesh(ne))
         do i=1,ne/2
             read(54,*) tmp,tmp2,tmp3 
             emesh(ne/2+1-i)=cmplx(-tmp,eim,kind=16)
             emesh(ne/2+i)=cmplx(tmp,eim,kind=16)
             fexact(ne/2+1-i)=cmplx(tmp2,-tmp3,kind=16)
             fexact(ne/2+i)=cmplx(tmp2,tmp3,kind=16)
         enddo
      endif
    else
      stop "what the crack?"
    endif
    close(54)
    estart=real(emesh(1))
    estop=real(emesh(ne))
else
   if(estart>estop) stop "estart has to be smaller than estop"
   if(ne<=0) stop "positive number of points needed in real mesh"
    ! create mesh from estart and estop and ne
    allocate(emesh(ne))
    do i=1,ne       !Evaluate distance eim above real axis
      emesh(i) = cmplx(estart+(i-1)*(estop-estart)/(ne-1),eim,16)
    enddo
endif
end subroutine

subroutine getmmesh(filen,wmax,startindex,ind,wn,Gn)
! Get the Matsubara frequencies below wmax and corresponding values.
character(len=800),intent(in) :: filen
real(kind=16),intent(in) :: wmax
integer :: startindex
integer,allocatable,intent(out) :: ind(:)
complex(kind=16),allocatable,intent(out) :: wn(:),Gn(:)

integer :: i,j,ok,n
real(kind=16) :: tmp,tmp2,tmp3

!calculate how many points below wmax
OPEN(11,FILE=trim(filen))
i=1
j=1
read(11,*,iostat=ok) tmp
if(ok>0) stop "read error"
tmp2=abs(tmp-wmax)
do
    read(11,*,iostat=ok) tmp
    if(ok>0) stop "read error"
    if(ok<0) exit !no more lines to read
    i=i+1
    tmp3=abs(wmax-tmp)
    if(tmp3<tmp2) then
        j=i
    endif
    tmp2=tmp3
end do
close(11)
n=j !nbr of positive points below wmax

!allocate ind,wn,Gn
if(allocated(ind)) deallocate(ind)
if(allocated(wn)) deallocate(wn)
if(allocated(Gn)) deallocate(Gn)
if(startindex<=0) then
   allocate(ind(n+abs(startindex)+1))
   allocate(wn(n+abs(startindex)+1))
   allocate(Gn(n+abs(startindex)+1))
else
   allocate(ind(n))
   allocate(wn(n))
   allocate(Gn(n))
endif

!get indices and corresponding frequency and function values of these points
call getFirstIndices(filen,startindex,ind,wn,Gn)

end subroutine

end module 
