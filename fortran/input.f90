module input
implicit none

interface getIndices
   module procedure getFirstIndices,getFirstplusuniformIndices, getParticularIndices
end interface getIndices

contains
    
subroutine getFirstIndices(inputfile,startindex,ind,wn,Gn)
     implicit none
     character(len=800),intent(in) :: inputfile
     !integer,intent(in) :: n
     integer,intent(in) :: startindex
     integer,intent(out) :: ind(:)
     complex(kind=16),intent(out) :: wn(:),Gn(:) !matsubara points and Green's function
   
    real(kind=16),allocatable :: v(:,:) !stores all ntotal matsubara points
    integer :: i,ok !dummy index 
    integer :: ntotal  !nbr of elements in input file
    real(kind=16) :: tmp,tmp2,tmp3
    integer :: n
    
    n=size(ind,1)
    open(77,file=inputfile)
    i=0
    do 
      read(77,*,iostat=ok) tmp
      if(ok>0) stop "error in reading"
      if(ok<0) exit
      i=i+1
    enddo
    ntotal=i
    if(ntotal<n) stop "Too many wanted input points"
    !reading data to vector
    rewind(77)
    allocate(v(ntotal,3)) 
    do i=1,ntotal
       read(77,*) tmp,tmp2,tmp3
       v(i,:)=[ tmp,tmp2,tmp3 ]
    enddo
    close(77)

    do i=1,n
       ind(i)=startindex+i-1
    enddo

    do i=1,n
      if(ind(i)<=0) then
         wn(i)=cmplx(0q0,-v(abs(ind(i))+1,1),kind=16)
         Gn(i)=cmplx(v(abs(ind(i))+1,2),-v(abs(ind(i))+1,3),kind=16)
      else
         wn(i)= cmplx(0q0,v(ind(i),1),kind=16) 
         Gn(i)= cmplx(v(ind(i),2),v(ind(i),3),kind=16) 
      endif
    enddo

    deallocate(v)
end subroutine getFirstIndices

subroutine getFirstplusuniformIndices(inputfile,startindex,q,ind,wn,Gn)
! Picks ns=floor(n*q) subsequent points, starting at startindex.
! So these indices are:
! startindex,startindex+1,startindex,...,startindex+ns-1
! The remaining nu=n-ns points are distributed uniformly between
! indices: i=startindex+ns and f=i+ns-1
! There indices are:
! x_j = floor(i+(f-i)/(nu+1)*j) , with j \in {1,2,...,nu}
implicit none
character(len=800),intent(in) :: inputfile
!integer,intent(in) :: n  !nbr of points to pick
integer,intent(in) :: startindex ! first index in file to use
real(kind=16),intent(in) :: q ! floor(n*q) is the number subsequent points
integer,intent(out) :: ind(:) ! list with indices to use 
complex(kind=16),intent(out) :: wn(:),Gn(:) !matsubara points and Green's function

real(kind=16),allocatable :: v(:,:) !stores all ntotal matsubara points
real(kind=16) :: tmp,tmp2,tmp3
integer :: ns,nu,i,k,f,ok !dummy index 
integer :: ntotal  !nbr of elements in input file
integer :: n  !nbr of points to pick

n=size(ind,1)
!count number of lines in input file
open(77,file=inputfile)
i=0
do 
   read(77,*,iostat=ok) tmp
   if(ok>0) stop "error in reading"
   if(ok<0) exit
   i=i+1
enddo
ntotal=i
close(77)
if(ntotal<n) stop "Too many wanted input points"
!reading data to vector
allocate(v(ntotal,3)) 
open(unit=78,file=inputfile)
do k=1,ntotal
   READ(78,*) tmp,tmp2,tmp3
   v(k,:)=[ tmp,tmp2,tmp3 ]
enddo
close(78)
!get first indices 
ns=floor(n*q)
do k=1,ns
   ind(k) = startindex+k-1
enddo
nu=n-ns
i=startindex+ns  ! initial index, lower bound
f=i+ns-1  ! final index, upper bound
!get the remaining indices
do k=1,nu
   ind(ns+k)=floor(i+(1d0*(f-i))/(nu+1)*k)
enddo
!save wn and Gn
do k=1,n
   if(ind(k)<=0) then
      wn(k)=cmplx(0q0,-v(abs(ind(k))+1,1),kind=16)
      Gn(k)=cmplx(v(abs(ind(k))+1,2),-v(abs(ind(k))+1,3),kind=16)
   else
      wn(k)=cmplx(0q0,v(ind(k),1),kind=16)
      Gn(k)=cmplx(v(ind(k),2),v(ind(k),3),kind=16)
   endif
enddo
deallocate(v)
end subroutine getFirstplusuniformIndices

subroutine getParticularIndices(inputfile,startindex,q,wsplit,wend,p,ind,wn,Gn)     
    implicit none
    character(len=800),intent(in) :: inputfile !The matsubara input file
    !integer,intent(in) :: n  !nbr of matsubara points to find
    integer,intent(in) :: startindex !the first index to use in pade
    real(kind=16),intent(in) :: q,wsplit,wend !procentage of points below wsplit. wend is highest considered freq
    integer,intent(in) :: p  !polynomial order of matsubara point distribution.
    integer,intent(out) :: ind(:)  !indices of the matsubara points
    complex(kind=16),intent(out) :: wn(:),Gn(:) !matsubara points and Green's function
    
    integer :: ntotal  !nbr of elements in input file
    integer :: n  !nbr of matsubara points to find
    real(kind=16),allocatable :: v(:) !stores all ntotal matsubara points
    real(kind=16),allocatable :: vs(:) 
    real(kind=16) :: vstart
    integer :: ncut,nncut !ncut is a index, correspoinding to wenda. nncut is a index, corresponding to wsplit
    integer :: ned,nlength !nbr of uniform and non-uniformly spaced points

    integer :: ok,i,j !dummy indices
    real(kind=16) :: tmp,tmp2,tmp3
    real(kind=16) :: delta
      
    n=size(ind)
    open(77,file=inputfile)
    i=0
    do 
      read(77,*,iostat=ok) tmp
      if(ok>0) stop "error in reading"
      if(ok<0) exit
      i=i+1
    enddo
    ntotal=i
    close(77)

    if(ntotal<n) stop "Too many wanted input points"

    allocate(v(ntotal)) 
    !reading data to vector
    open(unit=78,file=inputfile)
    do i=1,ntotal
      READ(78,*) v(i)
    enddo
    close(78)

    nncut = minloc(abs(v-wsplit),1)
    ncut = minloc(abs(v-wend),1)
    
    if( v(nncut) < 0 .or. v(ncut) < 0) stop "Can not select good matsubara mesh"

    ned = floor(q*n) !Number of evenly spaced points
    nlength = n-ned !Number of unevenly distributed points

    if(ned+1 >= nncut) then
        nncut=ned+1
    endif

    ! get indices of matsubara points to be used in pade approximants 
    allocate(vs(ned))
    vstart=v(startindex)
    !delta = (v(nncut)-v(1))/(ned-1) 
    delta = (v(nncut)-vstart)/(ned-1) 
    do i=1,ned
        !vs(i) = v(1)+(i-1)*delta
        vs(i) = vstart+(i-1)*delta
        if(i==1) then
            ind(i) = minloc(abs(vs(i)-v),1)
            if(ind(i) /= startindex ) stop "something wrong"
        else
            ind(i) = max(minloc(abs(vs(i)-v),1),ind(i-1)+1)
        endif
    enddo 
    deallocate(vs)
    allocate(vs(nlength))
    delta = (v(ncut)-v(nncut+1))/(nlength-1)**p
    do i=1,nlength
        vs(i) = v(nncut+1)+(i-1)**p*delta
        !ind(ned+i) = minloc(abs(vs(i)-v),1)
        !ind(ned+i) = ind(ned+i-1)+1
        ind(ned+i) = max(minloc(abs(vs(i)-v),1),ind(ned+i-1)+1)
    enddo
    deallocate(vs)

    !checks
    if(ind(1) .ne. startindex) stop "first index should not correct"
    if(ind(ned) /= nncut) stop "Too many points in uniform region"
    if(ind(ned) .ne. nncut) stop "middle index wrong"
    if(ind(n) .ne. ncut) stop "last index wrong"
    do i=1,n
        if(ind(i) .eq. 0) stop "zeros indices not allowed"
    enddo
    do i=1,n-1
        if(ind(i+1) .eq. ind(i)) stop "two input points are the same"
    enddo
    ! return also wn and G(wn)
    open(unit=79,file=inputfile)
    j=1
    do i=1,ntotal
        read(79,*) tmp,tmp2,tmp3
        if(i .eq. ind(j)) then
            wn(j) = cmplx(0q0,tmp,16)
            Gn(j) = cmplx(tmp2,tmp3,16)
            if(j .eq. n) exit
            j=j+1
        endif
    enddo
    close(79)
   deallocate(v)
end subroutine getParticularIndices

end module
