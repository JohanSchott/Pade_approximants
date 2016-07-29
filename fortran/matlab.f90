module matlab
implicit none

real(kind=16),parameter  :: pi=3.14159265359

contains

integer function filelen(filen) !Counts number of rows in file filen
implicit none
character(len=800),intent(in) :: filen
integer :: i,ok
character(len=800) :: tmp
open(77,file=trim(filen))
!write(*,*) trim(filen)
i=0
do
   read(77,*,iostat=ok) tmp
   if(ok>0) stop "read error"
   if(ok<0) exit !no more lines to read
   i=i+1
end do
close(77)
filelen = i
end function

character(len=800) function int2str(j)
implicit none
character(len=800) :: jstr
integer,intent(in) :: j
if(j<10) then
    write(jstr,'(I1)') j
elseif(j<100) then
    write(jstr,'(I2)') j
elseif(j<1000) then
    write(jstr,'(I3)') j
elseif(j<10000) then
    write(jstr,'(I4)') j
else
    write(jstr,'(I5)') j
endif
int2str=jstr
end function

real(kind=16) function integ(x,y)
!integrate y as function of x
implicit none
real(kind=16),intent(in) :: x(:),y(:)
real(kind=16),allocatable :: f(:)
integer :: n
n=size(x,1)
allocate(f(n))
f(1)=(x(2)-x(1))/2
f(n)=(x(n)-x(n-1))/2
f(2:n-1)=(x(3:n)-x(1:n-2))/2
integ=sum(f*y)
deallocate(f)
end function

subroutine interp(w,f,e,g)
!Linear interpolate values (w,f) to points e and return in g
!If a point e(i) is outside window of w, put g(i) to endpoint value of f 
!The two points closest to e(i) in w, w(j),w(j+1), are used to calculate g(i) according to:
!g(i) = f(j)+(f(j+1)-f(j))*(e(i)-w(j))/(w(j+1)-w(j))
implicit none
real(kind=16),intent(in) :: w(:),f(:),e(:)
real(kind=16),allocatable :: g(:)

integer :: n,m,i,j
n=size(w,1)
m=size(e,1)
if(allocated(g)) deallocate(g)
allocate(g(m))
!interpolate one point at a time: e(i) -> g(i)
do i=1,m
   if(e(i)<w(1)) then
      g(i)=f(1)
   elseif( e(i)>w(n) ) then
      g(i)=f(n)
   else
      !find w(j) and w(j+1), thus find index j
      j=minloc(abs(w-e(i)),1)
      if(w(j)>e(i) .or. j==n) j=j-1
      g(i)=f(j)+(f(j+1)-f(j))*(e(i)-w(j))/(w(j+1)-w(j))
   endif
enddo
end subroutine

function sort(A) !bubble sort algorithm, smallest first.
implicit none
real(kind=16),allocatable :: sort(:)
real(kind=16),intent(in) :: A(:)

real(kind=16) :: b
integer :: i,n,newn
n=size(A)
allocate(sort(n))
sort=A
do while(n>1) 
    newn=0
    do i=1,n-1 
        if(sort(i)>sort(i+1)) then
            b=sort(i)
            sort(i) = sort(i+1)
            sort(i+1) = b
            newn=i
        endif
    enddo
    n=newn
enddo
end function

function sortp(A) !bubble sort algorithm, sort with smallest first. Returns permutation array
! Eg: A=[ 3.1  1.3  2.5  1.3 ]    ->    A=[1.3 1.3 2.5 3.1] and sortp=[2 4 3 1]  so A=A_old(c)
implicit none
integer,allocatable :: sortp(:)
real(kind=16),intent(in) :: A(:)

integer,allocatable :: d(:),e(:)
real(kind=16),allocatable :: x(:),y(:)
integer :: i,j,k,n,m
n=size(A,1)
allocate(x(n),y(n),d(n),sortp(n))
x=A
y=sort(x) !sorted list
d=[ (i,i=1,n )  ]
do i=1,n !look which indices connect unsorted and sorted array
   j=minloc(abs(y(i)-x(d)),1)
   sortp(i)=d(j)
   !remove d(j) from d
   if(size(d,1)==1) then
      deallocate(d)
   else
      allocate(e(size(d,1)))
      e=d
      deallocate(d)
      allocate(d(size(e,1)-1))
      m=0
      do k=1,size(e,1)
         if(k/=j) then
            m=m+1
            d(m)=e(k)
         endif
      enddo
      deallocate(e)
   endif
enddo
deallocate(x,y)
end function

end module
