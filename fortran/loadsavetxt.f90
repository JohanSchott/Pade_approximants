module loadsavetxt
! This module contains interfaces openf and save2f.
! openf open file and save data to real variables of quadruple precision.
! If the header rows start with "#", they are treated as comments. 

implicit none

interface openf
    module procedure open1,open2,open3,open4
end interface

interface save2f
    module procedure save1r2f,save1c2f,save2r8_2f,save2r2f,save2c2f,save3r8_2f,save3r2f,save4r2f,save5r2f
end interface

contains

subroutine savefout(zout,fout)
    implicit none
    complex(kind=16),intent(in) :: zout(:)
    complex(kind=8),intent(in)  :: fout(:)
    real(kind=8) :: pi = 3.1415
    call save2f('pade_fout',real(zout,kind=8),real(fout),aimag(fout))
    call save2f('pade_A',real(zout,kind=8),-1d0/pi*aimag(fout))
end subroutine

subroutine readzoutdata(zout)
    implicit none
    complex(kind=16),allocatable,intent(out) :: zout(:)
    real(kind=16) :: ws,wf,eim
    integer :: N
    real(kind=16),allocatable :: w(:)
    integer :: i
    open(80,file='pade.par')
    read(80,*) ! read the header 
    read(80,*) ws,wf,N
    read(80,*) eim
    close(80)
    allocate(w(N),zout(N))
    w = (/ (ws+(wf-ws)/(N-1)*(i-1),i=1,N) /)
    zout = cmplx(w,eim,kind=16)
    deallocate(w)
end subroutine

subroutine readinputdata(zin,fin)
    implicit none
    complex(kind=16),allocatable,intent(out) :: zin(:),fin(:)
    real(kind=16),allocatable :: z_re(:),z_im(:)
    real(kind=16),allocatable :: re(:),im(:)
    call openf('pade.in',z_re,z_im,re,im)
    allocate(zin(size(z_re)),fin(size(z_re)))
    zin = cmplx(z_re,z_im,kind=16)
    fin = cmplx(re,im,kind=16)
end subroutine

integer function filelen(filen) !Counts number of rows in file filen
    implicit none
    character(len=*),intent(in) :: filen
    integer :: i,ok
    character(len=800) :: tmp
    open(70,file=trim(filen))
    i=0
    do
       read(70,*,iostat=ok) tmp
       if(ok>0) stop "read error"
       if(ok<0) exit !no more lines to read
       i=i+1
    end do
    close(70)
    filelen = i
end function

integer function readstart(filen)
    implicit none
    character(len=*),intent(in) :: filen
    integer :: i,n
    character(len=800)  :: s
    n=filelen(filen)
    open(70,file=trim(filen))
    do i=1,n
        read(70,'(a)') s
        if ( s(1:1) /= "#" ) then
            close(70)
            readstart = i
            return  
        endif
    enddo
end function

subroutine open1(filen,x)
    implicit none
    character(len=*),intent(in) :: filen
    real(kind=16),allocatable,intent(out) :: x(:)
    integer :: k,i,j,ok
    k=readstart(filen)
    i=filelen(filen)
    if(allocated(x)) deallocate(x)
    allocate(x(i-k+1))
    open(70,file=trim(filen),iostat=ok)
    do j=1,k-1
        read(70,*)
    enddo
    do j=1,i-k+1
        read(70,*) x(j)
    enddo
    close(70)
end subroutine

subroutine open2(filen,x,y)
    implicit none
    character(len=*),intent(in) :: filen
    real(kind=16),allocatable,intent(out) :: x(:),y(:)
    integer :: k,i,j,ok
    k=readstart(filen)
    i=filelen(filen)
    if(allocated(x)) deallocate(x)
    if(allocated(y)) deallocate(y)
    allocate(x(i-k+1),y(i-k+1))
    open(70,file=trim(filen),iostat=ok)
    do j=1,k-1
        read(70,*)
    enddo
    do j=1,i-k+1
        read(70,*) x(j),y(j)
    enddo
    close(70)
end subroutine

subroutine open3(filen,x,y,z)
    implicit none
    character(len=*),intent(in) :: filen
    real(kind=16),allocatable,intent(out) :: x(:),y(:),z(:)
    integer :: k,i,j,ok
    k=readstart(filen)
    i=filelen(filen)
    if(allocated(x)) deallocate(x)
    if(allocated(y)) deallocate(y)
    if(allocated(z)) deallocate(z)
    allocate(x(i-k+1),y(i-k+1),z(i-k+1))
    open(70,file=trim(filen),iostat=ok)
    do j=1,k-1
        read(70,*)
    enddo
    do j=1,i-k+1
        read(70,*) x(j),y(j),z(j)
    enddo
    close(70)
end subroutine

subroutine open4(filen,a,b,c,d)
    implicit none
    character(len=*),intent(in) :: filen
    real(kind=16),allocatable,intent(out) :: a(:),b(:),c(:),d(:)
    integer :: k,i,j,n,ok
    k=readstart(filen)
    i=filelen(filen)
    if(allocated(a)) deallocate(a)
    if(allocated(b)) deallocate(b)
    if(allocated(c)) deallocate(c)
    if(allocated(d)) deallocate(d)
    n = i-k+1
    allocate(a(n),b(n),c(n),d(n))
    open(70,file=trim(filen),iostat=ok)
    do j=1,k-1
        read(70,*)
    enddo
    do j=1,n
        read(70,*) a(j),b(j),c(j),d(j)
    enddo
    close(70)
end subroutine

subroutine save1r2f(filen,x)
    implicit none
    character(len=*),intent(in) :: filen
    real(kind=16),intent(in) :: x(:)
    integer :: i,ok
    open(70,file=trim(filen),iostat=ok)
    do i=1,size(x,1)
       write(70,'(E13.4)') x(i)
    enddo
    close(70)
end subroutine

subroutine save1c2f(filen,x)
    implicit none
    character(len=*),intent(in) :: filen
    complex(kind=16),intent(in) :: x(:)
    integer :: i,ok
    open(70,file=trim(filen),iostat=ok)
    do i=1,size(x,1)
       write(70,'(2E13.4)') real(x(i)),aimag(x(i))
    enddo
    close(70)
end subroutine

subroutine save2r8_2f(filen,x,y)
    implicit none
    character(len=*),intent(in) :: filen
    real(kind=8),intent(in) :: x(:),y(:)
    integer :: i,ok
    open(70,file=trim(filen),iostat=ok)
    do i=1,size(x,1)
       write(70,'(2E13.4)') x(i),y(i)
    enddo
    close(70)
end subroutine

subroutine save2r2f(filen,x,y)
    implicit none
    character(len=*),intent(in) :: filen
    real(kind=16),intent(in) :: x(:),y(:)
    integer :: i,ok
    open(70,file=trim(filen),iostat=ok)
    do i=1,size(x,1)
       write(70,'(2E13.4)') x(i),y(i)
    enddo
    close(70)
end subroutine

subroutine save2c2f(filen,x,y)
    implicit none
    character(len=*),intent(in) :: filen
    complex(kind=16),intent(in) :: x(:),y(:)
    integer :: i,ok
    open(70,file=trim(filen),iostat=ok)
    do i=1,size(x,1)
       write(70,'(4E13.4)') real(x(i)),aimag(x(i)),real(y(i)),aimag(y(i))
    enddo
    close(70)
end subroutine

subroutine save3r8_2f(filen,x,y,z)
    implicit none
    character(len=*),intent(in) :: filen
    real(kind=8),intent(in) :: x(:),y(:),z(:)
    integer :: i,ok
    open(70,file=trim(filen),iostat=ok)
    do i=1,size(x,1)
       write(70,'(3E13.4)') x(i),y(i),z(i)
    enddo
    close(70)
end subroutine

subroutine save3r2f(filen,x,y,z)
    implicit none
    character(len=*),intent(in) :: filen
    real(kind=16),intent(in) :: x(:),y(:),z(:)
    integer :: i,ok
    open(70,file=trim(filen),iostat=ok)
    do i=1,size(x,1)
       write(70,'(3E13.4)') x(i),y(i),z(i)
    enddo
    close(70)
end subroutine

subroutine save4r2f(filen,x,y,z,g)
    implicit none
    character(len=*),intent(in) :: filen
    real(kind=16),intent(in) :: x(:),y(:),z(:),g(:)
    integer :: i,ok
    open(70,file=trim(filen),iostat=ok)
    do i=1,size(x,1)
       write(70,'(4E13.4)') x(i),y(i),z(i),g(i)
    enddo
    close(70)
end subroutine

subroutine save5r2f(filen,x,y,z,g,h)
    implicit none
    character(len=*),intent(in) :: filen
    real(kind=16),intent(in) :: x(:),y(:),z(:),g(:),h(:)
    integer :: i,ok
    open(70,file=trim(filen),iostat=ok)
    do i=1,size(x,1)
       write(70,'(5E13.4)') x(i),y(i),z(i),g(i),h(i)
    enddo
    close(70)
end subroutine

end module 
