module padem

use pickpoints
use leastsquare
use poles
implicit none
include 'gitversion'

contains

subroutine pade(zin,fin,zout,fout)
    implicit none
    ! input arguments to the subroutine:
    complex(kind=16),intent(in) :: zin(:),fin(:)        ! Input points and values 
    complex(kind=16),intent(in) :: zout(:)              ! Output points 
    complex(kind=8),intent(out),allocatable :: fout(:)  ! Output values (at points zout)
    ! variables read from the parameter file pade.par:
    integer :: nminstart,nminfinish,nminstep
    integer :: Mmin,Mmax,Mstep
    integer :: Nmin,Nmax,Nstep
    logical :: diagonalPade
    integer :: solver 
    logical :: calcpoles   
    integer :: nprecision
    integer :: mtr
    real(kind=8) :: c1v,c2v
    ! help variables
    integer,allocatable :: nmins(:),Ms(:),Ns(:)         ! the arrays to loop over   
    complex(kind=16),allocatable :: zin_p(:),fin_p(:)   ! Picked points used for the analytical continuation
    type conf
       integer :: nmin                             ! index of first matsubara point.
       integer :: M                                ! number of input points used in analytical continuation.
       integer :: N                                ! number of Pade coefficients.
       integer :: w                                ! weight in average (either 0 or 1)
       complex(kind=16),allocatable :: x(:)        ! Pade approximant coefficients.
       complex(kind=8),allocatable :: P(:)         ! Pade approximant evaluated at the points zout.
    end type
    complex(kind=8),allocatable :: P(:,:)          ! For all continuations: Pade approximant evaluated at the points zout.
    type(conf),allocatable :: c(:)                 ! configuration data for each continuation is stored here.
    complex(kind=16),allocatable :: pole(:),res(:) ! poles and residues for one continuation.
    ! dummy variables
    integer :: i,j,k,b,counter

    ! reading the parameter file
    open(80,file='pade.par')
    read(80,*) ! read header
    read(80,*) ! read line about real axis mesh
    read(80,*) ! read line about distance above real axis
    read(80,*) nminstart,nminfinish,nminstep
    read(80,*) Mmin,Mmax,Mstep
    read(80,*) Nmin,Nmax,Nstep
    read(80,*) diagonalPade
    read(80,*) solver
    read(80,*) calcpoles
    read(80,*) nprecision
    read(80,*) mtr
    read(80,*) c1v,c2v
    close(80)
    
    ! write git version number to the information file 
    open(90,file='pade_info')
    write(90,'(a,a)') " git hash number: ", trim(git_version)

    ! construct nmins, Ms, Ns
    i = floor((nminfinish-nminstart)/(1d0*nminstep))
    allocate(nmins(i+1))
    nmins = (/ (nminstart+j*nminstep,j=0,i) /)
    i = floor((Mmax-Mmin)/(1d0*Mstep))
    allocate(Ms(i+1))
    Ms = (/ (Mmin+j*Mstep,j=0,i) /)
    i = floor((Nmax-Nmin)/(1d0*Nstep))
    allocate(Ns(i+1))
    Ns = (/ (Nmin+j*Nstep,j=0,i) /)

    write(90,*) 'nmins = ',nmins
    write(90,*) 'Ms = ',Ms
    write(90,*) 'Ns = ',Ns
    write(90,*)
    
    allocate(c(size(nmins)*size(Ms)*size(Ns))) ! allocate more memory than necessary, since size initially is unknown.
    counter = 0
    do i=1,size(nmins)
        do j=1,size(Ms)
            call pick_points(zin,fin,nmins(i),Ms(j),zin_p,fin_p)
            if (.not. diagonalPade) then
                do k=1,size(Ns)
                    if (Ns(k) <= Ms(j)) then
                        counter = counter + 1
                        write(90,'(a,3I4)') "(nmin,M,N)=",nmins(i),Ms(j),Ns(k)
                        c(counter)%nmin = nmins(i)
                        c(counter)%M = Ms(j)
                        c(counter)%N = Ns(k) 
                        allocate(c(counter)%x(Ns(k)),c(counter)%P(size(zout)))
                        call acPade(zin_p,fin_p,Ns(k),zout,solver,nprecision,mtr,c(counter)%x,c(counter)%P)
                        if (calcpoles) then
                            write(90,'(a)') "Calculate poles and residues"
                            write(90,'(a)') "Re[pole], Im[pole], Re[res], Im[res],nmin,M,N"
                            allocate(pole(Ns(k)/2),res(Ns(k)/2))
                            call getpoles(c(counter)%x,pole,res)
                            do b=1,size(pole)
                                write(90,'(4E15.5,3I4)') real(pole(b)),aimag(pole(b)),real(res(b)),aimag(res(b)), &
                                                          nmins(i),Ms(j),Ns(k)
                            enddo
                            deallocate(pole,res)
                        endif
                    endif
                enddo
            else
                counter = counter + 1
                write(90,'(a,3I4)') "(nmin,M,N)=",nmins(i),Ms(j),Ms(j) 
                c(counter)%nmin = nmins(i)
                c(counter)%M = Ms(j)
                c(counter)%N = Ms(j) 
                allocate(c(counter)%x(Ms(j)),c(counter)%P(size(zout)))
                call acPade(zin_p,fin_p,Ms(j),zout,solver,nprecision,mtr,c(counter)%x,c(counter)%P)
                if (calcpoles) then
                    write(90,'(a)') "Calculate poles and residues"
                    write(90,'(a)') "Re[pole], Im[pole], Re[res], Im[res],nmin,M,N"
                    allocate(pole(Ms(j)/2),res(Ms(j)/2))
                    call getpoles(c(counter)%x,pole,res)
                    do b=1,size(pole)
                        write(90,'(4E15.5,3I4)') real(pole(b)),aimag(pole(b)),real(res(b)),aimag(res(b)),nmins(i),Ms(j),Ms(j)
                    enddo
                    deallocate(pole,res)
                endif
            endif
        enddo
    enddo
    write(90,'(a,I4)') "counter = ", counter

    allocate(P(size(zout),counter))
    do i=1,counter
        P(:,i) = c(i)%P
    enddo
    call getAveragePade(P,c1v,c2v,c(1:counter)%w,fout)     ! create average 
    write(90,'(a,I4)') 'number of continuations with weight:',sum(c(1:counter)%w)
    
    ! write all spectra to file 
    open(71,file='pade_fout_all') 
    write(71,'(a)') "# nmin, M, N, w, Re[zout], Re[f(zout)], Im[f(zout)]"
    do i=1,counter
        do j=1,size(zout)
            write(71,'(4I4,3E13.3)') c(i)%nmin,c(i)%M,c(i)%N,c(i)%w,real(zout(j)),real(c(i)%P(j)),aimag(c(i)%P(j))
        enddo
        write(71,*)
    enddo
    close(71)

    close(90) ! close the file: pade_info
    do i=1,counter
        deallocate(c(i)%x,c(i)%P)
    enddo
    deallocate(c)
    deallocate(P)
    deallocate(nmins)
    deallocate(Ms)
    deallocate(Ns)

end subroutine

subroutine acPade(z,f,N,zout,solver,nprecision,mtr,x,fout)
    implicit none
    complex(kind=16),intent(in) :: z(:),f(:)
    integer :: N                                 ! number of Pade approximant coefficients 
    complex(kind=16),intent(in) :: zout(:)
    integer :: solver
    integer :: nprecision
    integer :: mtr
    complex(kind=16),intent(inout) :: x(:)
    complex(kind=8),intent(out) :: fout(:)
    ! help variables
    complex(kind=16),allocatable  :: fq(:)

    if ( .not. (mod(N,2) == 0 )) then
        stop "The number of Pade coefficients in the acPade routine should be even." 
    endif    
    ! Beach's matrix formulation
    if (solver == 0) then
        ! construct matrix and rhs in linear system of equations and seek a LS solution. 
        call padeMatrix(z,f,N,nprecision,mtr,x)
        allocate(fq(size(zout)))
        call epade(zout,x,fq)
        ! convert from quadruple to double precision
        fout = fq
        deallocate(fq)
    ! Thiele's algorithm 
    elseif ( solver == 1) then
        stop "Thiele's algorithm is not implemented yet."
    else
        stop "Solver value not valid."
    endif
end subroutine

subroutine padeMatrix(z,f,N,nprecision,mtr,x) 
    ! Returns the obtained Padé coefficients.
    implicit none
    complex(kind=16),intent(in) :: z(:),f(:)
    integer :: N
    integer :: nprecision
    integer :: mtr
    complex(kind=16),intent(out) :: x(:) ! Pade coefficients
    ! help variables
    integer :: M ! number of points
    integer :: r,i,j  
    complex(kind=16),allocatable :: A(:,:),y(:)
    complex(kind=16),allocatable :: errv(:)
    complex(kind=16),allocatable :: fp(:)
    ! double precision variables
    complex(kind=8),allocatable :: zd(:),fd(:)
    complex(kind=8),allocatable :: Ad(:,:),yd(:)
    
    M = size(z)
    r = N/2

    allocate(A(M,N),y(M))

    if (nprecision == 64 ) then ! double precision
        allocate(zd(M),fd(M))
        allocate(Ad(M,N),yd(M))
        zd = z
        fd = f
        yd = fd*zd**r 
        do i=1,M
            Ad(i,:r) = zd(i)**( (/ (j,j=0,r-1)   /) )
            Ad(i,r+1:) = -fd(i)*zd(i)**( (/ (j,j=0,r-1)   /) )
        enddo
        deallocate(zd)
        deallocate(fd)
        y = yd
        A = Ad
        deallocate(yd) 
        deallocate(Ad) 
    elseif (nprecision == 128 ) then ! quadruple precision
        y = f*z**r 
        do i=1,M
            A(i,:r) = z(i)**( (/ (j,j=0,r-1)   /) )
            A(i,r+1:) = -f(i)*z(i)**( (/ (j,j=0,r-1)   /) )
        enddo
    endif 
    call zls(A,y,mtr,nprecision,x)
    ! print two error checks: 
    allocate(errv(M))
    errv=matmul(A,x)-y
    write(90,'(a,E10.4)') '|k*x-y| = ',sqrt(sum(real(errv)**2+aimag(errv)**2))
    allocate(fp(M))
    call epade(z,x,fp)
    write(90,'(a,E10.4)') '|f(zin)-Pade(zin)| = ',sqrt(sum(abs(f-fp)**2))
    deallocate(errv)
    deallocate(fp)
    deallocate(A) 
    deallocate(y) 
end subroutine

subroutine epade(z,x,f)
    ! Returns the value of the Padé approximant at the points z.
    implicit none
    complex(kind=16),intent(in) :: z(:)  ! points in the complex plane
    complex(kind=16),intent(in) :: x(:)  ! Pade coefficients
    complex(kind=16),intent(out) :: f(:) ! values of the Pade approximant at the points z.
    integer :: i,r
    complex(kind=16),allocatable :: num(:),denum(:)
    r = size(x)/2
    allocate(num(size(z)),denum(size(z)))
    num = 0q0
    denum = 0q0
    do i=1,r
        num = num + x(i)*z**(i-1)
        denum = denum + x(r+i)*z**(i-1)
    enddo  
    denum = denum + z**r
    f = num/denum
    deallocate(num,denum)
end subroutine

subroutine getAveragePade(P,c1v,c2v,w,fout)
    implicit none
    complex(kind=8),intent(in)                :: P(:,:)  ! Pade approximants evaluated at points zout
    real(kind=8),intent(in)                   :: c1v,c2v ! average criteria parameters
    integer,intent(out)                       :: w(:)    ! averaging weights
    complex(kind=8),allocatable,intent(inout) :: fout(:) ! Average of Pade approximants P weighted with weightt w.
    ! help variables
    integer,allocatable      :: g(:)       ! index list with Pade approximants in P fullfilling Im(P)<=0      
    real(kind=8),allocatable :: A(:,:)     ! Imaginary part of Pade approximants, fullfilling Im(P)<=0 
    real(kind=8),allocatable :: delta(:,:) ! store distances between the continuations in A
    real(kind=8),allocatable :: d(:)       ! store accumulative distances to the other continuations in A
    integer,allocatable      :: ind(:)     
    logical,allocatable      :: c1(:),c2(:),mask(:)  
    ! dummy variables
    integer :: i,j,k,n,c

    n = size(P,1) ! nbr of points zout
    c = size(P,2) ! nbr of continuations

    ! count the number of "physical" continuations, meaning those who fulfill Im(P) <= 0
    k=0
    do i=1,c
        if (all(aimag(P(:,i))<=0)) then
            k = k+1    
        endif
    enddo 
    if(k==0) stop "No 'physical' continuations found"
    write(90,'(a,I4)') 'number of continuations with Im(f(zout))<=0 is:',k  
    allocate(g(k),A(size(P,1),k))
    k=0
    do i=1,c
        if (all(aimag(P(:,i))<=0)) then
            k = k+1
            g(k) = i
            A(:,k) = aimag(P(:,i))    
        endif
    enddo
    w = 0
    if(allocated(fout)) deallocate(fout)
    allocate(fout(size(P,1)))
    fout = 0
    if(k==1) then
        w(g(1)) = 1
        fout = P(:,g(1)) 
    else
        allocate(delta(k,k),d(k))
        delta = 0
        d = 0
        do i=1,k
            do j=i+1,k
                delta(i,j) = sum(abs(A(:,i)-A(:,j)))
            enddo 
        enddo
        do i=1,k
            d(i) = sum(delta(:i-1,i))+sum(delta(i,i+1:)) 
        enddo 
        ! calculate if to include or reject continuations based on two criteria:
        allocate(c1(k),c2(k),mask(k))
        c1 = ( d <= c1v*sum(d)/k ) 
        ind = sortp(d)
        c2 = .false.
        c2(ind(1:int(c2v*k))) = .true. ! the c2v % lowest continuations fulfill criterion 2.
        mask = c1 .and. c2
        j = 0
        do i=1,k
            if( mask(i) ) then
                j = j+1
                w(g(i)) = 1
                fout = fout + P(:,g(i))
            endif
        enddo
        fout = fout/j ! divide with the number of continuations contributing in the average
        deallocate(delta,d,c1,c2,mask)
    endif
    deallocate(g,A)
end subroutine

function sort(A) !bubble sort algorithm, smallest first.
implicit none
real(kind=8),allocatable :: sort(:)
real(kind=8),intent(in) :: A(:)

real(kind=8) :: b
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
! Eg: A=[ 3.1  1.3  2.5  1.3 ]    ->    A_sorted=[1.3 1.3 2.5 3.1] and sortp=[2 4 3 1]  so A=A_old(c)
implicit none
integer,allocatable :: sortp(:)
real(kind=8),intent(in) :: A(:)

integer,allocatable :: d(:),e(:)
real(kind=8),allocatable :: x(:),y(:)
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
