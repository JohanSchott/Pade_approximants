module leastsquare
implicit none

contains

subroutine zls(Ao,bo,mtr,pr,x)
    !Solves equation A*x=b in Least Square sense, using complex variables
    !mtr=0 sends A and b to lapack. This system is solved in Least Square sense. 
    !mtr=1 sends A^dagger*A and A^dagger*b to lapack. This system is solved exact.
    !mtr=2 sends A to lapack for SVD. Construct explicitly Least Square solution.
    !mtr=3 sends A to lapack for SVD. Construct explicitly Least Square solution, but truncate modes with small singular values.
    implicit none
    complex(kind=16),intent(in) :: Ao(:,:),bo(:)
    integer,intent(in) :: mtr
    integer,intent(in) :: pr
    complex(kind=16)   :: x(:)
    !internal
    integer :: i,N,M,p,r,mo
    complex(kind=16),allocatable :: A(:,:),b(:),k(:,:),g(:) 
    complex(kind=16),allocatable :: u(:,:),v(:,:)
    real(kind=16),allocatable :: s(:)
    real(kind=16) :: alpha

    N=size(Ao,1)
    M=size(Ao,2)
    p=min(N,M)
    allocate(A(N,M),b(N))
    b=bo !copy
    A=Ao !copy
    ! Let Lapack solve equation: A*x=b in LS sense
    if(mtr==0) then
       ! solves A*x=b using QR or LQ factorization. Full rank is assumed.
       call zgelswrapper(A,b,pr)
       ! solves A*x=b using Singular Value Decomposition SVD (only double precision implemented).  
       !call zgelsswrapper(A,b)  
       ! solves A*x=b using Singular Value Decomposition SVD with Divide and Conquer (only double precision implemented) 
       !call zgelsdwrapper(A,b)
    ! Let Lapack exactly solve normal equation: A^dagger*A*x=A^dagger*b
    elseif(mtr==1) then
       allocate(k(M,M),g(M))
       k=matmul(conjg(transpose(A)),A)
       g=matmul(conjg(transpose(A)),b)
       ! solves k*x=g using QR or LQ factorization. Full rank is assumed.
       call zgelswrapper(k,g,pr)
       ! solves A*x=b using Singular Value Decomposition SVD (only double precision implemented).  
       !call zgelsswrapper(k,g) 
       ! solves A*x=b using Singular Value Decomposition SVD with Divide and Conquer (only double precision implemented) 
       !call zgelsdwrapper(k,g) 
       b(1:M)=g
       deallocate(k,g)
    ! Explicit SVD
    elseif(mtr==2 .or. mtr==3) then
       allocate(u(N,N),s(p),v(M,M))
       ! computes the singular value decomposition of a complex matrix using Divide and Conquer. (only double precision implemented)
       call zgesddwrapper(A,u,s,v,pr)
       write(90,'(a)') "SVD's singular values: s_i and s_i/s_1"
       do i=1,p
          write(90,'(2E13.5)') s(i),s(i)/s(1)
       enddo
       write(*,*)
       if(pr==64) then
          r=minloc(abs(s-s(1)*0.5q0*10**(-15q0)),1)
       elseif(pr==128) then
          r=minloc(abs(s-s(1)*0.5q0*10**(-30q0)),1)
       else
          stop "precision not supported"
       endif
       write(90,'(I4,a)') r," is estimated SVD rank."
       write(90,'(E15.5,a)') s(1)/s(r)," is estimated condition number"
       if(mtr==2) then
          mo=r
       ! Explicit SVD and use singular values s(i)/s(1)>alpha
       elseif(mtr==3) then
          !This is just an example value.
          alpha=10**(-4q0)
          mo=minloc(abs(s-s(1)*alpha),1) 
       endif
       write(90,'(I4,a,I4)') mo," SVD modes used for solution, out of",min(M,N)
       ! construct solution here
       allocate(g(M))
       g=cmplx(0q0,0q0,kind=16)
       do i=1,mo
          g=g+1q0/s(i)*dot_product(u(:,i),b)*v(:,i)
       enddo
       b(1:M)=g
       deallocate(g)
       deallocate(u,s,v)
    else
       stop "wrong mtr value"
    endif
    x=b(1:M)
    deallocate(A,b)
end subroutine

subroutine zgesddwrapper(A,u,s,v,pr)
    ! Preforms SVD factorization of matrix A in precision pr
    implicit none
    complex(kind=16),intent(in) :: A(:,:)
    complex(kind=16),intent(out) :: u(:,:),v(:,:)
    real(kind=16),intent(out) :: s(:)
    integer,intent(in) :: pr
    !internal
    integer :: N,M,p
    integer :: lwork,lrwork,info
    integer,allocatable :: iwork(:)
    !double precison
    complex(kind=8),allocatable :: Ad(:,:)
    complex(kind=8),allocatable :: ud(:,:),vd(:,:)
    real(kind=8),allocatable :: sd(:)
    complex(kind=8),allocatable :: workd(:)
    real(kind=8),allocatable :: rworkd(:)

    N=size(A,1)
    M=size(A,2)
    p=min(N,M)

    if(pr==64) then  !double precision
       allocate(Ad(N,M),ud(N,N),sd(p),vd(M,M))
       Ad=A
       !call lapack routine
       lwork=2*min(N,M)*min(N,M)+2*min(N,M)+max(N,M)
       lrwork=min(N,M)*max(5*min(N,M)+7,2*max(N,M)+2*min(N,M)+1)
       allocate(workd(lwork),rworkd(lrwork),iwork(8*p))
       lwork=-1
       call zgesdd('A',N,M,Ad,N,sd,ud,N,vd,M,workd,lwork,rworkd,iwork,info)
       if(info<0) then
          stop "SVD factorization had wrong input"
       elseif(info>0) then
          stop "SVD factorization did not converge"
       endif
       lwork=int(real(workd(1)))  !optimal lwork
       deallocate(workd)
       allocate(workd(lwork))
       call zgesdd('A',N,M,Ad,N,sd,ud,N,vd,M,workd,lwork,rworkd,iwork,info)
       deallocate(workd,rworkd,iwork)
       vd=transpose(conjg(vd)) !transform v so not equal to v^dagger
       u=ud
       s=sd
       v=vd
       deallocate(Ad,ud,sd,vd) 
    elseif(pr==128) then !quad precision
       stop "Precison not supported yet."
    else
       stop "Precision not supported"
    endif
end subroutine

subroutine zgelswrapper(A,b,pr)
    !solves overdetermined or underdetermined complex linear systems, assuming full rank of matrix A.
    !Accepts double (pr=64) and quad (pr=128)  precision.
    implicit none
    complex(kind=16),intent(in) :: A(:,:)
    complex(kind=16),intent(inout) :: b(:)
    integer,intent(in) :: pr
    integer :: N,M
    !lapack
    integer :: lwork
    complex(kind=16),allocatable :: work(:)
    integer :: info
    !double precision variables
    complex(kind=8),allocatable :: workd(:)
    complex(kind=8),allocatable :: Ad(:,:),bd(:)
    N=size(A,1)
    M=size(A,2)
    if(N<M) stop "N should be bigger or equal to M"
    allocate(work(2*M))
    lwork=-1
    if(pr==64) then
       allocate(Ad(N,M),bd(N),workd(2*M))
       Ad=A
       bd=b
       call zgels('N',N,M,1,Ad,N,bd,N,workd,lwork,info)
       if(info .ne. 0) stop "least square routine do not work properly"
       lwork=ceiling(real(workd(1))) !get optimal lwork
       deallocate(workd)
       if(lwork < 2*M ) stop "lwork is simply too low. Increase lwork in lapace routine."
       allocate(workd(lwork))
       call zgels('N',N,M,1,Ad,N,bd,N,workd,lwork,info)
       b=bd ! convert solution to quad precision variable
       deallocate(Ad,bd,workd)
    elseif(pr==128) then
       call zgelsquad('N',N,M,1,A,N,b,N,work,lwork,info)
       if(info .ne. 0) stop "least square routine do not work properly"
       lwork=ceiling(real(work(1))) !get optimal lwork
       deallocate(work)
       if(lwork < 2*M ) stop "lwork is simply too low. Increase lwork in lapace routine."
       allocate(work(lwork))
       call zgelsquad('N',N,M,1,A,N,b,N,work,lwork,info)
    else
       stop "Lapack's zgels can only handle dbl or quad"
    endif
    if(info<0) then
        stop "an argument with illegal value."
    elseif(info>0) then
       write(90,*) "Matrix A does not have full rank, which is assumed. LS-solution could not be computed"
    endif
    deallocate(work)
    !b(1:M) is the solution
end subroutine

subroutine zgelsswrapper(A,b)
    !compute the minimum norm solution to complex linear least squares problem
    implicit none
    complex(kind=16),intent(inout) :: A(:,:),b(:)
    integer :: N,M
    complex(kind=8),allocatable :: Ad(:,:),bd(:)
    !lapack
    integer :: rank !effective rank of matrix A
    integer :: lwork,lrwork
    integer :: info
    complex(kind=8),allocatable :: work(:)
    real(kind=8),allocatable :: s(:)
    real(kind=8) :: rcond
    real(kind=8),allocatable :: rwork(:)
    N=size(A,1)
    M=size(A,2)
    allocate(Ad(N,M),bd(N))
    Ad=A
    bd=b
    if(N<M) stop "N should be bigger than or equal to M"
    lwork=-1 
    lrwork=max(5*M-4,1)
    allocate(s(M))
    allocate(rwork(lrwork)) 
    allocate(work(2*(2*N+M))) !twice times minimum, unnecessary big ??
    ! parameter to determine effective rank of A.
    ! singular values s_i <= rcond*s_1 are treated as zero.
    ! If rcond < 0, machine precision is used instead.
    rcond=-1d0 
    call ZGELSS(N,M,1,Ad,N,bd,N,s,rcond,rank,work,lwork,rwork,info)
    if(info .ne. 0) stop "least square routine do not work properly"
    lwork=ceiling(real(work(1))) !get optimal lwork
    deallocate(work)
    if(lwork < 2*N+M ) stop "lwork is simply too low. Increase lwork in lapack routine."
    allocate(work(lwork))
    call ZGELSS(N,M,1,Ad,N,bd,N,s,rcond,rank,work,lwork,rwork,info) !bd(1:M) contains solution
    if(info<0) then
        stop "an argument with illegal value."
    elseif(info>0) then
        write(90,*) "least square algorithm failed to converge."
    endif
    write(90,*) " "
    write(90,'(A,E12.4)') "Condition number: ",s(1)/s(M) 
    write(90,'(A,I3,A,I3,A)') "Effective rank of matrix:",rank," (",min(M,N)," for full rank)"
    deallocate(s,work,rwork)
    b=bd !convert from double to quad precision
    deallocate(Ad,bd)
end subroutine

subroutine zgelsdwrapper(A,b)
    !compute the minimum norm solution to complex linear least squares problem
    implicit none
    complex(kind=16),intent(inout) :: A(:,:),b(:)
    integer :: N,M
    complex(kind=8),allocatable :: Ad(:,:),bd(:)
    !lapack
    integer :: rank !effective rank of matrix A
    integer :: lwork,lrwork,liwork
    integer :: info
    complex(kind=8),allocatable :: work(:)
    real(kind=8),allocatable :: s(:)
    real(kind=8) :: rcond
    real(kind=8),allocatable :: rwork(:)
    integer,allocatable :: iwork(:)
    integer :: NLVL,SMLSIZ

    N=size(A,1)
    M=size(A,2)
    allocate(Ad(N,M),bd(N))
    Ad=A
    bd=b
    if(N<M) stop "N should be bigger than or equal to M"
    SMLSIZ=25 !normal value according to netlib.org
    NLVL = MAX( 0, INT( log(1d0*M/(SMLSIZ+1))/log(2d0) ) + 1 )
    lwork=-1 
    lrwork=max(5*M-4,1)
    liwork = max(1, 3*M*NLVL + 11*M)
    allocate(s(M))
    allocate(work(2*(2*N+M))) !twice times minimum, unnecessary big ??
    allocate(rwork(lrwork)) 
    allocate(iwork(liwork))
    rcond=-1  !parameter to determine effective rank of A
    call ZGELSD(N,M,1,Ad,N,bd,N,s,rcond,rank,work,lwork,rwork,iwork,info)
    if(info .ne. 0) stop "least square routine do not work properly"
    lwork=ceiling(real(work(1))) !get optimal lwork
    lrwork=ceiling(real(rwork(1)))
    liwork=iwork(1)
    deallocate(work,rwork,iwork)

    if(lwork < 2*M+M ) stop "lwork is simply too low. Increase lwork in lapack routine."
    if(lrwork <= 1 ) stop "lrwork is simply too low. Increase lrwork in lapack routine."
    if(liwork <= 1 ) stop "liwork is simply too low. Increase liwork in lapack routine."
    allocate(work(lwork),rwork(lrwork),iwork(liwork))
    call ZGELSD(N,M,1,Ad,N,bd,N,s,rcond,rank,work,lwork,rwork,iwork,info) !bd(1:M) contains solution
    if(info<0) then
        stop "an argument with illegal value."
    elseif(info>0) then
        write(90,*) "least square algorithm failed to converge."
    endif
    write(90,*) " "
    write(90,'(A,E12.4)') "Condition number: ",s(1)/s(M) 
    write(90,'(A,I3,A,I3,A)') "Effective rank of matrix:",rank," (",min(M,N)," for full rank)"
    deallocate(s,work,rwork,iwork)
    b=bd !convert from double to quad precision
    deallocate(Ad,bd)
end subroutine

end module 
