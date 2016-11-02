module poles

implicit none

contains

subroutine getpoles(p,q,res)
    IMPLICIT NONE
    complex(kind=16),intent(inout) :: p(:) !cofficients
    complex(kind=16),intent(inout) :: q(:),res(:) !pole position and residue
    ! help variables
    complex(kind=16),allocatable :: z(:) !zeros in pade
    complex(kind=16) :: c !pre factor
    integer :: M !nbr of coefficients

    if( sum(p) == cmplx(0q0,0q0,kind=16)) then
        write(*,*) "All coefficients are 0" 
        q=0
        res=1
        return
    endif
    M=size(p)
    if(mod(M,2)==0) then
        allocate(z(M/2-1))  !M/2-1 roots in nominator
    else
        allocate(z((M-1)/2)) !(M-1)/2 roots in nominator
    endif
    call p_factorize(p,q,z,c) !calculate pole positons q and zero positions z
    call p_residues(z,q,c,res) !calculate residues res
    deallocate(z)
end subroutine

SUBROUTINE p_factorize(pq,qq,zq,cq)
!**********************************************************
!
!       p_factorize
!
!>      Gives the poles and zeros of Pade approximant
!             f(z)=C*(Product_i (z-z(i)))/(Product_i (z-q(i))) 
!      
!>   Input:   p(nom)  coefficients of pade approximant.
!                     Pade approximant has shape either
!
!>                   f(z)=(p(1)+p(2)z+...+p(r)z**(r-1))/
!                        (p(r+1)+...+p(2r)z**(r-1)+z**r)
!                    if nom is even or
!>                   f(z)=(p(1)+p(2)z+...+p(r)z**(r-1)+p(r+1)z**r)/
!                        (p(r+2)+...+p(2r+1)z**(r-1)+z**r)
!                    if nom is odd. 
! 
!>            nom  number of Matsubara points
!>            r    degree of denominator
!
!>    Output: z(1:r-1) or z(1:r) , zeros of approximant. 1:r-1 if nom is even, 1:r if nom is odd
!>            q(1:r)      poles of approximant
!>            C         constant of factorisation
!
!**********************************************************
    IMPLICIT NONE
    complex(kind=16),intent(in)   :: pq(:)
    complex(kind=16),intent(out)  :: qq(:),zq(:)
    COMPLEX(KIND=16),intent(out)  :: cq

    complex(KIND=8),allocatable   :: A(:,:),B(:,:)
    complex(kind=8),allocatable   :: p(:),q(:),z(:)
    !INTEGER,PARAMETER             :: re=1d0
    INTEGER                       :: nom,r
    !dummy
    INTEGER                      ::      i,j
    nom=size(pq) !nbr of Pade coeff
    if(mod(nom,2)==0) then
       r=nom/2
       cq = pq(r)
       allocate(z(r-1))
       allocate(A(r-1,r-1),B(r,r))
    else
       r=(nom-1)/2
       cq = pq(r+1)
       allocate(z(r))
       allocate(A(r,r),B(r,r))
    endif
    allocate(p(nom),q(r))
    p=pq !convert pade coeff from quad to double
    A = (0d0,0d0)
    B = (0d0,0d0)

    ! Computes Companion Matrix of nominator and denominator
    ! polynomial
    if (r .eq. 1) then
       write(*,11) p(2)
    11   format("No zeros, pole at z=",F8.4,1x,E12.5)
    else 
       !construct companion matrix B (denominator) 
       do j=1,r-1
          B(j+1,j)= 1d0
       enddo
       if(mod(nom,2)==0) then
          do i=1,r
             B(i,r)=-p(r+i)
          enddo
       else
          do i=1,r
             B(i,r)=-p(r+1+i)
          enddo
       endif
       !construct companion matrix A (nominator) 
       if(mod(nom,2)==0) then
          do j=1,r-2
                A(j+1,j)= 1d0
          enddo
          if (p(r) .ne. 0) then
            do i=1,r-1
                A(i,r-1)=-p(i)/p(r)
            enddo
          else
            write(1,*) 'stop'         
            stop
          end if
       else
          do j=1,r-1
             A(j+1,j)= 1d0
          enddo
          if (p(r+1) .ne. 0) then
            do i=1,r
                A(i,r)=-p(i)/p(r+1)
            enddo
          else
                write(1,*) 'Sigma0 is calculated to be zero, so probably better to run with even Pade coefficients'         
                stop
          end if
       endif
       ! Finds eigenvalues of Companion Matrices
       call zgeevwrapper(A,z)  !r-1
       call zgeevwrapper(B,q)  !r 

    end if

    !convert back to quad variables
    qq=q
    zq=z
    deallocate(q,z,p)
    deallocate(A,B)
END SUBROUTINE

SUBROUTINE p_residues(z,q,C,w)
    !*******************************************************************
    !
    !       p_residues
    !
    !>      Can be used after Pade approximation and factorisation of
    !>      Pade polynomials to calculate the residues of the approximant
    !
    !>      Input:  z(1:r-1) or z(1:r)  zeros of approximant. size depending on if even or odd number of Pade coefficients
    !>              q(r)    poles of approximant
    !>              C       constant factor
    !>              r       degree of denominator (half of nom)
    !
    !       Output: w(r)    residues of poles q(r)
    !
    !*******************************************************************
    IMPLICIT NONE
    COMPLEX(KIND=16),intent(in)       ::  z(:),q(:),C
    COMPLEX(KIND=16)                  ::  w(:)
    INTEGER                           ::  i,j,r !dummy
    r=size(q)
    w(:)=C
    !Calculate residues
    do j=1,r
        do i=1,j-1
            w(j)=w(j)*(q(j)-z(i))/(q(j)-q(i))
        enddo
        do i=j+1,r
            w(j)=w(j)*(q(j)-z(i-1))/(q(j)-q(i))
        enddo
        if(size(z,1)==r) then
            w(j)=w(j)*(q(j)-z(r))
        endif
    enddo
end subroutine

!
!>    diagonalizes a general complex(kind=8) matrix and 
!!     returns the eigenvalues and eigenvectors
!!
! *******************************************************
subroutine zgeevwrapper(a,e)
    ! A wrapper around zgeev
    implicit none
    complex(KIND=8),intent(inout) :: a(:,:) !< input: Matrix, output: Eigenvectors
    complex(KIND=8),intent(out)   :: e(:)   !< Eigenvalues
    ! See lapack zgeev on netlib for more information 
    integer                       :: lwork,info,er
    complex(KIND=8),allocatable   :: work(:)
    real(KIND=8),allocatable      :: rwork(:)
    real(KIND=8)                  :: rwtmp(1) !< dummy
    complex(KIND=8)               :: v1tmp(1),v2tmp(2),wtmp(1) !< dummy,dummy, used to get proper work
    integer :: n      !< size(a,1)
    n=size(a,1)
    lwork = -1
    call zgeev('N','N',n,a,n,e,v1tmp,1,v2tmp,1,wtmp,lwork,rwtmp,info)
    if (info .ne. 0) stop
    lwork = max(6*n,int(wtmp(1)))
    allocate(work(lwork),rwork(2*n),stat=er)
    if (er .ne. 0) stop
    call zgeev('N','N',n,a,n,e,v1tmp,1,v2tmp,1,work,lwork,rwork,info)
    if (info .ne. 0) stop
    deallocate(work,rwork)
end subroutine

end module
