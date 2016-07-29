!*************************************************************
!  lapack_wrapper
!
!>    It contains a bunch of routines to execute general
!!       matrix operations.
!!
!!    These routines are based on LAPACK modules.
!!
!*************************************************************
module lapack_wrapper
   implicit none
   integer,save,private          :: zge_lwork = 33
   integer,save,private          :: zheev_lwork = 33
   integer,parameter,private     :: maxquad = 9

contains

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

      return
   end subroutine


! *******************************************************
!
!  zhinvert:
!
!>     inverts an hermitian complex(kind=8) matrix
!
! *******************************************************
   subroutine zhinvert(a,uplo)
   implicit none
   complex(KIND=8),intent(inout)       :: a(:,:) !< Input: matrix, Output: Inverse matrix
   character(len=1),optional,intent(in):: uplo   !< Optional character that tells lapack where the data is stored.
   integer                             :: lwork,info
   complex(KIND=8)                     :: dummy(1)
   complex(KIND=8),allocatable         :: work(:)
   integer                             :: ipiv(size(a,1)),n,i,j,er
   character(len=1)                    :: uplotmp

      n = size(a,1)
      if (present(uplo)) then
         uplotmp = uplo
      else
         uplotmp = "L"
      endif   

      lwork = -1
      ! Get lwork (See zhetrf documentation)
      call zhetrf(uplotmp,n,a,n,ipiv,dummy,lwork,info)
      if (info .ne. 0) stop
      lwork = nint(real(dummy(1)))
      allocate(work(lwork),stat=er)
      if (er .ne. 0) stop
      ! Make a LU decomposition of a
      call zhetrf(uplotmp,n,a,n,ipiv,work,lwork,info)
     if (info .ne. 0) stop

      ! Invert a
      call zhetri(uplotmp,n,a,n,ipiv,work,info)
      if (info /=0 ) stop
      deallocate(work)

      ! Copy to the upper/lower triangle
      if (uplotmp .eq. "L" .or. uplotmp .eq. "l") then
         do j= 1,n-1
            do i=j+1,n
               a(j,i) = dconjg(a(i,j))
            enddo
         enddo
      else
         do j= 1,n-1
            do i=j+1,n
               a(i,j) = dconjg(a(j,i))
            enddo
         enddo
      endif

      return
   end subroutine

! *******************************************************
!
!  zginvert:
!
!>    inverts a general complex(kind=8) matrix
!!
! *******************************************************
   subroutine zginvert(a)
   implicit none
   complex(KIND=8),intent(inout) :: a(:,:)   !< Input: matrix, Output: Inverse matrix
   integer                       :: lwork    !< error flag
   integer                       :: info     !< error flag
   complex(KIND=8)               :: work(size(a,1)*zge_lwork)  !< Temporary working array. zge_lwork is a module variable.
   integer                       :: ipiv(size(a,1))            !< Permutation indices
   integer                       :: n                          !< Dimension

      n = size(a,1)
      lwork = zge_lwork * n

      ! LU decomposition
      call zgetrf(n,n,a,n,ipiv,info)
      if (info .ne. 0) stop

      ! Inversion
      call zgetri(n,a,n,ipiv,work,lwork,info)
      if (info /=0 ) stop 
      zge_lwork = ceiling(dble(work(1))/n)
      if (zge_lwork .lt. 1 ) zge_lwork=1

      return
   end subroutine zginvert


!*******************************************************
!
!  zgetrf:
!  zgetrs:
!>      Solves system Ax=b by LU decompositon and
!       back substitution(for square matrix), 
!       overwrites b by x
!
!*******************************************************
   subroutine zlusolve(a,b,n)
   implicit none
   complex(KIND=8), intent(inout)  :: a(n,n)
!trix
   integer                       :: info,info2     !< error flag
   integer                       :: ipiv(size(a,1))            !< Permutation
!indices
   integer                       :: n,i                       !< Dimension
   complex(KIND=8)               :: b(n),b1(n,1)              !< Right handside
!   complex(KIND=8)               :: x(n)
        integer                       :: nrhs                 !< nr righthands
        character*1                   :: trans='N'            !< form:no transp

  
      ! LU decomposition
      call zgetrf(n,n,a,n,ipiv,info)
      if (info /= 0) then
        write(*,*) info, "Warning: matrix inversion cannot be done"
        stop
      end if

      nrhs = 1
      
      do i =1,n
        b1(i,1)=b(i)
      enddo

      ! Solving
      call zgetrs(trans,n,nrhs,a,n,ipiv,b1,n,info2)
      if (info2 /= 0) then
        write(*,*) info, "Warning: matrix inversion cannot be done"
        stop
      end if
      
      do i =1,n
        b(i) = b1(i,1)
      enddo
        
      return
   end subroutine zlusolve

!*************************************************************************
!
!       z16solve       
!
!       wrapper to zcgesv
!
!>      (Input)        a       matrix of system ax = b to be solved
!>      (Input)        b       righthandside 
!       (Output)       x       solution to system ax=b
!>      (Input)        n       order of square matrix
!
!*************************************************************************
subroutine z16solve(a,b,n)
IMPLICIT NONE

INTEGER                 ::      n, nrhs, lda, ldb,i, ipiv(n), ldx, iter, info
COMPLEX(KIND=8)         ::      a(n,n), b(n),xp(n,1),bp(n)
COMPLEX(KIND=8)         ::      work(n,n),rwork(n)
COMPLEX                 ::      swork(n*(n+1))

nrhs = 1
lda = n
ldb = n
ldx = n

do i=1,n
  bp(i)=b(i)
enddo

call zcgesv(n,nrhs,a,lda,ipiv,b,ldb,xp,ldx,work,swork,rwork,iter,info)

if (info .ne. 0) then
   write(*,*) info, "Error: matrix system cannot be solved"
   stop
end if

do i=1,n
  b(i)=xp(i,1)
enddo

end subroutine

end module lapack_wrapper
