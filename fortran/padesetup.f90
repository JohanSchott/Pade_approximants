module padeSetup
use tikhonovm
use config
implicit none

contains

subroutine pade(wn,Gn,mtr,pr,rep,solver,cof,err)
!Tries to solve A*x=b using Lapack routines.
!A_{n,m},x_{m,1},b_{n,1}.
!r=rank(A),p=min(n,m),r<=p
!Problem can be formulated either as real value or complex value problem. 
!Residue res=Ax-b
!Returns solution in cof and ||res||_2
!Solution either:
!1)Least Square solution:
!  Implemented for double and quad precision complex numbers. Not for real numbers.
!  Means minimizing ||A*x-b||_2, with extra constrain of minimizing ||x||_2 if needed.
!  Equivalent to solving normal equation A^dagger*A*x=A^dagger*y. This puts null space components (free modes) to zero. Many solvers: QR,SVD,...
!  Using SVD, one can also calculate LS solution as x=sum_i=1^r 1/s_i*(u_i^dagger b) v_i
!  OBS: Modified SVD solution also included here. Truncate singular values smaller than a certain number: x=sum_i=1^j 1/s_i*(u_i^dagger b) v_i , with j<r
!2)Tikhonov SVD solution:
!  Not implemented.
!  Filter away singular value components, using parameter alpha.
!  x=sum_i=1^r s_i/(s_i^2+alpha^2) (u_i^dagger*b) v_i
!  x is actually just a LS solution to slightly modified problem:
!  |A      |     |b|
!  |       |*x = | |=Am*x=bm, namely solves normal equation: Am^dagger*Am*x=Am^dagger*bm. So instead of LS solution for A and b, find LS solution for Am and bm.
!  |alpha*I|     |0|
!  Optimal parameter alpha value can be obtained by minimizing ||res||_2+||x||_2 
!3)Non-Negative LS solution:
!  Not implemented but don't think make sense to enforce non-negativity on Pade coefficients
!4)Non-Negative Tikhonov solution:
!  Not implemented but don't think make sense to enforce non-negativity on Pade coefficients
!  Same algorithm as Non-Negative LS but for matrix Am instead of for A.
   implicit none
   complex(kind=16),intent(in) :: wn(:),Gn(:)
   integer,intent(in) :: mtr
   integer*8,intent(in) :: pr
   integer,intent(in) :: rep
   integer,intent(in) :: solver  ! Which solver algorithm to use
   complex(kind=16),intent(inout) :: cof(:)
   real(kind=16),intent(out) :: err(:)
   !internal
   integer :: pre
   pre=pr

   if(rep==1) then !complex
      call pade_comp(wn,Gn,mtr,pre,solver,cof,err)
   else !real 
      call pade_real(wn,Gn,mtr,pre,solver,cof,err,rep)
   endif
end subroutine

subroutine pade_real(wn,Gn,mtr,pr,solver,cof,err,rep)
!Using real variables  
   implicit none
   complex(kind=16),intent(in) :: wn(:),Gn(:)
   integer,intent(in) :: mtr
   integer,intent(in) :: pr
   integer,intent(in) :: solver  ! Which solver algorithm to use
   complex(kind=16),intent(inout) :: cof(:)
   real(kind=16),intent(out) :: err(:)
   integer,intent(in) :: rep

   ! internal
   integer :: N,M,r
   integer :: i,k 
   real(kind=16),allocatable :: Gr(:),Gi(:),w(:)
   real(kind=16),allocatable :: A(:,:),b(:),Ar(:,:),Ai(:,:),br(:),bi(:)
   real(kind=16),allocatable :: sol(:)
   !real(kind=16),allocatable :: Ao(:,:),bo(:)
   real(kind=16),allocatable :: errv(:)
   complex(kind=16),allocatable :: Gp(:)

   N=size(wn,1)
   M=size(cof,1)
   if(N .lt. M) stop "N>=M is required"
   if(mod(M,2)==1 ) then
      stop "Pade with constant asymptotic when |z| big is not implemented yet!" 
   else
      r=M/2
      !help variables
      allocate(Gr(N),Gi(N),w(N))
      w=aimag(wn)
      Gr=real(Gn)
      Gi=aimag(Gn)

      !build br and bi from w,Gr and Gi
      allocate(Ar(N,M),Ai(N,M),br(N),bi(N))
      if(mod(r,4)==0) then
         br=Gr*w**r
         bi=Gi*w**r
      elseif(mod(r,4)==1) then
         br=-Gi*w**r
         bi=Gr*w**r
      elseif(mod(r,4)==2) then
         br=-Gr*w**r
         bi=-Gi*w**r
      elseif(mod(r,4)==3) then
         br=Gi*w**r
         bi=-Gr*w**r
      else
         stop "What the crack?"
      endif

      !build Ar and Ai from w,Gr and Gi
      do k=0,r-1
         if(mod(k,4)==0) then
            Ar(:,k+1)=w**k
            Ar(:,k+1+r)=-Gr*w**k
            Ai(:,k+1)=0q0
            Ai(:,k+1+r)=-Gi*w**k
         elseif(mod(k,4)==1) then
            Ar(:,k+1)=0q0
            Ar(:,k+1+r)=Gi*w**k
            Ai(:,k+1)=w**k
            Ai(:,k+1+r)=-Gr*w**k
         elseif(mod(k,4)==2) then
            Ar(:,k+1)=-w**k
            Ar(:,k+1+r)=Gr*w**k
            Ai(:,k+1)=0q0
            Ai(:,k+1+r)=Gi*w**k
         elseif(mod(k,4)==3) then
            Ar(:,k+1)=0q0
            Ar(:,k+1+r)=-Gi*w**k
            Ai(:,k+1)=-w**k
            Ai(:,k+1+r)=Gr*w**k
         endif
      enddo
      deallocate(Gr,Gi,w)
     
      !build A from Ar and Ai
      if(rep==2) then
         allocate(A(2*N,2*M))
         A(1:N,1:M)=Ar
         A(1:N,M+1:2*M)=-Ai
         A(N+1:2*N,1:M)=Ai
         A(N+1:2*N,M+1:2*M)=Ar
      elseif(rep==3) then
         allocate(A(2*N,M))
         A(1:N,:)=Ar
         A(N+1:2*N,:)=Ai
      else
         stop "Only two real pade representations are implemented"
      endif
      deallocate(Ar,Ai)

      !build b from br and bi
      allocate(b(2*N))
      b(1:N)=br
      b(N+1:2*N)=bi
      deallocate(br,bi)

      !allocate(Ao(size(A,1),size(A,2)),bo(size(b,1))) !copies 
      !Ao=A
      !bo=b
     
      if(solver==1) then
         call ls(A,b,mtr,pr,sol)
      elseif(solver==2) then
         call tikhonov(A,b,mtr,pr,.false.,sol)
      else
         stop "This solver algorithm is not implemented"
      endif
      
      !copy solution and calculate error
      allocate(errv(2*N))
      if(rep==2) then
         !errv=matmul(Ao,b(1:2*M))-bo
         errv=matmul(A,sol)-b
         do i=1,M
            cof(i)=cmplx(sol(i),sol(M+i),kind=16)
         enddo
      elseif(rep==3) then
         !errv=matmul(Ao,b(1:M))-bo
         errv=matmul(A,sol)-b
         do i=1,M
            cof(i)=cmplx(sol(i),0q0,kind=16)
         enddo
      endif
      err(1)=sqrt(sum(errv**2))  ! |k*x-b|
      allocate(Gp(N))
      call eval_pade(cof,wn,Gp)
      err(2)=sqrt(sum(abs(Gn-Gp)**2)) ! |G-Gpade|

      deallocate(A,b)
      deallocate(sol,errv,Gp)
   endif
end subroutine

subroutine pade_comp(wn,Gn,mtr,pr,solver,cof,err)
!Using complex variables  
   implicit none
   complex(kind=16),intent(in) :: wn(:),Gn(:)
   integer,intent(in) :: mtr
   integer,intent(in) :: pr
   integer,intent(in) :: solver  ! Which solver algorithm to use
   complex(kind=16),intent(inout) :: cof(:)
   real(kind=16),intent(out) :: err(:)

   ! internal
   integer :: N,M,r
   integer :: i 
   complex(kind=16),allocatable :: A(:,:),b(:)
   complex(kind=16),allocatable :: sol(:)
   complex(kind=16),allocatable :: errv(:),Gp(:)

   N=size(wn,1)
   M=size(cof,1)
   if(N .lt. M) stop "N>=M is required"
   if(mod(M,2)==1 ) then
      r=(M-1)/2 
      ! construct matrix 
      allocate(A(N,M),b(N))
      do i=1,r
          A(:,i) = wn(:)**(i-1)
          A(:,r+1+i) = -Gn(:)*(wn(:)**(i-1))
      enddo
      A(:,r+1)=wn(:)**r
      b(:) = Gn(:)*(wn(:)**r) !construct vector
   else
      r=M/2
      ! construct matrix 
      allocate(A(N,M),b(N))
      do i=1,r
          A(:,i) = wn(:)**(i-1)
          A(:,r+i) = -Gn(:)*(wn(:)**(i-1))
      enddo
      b(:) = Gn(:)*(wn(:)**r) !construct vector
   endif
   if(solver==1) then
      call ls(A,b,mtr,pr,sol)
   elseif(solver==2) then
      call tikhonov(A,b,mtr,pr,sol)
   else
      stop "This solver algorithm is not implemented"
   endif
   allocate(errv(N))
   !errv=matmul(Ao,cof)-bo
   errv=matmul(A,sol)-b
   err(1)=sqrt(sum(real(errv)**2+aimag(errv)**2)) ! |k*x-b|
   allocate(Gp(N))
   call eval_pade(sol,wn,Gp)
   err(2)=sqrt(sum(abs(Gn-Gp)**2)) ! |G-Gpade|
   deallocate(A,b)
   deallocate(errv,Gp)
   !deallocate(Ao,bo)
   cof=sol
   deallocate(sol)
end subroutine

end module 
