module leastsquare
use matlab
implicit none

interface ls
   module procedure dls,zls
end interface ls

contains

subroutine nnls(A,b,mtr,pr,x)
!Non Negative LS routine. Solve with LS A*x=b with contrain x>0
!Algorithm taken from book "Solving Least Squares Problems" by Charles Lawson and Richard Hanson.
!Implemented by Johan Schott at Uppsala University 2015
!Iterativly solve reduced least squares problems.
implicit none
real(kind=16),intent(in) :: A(:,:),b(:)
integer,intent(in) :: pr,mtr
real(kind=16),allocatable,intent(inout) :: x(:)

integer :: m,n,t,q
integer,allocatable :: P(:),R(:)
integer,allocatable :: C(:),D(:),F(:)
real(kind=16),allocatable ::  E(:,:),w(:),z(:),zs(:)
real(kind=16) :: alpha

integer :: i,j,k,g,h,ia,ib

m=size(A,1)
n=size(A,2)
allocate(z(n),w(n))

!Step 1
allocate(R(n))
R=[ (i,i=1,n) ]
if(allocated(x)) deallocate(x)
allocate(x(n))
x=0q0
ia=0
do !Loop A
   !Step 2
   w=matmul(transpose(A),b-matmul(A,x))
   if(allocated(F)) then
      w(F)=0q0
   endif
   !Step 3
   if(.not. allocated(R)) exit
   j=0
   do i=1,size(R,1)
      if(w(R(i)) .le. 0) j=j+1
   enddo
   if(j .eq. size(R,1)) exit
   !Step 4
   t=maxloc(w(R),1)
   t=R(t)
   !Step 5
   if( .not. allocated(P)) then
      !Create P and add t to P
      allocate(P(1))
      P(1)=t
      !Remove t from R. 
      allocate(C(size(R,1)))
      C=R
      deallocate(R)
      allocate(R(size(C,1)-1))
      k=1
      do i=1,size(C,1)
         if(C(i) .ne. t  ) then
            R(k)=C(i)
            k=k+1
         endif
      enddo
      deallocate(C)
   else
      !Add t to P
      allocate(C(size(P,1)+1))
      C(1:size(P,1))=P
      C(size(P,1)+1)=t
      deallocate(P)
      allocate(P(size(C,1)))
      P=C
      deallocate(C)
      !Remove t from R. 
      if(size(R,1)==1) then
         deallocate(R)
      else
         allocate(C(size(R,1)))
         C=R
         deallocate(R)
         allocate(R(size(C,1)-1))
         k=1
         do i=1,size(C,1)
            if(C(i).ne. t  ) then
               R(k)=C(i)
               k=k+1
            endif
         enddo
         deallocate(C)
      endif
   endif
   ib=0
   do !Loop B
      !Step 6
      allocate(E(m,size(P,1)),zs(size(P,1)))
      do i=1,size(P,1)
         E(:,i)=A(:,P(i))
      enddo
      call ls(E,b,mtr,pr,zs)
      z=0q0
      do i=1,size(P,1)
         z(P(i))=zs(i)
      enddo
      deallocate(zs,E)

      if(ib==0 .and. z(t)<=0q0) then
         !w(t) is probably very small positive.
         !should have been zero but due to round off..
         !solve by moving index t back to R from P and set w(t)=0q0
         !Can only exit loop and hope round offs will be ok next time??

         !write(*,*) "P=",P
         !write(*,*)
         !write(*,*) "z(P)=",z(P)
         !write(*,*)
         !write(*,*) "t=",t
         !write(*,*) "z(t)=",z(t)
         !write(*,*) "w(t)=",w(t)
         !write(*,*) "x(t)=",x(t)

         if(.not. allocated(R)) then
            !Create R and add t to R
            allocate(R(1))
            R(1)=t
            !Remove t from P. 
            allocate(C(size(P,1)))
            C=P
            deallocate(P)
            allocate(P(size(C,1)-1))
            P=C(1:size(C,1)-1)
            deallocate(C)
         else
            !Add t to R
            allocate(C(size(R,1)+1))
            C(1:size(R,1))=R
            C(size(R,1)+1)=t
            deallocate(R)
            allocate(R(size(C,1)))
            R=C
            deallocate(C)
            !Remove t from P. 
            if(size(P,1)==1) then
               deallocate(P)
            else
               allocate(C(size(P,1)))
               C=P
               deallocate(P)
               allocate(P(size(C,1)-1))
               P=C(1:size(C,1)-1)
               deallocate(C)
            endif
         endif
         if(allocated(F)) then
            !Add t to F
            allocate(C(size(F,1)+1))
            C(1:size(F,1))=F
            C(size(C,1))=t
            deallocate(F)
            allocate(F(size(C,1)))
            F=C
            deallocate(C)
         else
            allocate(F(1))
            F(1)=t
         endif
         write(*,*) "t=",t
         write(*,'(a,E10.3,a)') "Warning: w(t)=",w(t),", is not zero due to round off errors."
         !stop "fuckit"
         exit
      else
         if(allocated(F)) deallocate(F)
      endif

      !Step 7
      k=0
      do i=1,size(P,1)
         if(z(P(i))>0) k=k+1
      enddo
      if(k .eq. size(P,1)) then
         x=z
         exit
      endif
      !Step 8
      if(size(P,1)==1) then
         q=1
      else
         q=minloc(x(P)/(x(P)-z(P)),1,z(P)<=0q0)
      endif
      q=P(q)
      !Step 9
      alpha=x(q)/(x(q)-z(q))
      if(alpha<0) then
         stop "alpha<0, not ok"
      elseif(alpha==0) then
         write(*,*) "Warning message from NNLS:"
         write(*,*) "alpha is zero due to round off errors"
         write(*,*) "It's probably alright since modified step 10 in the algorithm" 
         write(*,*) "Printouts below is just help in case one anyhow suspects something went wrong"
         write(*,*)
         write(*,*) "t=",t
         write(*,*) "w(t)=",w(t)
         write(*,*) "z(t)=",z(t)
         write(*,*)
         write(*,*) "P=",P
         write(*,*)
         write(*,*) "z(P)=",z(P)
         write(*,*)
         write(*,*) "x(P)=",x(P)
         write(*,*)
         write(*,*) "x(P)/(x(P)-z(P))=",x(P)/(x(P)-z(P))
         write(*,*)
         write(*,*) "q=",q
         write(*,*) "z(q)=",z(q)
         write(*,*) "x(q)=",x(q)
         write(*,*) "x(q)/(x(q)-z(q)=",x(q)/(x(q)-z(q))
         write(*,*) "tiny(z(q))=",tiny(z(q))
      endif
      if(alpha>1 ) stop "alpha>1, not ok"
      !Step 10
      where(x==x(q)) 
         x=0q0
      elsewhere
         x=x+alpha*(z-x)
      endwhere
     
      if(alpha<=0 .or. alpha>1 ) then
         write(*,*) "x(q)=",x(q) ," again..."
      endif
      
      !Step 11
      !Identify all indices j in P fulfilling x(P(j))==0
      !Store these j values in D
      allocate(C(size(P,1)))
      k=0 !variable storing number of indices to move
      do i=1,size(P,1)
         if( x(P(i)) <= 0q0 ) then !this should be .eq. but should not have negative so ok. Could also have x(P(i))<tiny(x(1)) instead.
            k=k+1
            C(k)=i
         endif
      enddo
      allocate(D(k))
      D=C(1:k)
      deallocate(C)
      !Add P(D) to R
      allocate(C(size(R,1)+k))
      C(1:size(R,1))=R
      C(size(R,1)+1:size(R,1)+k)=P(D)
      deallocate(R)
      allocate(R(size(C,1)))
      R=C
      deallocate(C)
      !Remove P(D) from P
      if(size(P,1)-k==0) then
         !nothing left
         deallocate(P)
      else
         allocate(C(size(P,1)-k))
         h=0
         do i=1,size(P,1)
            g=0
            do j=1,size(D,1)
               if(i==D(j)) then
                  g=g+1
               endif
            enddo
            if(g==0) then
               h=h+1
               C(h)=P(i) 
            endif
         enddo
         deallocate(P)
         allocate(P(size(C,1)))
         P=C
         deallocate(C)
      endif
      deallocate(D)

      ib=ib+1
      if(ib==2*n) then
         write(*,*) "Loop B terminated without fulfilling conditions"
         exit
      endif
   enddo
   ia=ia+1
   if(ia==3*n) exit
enddo
if(allocated(R)) deallocate(R)
if(allocated(F)) deallocate(F)
deallocate(z,w,P)

!write(*,'(a,I5,a)') "Did ",ib,"loop B iterations in total"
write(*,'(a,I5,a)') "Did ",ia," loop A iterations"
write(*,'(a,I5,a)') "Did ",ia+ib," LS calls in NNLS routine"
write(*,*)
end subroutine

!Not implemented completly!!!!
subroutine dls(Ao,bo,mtr,pr,x)
!Solves equation A*x=b in Least Square sense, using real variables.
!mtr=0 sends A and b to lapack. This system is solved in Least Square sense. 
!mtr=1 sends A^T*A and A^T*b to lapack. This system is solved exact.
!mtr=2 sends A to lapack for SVD. Construct explicitly Least Square solution.
!mtr=3 sends A to lapack for SVD. Construct explicitly Least Square solution, but truncate modes with small singular values.
implicit none
real(kind=16),intent(in) :: Ao(:,:),bo(:)
integer,intent(in) :: mtr
integer,intent(in) :: pr
real(kind=16),allocatable,intent(inout) :: x(:)
!internal
integer :: i,N,M,p,r,mo
real(kind=16),allocatable :: A(:,:),b(:),k(:,:),g(:) 
real(kind=16),allocatable :: u(:,:),v(:,:)
real(kind=16),allocatable :: s(:)
real(kind=16) :: alpha

N=size(Ao,1)
M=size(Ao,2)
p=min(N,M)
allocate(A(N,M),b(N))
b=bo !copy
A=Ao !copy

if(mtr==0) then ! Let Lapack solve equation: A*x=b in LS sense
   !call dgelswrapper(A,b,pr) !solves A*x=b using QR or LQ factorization
   !call dgelsswrapper(A,b)  !solves A*x=b using Singular Value Decomposition SVD (double precision only)
   call dgelsdwrapper(A,b,pr)  !solves A*x=b using Singular Value Decomposition SVD with Divide and Conquer
elseif(mtr==1) then ! Let Lapack exactly solve normal equation: A^dagger*A*x=A^dagger*b
   allocate(k(M,M),g(M))
   k=matmul(transpose(A),A)
   g=matmul(transpose(A),b)
   call dgelswrapper(k,g,pr) !solves k*x=g using QR or LQ factorization
   !call dgelsswrapper(k,g) !solvesm k*x=g using Singular Value Decomposition SVD (double precision only)
   !call dgelsdwrapper(k,g,pr) !solves k*x=g using Singular Value Decomposition SVD with Divide and Conquer
   b(1:M)=g
   deallocate(k,g)
elseif(mtr==2 .or. mtr==3) then !Explicitly do SVD and use singular values bigger than machine precision
   stop "SVD for double not yet implemented!"
   allocate(u(N,N),s(p),v(M,M))
   !call dgesddwrapper(A,u,s,v,pr)
   write(*,'(a)') "SVD's singular values="
   do i=1,p
      write(*,'(2E13.5)') s(i),s(i)/s(1)
   enddo
   write(*,*)
   if(pr==64) then
      r=minloc(abs(s-s(1)*0.5q0*10**(-15q0)),1)
   elseif(pr==128) then
      r=minloc(abs(s-s(1)*0.5q0*10**(-30q0)),1)
   else
      stop "precision not supported"
   endif
   write(*,'(I4,a)') r," is estimated SVD rank."
   write(*,'(E15.5,a)') s(1)/s(r)," is estimated condition number"
   if(mtr==2) then
      mo=r
   elseif(mtr==3) then !Explicitly do SVD and use singular values s(i)/s(1)>alpha
      alpha=10**(-4q0) !This is just an example
      mo=minloc(abs(s-s(1)*alpha),1) 
   endif
   write(*,'(I4,a,I4)') mo," SVD modes used for solution, out of",M
   !construct solution here
   allocate(g(M))
   g=0q0
   do i=1,mo
      g=g+1q0/s(i)*dot_product(u(:,i),b)*v(:,i)
   enddo
   b(1:M)=g
   deallocate(g)
   deallocate(u,s,v)
else
   stop "wrong mtr value"
endif
if(allocated(x)) deallocate(x)
allocate(x(M))
x=b(1:M)
deallocate(A,b)
end subroutine

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
complex(kind=16),allocatable :: x(:)
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

if(mtr==0) then ! Let Lapack solve equation: A*x=b in LS sense
   call zgelswrapper(A,b,pr)  !solves A*x=b using QR or LQ factorization
   !call zgelsswrapper(A,b) !solves A*x=b using Singular Value Decomposition SVD (double precision only)
   !call zgelsdwrapper(A,b)  !solves A*x=b using Singular Value Decomposition SVD with Divide and Conquer (double precision only) 
elseif(mtr==1) then ! Let Lapack exactly solve normal equation: A^dagger*A*x=A^dagger*b
   allocate(k(M,M),g(M))
   k=matmul(conjg(transpose(A)),A)
   g=matmul(conjg(transpose(A)),b)
   call zgelswrapper(k,g,pr) !solves k*x=g using QR or LQ factorization
   !call zgelsswrapper(k,g) !solvesm k*x=g using Singular Value Decomposition SVD (double precision only)
   !call zgelsdwrapper(k,g) !solves k*x=g using Singular Value Decomposition SVD with Divide and Conquer (double precision only)
   b(1:M)=g
   deallocate(k,g)
elseif(mtr==2 .or. mtr==3) then !Explicitly do SVD and use singular values bigger than machine precision
   allocate(u(N,N),s(p),v(M,M))
   call zgesddwrapper(A,u,s,v,pr)
   write(*,'(a)') "SVD's singular values="
   do i=1,p
      write(*,'(2E13.5)') s(i),s(i)/s(1)
   enddo
   write(*,*)
   if(pr==64) then
      r=minloc(abs(s-s(1)*0.5q0*10**(-15q0)),1)
   elseif(pr==128) then
      r=minloc(abs(s-s(1)*0.5q0*10**(-30q0)),1)
   else
      stop "precision not supported"
   endif
   write(*,'(I4,a)') r," is estimated SVD rank."
   write(*,'(E15.5,a)') s(1)/s(r)," is estimated condition number"
   if(mtr==2) then
      mo=r
   elseif(mtr==3) then !Explicitly do SVD and use singular values s(i)/s(1)>alpha
      alpha=10**(-4q0) !This is just an example
      mo=minloc(abs(s-s(1)*alpha),1) 
   endif
   write(*,'(I4,a,I4)') mo," SVD modes used for solution, out of",M
   !construct solution here
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
if(allocated(x)) deallocate(x)
allocate(x(M))
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
   stop "Precison not supported at the moment"
else
   stop "Precision not supported"
endif
end subroutine

subroutine dgelsdwrapper(A,b,pr)
!Computes the minimum-norm solution to a real linear least square problem: 
! minimize 2-norm( |A*x-b|)
!using the singular value decomposition (SVD) of A. A is an M-by-N matrix which
!may be rank-deficient.
!Accepts double (pr=64) and quad (pr=128)  precision.
implicit none
real(kind=16),intent(inout) :: A(:,:),b(:)
integer,intent(in) :: pr
!lapack
integer :: M,N,nrhs,lda,ldb,rank,lwork,info,liwork
integer,allocatable :: iwork(:)
real(kind=16),allocatable :: s(:),work(:)
real(kind=16) :: rcond
!double precision variables
real(kind=8),allocatable :: sd(:),workd(:)
real(kind=8) :: rcondd
real(kind=8),allocatable :: Ad(:,:),bd(:)

M=size(A,1)
N=size(A,2)
if(M<N) stop "N should be bigger or equal to M"
nrhs=1 ! assume only one right-hand side
lda=max(N,M)
ldb=max(M,nrhs)
lwork=-1 !to determine optimal work size
rcond=-1 !use to determine the effective rank of A. Singular values s(i) <=rcond*s(1) are treated as zero. If rcond<0, machine precision is used instead.
if(pr==64) then
   allocate(Ad(M,N),bd(M))
   Ad=A
   bd=b
   rcondd=rcond
   allocate(iwork(1),workd(1))
   allocate(sd(min(N,M)))
   call dgelsd(M,N,nrhs,Ad,lda,Bd,ldb,sd,rcondd,rank,workd,lwork,iwork,info) !find optimal lwork and liwork
   if(info<0) then
      stop "illegal argument"
   elseif(info>0) then
      stop "SVD failed to converge"
   endif
   lwork=workd(1)
   liwork=iwork(1)
   deallocate(workd,iwork)
   allocate(workd(lwork),iwork(liwork))
   call dgelsd(M,N,nrhs,Ad,lda,Bd,ldb,sd,rcondd,rank,workd,lwork,iwork,info)
   deallocate(workd,iwork)
   b=bd !transfer solution to quad variable
   deallocate(Ad,bd)
   deallocate(sd)
elseif(pr==128) then
   allocate(iwork(1),work(1))
   allocate(s(min(N,M)))
   call dgelsdquad(M,N,nrhs,A,lda,B,ldb,s,rcond,rank,work,lwork,iwork,info) !find optimal lwork and liwork
   if(info<0) then
      stop "illegal argument"
   elseif(info>0) then
      stop "SVD failed to converge"
   endif
   lwork=work(1)
   liwork=iwork(1)
   deallocate(work,iwork)
   allocate(work(lwork),iwork(liwork))
   call dgelsdquad(M,N,nrhs,A,lda,B,ldb,s,rcond,rank,work,lwork,iwork,info)
   deallocate(work,iwork)
   deallocate(s)
else
   stop "dgelsd only implemented for double and quad precsion"
endif
if(info .ne. 0 ) stop "Something went wrong in dgelsdwrapper"
!b(1:N) is the solution
end subroutine

subroutine dgelswrapper(A,b,pr)
!solves overdetermined or underdetermined real linear systems, assuming full rank of matrix A.
!Accepts double (pr=64) and quad (pr=128)  precision.
implicit none
real(kind=16),intent(in) :: A(:,:)
real(kind=16),intent(inout) :: b(:)
integer,intent(in) :: pr
integer :: N,M
!lapack
integer :: lwork
real(kind=16),allocatable :: work(:)
integer :: info
!double precision variables
real(kind=8),allocatable :: workd(:)
real(kind=8),allocatable :: Ad(:,:),bd(:)
N=size(A,1)
M=size(A,2)
if(N<M) stop "N should be bigger or equal to M"
allocate(work(2*M))
lwork=-1
if(pr==64) then
   allocate(Ad(N,M),bd(N),workd(2*M))
   Ad=A
   bd=b
   call dgels('N',N,M,1,Ad,N,bd,N,workd,lwork,info)
   if(info .ne. 0) stop "least square routine do not work properly"
   lwork=ceiling(workd(1)) !get optimal lwork
   deallocate(workd)
   if(lwork < 2*M ) stop "lwork is simply too low. Increase lwork in lapace routine."
   allocate(workd(lwork))
   call dgels('N',N,M,1,Ad,N,bd,N,workd,lwork,info)
   b=bd ! convert solution to quad precision variable
   deallocate(Ad,bd,workd)
elseif(pr==128) then
   stop "quad precision for dgels is yet not implemented!"
   !call dgelsquad('N',N,M,1,A,N,b,N,work,lwork,info)
   if(info .ne. 0) stop "least square routine do not work properly"
   lwork=ceiling(work(1)) !get optimal lwork
   deallocate(work)
   if(lwork < 2*M ) stop "lwork is simply too low. Increase lwork in lapace routine."
   allocate(work(lwork))
   !call dgelsquad('N',N,M,1,A,N,b,N,work,lwork,info)
else
   stop "Lapack's zgels can only handle dbl or quad"
endif
if(info<0) then
    stop "an argument with illegal value."
elseif(info>0) then
   write(*,*) "Matrix A does not have full rank, which is assumed. LS-solution could not be computed"
endif
deallocate(work)
!b(1:M) is the solution
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
   write(*,*) "Matrix A does not have full rank, which is assumed. LS-solution could not be computed"
endif
deallocate(work)
!b(1:M) is the solution
end subroutine

subroutine zgelsquadwrapper(A,b)
!solves overdetermined or underdetermined complex linear systems, assuming full rank of matrix A
!Uses quadruple precision
implicit none
complex(kind=16),intent(inout) :: A(:,:),b(:)
integer :: N,M
!lapack
integer :: lwork
complex(kind=16),allocatable :: work(:)
integer :: info
N=size(A,1)
M=size(A,2)
if(N<M) stop "N should be bigger or equal to M"
allocate(work(2*M))
lwork=-1
call zgelsquad('N',N,M,1,A,N,b,N,work,lwork,info)
if(info .ne. 0) stop "least square routine do not work properly"
lwork=ceiling(real(work(1))) !get optimal lwork
deallocate(work)
if(lwork < 2*M ) stop "lwork is simply too low. Increase lwork in lapace routine."
allocate(work(lwork))
call zgelsquad('N',N,M,1,A,N,b,N,work,lwork,info)
if(info<0) then
    stop "an argument with illegal value."
elseif(info>0) then
write(*,*) "Matrix A does not have full rank, which is assumed. LS-solution could not be computed"
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
rcond=-1d0  !parameter to determine effective rank of A
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
    write(*,*) "least square algorithm failed to converge."
endif
write(*,*) " "
write(*,'(A,E12.4)') "Condition number: ",s(1)/s(M) 
write(*,'(A,I3,A,I3,A)') "Effective rank of matrix:",rank," (compare with M=",M,")"
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
rcond=-1d0  !parameter to determine effective rank of A
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
write(*,*) "liwork=",liwork
write(*,*) "lwork=",lwork
call ZGELSD(N,M,1,Ad,N,bd,N,s,rcond,rank,work,lwork,rwork,iwork,info) !bd(1:M) contains solution
if(info<0) then
    stop "an argument with illegal value."
elseif(info>0) then
    write(*,*) "least square algorithm failed to converge."
endif
write(*,*) " "
write(*,'(A,E12.4)') "Condition number: ",s(1)/s(M) 
write(*,'(A,I3,A,I3,A)') "Effective rank of matrix:",rank," (compare with M=",M,")"
deallocate(s,work,rwork,iwork)
b=bd !convert from double to quad precision
deallocate(Ad,bd)
end subroutine

end module 
