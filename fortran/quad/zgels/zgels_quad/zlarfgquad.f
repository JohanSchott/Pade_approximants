*> \brief \b zlarfgquad generates an elementary reflector (Householder matrix).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download zlarfgquad + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarfg.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarfg.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarfg.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE zlarfgquad( N, ALPHA, X, INCX, TAU )
* 
*       .. Scalar Arguments ..
*       INTEGER            INCX, N
*       complex*32         ALPHA, TAU
*       ..
*       .. Array Arguments ..
*       complex*32         X( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> zlarfgquad generates a complex elementary reflector H of order n, such
*> that
*>
*>       H**H * ( alpha ) = ( beta ),   H**H * H = I.
*>              (   x   )   (   0  )
*>
*> where alpha and beta are scalars, with beta real, and x is an
*> (n-1)-element complex vector. H is represented in the form
*>
*>       H = I - tau * ( 1 ) * ( 1 v**H ) ,
*>                     ( v )
*>
*> where tau is a complex scalar and v is a complex (n-1)-element
*> vector. Note that H is not hermitian.
*>
*> If the elements of x are all zero and alpha is real, then tau = 0
*> and H is taken to be the unit matrix.
*>
*> Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1 .
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the elementary reflector.
*> \endverbatim
*>
*> \param[in,out] ALPHA
*> \verbatim
*>          ALPHA is complex*32
*>          On entry, the value alpha.
*>          On exit, it is overwritten with the value beta.
*> \endverbatim
*>
*> \param[in,out] X
*> \verbatim
*>          X is complex*32 array, dimension
*>                         (1+(N-2)*abs(INCX))
*>          On entry, the vector x.
*>          On exit, it is overwritten with the vector v.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>          The increment between elements of X. INCX > 0.
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is complex*32
*>          The value tau.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date September 2012
*
*> \ingroup complex16OTHERauxiliary
*
*  =====================================================================
      SUBROUTINE zlarfgquad( N, ALPHA, X, INCX, TAU )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      complex*32         ALPHA, TAU
*     ..
*     .. Array Arguments ..
      complex*32         X( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      real*16   ONE, ZERO
      PARAMETER          ( ONE = 1.0q0, ZERO = 0.0q0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, KNT
      real*16   ALPHI, ALPHR, BETA, RSAFMN, SAFMIN, XNORM
*     ..
*     .. External Functions ..
      real*16   dlamchquad, dlapy3quad, dznrm2quad
      complex*32         zladivquad
      EXTERNAL           dlamchquad, dlapy3quad, dznrm2quad, zladivquad
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, real, CMPLX, AIMAG, SIGN
*     ..
*     .. External Subroutines ..
      EXTERNAL           zdscalquad, zscalquad
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 ) THEN
         TAU = ZERO
         RETURN
      END IF
*
      XNORM = dznrm2quad( N-1, X, INCX )
      ALPHR = real( ALPHA )
      ALPHI = AIMAG( ALPHA )
*
      IF( XNORM.EQ.ZERO .AND. ALPHI.EQ.ZERO ) THEN
*
*        H  =  I
*
         TAU = ZERO
      ELSE
*
*        general case
*
         BETA = -SIGN( dlapy3quad( ALPHR, ALPHI, XNORM ), ALPHR )
         SAFMIN = dlamchquad( 'S' ) / dlamchquad( 'E' )
         RSAFMN = ONE / SAFMIN
*
         KNT = 0
         IF( ABS( BETA ).LT.SAFMIN ) THEN
*
*           XNORM, BETA may be inaccurate; scale X and recompute them
*
   10       CONTINUE
            KNT = KNT + 1
            CALL zdscalquad( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHI = ALPHI*RSAFMN
            ALPHR = ALPHR*RSAFMN
            IF( ABS( BETA ).LT.SAFMIN )
     $         GO TO 10
*
*           New BETA is at most 1, at least SAFMIN
*
            XNORM = dznrm2quad( N-1, X, INCX )
            ALPHA = CMPLX( ALPHR, ALPHI,16 )
            BETA = -SIGN( dlapy3quad( ALPHR, ALPHI, XNORM ), ALPHR )
         END IF
         TAU = CMPLX( ( BETA-ALPHR ) / BETA, -ALPHI / BETA,16 )
         ALPHA = zladivquad( CMPLX( ONE,0q0,16 ), ALPHA-BETA )
         CALL zscalquad( N-1, ALPHA, X, INCX )
*
*        If ALPHA is subnormal, it may lose relative accuracy
*
         DO 20 J = 1, KNT
            BETA = BETA*SAFMIN
 20      CONTINUE
         ALPHA = BETA
      END IF
*
      RETURN
*
*     End of zlarfgquad
*
      END
