*> \brief \b dladivquad performs complex division in real arithmetic, avoiding unnecessary overflow.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download dladivquad + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dladiv.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dladiv.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dladiv.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE dladivquad( A, B, C, D, P, Q )
* 
*       .. Scalar Arguments ..
*       real*16   A, B, C, D, P, Q
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> dladivquad performs complex division in  real arithmetic
*>
*>                       a + i*b
*>            p + i*q = ---------
*>                       c + i*d
*>
*> The algorithm is due to Michael Baudin and Robert L. Smith
*> and can be found in the paper
*> "A Robust Complex Division in Scilab"
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] A
*> \verbatim
*>          A is real*16
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is real*16
*> \endverbatim
*>
*> \param[in] C
*> \verbatim
*>          C is real*16
*> \endverbatim
*>
*> \param[in] D
*> \verbatim
*>          D is real*16
*>          The scalars a, b, c, and d in the above expression.
*> \endverbatim
*>
*> \param[out] P
*> \verbatim
*>          P is real*16
*> \endverbatim
*>
*> \param[out] Q
*> \verbatim
*>          Q is real*16
*>          The scalars p and q in the above expression.
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
*> \date January 2013
*
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      SUBROUTINE dladivquad( A, B, C, D, P, Q )
*
*  -- LAPACK auxiliary routine (version 3.5.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     January 2013
*
*     .. Scalar Arguments ..
      real*16   A, B, C, D, P, Q
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      real*16   BS
      PARAMETER          ( BS = 2.0q0 )
      real*16   HALF
      PARAMETER          ( HALF = 0.5q0 )
      real*16   TWO
      PARAMETER          ( TWO = 2.0q0 )
*
*     .. Local Scalars ..
      real*16   AA, BB, CC, DD, AB, CD, S, OV, UN, BE, EPS
*     ..
*     .. External Functions ..
      real*16   dlamchquad
      EXTERNAL           dlamchquad
*     ..
*     .. External Subroutines ..
      EXTERNAL           dladivquad1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
      AA = A
      BB = B
      CC = C
      DD = D
      AB = MAX( ABS(A), ABS(B) )
      CD = MAX( ABS(C), ABS(D) )
      S = 1.0q0
      
      OV = dlamchquad( 'Overflow threshold' )
      UN = dlamchquad( 'Safe minimum' )
      EPS = dlamchquad( 'Epsilon' )
      BE = BS / (EPS*EPS)
      
      IF( AB >= HALF*OV ) THEN
         AA = HALF * AA
         BB = HALF * BB
         S  = TWO * S
      END IF
      IF( CD >= HALF*OV ) THEN
         CC = HALF * CC
         DD = HALF * DD
         S  = HALF * S
      END IF
      IF( AB <= UN*BS/EPS ) THEN
         AA = AA * BE
         BB = BB * BE
         S  = S / BE
      END IF
      IF( CD <= UN*BS/EPS ) THEN
         CC = CC * BE
         DD = DD * BE
         S  = S * BE
      END IF
      IF( ABS( D ).LE.ABS( C ) ) THEN
         CALL dladivquad1(AA, BB, CC, DD, P, Q)
      ELSE
         CALL dladivquad1(BB, AA, DD, CC, P, Q)
         Q = -Q
      END IF
      P = P * S
      Q = Q * S
*
      RETURN
*
*     End of dladivquad
*
      END

      

      SUBROUTINE dladivquad1( A, B, C, D, P, Q )
*
*  -- LAPACK auxiliary routine (version 3.5.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     January 2013
*
*     .. Scalar Arguments ..
      real*16   A, B, C, D, P, Q
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      real*16   ONE
      PARAMETER          ( ONE = 1.0q0 )
*
*     .. Local Scalars ..
      real*16   R, T
*     ..
*     .. External Functions ..
      real*16   dladivquad2
      EXTERNAL           dladivquad2
*     ..
*     .. Executable Statements ..
*
      R = D / C
      T = ONE / (C + D * R)
      P = dladivquad2(A, B, C, D, R, T)
      A = -A
      Q = dladivquad2(B, A, C, D, R, T)
*
      RETURN
*
*     End of dladivquad1
*
      END

      real*16 FUNCTION dladivquad2( A, B, C, D, R, T )
*
*  -- LAPACK auxiliary routine (version 3.5.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     January 2013
*
*     .. Scalar Arguments ..
      real*16   A, B, C, D, R, T
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      real*16   ZERO
      PARAMETER          ( ZERO = 0.0q0 )
*
*     .. Local Scalars ..
      real*16   BR
*     ..
*     .. Executable Statements ..
*
      IF( R.NE.ZERO ) THEN
         BR = B * R
         if( BR.NE.ZERO ) THEN
            dladivquad2 = (A + BR) * T
         ELSE
            dladivquad2 = A * T + (B * T) * R
         END IF
      ELSE
         dladivquad2 = (A + D * (B / C)) * T
      END IF
*
      RETURN
*
*     End of dladivquad12
*
      END
