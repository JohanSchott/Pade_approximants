*> \brief \b dlapy3quad returns sqrt(x2+y2+z2).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download dlapy3quad + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlapy3.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlapy3.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlapy3.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       real*16 FUNCTION dlapy3quad( X, Y, Z )
* 
*       .. Scalar Arguments ..
*       real*16   X, Y, Z
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> dlapy3quad returns sqrt(x**2+y**2+z**2), taking care not to cause
*> unnecessary overflow.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] X
*> \verbatim
*>          X is real*16
*> \endverbatim
*>
*> \param[in] Y
*> \verbatim
*>          Y is real*16
*> \endverbatim
*>
*> \param[in] Z
*> \verbatim
*>          Z is real*16
*>          X, Y and Z specify the values x, y and z.
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
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      real*16 FUNCTION dlapy3quad( X, Y, Z )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      real*16   X, Y, Z
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      real*16   ZERO
      PARAMETER          ( ZERO = 0.0q0 )
*     ..
*     .. Local Scalars ..
      real*16   W, XABS, YABS, ZABS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
      XABS = ABS( X )
      YABS = ABS( Y )
      ZABS = ABS( Z )
      W = MAX( XABS, YABS, ZABS )
      IF( W.EQ.ZERO ) THEN
*     W can be zero for max(0,nan,0)
*     adding all three entries together will make sure
*     NaN will not disappear.
         dlapy3quad =  XABS + YABS + ZABS
      ELSE
         dlapy3quad = W*SQRT( ( XABS / W )**2+( YABS / W )**2+
     $            ( ZABS / W )**2 )
      END IF
      RETURN
*
*     End of dlapy3quad
*
      END
