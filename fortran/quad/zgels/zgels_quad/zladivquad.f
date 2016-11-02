*> \brief \b zladivquad performs complex division in real arithmetic, avoiding unnecessary overflow.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download zladivquad + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zladiv.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zladiv.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zladiv.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       complex*32     FUNCTION zladivquad( X, Y )
* 
*       .. Scalar Arguments ..
*       complex*32         X, Y
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> zladivquad := X / Y, where X and Y are complex.  The computation of X / Y
*> will not overflow on an intermediary step unless the results
*> overflows.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] X
*> \verbatim
*>          X is complex*32
*> \endverbatim
*>
*> \param[in] Y
*> \verbatim
*>          Y is complex*32
*>          The complex scalars X and Y.
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
      complex*32     FUNCTION zladivquad( X, Y )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      complex*32         X, Y
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      real*16   ZI, ZR
*     ..
*     .. External Subroutines ..
      EXTERNAL           dladivquad
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          real, CMPLX, AIMAG
*     ..
*     .. Executable Statements ..
*
      CALL dladivquad( real( X ), AIMAG( X ), real( Y ), AIMAG( Y ),
     $             ZR,ZI )
      zladivquad = CMPLX( ZR, ZI,16 )
*
      RETURN
*
*     End of zladivquad
*
      END
