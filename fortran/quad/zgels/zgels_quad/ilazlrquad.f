*> \brief \b ilazlrquad scans a matrix for its last non-zero row.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ilazlrquad + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilazlr.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilazlr.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilazlr.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       INTEGER FUNCTION ilazlrquad( M, N, A, LDA )
* 
*       .. Scalar Arguments ..
*       INTEGER            M, N, LDA
*       ..
*       .. Array Arguments ..
*       complex*32         A( LDA, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ilazlrquad scans A for its last non-zero row.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is complex*32 array, dimension (LDA,N)
*>          The m by n matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A. LDA >= max(1,M).
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
      INTEGER FUNCTION ilazlrquad( M, N, A, LDA )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      INTEGER            M, N, LDA
*     ..
*     .. Array Arguments ..
      complex*32         A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      complex*32       ZERO
      PARAMETER ( ZERO = (0.0q0, 0.0q0) )
*     ..
*     .. Local Scalars ..
      INTEGER I, J
*     ..
*     .. Executable Statements ..
*
*     Quick test for the common case where one corner is non-zero.
      IF( M.EQ.0 ) THEN
         ilazlrquad = M
      ELSE IF( A(M, 1).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
         ilazlrquad = M
      ELSE
*     Scan up each column tracking the last zero row seen.
         ilazlrquad = 0
         DO J = 1, N
            I=M
            DO WHILE((A(MAX(I,1),J).EQ.ZERO).AND.(I.GE.1))
               I=I-1
            ENDDO
            ilazlrquad = MAX( ilazlrquad, I )
         END DO
      END IF
      RETURN
      END
