*> \brief \b ztrtrsquad
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ztrtrsquad + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrtrs.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrtrs.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrtrs.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ztrtrsquad( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB,
*                          INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          DIAG, TRANS, UPLO
*       INTEGER            INFO, LDA, LDB, N, NRHS
*       ..
*       .. Array Arguments ..
*       complex*32         A( LDA, * ), B( LDB, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ztrtrsquad solves a triangular system of the form
*>
*>    A * X = B,  A**T * X = B,  or  A**H * X = B,
*>
*> where A is a triangular matrix of order N, and B is an N-by-NRHS
*> matrix.  A check is made to verify that A is nonsingular.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          = 'U':  A is upper triangular;
*>          = 'L':  A is lower triangular.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          Specifies the form of the system of equations:
*>          = 'N':  A * X = B     (No transpose)
*>          = 'T':  A**T * X = B  (Transpose)
*>          = 'C':  A**H * X = B  (Conjugate transpose)
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>          = 'N':  A is non-unit triangular;
*>          = 'U':  A is unit triangular.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>          The number of right hand sides, i.e., the number of columns
*>          of the matrix B.  NRHS >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is complex*32 array, dimension (LDA,N)
*>          The triangular matrix A.  If UPLO = 'U', the leading N-by-N
*>          upper triangular part of the array A contains the upper
*>          triangular matrix, and the strictly lower triangular part of
*>          A is not referenced.  If UPLO = 'L', the leading N-by-N lower
*>          triangular part of the array A contains the lower triangular
*>          matrix, and the strictly upper triangular part of A is not
*>          referenced.  If DIAG = 'U', the diagonal elements of A are
*>          also not referenced and are assumed to be 1.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is complex*32 array, dimension (LDB,NRHS)
*>          On entry, the right hand side matrix B.
*>          On exit, if INFO = 0, the solution matrix X.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
*>          > 0: if INFO = i, the i-th diagonal element of A is zero,
*>               indicating that the matrix is singular and the solutions
*>               X have not been computed.
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
*> \date November 2011
*
*> \ingroup complex16OTHERcomputational
*
*  =====================================================================
      SUBROUTINE ztrtrsquad( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB,
     $                   INFO )
*
*  -- LAPACK computational routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      complex*32         A( LDA, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      complex*32         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0q0, 0.0q0 ),
     $                   ONE = ( 1.0q0, 0.0q0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT
*     ..
*     .. External Functions ..
      LOGICAL            lsamequad
      EXTERNAL           lsamequad
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerblaquad, ztrsmquad
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NOUNIT = lsamequad( DIAG, 'N' )
      IF( .NOT.lsamequad( UPLO, 'U' ) .AND.
     $   .NOT.lsamequad( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.lsamequad( TRANS, 'N' ) .AND. .NOT.
     $ lsamequad( TRANS, 'T' ) .AND. .NOT.lsamequad( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOUNIT .AND. .NOT.lsamequad( DIAG, 'U' ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL xerblaquad( 'ztrtrsquad', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Check for singularity.
*
      IF( NOUNIT ) THEN
         DO 10 INFO = 1, N
            IF( A( INFO, INFO ).EQ.ZERO )
     $         RETURN
   10    CONTINUE
      END IF
      INFO = 0
*
*     Solve A * x = b,  A**T * x = b,  or  A**H * x = b.
*
      CALL ztrsmquad( 'Left', UPLO, TRANS, DIAG, N, NRHS, ONE,
     $     A, LDA, B,    LDB )
*
      RETURN
*
*     End of ztrtrsquad
*
      END
