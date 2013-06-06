      SUBROUTINE ZGGEVR( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA, BETA,
     $                  VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )
*
*  -- LAPACK driver routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     MODIFIED BY R.G.V. TO USE MY VERSION OF ZHGEQZ.f (ZHGERV)
*
*     .. Scalar Arguments ..
      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ),
     $                   BETA( * ), VL( LDVL, * ), VR( LDVR, * ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZGGEV computes for a pair of N-by-N complex nonsymmetric matrices
*  (A,B), the generalized eigenvalues, and optionally, the left and/or
*  right generalized eigenvectors.
*
*  A generalized eigenvalue for a pair of matrices (A,B) is a scalar
*  lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
*  singular. It is usually represented as the pair (alpha,beta), as
*  there is a reasonable interpretation for beta=0, and even for both
*  being zero.
*
*  The right generalized eigenvector v(j) corresponding to the
*  generalized eigenvalue lambda(j) of (A,B) satisfies
*
*               A * v(j) = lambda(j) * B * v(j).
*
*  The left generalized eigenvector u(j) corresponding to the
*  generalized eigenvalues lambda(j) of (A,B) satisfies
*
*               u(j)**H * A = lambda(j) * u(j)**H * B
*
*  where u(j)**H is the conjugate-transpose of u(j).
*
*  Arguments
*  =========
*
*  JOBVL   (input) CHARACTER*1
*          = 'N':  do not compute the left generalized eigenvectors;
*          = 'V':  compute the left generalized eigenvectors.
*
*  JOBVR   (input) CHARACTER*1
*          = 'N':  do not compute the right generalized eigenvectors;
*          = 'V':  compute the right generalized eigenvectors.
*
*  N       (input) INTEGER
*          The order of the matrices A, B, VL, and VR.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA, N)
*          On entry, the matrix A in the pair (A,B).
*          On exit, A has been overwritten.
*
*  LDA     (input) INTEGER
*          The leading dimension of A.  LDA >= max(1,N).
*
*  B       (input/output) COMPLEX*16 array, dimension (LDB, N)
*          On entry, the matrix B in the pair (A,B).
*          On exit, B has been overwritten.
*
*  LDB     (input) INTEGER
*          The leading dimension of B.  LDB >= max(1,N).
*
*  ALPHA   (output) COMPLEX*16 array, dimension (N)
*  BETA    (output) COMPLEX*16 array, dimension (N)
*          On exit, ALPHA(j)/BETA(j), j=1,...,N, will be the
*          generalized eigenvalues.
*
*          Note: the quotients ALPHA(j)/BETA(j) may easily over- or
*          underflow, and BETA(j) may even be zero.  Thus, the user
*          should avoid naively computing the ratio alpha/beta.
*          However, ALPHA will be always less than and usually
*          comparable with norm(A) in magnitude, and BETA always less
*          than and usually comparable with norm(B).
*
*  VL      (output) COMPLEX*16 array, dimension (LDVL,N)
*          If JOBVL = 'V', the left generalized eigenvectors u(j) are
*          stored one after another in the columns of VL, in the same
*          order as their eigenvalues.
*          Each eigenvector is scaled so the largest component has
*          abs(real part) + abs(imag. part) = 1.
*          Not referenced if JOBVL = 'N'.
*
*  LDVL    (input) INTEGER
*          The leading dimension of the matrix VL. LDVL >= 1, and
*          if JOBVL = 'V', LDVL >= N.
*
*  VR      (output) COMPLEX*16 array, dimension (LDVR,N)
*          If JOBVR = 'V', the right generalized eigenvectors v(j) are
*          stored one after another in the columns of VR, in the same
*          order as their eigenvalues.
*          Each eigenvector is scaled so the largest component has
*          abs(real part) + abs(imag. part) = 1.
*          Not referenced if JOBVR = 'N'.
*
*  LDVR    (input) INTEGER
*          The leading dimension of the matrix VR. LDVR >= 1, and
*          if JOBVR = 'V', LDVR >= N.
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,2*N).
*          For good performance, LWORK must generally be larger.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  RWORK   (workspace/output) DOUBLE PRECISION array, dimension (8*N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          =1,...,N:
*                The QZ iteration failed.  No eigenvectors have been
*                calculated, but ALPHA(j) and BETA(j) should be
*                correct for j=INFO+1,...,N.
*          > N:  =N+1: other then QZ iteration failed in DHGEQZ,
*                =N+2: error return from DTGEVC.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D0, 0.0D0 ),
     $                   CONE = ( 1.0D0, 0.0D0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILASCL, ILBSCL, ILV, ILVL, ILVR, LQUERY
      CHARACTER          CHTEMP
      INTEGER            ICOLS, IERR, IHI, IJOBVL, IJOBVR, ILEFT, ILO,
     $                   IN, IRIGHT, IROWS, IRWRK, ITAU, IWRK, JC, JR,
     $                   LWKMIN, LWKOPT
      DOUBLE PRECISION   ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS,
     $                   SMLNUM, TEMP
      COMPLEX*16         X
*     ..
*     .. Local Arrays ..
      LOGICAL            LDUMMA( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLABAD, XERBLA, ZGEQRF, ZGGBAK, ZGGBAL, ZGGHRD,
     $                   ZHGERV, ZLACPY, ZLASCL, ZLASET, ZTGEVC, ZUNGQR,
     $                   ZUNMQR
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, ZLANGE
      EXTERNAL           LSAME, ILAENV, DLAMCH, ZLANGE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG, MAX, SQRT
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   ABS1
*     ..
*     .. Statement Function definitions ..
      ABS1( X ) = ABS( DBLE( X ) ) + ABS( DIMAG( X ) )
*     ..
*     .. Executable Statements ..
*
*     Decode the input arguments
*
      IF( LSAME( JOBVL, 'N' ) ) THEN
         IJOBVL = 1
         ILVL = .FALSE.
      ELSE IF( LSAME( JOBVL, 'V' ) ) THEN
         IJOBVL = 2
         ILVL = .TRUE.
      ELSE
         IJOBVL = -1
         ILVL = .FALSE.
      END IF
*
      IF( LSAME( JOBVR, 'N' ) ) THEN
         IJOBVR = 1
         ILVR = .FALSE.
      ELSE IF( LSAME( JOBVR, 'V' ) ) THEN
         IJOBVR = 2
         ILVR = .TRUE.
      ELSE
         IJOBVR = -1
         ILVR = .FALSE.
      END IF
      ILV = ILVL .OR. ILVR
*
*     Test the input arguments
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      IF( IJOBVL.LE.0 ) THEN
         INFO = -1
      ELSE IF( IJOBVR.LE.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDVL.LT.1 .OR. ( ILVL .AND. LDVL.LT.N ) ) THEN
         INFO = -11
      ELSE IF( LDVR.LT.1 .OR. ( ILVR .AND. LDVR.LT.N ) ) THEN
         INFO = -13
      END IF
*
*     Compute workspace
*      (Note: Comments in the code beginning "Workspace:" describe the
*       minimal amount of workspace needed at that point in the code,
*       as well as the preferred amount for good performance.
*       NB refers to the optimal block size for the immediately
*       following subroutine, as returned by ILAENV. The workspace is
*       computed assuming ILO = 1 and IHI = N, the worst case.)
*
      IF( INFO.EQ.0 ) THEN
         LWKMIN = MAX( 1, 2*N )
         LWKOPT = MAX( 1, N + N*ILAENV( 1, 'ZGEQRF', ' ', N, 1, N, 0 ) )
         LWKOPT = MAX( LWKOPT, N +
     $                 N*ILAENV( 1, 'ZUNMQR', ' ', N, 1, N, 0 ) )
         IF( ILVL ) THEN
            LWKOPT = MAX( LWKOPT, N +
     $                    N*ILAENV( 1, 'ZUNGQR', ' ', N, 1, N, -1 ) )
         END IF
         WORK( 1 ) = LWKOPT
*
         IF( LWORK.LT.LWKMIN .AND. .NOT.LQUERY )
     $      INFO = -15
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGGEV ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Get machine constants
*
      EPS = DLAMCH( 'E' )*DLAMCH( 'B' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM
*
*     Scale A if max element outside range [SMLNUM,BIGNUM]
*
      ANRM = ZLANGE( 'M', N, N, A, LDA, RWORK )
      ILASCL = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         ANRMTO = SMLNUM
         ILASCL = .TRUE.
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         ANRMTO = BIGNUM
         ILASCL = .TRUE.
      END IF
      IF( ILASCL )
     $   CALL ZLASCL( 'G', 0, 0, ANRM, ANRMTO, N, N, A, LDA, IERR )
*
*     Scale B if max element outside range [SMLNUM,BIGNUM]
*
      BNRM = ZLANGE( 'M', N, N, B, LDB, RWORK )
      ILBSCL = .FALSE.
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
         BNRMTO = SMLNUM
         ILBSCL = .TRUE.
      ELSE IF( BNRM.GT.BIGNUM ) THEN
         BNRMTO = BIGNUM
         ILBSCL = .TRUE.
      END IF
      IF( ILBSCL )
     $   CALL ZLASCL( 'G', 0, 0, BNRM, BNRMTO, N, N, B, LDB, IERR )
*
*     Permute the matrices A, B to isolate eigenvalues if possible
*     (Real Workspace: need 6*N)
*
      ILEFT = 1
      IRIGHT = N + 1
      IRWRK = IRIGHT + N
      CALL ZGGBAL( 'P', N, A, LDA, B, LDB, ILO, IHI, RWORK( ILEFT ),
     $             RWORK( IRIGHT ), RWORK( IRWRK ), IERR )
*
*     Reduce B to triangular form (QR decomposition of B)
*     (Complex Workspace: need N, prefer N*NB)
*
      IROWS = IHI + 1 - ILO
      IF( ILV ) THEN
         ICOLS = N + 1 - ILO
      ELSE
         ICOLS = IROWS
      END IF
      ITAU = 1
      IWRK = ITAU + IROWS
      CALL ZGEQRF( IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ),
     $             WORK( IWRK ), LWORK+1-IWRK, IERR )
*
*     Apply the orthogonal transformation to matrix A
*     (Complex Workspace: need N, prefer N*NB)
*
      CALL ZUNMQR( 'L', 'C', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB,
     $             WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ),
     $             LWORK+1-IWRK, IERR )
*
*     Initialize VL
*     (Complex Workspace: need N, prefer N*NB)
*
      IF( ILVL ) THEN
         CALL ZLASET( 'Full', N, N, CZERO, CONE, VL, LDVL )
         IF( IROWS.GT.1 ) THEN
            CALL ZLACPY( 'L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB,
     $                   VL( ILO+1, ILO ), LDVL )
         END IF
         CALL ZUNGQR( IROWS, IROWS, IROWS, VL( ILO, ILO ), LDVL,
     $                WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR )
      END IF
*
*     Initialize VR
*
      IF( ILVR )
     $   CALL ZLASET( 'Full', N, N, CZERO, CONE, VR, LDVR )
*
*     Reduce to generalized Hessenberg form
*
      IF( ILV ) THEN
*
*        Eigenvectors requested -- work on whole matrix.
*
         CALL ZGGHRD( JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, VL,
     $                LDVL, VR, LDVR, IERR )
      ELSE
         CALL ZGGHRD( 'N', 'N', IROWS, 1, IROWS, A( ILO, ILO ), LDA,
     $                B( ILO, ILO ), LDB, VL, LDVL, VR, LDVR, IERR )
      END IF
*
*     Perform QZ algorithm (Compute eigenvalues, and optionally, the
*     Schur form and Schur vectors)
*     (Complex Workspace: need N)
*     (Real Workspace: need N)
*
      IWRK = ITAU
      IF( ILV ) THEN
         CHTEMP = 'S'
      ELSE
         CHTEMP = 'E'
      END IF
      CALL ZHGERV( CHTEMP, JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB,
     $             ALPHA, BETA, VL, LDVL, VR, LDVR, WORK( IWRK ),
     $             LWORK+1-IWRK, RWORK( IRWRK ), IERR )
      IF( IERR.NE.0 ) THEN
         IF( IERR.GT.0 .AND. IERR.LE.N ) THEN
            INFO = IERR
         ELSE IF( IERR.GT.N .AND. IERR.LE.2*N ) THEN
            INFO = IERR - N
         ELSE
            INFO = N + 1
         END IF
         GO TO 70
      END IF
*
*     Compute Eigenvectors
*     (Real Workspace: need 2*N)
*     (Complex Workspace: need 2*N)
*
      IF( ILV ) THEN
         IF( ILVL ) THEN
            IF( ILVR ) THEN
               CHTEMP = 'B'
            ELSE
               CHTEMP = 'L'
            END IF
         ELSE
            CHTEMP = 'R'
         END IF
*
         CALL ZTGEVC( CHTEMP, 'B', LDUMMA, N, A, LDA, B, LDB, VL, LDVL,
     $                VR, LDVR, N, IN, WORK( IWRK ), RWORK( IRWRK ),
     $                IERR )
         IF( IERR.NE.0 ) THEN
            INFO = N + 2
            GO TO 70
         END IF
*
*        Undo balancing on VL and VR and normalization
*        (Workspace: none needed)
*
         IF( ILVL ) THEN
            CALL ZGGBAK( 'P', 'L', N, ILO, IHI, RWORK( ILEFT ),
     $                   RWORK( IRIGHT ), N, VL, LDVL, IERR )
            DO 30 JC = 1, N
               TEMP = ZERO
               DO 10 JR = 1, N
                  TEMP = MAX( TEMP, ABS1( VL( JR, JC ) ) )
   10          CONTINUE
               IF( TEMP.LT.SMLNUM )
     $            GO TO 30
               TEMP = ONE / TEMP
               DO 20 JR = 1, N
                  VL( JR, JC ) = VL( JR, JC )*TEMP
   20          CONTINUE
   30       CONTINUE
         END IF
         IF( ILVR ) THEN
            CALL ZGGBAK( 'P', 'R', N, ILO, IHI, RWORK( ILEFT ),
     $                   RWORK( IRIGHT ), N, VR, LDVR, IERR )
            DO 60 JC = 1, N
               TEMP = ZERO
               DO 40 JR = 1, N
                  TEMP = MAX( TEMP, ABS1( VR( JR, JC ) ) )
   40          CONTINUE
               IF( TEMP.LT.SMLNUM )
     $            GO TO 60
               TEMP = ONE / TEMP
               DO 50 JR = 1, N
                  VR( JR, JC ) = VR( JR, JC )*TEMP
   50          CONTINUE
   60       CONTINUE
         END IF
      END IF
*
*     Undo scaling if necessary
*
      IF( ILASCL )
     $   CALL ZLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHA, N, IERR )
*
      IF( ILBSCL )
     $   CALL ZLASCL( 'G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR )
*
   70 CONTINUE
      WORK( 1 ) = LWKOPT
*
      RETURN
*
*     End of ZGGEV
*
      END


      SUBROUTINE ZGGEVV( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, B, LDB,
     $                   ALPHA, BETA, VL, LDVL, VR, LDVR, ILO, IHI,
     $                   LSCALE, RSCALE, ABNRM, BBNRM, RCONDE, RCONDV,
     $                   WORK, LWORK, RWORK, IWORK, BWORK, INFO )
*
*  -- LAPACK driver routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          BALANC, JOBVL, JOBVR, SENSE
      INTEGER            IHI, ILO, INFO, LDA, LDB, LDVL, LDVR, LWORK, N
      DOUBLE PRECISION   ABNRM, BBNRM
*     ..
*     .. Array Arguments ..
      LOGICAL            BWORK( * )
      INTEGER            IWORK( * )
      DOUBLE PRECISION   LSCALE( * ), RCONDE( * ), RCONDV( * ),
     $                   RSCALE( * ), RWORK( * )
      COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ),
     $                   BETA( * ), VL( LDVL, * ), VR( LDVR, * ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZGGEVX computes for a pair of N-by-N complex nonsymmetric matrices
*  (A,B) the generalized eigenvalues, and optionally, the left and/or
*  right generalized eigenvectors.
*
*  Optionally, it also computes a balancing transformation to improve
*  the conditioning of the eigenvalues and eigenvectors (ILO, IHI,
*  LSCALE, RSCALE, ABNRM, and BBNRM), reciprocal condition numbers for
*  the eigenvalues (RCONDE), and reciprocal condition numbers for the
*  right eigenvectors (RCONDV).
*
*  A generalized eigenvalue for a pair of matrices (A,B) is a scalar
*  lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
*  singular. It is usually represented as the pair (alpha,beta), as
*  there is a reasonable interpretation for beta=0, and even for both
*  being zero.
*
*  The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
*  of (A,B) satisfies
*                   A * v(j) = lambda(j) * B * v(j) .
*  The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
*  of (A,B) satisfies
*                   u(j)**H * A  = lambda(j) * u(j)**H * B.
*  where u(j)**H is the conjugate-transpose of u(j).
*
*
*  Arguments
*  =========
*
*  BALANC  (input) CHARACTER*1
*          Specifies the balance option to be performed:
*          = 'N':  do not diagonally scale or permute;
*          = 'P':  permute only;
*          = 'S':  scale only;
*          = 'B':  both permute and scale.
*          Computed reciprocal condition numbers will be for the
*          matrices after permuting and/or balancing. Permuting does
*          not change condition numbers (in exact arithmetic), but
*          balancing does.
*
*  JOBVL   (input) CHARACTER*1
*          = 'N':  do not compute the left generalized eigenvectors;
*          = 'V':  compute the left generalized eigenvectors.
*
*  JOBVR   (input) CHARACTER*1
*          = 'N':  do not compute the right generalized eigenvectors;
*          = 'V':  compute the right generalized eigenvectors.
*
*  SENSE   (input) CHARACTER*1
*          Determines which reciprocal condition numbers are computed.
*          = 'N': none are computed;
*          = 'E': computed for eigenvalues only;
*          = 'V': computed for eigenvectors only;
*          = 'B': computed for eigenvalues and eigenvectors.
*
*  N       (input) INTEGER
*          The order of the matrices A, B, VL, and VR.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA, N)
*          On entry, the matrix A in the pair (A,B).
*          On exit, A has been overwritten. If JOBVL='V' or JOBVR='V'
*          or both, then A contains the first part of the complex Schur
*          form of the "balanced" versions of the input A and B.
*
*  LDA     (input) INTEGER
*          The leading dimension of A.  LDA >= max(1,N).
*
*  B       (input/output) COMPLEX*16 array, dimension (LDB, N)
*          On entry, the matrix B in the pair (A,B).
*          On exit, B has been overwritten. If JOBVL='V' or JOBVR='V'
*          or both, then B contains the second part of the complex
*          Schur form of the "balanced" versions of the input A and B.
*
*  LDB     (input) INTEGER
*          The leading dimension of B.  LDB >= max(1,N).
*
*  ALPHA   (output) COMPLEX*16 array, dimension (N)
*  BETA    (output) COMPLEX*16 array, dimension (N)
*          On exit, ALPHA(j)/BETA(j), j=1,...,N, will be the generalized
*          eigenvalues.
*
*          Note: the quotient ALPHA(j)/BETA(j) ) may easily over- or
*          underflow, and BETA(j) may even be zero.  Thus, the user
*          should avoid naively computing the ratio ALPHA/BETA.
*          However, ALPHA will be always less than and usually
*          comparable with norm(A) in magnitude, and BETA always less
*          than and usually comparable with norm(B).
*
*  VL      (output) COMPLEX*16 array, dimension (LDVL,N)
*          If JOBVL = 'V', the left generalized eigenvectors u(j) are
*          stored one after another in the columns of VL, in the same
*          order as their eigenvalues.
*          Each eigenvector will be scaled so the largest component
*          will have abs(real part) + abs(imag. part) = 1.
*          Not referenced if JOBVL = 'N'.
*
*  LDVL    (input) INTEGER
*          The leading dimension of the matrix VL. LDVL >= 1, and
*          if JOBVL = 'V', LDVL >= N.
*
*  VR      (output) COMPLEX*16 array, dimension (LDVR,N)
*          If JOBVR = 'V', the right generalized eigenvectors v(j) are
*          stored one after another in the columns of VR, in the same
*          order as their eigenvalues.
*          Each eigenvector will be scaled so the largest component
*          will have abs(real part) + abs(imag. part) = 1.
*          Not referenced if JOBVR = 'N'.
*
*  LDVR    (input) INTEGER
*          The leading dimension of the matrix VR. LDVR >= 1, and
*          if JOBVR = 'V', LDVR >= N.
*
*  ILO     (output) INTEGER
*  IHI     (output) INTEGER
*          ILO and IHI are integer values such that on exit
*          A(i,j) = 0 and B(i,j) = 0 if i > j and
*          j = 1,...,ILO-1 or i = IHI+1,...,N.
*          If BALANC = 'N' or 'S', ILO = 1 and IHI = N.
*
*  LSCALE  (output) DOUBLE PRECISION array, dimension (N)
*          Details of the permutations and scaling factors applied
*          to the left side of A and B.  If PL(j) is the index of the
*          row interchanged with row j, and DL(j) is the scaling
*          factor applied to row j, then
*            LSCALE(j) = PL(j)  for j = 1,...,ILO-1
*                      = DL(j)  for j = ILO,...,IHI
*                      = PL(j)  for j = IHI+1,...,N.
*          The order in which the interchanges are made is N to IHI+1,
*          then 1 to ILO-1.
*
*  RSCALE  (output) DOUBLE PRECISION array, dimension (N)
*          Details of the permutations and scaling factors applied
*          to the right side of A and B.  If PR(j) is the index of the
*          column interchanged with column j, and DR(j) is the scaling
*          factor applied to column j, then
*            RSCALE(j) = PR(j)  for j = 1,...,ILO-1
*                      = DR(j)  for j = ILO,...,IHI
*                      = PR(j)  for j = IHI+1,...,N
*          The order in which the interchanges are made is N to IHI+1,
*          then 1 to ILO-1.
*
*  ABNRM   (output) DOUBLE PRECISION
*          The one-norm of the balanced matrix A.
*
*  BBNRM   (output) DOUBLE PRECISION
*          The one-norm of the balanced matrix B.
*
*  RCONDE  (output) DOUBLE PRECISION array, dimension (N)
*          If SENSE = 'E' or 'B', the reciprocal condition numbers of
*          the eigenvalues, stored in consecutive elements of the array.
*          If SENSE = 'N' or 'V', RCONDE is not referenced.
*
*  RCONDV  (output) DOUBLE PRECISION array, dimension (N)
*          If JOB = 'V' or 'B', the estimated reciprocal condition
*          numbers of the eigenvectors, stored in consecutive elements
*          of the array. If the eigenvalues cannot be reordered to
*          compute RCONDV(j), RCONDV(j) is set to 0; this can only occur
*          when the true value would be very small anyway.
*          If SENSE = 'N' or 'E', RCONDV is not referenced.
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= max(1,2*N).
*          If SENSE = 'E', LWORK >= max(1,4*N).
*          If SENSE = 'V' or 'B', LWORK >= max(1,2*N*N+2*N).
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  RWORK   (workspace) REAL array, dimension (lrwork)
*          lrwork must be at least max(1,6*N) if BALANC = 'S' or 'B',
*          and at least max(1,2*N) otherwise.
*          Real workspace.
*
*  IWORK   (workspace) INTEGER array, dimension (N+2)
*          If SENSE = 'E', IWORK is not referenced.
*
*  BWORK   (workspace) LOGICAL array, dimension (N)
*          If SENSE = 'N', BWORK is not referenced.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          = 1,...,N:
*                The QZ iteration failed.  No eigenvectors have been
*                calculated, but ALPHA(j) and BETA(j) should be correct
*                for j=INFO+1,...,N.
*          > N:  =N+1: other than QZ iteration failed in ZHGEQZ.
*                =N+2: error return from ZTGEVC.
*
*  Further Details
*  ===============
*
*  Balancing a matrix pair (A,B) includes, first, permuting rows and
*  columns to isolate eigenvalues, second, applying diagonal similarity
*  transformation to the rows and columns to make the rows and columns
*  as close in norm as possible. The computed reciprocal condition
*  numbers correspond to the balanced matrix. Permuting rows and columns
*  will not change the condition numbers (in exact arithmetic) but
*  diagonal scaling will.  For further explanation of balancing, see
*  section 4.11.1.2 of LAPACK Users' Guide.
*
*  An approximate error bound on the chordal distance between the i-th
*  computed generalized eigenvalue w and the corresponding exact
*  eigenvalue lambda is
*
*       chord(w, lambda) <= EPS * norm(ABNRM, BBNRM) / RCONDE(I)
*
*  An approximate error bound for the angle between the i-th computed
*  eigenvector VL(i) or VR(i) is given by
*
*       EPS * norm(ABNRM, BBNRM) / DIF(i).
*
*  For further explanation of the reciprocal condition numbers RCONDE
*  and RCONDV, see section 4.11 of LAPACK User's Guide.
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ),
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILASCL, ILBSCL, ILV, ILVL, ILVR, LQUERY, NOSCL,
     $                   WANTSB, WANTSE, WANTSN, WANTSV
      CHARACTER          CHTEMP
      INTEGER            I, ICOLS, IERR, IJOBVL, IJOBVR, IN, IROWS,
     $                   ITAU, IWRK, IWRK1, J, JC, JR, M, MAXWRK, MINWRK
      DOUBLE PRECISION   ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS,
     $                   SMLNUM, TEMP
      COMPLEX*16         X
*     ..
*     .. Local Arrays ..
      LOGICAL            LDUMMA( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLABAD, DLASCL, XERBLA, ZGEQRF, ZGGBAK, ZGGBAL,
     $                   ZGGHRD, ZHGERV, ZLACPY, ZLASCL, ZLASET, ZTGEVC,
     $                   ZTGSNA, ZUNGQR, ZUNMQR
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, ZLANGE
      EXTERNAL           LSAME, ILAENV, DLAMCH, ZLANGE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG, MAX, SQRT
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   ABS1
*     ..
*     .. Statement Function definitions ..
      ABS1( X ) = ABS( DBLE( X ) ) + ABS( DIMAG( X ) )
*     ..
*     .. Executable Statements ..
*
*     Decode the input arguments
*
      IF( LSAME( JOBVL, 'N' ) ) THEN
         IJOBVL = 1
         ILVL = .FALSE.
      ELSE IF( LSAME( JOBVL, 'V' ) ) THEN
         IJOBVL = 2
         ILVL = .TRUE.
      ELSE
         IJOBVL = -1
         ILVL = .FALSE.
      END IF
*
      IF( LSAME( JOBVR, 'N' ) ) THEN
         IJOBVR = 1
         ILVR = .FALSE.
      ELSE IF( LSAME( JOBVR, 'V' ) ) THEN
         IJOBVR = 2
         ILVR = .TRUE.
      ELSE
         IJOBVR = -1
         ILVR = .FALSE.
      END IF
      ILV = ILVL .OR. ILVR
*
      NOSCL  = LSAME( BALANC, 'N' ) .OR. LSAME( BALANC, 'P' )
      WANTSN = LSAME( SENSE, 'N' )
      WANTSE = LSAME( SENSE, 'E' )
      WANTSV = LSAME( SENSE, 'V' )
      WANTSB = LSAME( SENSE, 'B' )
*
*     Test the input arguments
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.( NOSCL .OR. LSAME( BALANC,'S' ) .OR.
     $    LSAME( BALANC, 'B' ) ) ) THEN
         INFO = -1
      ELSE IF( IJOBVL.LE.0 ) THEN
         INFO = -2
      ELSE IF( IJOBVR.LE.0 ) THEN
         INFO = -3
      ELSE IF( .NOT.( WANTSN .OR. WANTSE .OR. WANTSB .OR. WANTSV ) )
     $          THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDVL.LT.1 .OR. ( ILVL .AND. LDVL.LT.N ) ) THEN
         INFO = -13
      ELSE IF( LDVR.LT.1 .OR. ( ILVR .AND. LDVR.LT.N ) ) THEN
         INFO = -15
      END IF
*
*     Compute workspace
*      (Note: Comments in the code beginning "Workspace:" describe the
*       minimal amount of workspace needed at that point in the code,
*       as well as the preferred amount for good performance.
*       NB refers to the optimal block size for the immediately
*       following subroutine, as returned by ILAENV. The workspace is
*       computed assuming ILO = 1 and IHI = N, the worst case.)
*
      IF( INFO.EQ.0 ) THEN
         IF( N.EQ.0 ) THEN
            MINWRK = 1
            MAXWRK = 1
         ELSE
            MINWRK = 2*N
            IF( WANTSE ) THEN
               MINWRK = 4*N
            ELSE IF( WANTSV .OR. WANTSB ) THEN
               MINWRK = 2*N*( N + 1)
            END IF
            MAXWRK = MINWRK
            MAXWRK = MAX( MAXWRK,
     $                    N + N*ILAENV( 1, 'ZGEQRF', ' ', N, 1, N, 0 ) )
            MAXWRK = MAX( MAXWRK,
     $                    N + N*ILAENV( 1, 'ZUNMQR', ' ', N, 1, N, 0 ) )
            IF( ILVL ) THEN
               MAXWRK = MAX( MAXWRK, N +
     $                       N*ILAENV( 1, 'ZUNGQR', ' ', N, 1, N, 0 ) )
            END IF
         END IF
         WORK( 1 ) = MAXWRK
*
         IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
            INFO = -25
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGGEVX', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Get machine constants
*
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM
*
*     Scale A if max element outside range [SMLNUM,BIGNUM]
*
      ANRM = ZLANGE( 'M', N, N, A, LDA, RWORK )
      ILASCL = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         ANRMTO = SMLNUM
         ILASCL = .TRUE.
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         ANRMTO = BIGNUM
         ILASCL = .TRUE.
      END IF
      IF( ILASCL )
     $   CALL ZLASCL( 'G', 0, 0, ANRM, ANRMTO, N, N, A, LDA, IERR )
*
*     Scale B if max element outside range [SMLNUM,BIGNUM]
*
      BNRM = ZLANGE( 'M', N, N, B, LDB, RWORK )
      ILBSCL = .FALSE.
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
         BNRMTO = SMLNUM
         ILBSCL = .TRUE.
      ELSE IF( BNRM.GT.BIGNUM ) THEN
         BNRMTO = BIGNUM
         ILBSCL = .TRUE.
      END IF
      IF( ILBSCL )
     $   CALL ZLASCL( 'G', 0, 0, BNRM, BNRMTO, N, N, B, LDB, IERR )
*
*     Permute and/or balance the matrix pair (A,B)
*     (Real Workspace: need 6*N if BALANC = 'S' or 'B', 1 otherwise)
*
      CALL ZGGBAL( BALANC, N, A, LDA, B, LDB, ILO, IHI, LSCALE, RSCALE,
     $             RWORK, IERR )
*
*     Compute ABNRM and BBNRM
*
      ABNRM = ZLANGE( '1', N, N, A, LDA, RWORK( 1 ) )
      IF( ILASCL ) THEN
         RWORK( 1 ) = ABNRM
         CALL DLASCL( 'G', 0, 0, ANRMTO, ANRM, 1, 1, RWORK( 1 ), 1,
     $                IERR )
         ABNRM = RWORK( 1 )
      END IF
*
      BBNRM = ZLANGE( '1', N, N, B, LDB, RWORK( 1 ) )
      IF( ILBSCL ) THEN
         RWORK( 1 ) = BBNRM
         CALL DLASCL( 'G', 0, 0, BNRMTO, BNRM, 1, 1, RWORK( 1 ), 1,
     $                IERR )
         BBNRM = RWORK( 1 )
      END IF
*
*     Reduce B to triangular form (QR decomposition of B)
*     (Complex Workspace: need N, prefer N*NB )
*
      IROWS = IHI + 1 - ILO
      IF( ILV .OR. .NOT.WANTSN ) THEN
         ICOLS = N + 1 - ILO
      ELSE
         ICOLS = IROWS
      END IF
      ITAU = 1
      IWRK = ITAU + IROWS
      CALL ZGEQRF( IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ),
     $             WORK( IWRK ), LWORK+1-IWRK, IERR )
*
*     Apply the unitary transformation to A
*     (Complex Workspace: need N, prefer N*NB)
*
      CALL ZUNMQR( 'L', 'C', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB,
     $             WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ),
     $             LWORK+1-IWRK, IERR )
*
*     Initialize VL and/or VR
*     (Workspace: need N, prefer N*NB)
*
      IF( ILVL ) THEN
         CALL ZLASET( 'Full', N, N, CZERO, CONE, VL, LDVL )
         IF( IROWS.GT.1 ) THEN
            CALL ZLACPY( 'L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB,
     $                   VL( ILO+1, ILO ), LDVL )
         END IF
         CALL ZUNGQR( IROWS, IROWS, IROWS, VL( ILO, ILO ), LDVL,
     $                WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR )
      END IF
*
      IF( ILVR )
     $   CALL ZLASET( 'Full', N, N, CZERO, CONE, VR, LDVR )
*
*     Reduce to generalized Hessenberg form
*     (Workspace: none needed)
*
      IF( ILV .OR. .NOT.WANTSN ) THEN
*
*        Eigenvectors requested -- work on whole matrix.
*
         CALL ZGGHRD( JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, VL,
     $                LDVL, VR, LDVR, IERR )
      ELSE
         CALL ZGGHRD( 'N', 'N', IROWS, 1, IROWS, A( ILO, ILO ), LDA,
     $                B( ILO, ILO ), LDB, VL, LDVL, VR, LDVR, IERR )
      END IF
*
*     Perform QZ algorithm (Compute eigenvalues, and optionally, the
*     Schur forms and Schur vectors)
*     (Complex Workspace: need N)
*     (Real Workspace: need N)
*
      IWRK = ITAU
      IF( ILV .OR. .NOT.WANTSN ) THEN
         CHTEMP = 'S'
      ELSE
         CHTEMP = 'E'
      END IF
*
      CALL ZHGERV( CHTEMP, JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB,
     $             ALPHA, BETA, VL, LDVL, VR, LDVR, WORK( IWRK ),
     $             LWORK+1-IWRK, RWORK, IERR )
      IF( IERR.NE.0 ) THEN
         IF( IERR.GT.0 .AND. IERR.LE.N ) THEN
            INFO = IERR
         ELSE IF( IERR.GT.N .AND. IERR.LE.2*N ) THEN
            INFO = IERR - N
         ELSE
            INFO = N + 1
         END IF
         GO TO 90
      END IF
*
*     Compute Eigenvectors and estimate condition numbers if desired
*     ZTGEVC: (Complex Workspace: need 2*N )
*             (Real Workspace:    need 2*N )
*     ZTGSNA: (Complex Workspace: need 2*N*N if SENSE='V' or 'B')
*             (Integer Workspace: need N+2 )
*
      IF( ILV .OR. .NOT.WANTSN ) THEN
         IF( ILV ) THEN
            IF( ILVL ) THEN
               IF( ILVR ) THEN
                  CHTEMP = 'B'
               ELSE
                  CHTEMP = 'L'
               END IF
            ELSE
               CHTEMP = 'R'
            END IF
*
            CALL ZTGEVC( CHTEMP, 'B', LDUMMA, N, A, LDA, B, LDB, VL,
     $                   LDVL, VR, LDVR, N, IN, WORK( IWRK ), RWORK,
     $                   IERR )
            IF( IERR.NE.0 ) THEN
               INFO = N + 2
               GO TO 90
            END IF
         END IF
*
         IF( .NOT.WANTSN ) THEN
*
*           compute eigenvectors (DTGEVC) and estimate condition
*           numbers (DTGSNA). Note that the definition of the condition
*           number is not invariant under transformation (u,v) to
*           (Q*u, Z*v), where (u,v) are eigenvectors of the generalized
*           Schur form (S,T), Q and Z are orthogonal matrices. In order
*           to avoid using extra 2*N*N workspace, we have to
*           re-calculate eigenvectors and estimate the condition numbers
*           one at a time.
*
            DO 20 I = 1, N
*
               DO 10 J = 1, N
                  BWORK( J ) = .FALSE.
   10          CONTINUE
               BWORK( I ) = .TRUE.
*
               IWRK = N + 1
               IWRK1 = IWRK + N
*
               IF( WANTSE .OR. WANTSB ) THEN
                  CALL ZTGEVC( 'B', 'S', BWORK, N, A, LDA, B, LDB,
     $                         WORK( 1 ), N, WORK( IWRK ), N, 1, M,
     $                         WORK( IWRK1 ), RWORK, IERR )
                  IF( IERR.NE.0 ) THEN
                     INFO = N + 2
                     GO TO 90
                  END IF
               END IF
*
               CALL ZTGSNA( SENSE, 'S', BWORK, N, A, LDA, B, LDB,
     $                      WORK( 1 ), N, WORK( IWRK ), N, RCONDE( I ),
     $                      RCONDV( I ), 1, M, WORK( IWRK1 ),
     $                      LWORK-IWRK1+1, IWORK, IERR )
*
   20       CONTINUE
         END IF
      END IF
*
*     Undo balancing on VL and VR and normalization
*     (Workspace: none needed)
*
      IF( ILVL ) THEN
         CALL ZGGBAK( BALANC, 'L', N, ILO, IHI, LSCALE, RSCALE, N, VL,
     $                LDVL, IERR )
*
         DO 50 JC = 1, N
            TEMP = ZERO
            DO 30 JR = 1, N
               TEMP = MAX( TEMP, ABS1( VL( JR, JC ) ) )
   30       CONTINUE
            IF( TEMP.LT.SMLNUM )
     $         GO TO 50
            TEMP = ONE / TEMP
            DO 40 JR = 1, N
               VL( JR, JC ) = VL( JR, JC )*TEMP
   40       CONTINUE
   50    CONTINUE
      END IF
*
      IF( ILVR ) THEN
         CALL ZGGBAK( BALANC, 'R', N, ILO, IHI, LSCALE, RSCALE, N, VR,
     $                LDVR, IERR )
         DO 80 JC = 1, N
            TEMP = ZERO
            DO 60 JR = 1, N
               TEMP = MAX( TEMP, ABS1( VR( JR, JC ) ) )
   60       CONTINUE
            IF( TEMP.LT.SMLNUM )
     $         GO TO 80
            TEMP = ONE / TEMP
            DO 70 JR = 1, N
               VR( JR, JC ) = VR( JR, JC )*TEMP
   70       CONTINUE
   80    CONTINUE
      END IF
*
*     Undo scaling if necessary
*
      IF( ILASCL )
     $   CALL ZLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHA, N, IERR )
*
      IF( ILBSCL )
     $   CALL ZLASCL( 'G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR )
*
   90 CONTINUE
      WORK( 1 ) = MAXWRK
*
      RETURN
*
*     End of ZGGEVX
*
      END


      SUBROUTINE ZHGERV( JOB, COMPQ, COMPZ, N, ILO, IHI, H, LDH, T, LDT,
     $                   ALPHA, BETA, Q, LDQ, Z, LDZ, WORK, LWORK,
     $                   RWORK, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*
*  ROUTINE MODIFIED BY RV! SET ULP to 1E-15 TO ACHIEVE CONVERGENCE
*  FOR SIMULATIONS INLCUDING PMLS!
*
*     .. Scalar Arguments ..
      CHARACTER          COMPQ, COMPZ, JOB
      INTEGER            IHI, ILO, INFO, LDH, LDQ, LDT, LDZ, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         ALPHA( * ), BETA( * ), H( LDH, * ),
     $                   Q( LDQ, * ), T( LDT, * ), WORK( * ),
     $                   Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  ZHGEQZ computes the eigenvalues of a complex matrix pair (H,T),
*  where H is an upper Hessenberg matrix and T is upper triangular,
*  using the single-shift QZ method.
*  Matrix pairs of this type are produced by the reduction to
*  generalized upper Hessenberg form of a complex matrix pair (A,B):
*
*     A = Q1*H*Z1**H,  B = Q1*T*Z1**H,
*
*  as computed by ZGGHRD.
*
*  If JOB='S', then the Hessenberg-triangular pair (H,T) is
*  also reduced to generalized Schur form,
*
*     H = Q*S*Z**H,  T = Q*P*Z**H,
*
*  where Q and Z are unitary matrices and S and P are upper triangular.
*
*  Optionally, the unitary matrix Q from the generalized Schur
*  factorization may be postmultiplied into an input matrix Q1, and the
*  unitary matrix Z may be postmultiplied into an input matrix Z1.
*  If Q1 and Z1 are the unitary matrices from ZGGHRD that reduced
*  the matrix pair (A,B) to generalized Hessenberg form, then the output
*  matrices Q1*Q and Z1*Z are the unitary factors from the generalized
*  Schur factorization of (A,B):
*
*     A = (Q1*Q)*S*(Z1*Z)**H,  B = (Q1*Q)*P*(Z1*Z)**H.
*
*  To avoid overflow, eigenvalues of the matrix pair (H,T)
*  (equivalently, of (A,B)) are computed as a pair of complex values
*  (alpha,beta).  If beta is nonzero, lambda = alpha / beta is an
*  eigenvalue of the generalized nonsymmetric eigenvalue problem (GNEP)
*     A*x = lambda*B*x
*  and if alpha is nonzero, mu = beta / alpha is an eigenvalue of the
*  alternate form of the GNEP
*     mu*A*y = B*y.
*  The values of alpha and beta for the i-th eigenvalue can be read
*  directly from the generalized Schur form:  alpha = S(i,i),
*  beta = P(i,i).
*
*  Ref: C.B. Moler & G.W. Stewart, "An Algorithm for Generalized Matrix
*       Eigenvalue Problems", SIAM J. Numer. Anal., 10(1973),
*       pp. 241--256.
*
*  Arguments
*  =========
*
*  JOB     (input) CHARACTER*1
*          = 'E': Compute eigenvalues only;
*          = 'S': Computer eigenvalues and the Schur form.
*
*  COMPQ   (input) CHARACTER*1
*          = 'N': Left Schur vectors (Q) are not computed;
*          = 'I': Q is initialized to the unit matrix and the matrix Q
*                 of left Schur vectors of (H,T) is returned;
*          = 'V': Q must contain a unitary matrix Q1 on entry and
*                 the product Q1*Q is returned.
*
*  COMPZ   (input) CHARACTER*1
*          = 'N': Right Schur vectors (Z) are not computed;
*          = 'I': Q is initialized to the unit matrix and the matrix Z
*                 of right Schur vectors of (H,T) is returned;
*          = 'V': Z must contain a unitary matrix Z1 on entry and
*                 the product Z1*Z is returned.
*
*  N       (input) INTEGER
*          The order of the matrices H, T, Q, and Z.  N >= 0.
*
*  ILO     (input) INTEGER
*  IHI     (input) INTEGER
*          ILO and IHI mark the rows and columns of H which are in
*          Hessenberg form.  It is assumed that A is already upper
*          triangular in rows and columns 1:ILO-1 and IHI+1:N.
*          If N > 0, 1 <= ILO <= IHI <= N; if N = 0, ILO=1 and IHI=0.
*
*  H       (input/output) COMPLEX*16 array, dimension (LDH, N)
*          On entry, the N-by-N upper Hessenberg matrix H.
*          On exit, if JOB = 'S', H contains the upper triangular
*          matrix S from the generalized Schur factorization.
*          If JOB = 'E', the diagonal of H matches that of S, but
*          the rest of H is unspecified.
*
*  LDH     (input) INTEGER
*          The leading dimension of the array H.  LDH >= max( 1, N ).
*
*  T       (input/output) COMPLEX*16 array, dimension (LDT, N)
*          On entry, the N-by-N upper triangular matrix T.
*          On exit, if JOB = 'S', T contains the upper triangular
*          matrix P from the generalized Schur factorization.
*          If JOB = 'E', the diagonal of T matches that of P, but
*          the rest of T is unspecified.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T.  LDT >= max( 1, N ).
*
*  ALPHA   (output) COMPLEX*16 array, dimension (N)
*          The complex scalars alpha that define the eigenvalues of
*          GNEP.  ALPHA(i) = S(i,i) in the generalized Schur
*          factorization.
*
*  BETA    (output) COMPLEX*16 array, dimension (N)
*          The real non-negative scalars beta that define the
*          eigenvalues of GNEP.  BETA(i) = P(i,i) in the generalized
*          Schur factorization.
*
*          Together, the quantities alpha = ALPHA(j) and beta = BETA(j)
*          represent the j-th eigenvalue of the matrix pair (A,B), in
*          one of the forms lambda = alpha/beta or mu = beta/alpha.
*          Since either lambda or mu may overflow, they should not,
*          in general, be computed.
*
*  Q       (input/output) COMPLEX*16 array, dimension (LDQ, N)
*          On entry, if COMPZ = 'V', the unitary matrix Q1 used in the
*          reduction of (A,B) to generalized Hessenberg form.
*          On exit, if COMPZ = 'I', the unitary matrix of left Schur
*          vectors of (H,T), and if COMPZ = 'V', the unitary matrix of
*          left Schur vectors of (A,B).
*          Not referenced if COMPZ = 'N'.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.  LDQ >= 1.
*          If COMPQ='V' or 'I', then LDQ >= N.
*
*  Z       (input/output) COMPLEX*16 array, dimension (LDZ, N)
*          On entry, if COMPZ = 'V', the unitary matrix Z1 used in the
*          reduction of (A,B) to generalized Hessenberg form.
*          On exit, if COMPZ = 'I', the unitary matrix of right Schur
*          vectors of (H,T), and if COMPZ = 'V', the unitary matrix of
*          right Schur vectors of (A,B).
*          Not referenced if COMPZ = 'N'.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1.
*          If COMPZ='V' or 'I', then LDZ >= N.
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
*          On exit, if INFO >= 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,N).
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          = 1,...,N: the QZ iteration did not converge.  (H,T) is not
*                     in Schur form, but ALPHA(i) and BETA(i),
*                     i=INFO+1,...,N should be correct.
*          = N+1,...,2*N: the shift calculation failed.  (H,T) is not
*                     in Schur form, but ALPHA(i) and BETA(i),
*                     i=INFO-N+1,...,N should be correct.
*
*  Further Details
*  ===============
*
*  We assume that complex ABS works as long as its value is less than
*  overflow.
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ),
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILAZR2, ILAZRO, ILQ, ILSCHR, ILZ, LQUERY
      INTEGER            ICOMPQ, ICOMPZ, IFIRST, IFRSTM, IITER, ILAST,
     $                   ILASTM, IN, ISCHUR, ISTART, J, JC, JCH, JITER,
     $                   JR, MAXIT
      DOUBLE PRECISION   ABSB, ANORM, ASCALE, ATOL, BNORM, BSCALE, BTOL,
     $                   C, SAFMIN, TEMP, TEMP2, TEMPR, ULP
      COMPLEX*16         ABI22, AD11, AD12, AD21, AD22, CTEMP, CTEMP2,
     $                   CTEMP3, ESHIFT, RTDISC, S, SHIFT, SIGNBC, T1,
     $                   U12, X
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, ZLANHS
      EXTERNAL           LSAME, DLAMCH, ZLANHS
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARTG, ZLASET, ZROT, ZSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DCONJG, DIMAG, MAX, MIN,
     $                   SQRT
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   ABS1
*     ..
*     .. Statement Function definitions ..
      ABS1( X ) = ABS( DBLE( X ) ) + ABS( DIMAG( X ) )
*     ..
*     .. Executable Statements ..
*
*     Decode JOB, COMPQ, COMPZ
*
      IF( LSAME( JOB, 'E' ) ) THEN
         ILSCHR = .FALSE.
         ISCHUR = 1
      ELSE IF( LSAME( JOB, 'S' ) ) THEN
         ILSCHR = .TRUE.
         ISCHUR = 2
      ELSE
         ISCHUR = 0
      END IF
*
      IF( LSAME( COMPQ, 'N' ) ) THEN
         ILQ = .FALSE.
         ICOMPQ = 1
      ELSE IF( LSAME( COMPQ, 'V' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 2
      ELSE IF( LSAME( COMPQ, 'I' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 3
      ELSE
         ICOMPQ = 0
      END IF
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ILZ = .FALSE.
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 2
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 3
      ELSE
         ICOMPZ = 0
      END IF
*
*     Check Argument Values
*
      INFO = 0
      WORK( 1 ) = MAX( 1, N )
      LQUERY = ( LWORK.EQ.-1 )
      IF( ISCHUR.EQ.0 ) THEN
         INFO = -1
      ELSE IF( ICOMPQ.EQ.0 ) THEN
         INFO = -2
      ELSE IF( ICOMPZ.EQ.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( ILO.LT.1 ) THEN
         INFO = -5
      ELSE IF( IHI.GT.N .OR. IHI.LT.ILO-1 ) THEN
         INFO = -6
      ELSE IF( LDH.LT.N ) THEN
         INFO = -8
      ELSE IF( LDT.LT.N ) THEN
         INFO = -10
      ELSE IF( LDQ.LT.1 .OR. ( ILQ .AND. LDQ.LT.N ) ) THEN
         INFO = -14
      ELSE IF( LDZ.LT.1 .OR. ( ILZ .AND. LDZ.LT.N ) ) THEN
         INFO = -16
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -18
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHGERV', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
*     WORK( 1 ) = CMPLX( 1 )
      IF( N.LE.0 ) THEN
         WORK( 1 ) = DCMPLX( 1 )
         RETURN
      END IF
*
*     Initialize Q and Z
*
      IF( ICOMPQ.EQ.3 )
     $   CALL ZLASET( 'Full', N, N, CZERO, CONE, Q, LDQ )
      IF( ICOMPZ.EQ.3 )
     $   CALL ZLASET( 'Full', N, N, CZERO, CONE, Z, LDZ )
*
*     Machine Constants
*
      IN = IHI + 1 - ILO
      SAFMIN = DLAMCH( 'S' )
      ULP = DLAMCH( 'E' )*DLAMCH( 'B' )
*
*     MANUALLY OVERRIDE ULP, SET TO 1e-13
*
      ULP = DLAMCH( 'B' )*1e-13
      ANORM = ZLANHS( 'F', IN, H( ILO, ILO ), LDH, RWORK )
      BNORM = ZLANHS( 'F', IN, T( ILO, ILO ), LDT, RWORK )
      ATOL = MAX( SAFMIN, ULP*ANORM )
      BTOL = MAX( SAFMIN, ULP*BNORM )
      ASCALE = ONE / MAX( SAFMIN, ANORM )
      BSCALE = ONE / MAX( SAFMIN, BNORM )



*
*
*     Set Eigenvalues IHI+1:N
*
      DO 10 J = IHI + 1, N
         ABSB = ABS( T( J, J ) )
         IF( ABSB.GT.SAFMIN ) THEN
            SIGNBC = DCONJG( T( J, J ) / ABSB )
            T( J, J ) = ABSB
            IF( ILSCHR ) THEN
               CALL ZSCAL( J-1, SIGNBC, T( 1, J ), 1 )
               CALL ZSCAL( J, SIGNBC, H( 1, J ), 1 )
            ELSE
               H( J, J ) = H( J, J )*SIGNBC
            END IF
            IF( ILZ )
     $         CALL ZSCAL( N, SIGNBC, Z( 1, J ), 1 )
         ELSE
            T( J, J ) = CZERO
         END IF
         ALPHA( J ) = H( J, J )
         BETA( J ) = T( J, J )
   10 CONTINUE
*
*     If IHI < ILO, skip QZ steps
*
      IF( IHI.LT.ILO )
     $   GO TO 190
*
*     MAIN QZ ITERATION LOOP
*
*     Initialize dynamic indices
*
*     Eigenvalues ILAST+1:N have been found.
*        Column operations modify rows IFRSTM:whatever
*        Row operations modify columns whatever:ILASTM
*
*     If only eigenvalues are being computed, then
*        IFRSTM is the row of the last splitting row above row ILAST;
*        this is always at least ILO.
*     IITER counts iterations since the last eigenvalue was found,
*        to tell when to use an extraordinary shift.
*     MAXIT is the maximum number of QZ sweeps allowed.
*

      ILAST = IHI
      IF( ILSCHR ) THEN
         IFRSTM = 1
         ILASTM = N
      ELSE
         IFRSTM = ILO
         ILASTM = IHI
      END IF
      IITER = 0
      ESHIFT = CZERO
*      MAXIT = 30*( IHI-ILO+1 )
      MAXIT = 60*( IHI-ILO+1 )
*
      DO 170 JITER = 1, MAXIT
*
*        Check for too many iterations.
*
         IF( JITER.GT.MAXIT )
     $      GO TO 180

*
*        Split the matrix if possible.
*
*        Two tests:
*           1: H(j,j-1)=0  or  j=ILO
*           2: T(j,j)=0
*
*        Special case: j=ILAST
*
         IF( ILAST.EQ.ILO ) THEN
            GO TO 60
         ELSE
            IF( ABS1( H( ILAST, ILAST-1 ) ).LE.ATOL ) THEN
               H( ILAST, ILAST-1 ) = CZERO
               GO TO 60
            END IF
         END IF
*
         IF( ABS( T( ILAST, ILAST ) ).LE.BTOL ) THEN
            T( ILAST, ILAST ) = CZERO
            GO TO 50
         END IF
*
*        General case: j<ILAST
*
         DO 40 J = ILAST - 1, ILO, -1
*
*           Test 1: for H(j,j-1)=0 or j=ILO
*
            IF( J.EQ.ILO ) THEN
               ILAZRO = .TRUE.
            ELSE
               IF( ABS1( H( J, J-1 ) ).LE.ATOL ) THEN
                  H( J, J-1 ) = CZERO
                  ILAZRO = .TRUE.
               ELSE
                  ILAZRO = .FALSE.
               END IF
            END IF
*
*           Test 2: for T(j,j)=0
*
            IF( ABS( T( J, J ) ).LT.BTOL ) THEN
               T( J, J ) = CZERO
*
*              Test 1a: Check for 2 consecutive small subdiagonals in A
*
               ILAZR2 = .FALSE.
               IF( .NOT.ILAZRO ) THEN
                  IF( ABS1( H( J, J-1 ) )*( ASCALE*ABS1( H( J+1,
     $                J ) ) ).LE.ABS1( H( J, J ) )*( ASCALE*ATOL ) )
     $                ILAZR2 = .TRUE.
               END IF
*
*              If both tests pass (1 & 2), i.e., the leading diagonal
*              element of B in the block is zero, split a 1x1 block off
*              at the top. (I.e., at the J-th row/column) The leading
*              diagonal element of the remainder can also be zero, so
*              this may have to be done repeatedly.
*
               IF( ILAZRO .OR. ILAZR2 ) THEN
                  DO 20 JCH = J, ILAST - 1
                     CTEMP = H( JCH, JCH )
                     CALL ZLARTG( CTEMP, H( JCH+1, JCH ), C, S,
     $                            H( JCH, JCH ) )
                     H( JCH+1, JCH ) = CZERO
                     CALL ZROT( ILASTM-JCH, H( JCH, JCH+1 ), LDH,
     $                          H( JCH+1, JCH+1 ), LDH, C, S )
                     CALL ZROT( ILASTM-JCH, T( JCH, JCH+1 ), LDT,
     $                          T( JCH+1, JCH+1 ), LDT, C, S )
                     IF( ILQ )
     $                  CALL ZROT( N, Q( 1, JCH ), 1, Q( 1, JCH+1 ), 1,
     $                             C, DCONJG( S ) )
                     IF( ILAZR2 )
     $                  H( JCH, JCH-1 ) = H( JCH, JCH-1 )*C
                     ILAZR2 = .FALSE.
                     IF( ABS1( T( JCH+1, JCH+1 ) ).GE.BTOL ) THEN
                        IF( JCH+1.GE.ILAST ) THEN
                           GO TO 60
                        ELSE
                           IFIRST = JCH + 1
                           GO TO 70
                        END IF
                     END IF
                     T( JCH+1, JCH+1 ) = CZERO
   20             CONTINUE
                  GO TO 50
               ELSE
*
*                 Only test 2 passed -- chase the zero to T(ILAST,ILAST)
*                 Then process as in the case T(ILAST,ILAST)=0
*
                  DO 30 JCH = J, ILAST - 1
                     CTEMP = T( JCH, JCH+1 )
                     CALL ZLARTG( CTEMP, T( JCH+1, JCH+1 ), C, S,
     $                            T( JCH, JCH+1 ) )
                     T( JCH+1, JCH+1 ) = CZERO
                     IF( JCH.LT.ILASTM-1 )
     $                  CALL ZROT( ILASTM-JCH-1, T( JCH, JCH+2 ), LDT,
     $                             T( JCH+1, JCH+2 ), LDT, C, S )
                     CALL ZROT( ILASTM-JCH+2, H( JCH, JCH-1 ), LDH,
     $                          H( JCH+1, JCH-1 ), LDH, C, S )
                     IF( ILQ )
     $                  CALL ZROT( N, Q( 1, JCH ), 1, Q( 1, JCH+1 ), 1,
     $                             C, DCONJG( S ) )
                     CTEMP = H( JCH+1, JCH )
                     CALL ZLARTG( CTEMP, H( JCH+1, JCH-1 ), C, S,
     $                            H( JCH+1, JCH ) )
                     H( JCH+1, JCH-1 ) = CZERO
                     CALL ZROT( JCH+1-IFRSTM, H( IFRSTM, JCH ), 1,
     $                          H( IFRSTM, JCH-1 ), 1, C, S )
                     CALL ZROT( JCH-IFRSTM, T( IFRSTM, JCH ), 1,
     $                          T( IFRSTM, JCH-1 ), 1, C, S )
                     IF( ILZ )
     $                  CALL ZROT( N, Z( 1, JCH ), 1, Z( 1, JCH-1 ), 1,
     $                             C, S )
   30             CONTINUE
                  GO TO 50
               END IF
            ELSE IF( ILAZRO ) THEN
*
*              Only test 1 passed -- work on J:ILAST
*
               IFIRST = J
               GO TO 70
            END IF
*
*           Neither test passed -- try next J
*
   40    CONTINUE
*
*        (Drop-through is "impossible")
*
         INFO = 2*N + 1
         GO TO 210
*
*        T(ILAST,ILAST)=0 -- clear H(ILAST,ILAST-1) to split off a
*        1x1 block.
*
   50    CONTINUE
         CTEMP = H( ILAST, ILAST )
         CALL ZLARTG( CTEMP, H( ILAST, ILAST-1 ), C, S,
     $                H( ILAST, ILAST ) )
         H( ILAST, ILAST-1 ) = CZERO
         CALL ZROT( ILAST-IFRSTM, H( IFRSTM, ILAST ), 1,
     $              H( IFRSTM, ILAST-1 ), 1, C, S )
         CALL ZROT( ILAST-IFRSTM, T( IFRSTM, ILAST ), 1,
     $              T( IFRSTM, ILAST-1 ), 1, C, S )
         IF( ILZ )
     $      CALL ZROT( N, Z( 1, ILAST ), 1, Z( 1, ILAST-1 ), 1, C, S )
*
*        H(ILAST,ILAST-1)=0 -- Standardize B, set ALPHA and BETA
*
   60    CONTINUE
         ABSB = ABS( T( ILAST, ILAST ) )
         IF( ABSB.GT.SAFMIN ) THEN
            SIGNBC = DCONJG( T( ILAST, ILAST ) / ABSB )
            T( ILAST, ILAST ) = ABSB
            IF( ILSCHR ) THEN
               CALL ZSCAL( ILAST-IFRSTM, SIGNBC, T( IFRSTM, ILAST ), 1 )
               CALL ZSCAL( ILAST+1-IFRSTM, SIGNBC, H( IFRSTM, ILAST ),
     $                     1 )
            ELSE
               H( ILAST, ILAST ) = H( ILAST, ILAST )*SIGNBC
            END IF
            IF( ILZ )
     $         CALL ZSCAL( N, SIGNBC, Z( 1, ILAST ), 1 )
         ELSE
            T( ILAST, ILAST ) = CZERO
         END IF
         ALPHA( ILAST ) = H( ILAST, ILAST )
         BETA( ILAST ) = T( ILAST, ILAST )
*
*        Go to next block -- exit if finished.
*
         ILAST = ILAST - 1
         IF( ILAST.LT.ILO )
     $      GO TO 190
*
*        Reset counters
*
         IITER = 0
         ESHIFT = CZERO
         IF( .NOT.ILSCHR ) THEN
            ILASTM = ILAST
            IF( IFRSTM.GT.ILAST )
     $         IFRSTM = ILO
         END IF
         GO TO 160
*
*        QZ step
*
*        This iteration only involves rows/columns IFIRST:ILAST.  We
*        assume IFIRST < ILAST, and that the diagonal of B is non-zero.
*
   70    CONTINUE
         IITER = IITER + 1
         IF( .NOT.ILSCHR ) THEN
            IFRSTM = IFIRST
         END IF
*
*        Compute the Shift.
*
*        At this point, IFIRST < ILAST, and the diagonal elements of
*        T(IFIRST:ILAST,IFIRST,ILAST) are larger than BTOL (in
*        magnitude)
*
         IF( ( IITER / 10 )*10.NE.IITER ) THEN
*
*           The Wilkinson shift (AEP p.512), i.e., the eigenvalue of
*           the bottom-right 2x2 block of A inv(B) which is nearest to
*           the bottom-right element.
*
*           We factor B as U*D, where U has unit diagonals, and
*           compute (A*inv(D))*inv(U).
*
            U12 = ( BSCALE*T( ILAST-1, ILAST ) ) /
     $            ( BSCALE*T( ILAST, ILAST ) )
            AD11 = ( ASCALE*H( ILAST-1, ILAST-1 ) ) /
     $             ( BSCALE*T( ILAST-1, ILAST-1 ) )
            AD21 = ( ASCALE*H( ILAST, ILAST-1 ) ) /
     $             ( BSCALE*T( ILAST-1, ILAST-1 ) )
            AD12 = ( ASCALE*H( ILAST-1, ILAST ) ) /
     $             ( BSCALE*T( ILAST, ILAST ) )
            AD22 = ( ASCALE*H( ILAST, ILAST ) ) /
     $             ( BSCALE*T( ILAST, ILAST ) )
            ABI22 = AD22 - U12*AD21
*
            T1 = HALF*( AD11+ABI22 )
            RTDISC = SQRT( T1**2+AD12*AD21-AD11*AD22 )
            TEMP = DBLE( T1-ABI22 )*DBLE( RTDISC ) +
     $             DIMAG( T1-ABI22 )*DIMAG( RTDISC )
            IF( TEMP.LE.ZERO ) THEN
               SHIFT = T1 + RTDISC
            ELSE
               SHIFT = T1 - RTDISC
            END IF
         ELSE
*
*           Exceptional shift.  Chosen for no particularly good reason.
*
            ESHIFT = ESHIFT + DCONJG( ( ASCALE*H( ILAST-1, ILAST ) ) /
     $               ( BSCALE*T( ILAST-1, ILAST-1 ) ) )
            SHIFT = ESHIFT
         END IF
*
*        Now check for two consecutive small subdiagonals.
*
         DO 80 J = ILAST - 1, IFIRST + 1, -1
            ISTART = J
            CTEMP = ASCALE*H( J, J ) - SHIFT*( BSCALE*T( J, J ) )
            TEMP = ABS1( CTEMP )
            TEMP2 = ASCALE*ABS1( H( J+1, J ) )
            TEMPR = MAX( TEMP, TEMP2 )
            IF( TEMPR.LT.ONE .AND. TEMPR.NE.ZERO ) THEN
               TEMP = TEMP / TEMPR
               TEMP2 = TEMP2 / TEMPR
            END IF
            IF( ABS1( H( J, J-1 ) )*TEMP2.LE.TEMP*ATOL )
     $         GO TO 90
   80    CONTINUE
*
         ISTART = IFIRST
         CTEMP = ASCALE*H( IFIRST, IFIRST ) -
     $           SHIFT*( BSCALE*T( IFIRST, IFIRST ) )
   90    CONTINUE
*
*        Do an implicit-shift QZ sweep.
*
*        Initial Q
*
         CTEMP2 = ASCALE*H( ISTART+1, ISTART )
         CALL ZLARTG( CTEMP, CTEMP2, C, S, CTEMP3 )
*
*        Sweep
*
         DO 150 J = ISTART, ILAST - 1
            IF( J.GT.ISTART ) THEN
               CTEMP = H( J, J-1 )
               CALL ZLARTG( CTEMP, H( J+1, J-1 ), C, S, H( J, J-1 ) )
               H( J+1, J-1 ) = CZERO
            END IF
*
            DO 100 JC = J, ILASTM
               CTEMP = C*H( J, JC ) + S*H( J+1, JC )
               H( J+1, JC ) = -DCONJG( S )*H( J, JC ) + C*H( J+1, JC )
               H( J, JC ) = CTEMP
               CTEMP2 = C*T( J, JC ) + S*T( J+1, JC )
               T( J+1, JC ) = -DCONJG( S )*T( J, JC ) + C*T( J+1, JC )
               T( J, JC ) = CTEMP2
  100       CONTINUE
            IF( ILQ ) THEN
               DO 110 JR = 1, N
                  CTEMP = C*Q( JR, J ) + DCONJG( S )*Q( JR, J+1 )
                  Q( JR, J+1 ) = -S*Q( JR, J ) + C*Q( JR, J+1 )
                  Q( JR, J ) = CTEMP
  110          CONTINUE
            END IF
*
            CTEMP = T( J+1, J+1 )
            CALL ZLARTG( CTEMP, T( J+1, J ), C, S, T( J+1, J+1 ) )
            T( J+1, J ) = CZERO
*
            DO 120 JR = IFRSTM, MIN( J+2, ILAST )
               CTEMP = C*H( JR, J+1 ) + S*H( JR, J )
               H( JR, J ) = -DCONJG( S )*H( JR, J+1 ) + C*H( JR, J )
               H( JR, J+1 ) = CTEMP
  120       CONTINUE
            DO 130 JR = IFRSTM, J
               CTEMP = C*T( JR, J+1 ) + S*T( JR, J )
               T( JR, J ) = -DCONJG( S )*T( JR, J+1 ) + C*T( JR, J )
               T( JR, J+1 ) = CTEMP
  130       CONTINUE
            IF( ILZ ) THEN
               DO 140 JR = 1, N
                  CTEMP = C*Z( JR, J+1 ) + S*Z( JR, J )
                  Z( JR, J ) = -DCONJG( S )*Z( JR, J+1 ) + C*Z( JR, J )
                  Z( JR, J+1 ) = CTEMP
  140          CONTINUE
            END IF
  150    CONTINUE
*
  160    CONTINUE
*
  170 CONTINUE

*
*     Drop-through = non-convergence
*
  180 CONTINUE
      INFO = ILAST
      GO TO 210
*
*     Successful completion of all QZ steps
*
  190 CONTINUE
*
*     Set Eigenvalues 1:ILO-1
*
      DO 200 J = 1, ILO - 1
         ABSB = ABS( T( J, J ) )
         IF( ABSB.GT.SAFMIN ) THEN
            SIGNBC = DCONJG( T( J, J ) / ABSB )
            T( J, J ) = ABSB
            IF( ILSCHR ) THEN
               CALL ZSCAL( J-1, SIGNBC, T( 1, J ), 1 )
               CALL ZSCAL( J, SIGNBC, H( 1, J ), 1 )
            ELSE
               H( J, J ) = H( J, J )*SIGNBC
            END IF
            IF( ILZ )
     $         CALL ZSCAL( N, SIGNBC, Z( 1, J ), 1 )
         ELSE
            T( J, J ) = CZERO
         END IF
         ALPHA( J ) = H( J, J )
         BETA( J ) = T( J, J )
  200 CONTINUE
*
*     Normal Termination
*
      INFO = 0
*
*     Exit (other than argument error) -- return optimal workspace size
*
  210 CONTINUE
      WORK( 1 ) = DCMPLX( N )
      RETURN
*
*     End of ZHGERV
*
      END
