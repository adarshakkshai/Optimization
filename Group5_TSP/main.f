C=================================================================
C BEGINNING OF FILE MAIN.FOR
C=================================================================
      PROGRAM MAIN
C
C     SAMPLE CALLING FOR SUBROUTINE CDT
C
C
C     .. Arrays in Common ..
      INTEGER X(10000)
C     ..
C     .. Local Scalars ..
      REAL ALPHA,AVSON,SPARS
      INTEGER ACTIVE,ERR,I,INF,LB0,LBC,LOPT,MAXND,N,NASS,NEXP,NPROBQ,
     +        ORDX,ZEUR,ZSTAR
C     ..
C     .. Local Arrays ..
      INTEGER FSTAR(100)
C     ..
C     .. External Subroutines ..
      EXTERNAL CDT,READA
C     ..
C     .. Common blocks ..
      COMMON X
C     ..
      ORDX = 10000
      INF = 99999999
      ALPHA = -1.
      MAXND = -1
      ZEUR = -1
C     READ N (WITH N .LE. 500) AND THE COST MATRIX
      OPEN (UNIT=1,FILE='INP.DAT',STATUS='OLD')
      READ (1,FMT=*) N
      CALL READA(N,X(1))
      CLOSE (1)
C
      CALL CDT(N,ORDX,X,MAXND,INF,ALPHA,ZEUR,ZSTAR,FSTAR,LB0,LBC,NEXP,
     +         NPROBQ,NASS,ACTIVE,LOPT,SPARS,AVSON,ERR)
C
      IF (ERR.EQ.0) GO TO 10
      WRITE (6,FMT='('' SOLUTION NOT OPTIMAL '')')
   10 WRITE (6,FMT=
     +'('' ZSTAR='',I8,'' LB0     ='',I8,'' LBC  ='',I8,      '' SPARS =
     +'',F8.4)') ZSTAR,LB0,LBC,SPARS
      WRITE (6,FMT=
     +'('' N.EXP='',I8,'' N.PROB.Q='',I8,'' N.ASS='',I8,     '' ACTIVE='
     +',I8,'' AV.SON='',F8.2)') NEXP,NPROBQ,NASS,ACTIVE,AVSON
      WRITE (6,FMT='(20I4)') (FSTAR(I),I=1,N)
      STOP
      END
C
      SUBROUTINE READA(N,A)
C
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      INTEGER A(N,N)
C     ..
C     .. Local Scalars ..
      INTEGER I,J
C     ..
      DO 10 I = 1,N
      READ (1,FMT=*) (A(I,J),J=1,N)
10    CONTINUE
      RETURN
      END
