      SUBROUTINE GET_DIST(XYZ,NUM_DIS,DIS_DAT,NUM_SVP,SVP_DAT,
     1                     DIS_RANG,SOUND_VEL)
      IMPLICIT REAL*8(A-H,O-Z)

       INTEGER(4) NUM_SVP,NUM_DIS,KFILE
       DIMENSION XYZ(4,5),SVP_DAT(3,50000),SVP_H(5000),SVP_V(5000)
       DIMENSION DIS_DAT(3,6000),DIS_RANG(3,6000)
       DIMENSION DMIN(6),DMAX(6),SOUND_VEL(6)
       DIMENSION DSVP(100),DH(100)
	 CHARACTER*256  FIL_TMP,FIL_NAM

       DMAXH = 0.0D0
       DMINH = 1.0D5
       SVP_H(1:500) = 0.0D0
       SVP_V(1:500) = 0.0D0
       DO I  = 1,3
        DO J  = 1,6000
           DIS_RANG(I,J) = 0.0D0
        ENDDO
       ENDDO

       NUM = 0
!
!C.....DETERMINE THE RANGE OF  HEIGHT DIFFERENCE OF EACH BASELINE
      DO I = 1, 3
      DO J = I+1,4
          NUM = NUM + 1
          IF(DABS(XYZ(I,4)).GT.DABS(XYZ(J,4))) THEN
                DMAX(NUM) =  DABS(XYZ(I,4)) + XYZ(I,5)
                DMIN(NUM) =  DABS(XYZ(J,4)) + XYZ(J,5)
             ELSE
                DMAX(NUM) =  DABS(XYZ(J,4)) + XYZ(J,5)
                DMIN(NUM) =  DABS(XYZ(I,4)) + XYZ(I,5)
          ENDIF

        ENDDO
       ENDDO
      DMAXH = MAX(DMAX(1),DMAX(2),DMAX(3),DMAX(4),DMAX(5)
     1           ,DMAX(6))
      DMINH = MIN(DMIN(1),DMIN(2),DMIN(3),DMIN(4),DMIN(5)
     1            ,DMIN(6))
C============================================================
C.....COMPUTE THE LINEAR FITTING COEFFICIENTS OF THE SOUND SPEED VARYING
        MM = 0

      DO J = 1,NUM_SVP
      IF(SVP_DAT(2,J).GE.DMINH-10.0D0
     1           .AND.SVP_DAT(2,J).LT.DMAXH+10.0D0) THEN
        MM = MM + 1
       SVP_H(MM) = SVP_DAT(2,J)
       SVP_V(MM) = SVP_DAT(3,J)
      ENDIF
      ENDDO
      CALL FIT(SVP_H,SVP_V,MM,A,B,SIGA,SIGB,CHI2)
C
C.....COMPUTE THE SOUND SPEED OF EACH BASELINE
      DO I = 1,6
         DHW = (DMIN(I)+DMAX(I))/2.0D0
         SOUND_VEL(I) = A + B*DHW
      ENDDO
C*******************************************************************
C.....COMPUTE RANGE OF EACH BASELINE
       DO I = 1,NUM_DIS
        IF(DIS_DAT(2,I).EQ.102.0D0) KK = 1
        IF(DIS_DAT(2,I).EQ.103.0D0) KK = 2
        IF(DIS_DAT(2,I).EQ.104.0D0) KK = 3
        IF(DIS_DAT(2,I).EQ.203.0D0) KK = 4
        IF(DIS_DAT(2,I).EQ.204.0D0) KK = 5
        IF(DIS_DAT(2,I).EQ.304.0D0) KK = 6
           DIS_RANG(1,I)= DIS_DAT(1,I)
           DIS_RANG(2,I)= DIS_DAT(2,I)
           DIS_RANG(3,I)= DIS_DAT(3,I)*SOUND_VEL(KK)
        ENDDO

        DC_SUM = 0.0D0


      DO I = 1, 3
          DO J = I+1,4

          KFILE = 20 + 1
          FIL_TMP = 'D:\OCEAN\NETSOL\DAT\0RESIDUAL.'
          WRITE(FIL_NAM,'(A31,2I2.2)')FIL_TMP,I,J
          FIL_NAM = TRIM(FIL_NAM)
          OPEN(KFILE,FILE=FIL_NAM,STATUS = 'UNKNOWN')

 	     DO IPRN = 1,NUM_DIS
               IK1 = INT(DIS_RANG(2,IPRN)/100)
               IK2 = DMOD(DIS_RANG(2,IPRN),100.0D0)
               DX1 =    XYZ(IK1,1) - XYZ(IK2,1)
               DY1 =    XYZ(IK1,2) - XYZ(IK2,2)
               DZ1 =    XYZ(IK1,3) - XYZ(IK2,3)
	         DR  = DSQRT(DX1**2+DY1**2+DZ1**2)
               DDR = DIS_RANG(3,IPRN)- DR

      IF(IK1.EQ.I.AND.IK2.EQ.J.OR.IK1.EQ.J.AND.IK2.EQ.I)
     1        WRITE(KFILE,'(4F12.3)')DIS_RANG(1:3,IPRN),DDR
      ENDDO
      ENDDO
      ENDDO

999   CONTINUE
        RETURN
	END
