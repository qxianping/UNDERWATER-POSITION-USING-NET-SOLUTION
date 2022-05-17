       PROGRAM NETSOL                                                   
C                                                                                                                                 
C...PURPOSE:  TO COMPUTE TRAS POSITION USING NET SOLUTION                                                                                          
       USE DFPORT
	 IMPLICIT REAL *8 (A-H,O-Z)
	 CHARACTER*166  EXEPATH,INPPATH
	 CHARACTER*166  FILE_POINT,FILE_DIST,FILE_SVP
	 CHARACTER*166  FILE_TEMP,FILE_RESULT
	 CHARACTER*166  FIL_TMP,FIL_NAM

       CHARACTER*180 BUFF180                                            
       INTEGER(4) IYEAR,IMN,IDAY,IHR,IMIN,NUM_SVP,NUM_DIS,KFILE
       INTEGER(4) NUM_DIS_USE
       DIMENSION XYZ(4,5),SVP_DAT(3,50000),SOUND_VEL(6)
       DIMENSION DIS_DAT(3,6000),DIS_RANG(3,6000),DL(6000)
       DIMENSION BLH(7,4),X_CHAMP(14),DN0(14)
       DIMENSION POS0(3),POS1(3),DNEU0(4,5),DNEU1(4,5),TEMP(3)
       DIMENSION DNEUW(4,5)

       PI=3.141592653589793D0
C.....
       X_CHAMP(1:14) = 0.0D0
       SOUND_VEL(1:6) = 15080
 	CALL GETARG(0,EXEPATH)
      CALL GETARG(1,INPPATH)
C.....
      OPEN (15,FILE='D:\OCEAN\NETSOL\FILENET.INP',STATUS = 'OLD')
      WRITE (FILE_POINT,'(256X)')
      WRITE (FILE_DIST,'(256X)')
      WRITE (FILE_SVP,'(256X)')
      WRITE (FILE_TEMP,'(256X)')
      WRITE (FILE_RESULT,'(256X)')

!C-------------------------------------------------------------
      CALL READ_XYZ(LFNAM,XYZ)
      READ (15,'(A)') FILE_DIST
      READ (15,'(A)') FILE_SVP
      READ (15,'(A)') FILE_TEMP
      READ (15,'(A)') FILE_RESULT
      CLOSE (15)

      OPEN (11,FILE=FILE_DIST,STATUS = 'OLD',ERR=999)
      OPEN (12,FILE=FILE_SVP,STATUS = 'OLD',ERR=999)
      OPEN (13,FILE=FILE_TEMP,STATUS = 'UNKNOWN')
      OPEN (14,FILE=FILE_RESULT,STATUS = 'UNKNOWN')

!C-------------------------------------------------------------
      CALL WRT_NEU(14,XYZ,DNEU0)
      DNEUW(1:4,1:5) = DNEU0(1:4,1:5)
C-------------------------------------------------------------
       CALL READ_DIST(11,NUM_DIS,DIS_DAT)
       CALL READ_SVP(12,NUM_SVP,SVP_DAT)
       CALL GET_DIST(XYZ,NUM_DIS,DIS_DAT,NUM_SVP,SVP_DAT
     1               ,DIS_RANG,SOUND_VEL)
C#############################################################

       NTER = 10
       ITER = 0
       DELT_RR = 1.0D5
       IDELT_R = 0
       DELT_R = 1.0D5
          NPAR = 14
          X_CHAMP(1:14) =  0.0D0
          DO IK = 1,4
          X_CHAMP((IK-1)*3+1) = DNEU0(IK,1)
          X_CHAMP((IK-1)*3+2) = DNEU0(IK,2)
          X_CHAMP((IK-1)*3+3) = DNEU0(IK,3)
          ENDDO

200   CONTINUE
!
!!!!C//////////2D ADJUSTMENT///////////////////////
!        WRITE(14,'(A40)') '==========4: 2D ADJUSTMENT============='
!       CALL NETADJUST(4,DNEU0,NUM_DIS,DIS_DAT,DIS_RANG,X_CHAMP,
!     1                DL,DN0,NPAR,SIGMA,DARC_INT,SOUND_VEL)
!         DT_CS = DABS(X_CHAMP(14))/(SIGMA*DSQRT(DN0(14)))
!        WRITE(14,'(F12.3,12F15.3)') SIGMA,X_CHAMP(13),X_CHAMP(14)
!     2         ,DT_CS
!        WRITE(14,*)
!C//////////3D ADJUSTMENT///////////////////////
        WRITE(14,'(A40)') '==========5: 3D ADJUSTMENT=============='
       CALL NETADJUST(5,DNEU0,NUM_DIS,DIS_DAT,DIS_RANG,X_CHAMP,
     1                DL,DN0,NPAR,SIGMA,DARC_INT,SOUND_VEL)
         DT_CS = DABS(X_CHAMP(14))/(SIGMA*DSQRT(DN0(14)))
        WRITE(14,'(F12.3,10F15.3)') SIGMA,X_CHAMP(13),X_CHAMP(14)
     2         ,DT_CS,SIGMA*DSQRT(DN0(13)),SIGMA*DSQRT(DN0(14))
        WRITE(14,*)
!
!C-------------------------------------------------------------
300   CONTINUE
          ITER =   ITER + 1
          DO IK = 1,4
          DNEU1(IK,1) = X_CHAMP((IK-1)*3+1)
          DNEU1(IK,2) = X_CHAMP((IK-1)*3+2) 
          DNEU1(IK,3) = X_CHAMP((IK-1)*3+3) 
          ENDDO
!C...
        WRITE(14,'(A60)') '===========NEU:AFTER ADJUSTMENT==========='
          DO I = 1,4
          WRITE(14,'(8X,3(F15.3,2X))') DNEU1(I,1:3)
          ENDDO
!        WRITE(14,*)
!C.....OUTPUT THE DIFFERENCE OF COORDINATE (AFTER-BEFORE ADJUSTMENT)
!C...
        WRITE(14,'(A40)') '=====THE DIFFERENCE OF COORDINATE========='
          DO MK = 1,4
         IF(IDELT_R.NE.1) WRITE(14,'(I4,F12.3,6F15.3)') ITER,SIGMA,
     1                DNEU1(MK,1:3) - DNEU0(MK,1:3),
     1                DNEU1(MK,1:3) - DNEUW(MK,1:3)
       
        ENDDO
        DELT_R = 0.0D0
        DO IK =1,4
        DELT_R = DELT_R + DSQRT((DNEU0(IK,1)-DNEU1(IK,1))**2
     1                +(DNEU0(IK,2)-DNEU1(IK,2))**2
     1                +(DNEU0(IK,3)-DNEU1(IK,3))**2)
        ENDDO

         DELT_R2 = DABS(DELT_R-DELT_RR)
         DELT_RR = DELT_R

          DO IK = 1,4
          DNEU0(IK,1) =  DNEU1(IK,1)
          DNEU0(IK,2) =  DNEU1(IK,2)
          DNEU0(IK,3) =  DNEU1(IK,3) 
          ENDDO
         IF(ITER.GE.NTER)GOTO 400
         IF(DELT_R.LT.1.0D-1)GOTO 400
         IF(DELT_R2.LT.5.0D-2)GOTO 400


          GOTO 200
C#############################################################
 400   CONTINUE
!
      WRITE(14,'(A40)') '=============LINE:SV/RANGE/DH=============='
          DO IMK = 16,21
                IF(IMK.EQ.16) IDSAT = 102
                IF(IMK.EQ.17) IDSAT = 103
                IF(IMK.EQ.18) IDSAT = 104
                IF(IMK.EQ.19) IDSAT = 203
                IF(IMK.EQ.20) IDSAT = 204
                IF(IMK.EQ.21) IDSAT = 304

               IK1 = INT(IDSAT/100)
               IK2 = DMOD(DBLE(IDSAT),100.0D0)
               DX1 =    X_CHAMP((IK1-1)*3+1) 
     1                - X_CHAMP((IK2-1)*3+1)
               DY1 =    X_CHAMP((IK1-1)*3+2) 
     1                - X_CHAMP((IK2-1)*3+2)
               DZ1 =    X_CHAMP((IK1-1)*3+3) 
     1                - X_CHAMP((IK2-1)*3+3)
	         DR = DSQRT(DX1**2+DY1**2+DZ1**2)
	         DDH = DZ1 

          WRITE(14,'(I4,F12.3,I6,6F15.3)') ITER,SIGMA,IDSAT,
     1        SOUND_VEL(IMK-15),X_CHAMP(13),DR,DDH
        ENDDO
        WRITE(14,*)
!!   
!
      WRITE(14,'(A40)') '==LINE/ALL_NUM./USED_NUM./MEAN/STD/RMS=='
!C.....OUTPUT THE SERIES OF RESIDUAL(ONE LINE ONE FILE)
      DO I = 1, 3
           DO J = I+1,4
          NUM_L = 0
          DLSUM = 0.0D0
          DLSUM2 = 0.0D0
          KFILE = 20 + 1
          FIL_TMP = 'D:\OCEAN\NETSOL\DAT\RESIDUAL.'
          WRITE(FIL_NAM,'(A30,2I2.2)')FIL_TMP,I,J
          FIL_NAM = TRIM(FIL_NAM)
          OPEN(KFILE,FILE=FIL_NAM,STATUS = 'UNKNOWN')

 	     DO IOBS = 1,NUM_DIS
               IK1 = INT(DIS_RANG(2,IOBS)/100)
               IK2 = DMOD(DIS_RANG(2,IOBS),100.0D0)
      IF(IK1.EQ.I.AND.IK2.EQ.J.AND.DABS(DL(IOBS)).LT.3.0*SIGMA) THEN
          NUM_L = NUM_L +1
          DLSUM = DLSUM + DL(IOBS)
          DLSUM2 = DLSUM2 + DL(IOBS)*DL(IOBS)

         WRITE(KFILE,'(F12.3,2F9.3)') DIS_RANG(1,IOBS)/3600.0D0,
     1                                 DIS_RANG(2,IOBS),DL(IOBS)
       ENDIF
      ENDDO 

C....OUTPUT THE STATISTICS OF RESIDUAL
          DMEAN = DLSUM/NUM_L
          DRMS = DSQRT(DLSUM2/NUM_L)

          DLSUM2 = 0.0D0
          NUM_IJ = 0
          NUM_IJ_USE = 0
 	   DO IOBS = 1,NUM_DIS
               IK1 = INT(DIS_RANG(2,IOBS)/100)
               IK2 = DMOD(DIS_RANG(2,IOBS),100.0D0)
          IF(IK1.EQ.I.AND.IK2.EQ.J) THEN
               NUM_IJ = NUM_IJ +1
            IF(DABS(DL(IOBS)).LT.3.0*SIGMA) THEN
              NUM_IJ_USE = NUM_IJ_USE + 1
              DLSUM2 = DLSUM2 + (DL(IOBS)-DMEAN)*(DL(IOBS)-DMEAN)
           ENDIF
          ENDIF
        ENDDO 
        DSTD = DSQRT(DLSUM2/NUM_L)

        WRITE(14,'(3X,I1,I2.2,2I6,5F15.3)') I,J,NUM_IJ,NUM_IJ_USE,
     1                                      DMEAN,DSTD,DRMS

      ENDDO
      ENDDO
        WRITE(14,'(A40)') '========================================'

C.....OUTPUT THE SERIES OF RESIDUAL
         NUM_OUTLIN = 0
 	   DO IOBS = 1,NUM_DIS
               IK1 = INT(DIS_RANG(2,IOBS)/100)
               IK2 = DMOD(DIS_RANG(2,IOBS),100.0D0)
          IF(IK1.LT.5.AND.IK2.LT.5.AND.DABS(DL(IOBS)).LT.3.0*SIGMA) THEN
          WRITE(13,'(F12.3,2F9.3)') DIS_RANG(1,IOBS)/3600.0D0,
     1                              DIS_RANG(2,IOBS),DL(IOBS)
          ENDIF
        ENDDO 


999    CONTINUE
      END




