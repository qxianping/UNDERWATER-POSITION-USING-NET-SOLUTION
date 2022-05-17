CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C NAME     : ADJUSTMENT                                             C
C AUTHOR   : QXP                                                    C 
C PURPOSE  : COMPUTE THE COORDINATE OF UNDERWATER OBJECT(X_CHAMP)  C
C            IN THIS SUBROUTINE, ALL CORRECTION HAVE BE DONE       C  
C USED SUBROUTINE  : DC_COM                                         C
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC           
      SUBROUTINE NETADJUST(ISYS,XYZ,NOBS,DIS_DAT,DIS_RANG,X_CHAMP0,
     1            DLW,DNW,NPAR,SIGMA,DARC_INT,SOUND_VEL)

       IMPLICIT REAL *8 (A-H,O-Z)

       DIMENSION XYZ(4,5),X_CHAMP0(14),DNW(14),SOUND_VEL(6)
       DIMENSION DIS_RANG(3,6000),DLW(6000),DIS_DAT(3,6000)
       DIMENSION DC_TEMP(14)
       INTEGER(4) NOBS,NPAR,IPRN,IOBS,ISYS

       DIMENSION DC(NPAR)

      REAL*8,  DIMENSION(:,:),ALLOCATABLE :: DA    
      REAL*8,  DIMENSION(:,:),ALLOCATABLE :: DN    
      REAL*8,  DIMENSION(:),ALLOCATABLE :: DL    
      REAL*8,  DIMENSION(:),ALLOCATABLE :: PP    
      REAL*8,  DIMENSION(:),ALLOCATABLE :: DP0    
      REAL*8,  DIMENSION(:),ALLOCATABLE :: X_CHAMP    
      REAL*8,  DIMENSION(:),ALLOCATABLE :: BLH    

                ALLOCATE(DA(NOBS,NPAR),   STAT=IRC)                                                                                                                                                                                                                         
                ALLOCATE(DN(NPAR,NPAR),   STAT=IRC)                                                                                                                                                                                                                         
                ALLOCATE(DL(NOBS),   STAT=IRC)                                                                                                                                                                                                                         
                ALLOCATE(PP(NOBS),   STAT=IRC)                                                                                                                                                                                                                         
                ALLOCATE(DP0(NOBS),   STAT=IRC)                                                                                                                                                                                                                         
                ALLOCATE(X_CHAMP(NPAR),   STAT=IRC)                                                                                                                                                                                                                         
                ALLOCATE(BLH(NPAR),   STAT=IRC)                                                                                                                                                                                                                         

          DELTR_0 = 1.0D9
          PI=3.141592653589793D0
          SIGMA = 1000000.0D0
          ITER = 0
          IK = 0
          DC(1:NPAR) = 0.0D0
          DNW(1:NPAR) = 0.0D0

          X_CHAMP(1:NPAR) =  X_CHAMP0(1:NPAR)

C
C....BEGIN COMPUTE THE COORDINATE
C....THIS IS THE HEAD OF CIRCLE
100     CONTINUE

C
C...INITIALIZE OBSERVATION	PP(IPRN) AND RESIDUAL DL(IPRN)
C...INITIALIZE THE OBSERVATION MATRIX DA(IPRN,NPAR)
	   IOBS = 0
 	     DO IPRN = 1,NOBS
	       PP(IPRN)=0.0D0
	       DL(IPRN)=0.0D0
	       DO J =1,NPAR
	         DA(IPRN,J)=0.0D0
	       ENDDO
	     ENDDO
C
C =======================================================                                                                                                                                                                                                                
C...BEGIN COMPUTE DA AND DL
               DSUM=0.0D0
               DALL_COSZ = 0.0D0
 	     DO IPRN = 1,NOBS
                IF(DIS_DAT(2,IPRN).EQ.102.0D0) KK = 1
                IF(DIS_DAT(2,IPRN).EQ.103.0D0) KK = 2
                IF(DIS_DAT(2,IPRN).EQ.104.0D0) KK = 3
                IF(DIS_DAT(2,IPRN).EQ.203.0D0) KK = 4
                IF(DIS_DAT(2,IPRN).EQ.204.0D0) KK = 5
                IF(DIS_DAT(2,IPRN).EQ.304.0D0) KK = 6

                IF(KK.EQ.3) THEN
                CONTINUE
                ENDIF
               IF(IPRN.GT.1) THEN
               DELARC = DIS_RANG(1,IPRN) - DIS_RANG(1,IPRN-1) 
               ENDIF
               PP(IPRN) = DIS_RANG(3,IPRN)
               IK1 = INT(DIS_RANG(2,IPRN)/100)
               IK2 = DMOD(DIS_RANG(2,IPRN),100.0D0)
               DX1 =    X_CHAMP((IK1-1)*3+1) 
     1                - X_CHAMP((IK2-1)*3+1)
               DY1 =    X_CHAMP((IK1-1)*3+2) 
     1                - X_CHAMP((IK2-1)*3+2)
               DZ1 =    X_CHAMP((IK1-1)*3+3) 
     1                - X_CHAMP((IK2-1)*3+3)
	         DR = DSQRT(DX1**2+DY1**2+DZ1**2)
	         DR_12 = DSQRT(DX1**2+DY1**2)

C ---------------------------------------------                                                                                                                                                                                                                 
C.....COMPUTE DA
C.....(K-1)*3+1: THE OBSERVATION COEFFICIENT OF N COORDINATE COMPONENT OF KTH STATION
C.....(K-1)*3+2: THE OBSERVATION COEFFICIENT OF E COORDINATE COMPONENT OF KTH STATION
C.....(K-1)*3+3: THE OBSERVATION COEFFICIENT OF U COORDINATE COMPONENT OF KTH STATION
C.....13:        THE OBSERVATION COEFFICIENT OF THE SOUND SPEED BIAS
C.....14:        THE OBSERVATION COEFFICIENT OF THE ACOUSTIC RAY BENDING
               DA(IPRN,(IK1-1)*3+1) =  DX1/DR
               DA(IPRN,(IK1-1)*3+2) =  DY1/DR
               DA(IPRN,(IK1-1)*3+3) =  DZ1/DR
               DA(IPRN,(IK2-1)*3+1) = -DX1/DR
               DA(IPRN,(IK2-1)*3+2) = -DY1/DR
               DA(IPRN,(IK2-1)*3+3) = -DZ1/DR
               DA(IPRN,13) =  DIS_DAT(3,IPRN)
               DA(IPRN,14) =  (DR_12*1.0D-3)**2
C.....COMPUTE THE WEIGHT DP0(IPRN)=2000.0/RANGE
                DP0(IPRN)= 1.0D0
!               IF(ISYS.GE.3)DP0(IPRN) = DCOS(DABS(DZ1)/DR)**2
!               DP0(IPRN)= (2000.0D0/DR)
               DL_SVP =   DA(IPRN,13)*X_CHAMP(13) 
     3                  + DA(IPRN,14)*X_CHAMP(14) 

!        WRITE(13,'(I4,2F12.3)')
!     1           ITER,DIS_RANG(1,IPRN)/3600.0D0,DA(IPRN,14)*X_CHAMP(14)
!
               DL(IPRN)=PP(IPRN)- DR - DL_SVP

         IF(DABS(DL(IPRN)).GT.3.0D0*SIGMA) DP0(IPRN)= 0.0D0
         IF(DABS(DP0(IPRN)).GT.1.0D-6)
     1         DSUM=DSUM+DL(IPRN)*DL(IPRN)
	  ENDDO
C...END OF COMPUTE DA AND DL
C =======================================================                                                                                                                                                                                                                
	       IF(NOBS.LT.NPAR) RETURN
             SIGMA = DSQRT(DSUM/(NOBS-NPAR))
!C =======================================================                                                                                                                                                                                                                
!C....CALL DC_COMP() FOR ADJUSTMENT
	       CALL DC_COMP(ISYS,DA,DL,DP0,DC,DN,NOBS,NPAR)
              DO IM = 1,NPAR
               X_CHAMP(IM) = X_CHAMP(IM) + DC(IM)
              ENDDO

              
              DELT = 0.0D0
             DO IM = 1,NPAR
              DELT = DELT + DC(IM)**2
             ENDDO
             DELTR= DSQRT(DELT)
	       ITER = ITER + 1
!
C     COMPUTE SIGMA
             DSUM=0.0D0
 
             DO IPRN=1,NOBS

           IF(DABS(DL(IPRN)).GT.3.0D0*SIGMA) DP0(IPRN)= 0.0D0

           IF(DABS(DP0(IPRN)).GT.1.0D-6)DSUM=DSUM+DL(IPRN)*DL(IPRN)
             ENDDO
             SIGMA = DSQRT(DSUM/(NOBS-NPAR))

       DO IM =1,NPAR
         X_CHAMP0(IM) = X_CHAMP(IM)
         DNW(IM) = DN(IM,IM)
       ENDDO

!          DO IM = 1,4
!          WRITE(13,'(I4,2F12.3,5F15.3)') ITER,SIGMA,DELTR,
!     1            X_CHAMP((IM-1)*3+1:(IM-1)*3+4)
!          ENDDO
!        WRITE(13,'(A40)') '========================================'
!        WRITE(13,*)
C 
C ---------------------------------------------                                                                                                                                                                                                                 
C.....DETERMINE WHETHER ITERATIVE OR NOT
        DELTR_12 = DABS(DELTR_0-DELTR)
	 IF(DELTR.LE.1.0D-2) GOTO 300 
       IF(DELTR_12.GT.1.0D-2) DELTR_0 = DELTR
	 IF(ITER.GT.30) GOTO 300 
	 IF(DELTR_12.GT.1.0D-2) GOTO 100 
C...END OF ONE TIMES COMPUTE THE COORDINATE
C   THIS IS THE END OF CIRCLE
300     CONTINUE

        DO IOBS =1,NOBS
           DLW(IOBS)=DL(IOBS)
        ENDDO

             DEALLOCATE(DA,STAT=IRC)   
             DEALLOCATE(DN,STAT=IRC)   
             DEALLOCATE(PP,STAT=IRC)   
             DEALLOCATE(DL,STAT=IRC)   
             DEALLOCATE(DP0,STAT=IRC)   
             DEALLOCATE(X_CHAMP,STAT=IRC)   
             DEALLOCATE(BLH,STAT=IRC)   

        RETURN
	END
