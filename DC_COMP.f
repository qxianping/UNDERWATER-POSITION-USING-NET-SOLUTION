CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C NAME     : DC_COM                                                 C
C AUTHOR   : QXP                                                    C 
C PURPOSE  : COMPUTE THE EQUATION:    DA0*DC=DL0                    C  
C INPUT  PARAMETERS:                                                C
C           DA0  --- THE MATRIX OF EQUATION                         C
C           DL0  --- THE REMAIN OF OBSERVATION                      C
C           NOBS --- THE NUMBER OF OBSERVATION                      C
C           NPAR --- THE NUMBER OF SOLVING PARAMETER                C
C OUTPUT PARAMETERS:                                                C
C           DC   --- THE RESULT OF SOLUTION                         C                                                     
C COMMON DATA      :                                                C
C USED SUBROUTINE  :   DC_CONTRAL(), MATRIX_REV(), MATRIX_MUL       C
C                      DCHOLS()                                     C
C                                                                   C
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC           
      SUBROUTINE DC_COMP(ISYS,DA0,DL0,DP0,DC,DN0,NOBS,NPAR)
       IMPLICIT REAL *8 (A-H,O-Z)


      DIMENSION DA0(NOBS,NPAR),DL0(NOBS),DC(NPAR),DP0(NOBS)
      DIMENSION DN0(NPAR,NPAR)

      REAL*8,  DIMENSION(:,:),ALLOCATABLE :: DN    
      REAL*8,  DIMENSION(:,:),ALLOCATABLE :: DA    
      REAL*8,  DIMENSION(:,:),ALLOCATABLE :: DAT    
      REAL*8,  DIMENSION(:),ALLOCATABLE :: DU    
      REAL*8,  DIMENSION(:),ALLOCATABLE :: DL    
      REAL*8,  DIMENSION(:),ALLOCATABLE :: DP    
      INTEGER(4) ISYS

                ALLOCATE(DA(NOBS,NPAR),   STAT=IRC)                                                                                                                                                                                                                         
                ALLOCATE(DAT(NPAR,NOBS),   STAT=IRC)                                                                                                                                                                                                                         
                ALLOCATE(DN(NPAR,NPAR),   STAT=IRC)                                                                                                                                                                                                                         
                ALLOCATE(DU(NPAR),   STAT=IRC)                                                                                                                                                                                                                         
                ALLOCATE(DL(NOBS),   STAT=IRC)                                                                                                                                                                                                                         
                ALLOCATE(DP(NOBS),   STAT=IRC)                                                                                                                                                                                                                         



           DN0(1:NPAR,1:NPAR) = 0.0D0
        DO I =1,NPAR
	         DU(I) = 0.0D0
	         DC(I) = 0.0D0
	      DO J = 1,NPAR
	         DN(I,J) = 0.0D0
	      ENDDO
	    ENDDO

		IOBS = 0
C
C...PUT DL0 INTO DL,DA0 INTO DA
	DO I = 1, NOBS
	   IF(DL0(I).NE.0.0)THEN
	        IOBS = IOBS + 1
	        DL(IOBS) = DL0(I)
	        DP(IOBS) = DP0(I)
	      DO J =1,NPAR
	        DA(IOBS,J) = DA0(I,J)
	      ENDDO
	   ENDIF
	ENDDO

	      CALL MATRIX_REV(DA,DAT,NOBS,NPAR)
C  COMPUTE ATP
          DO I = 1, NPAR
            DO J = 1, NOBS
              DAT(I,J) = DAT(I,J)*DP(J)
           ENDDO
          ENDDO
	      CALL MATRIX_MUL(DAT,DA,DN,NPAR,NOBS,NPAR)

C===========================================================
C.....ISYS=0:COMPUTE THE ORIGINAL RESIDUAL
C.....ISYS=1:COMPUTE THE SOUND SPEED BIAS
C.....ISYS=2:NOT USED
C.....ISYS=3:COMPUTE THE ACOUSTIC RAY BENDING
C.....ISYS=4:COMPUTE 2D COORDINATE 
C.....ISYS=5:COMPUTE 3D COORDINATE 
C....
         CALL DC_CONTRAL(ISYS,DN,NPAR)
C===========================================================

	      CALL MATRIX_MUL(DAT,DL,DU,NPAR,NOBS,1)

            CALL DCHOLS(DN,NPAR,IERR)
!            

            IF(IERR.NE.0.0) STOP 'FAIL TO INVERT MATRIS'
C     COMPUTE X (TDBIAS AND TBIAS)            
            DO I=1,NPAR
               DSUM=0.0D0
              DO J=1,NPAR
                  DSUM=DSUM+DN(I,J)*DU(J)
              ENDDO
                  DC(I)=DSUM
            ENDDO


            DO I=1,NPAR
              DO J=1,NPAR
                  DN0(I,J) = DN(I,J)
              ENDDO
            ENDDO

             DEALLOCATE(DA,STAT=IRC)   
             DEALLOCATE(DAT,STAT=IRC)   
             DEALLOCATE(DN,STAT=IRC)   
             DEALLOCATE(DU,STAT=IRC)   
             DEALLOCATE(DL,STAT=IRC)   
             DEALLOCATE(DP,STAT=IRC)   
	  RETURN
	END
