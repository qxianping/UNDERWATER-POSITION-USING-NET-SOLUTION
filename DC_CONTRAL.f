      SUBROUTINE  DC_CONTRAL(ISYS,DN,NPAR)
C...PURPOSE:  CONTRAL OF ADJUSTMENT

      IMPLICIT REAL *8 (A-H,O-Z)
      REAL*8 DN(NPAR,NPAR)
      INTEGER(4) ISYS
       
       DELTN = 1.0D0/(50.0D0*50.0D0)
C===========================================================
C.....ISYS=0:COMPUTE THE ORIGINAL RESIDUAL
        IF(ISYS.EQ.0) THEN
         DO I = 1,NPAR
          DN(I,I) =  DN(I,I) + 1.0D90
         ENDDO
            RETURN
        ENDIF
C===========================================================
C.....ISYS=1:COMPUTE THE SOUND SPEED BIAS
        IF(ISYS.EQ.1) THEN
         DO I = 1,NPAR-2
          DN(I,I) =  DN(I,I) + 1.0D90
         ENDDO
          DN(NPAR,NPAR) =  DN(NPAR,NPAR) + 1.0D90
          RETURN
        ENDIF

C===========================================================
C.....ISYS=3:COMPUTE THE ACOUSTIC RAY BENDING
        IF(ISYS.EQ.3) THEN
         DO I = 1,NPAR-1
          DN(I,I) =  DN(I,I) + 1.0D90
         ENDDO
         RETURN
        ENDIF
C===========================================================
C.....ISYS=4:COMPUTE 2D COORDINATE
        IF(ISYS.EQ.4) THEN
         DO I = 1,4
          DN((I-1)*3+1,(I-1)*3+1) =  DN((I-1)*3+1,(I-1)*3+1) + DELTN
          DN((I-1)*3+2,(I-1)*3+2) =  DN((I-1)*3+2,(I-1)*3+2) + DELTN
          DN((I-1)*3+3,(I-1)*3+3) =  DN((I-1)*3+3,(I-1)*3+3) +  1.0D90
         ENDDO
          DN(NPAR-1,NPAR-1) =  DN(NPAR-1,NPAR-1) + 1.0D0
          DN(NPAR,NPAR) =  DN(NPAR,NPAR) + 1.0D0
          RETURN
        ENDIF
C===========================================================
C.....ISYS=5:COMPUTE 3D COORDINATE
        IF(ISYS.EQ.5) THEN
         DO I = 1,4
          DN((I-1)*3+1,(I-1)*3+1) =  DN((I-1)*3+1,(I-1)*3+1) + DELTN
          DN((I-1)*3+2,(I-1)*3+2) =  DN((I-1)*3+2,(I-1)*3+2) + DELTN
          DN((I-1)*3+3,(I-1)*3+3) =  DN((I-1)*3+3,(I-1)*3+3) + DELTN
         ENDDO
          DN(NPAR-1,NPAR-1) =  DN(NPAR-1,NPAR-1) + 1.0D0
          DN(NPAR,NPAR) =  DN(NPAR,NPAR) + 0.0D0
          RETURN
        ENDIF

        RETURN

      END   
                                                                                                                                   
