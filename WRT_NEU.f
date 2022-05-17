      SUBROUTINE WRT_NEU(LFNAM,XYZ,DNEU)
C                                                                                                                                 
     
      IMPLICIT REAL*8(A-H,O-Z)

       INTEGER(4) LFNAM
       DIMENSION XYZ(4,5)
       DIMENSION DNEU(4,5)

!C-------------------------------------------------------------
       DO IK = 1, 4
           DNEU(IK,1:3) =  XYZ(IK,1:3)
       ENDDO
C.....OUTPUT THE COORDINATE BEFORE ADJUSTMENT
        WRITE(14,'(A60)') '==========NEU:BEFORE ADJUSTMENT============='
          DO I = 1,4
          WRITE(14,'(8X,3(F15.3,2X))') DNEU(I,1:3)
          ENDDO
        WRITE(LFNAM,*)

999   CONTINUE
       RETURN
	END
