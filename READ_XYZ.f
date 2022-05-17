      SUBROUTINE READ_XYZ(LFNAM,XYZ)
      IMPLICIT REAL*8(A-H,O-Z)

	 CHARACTER*166  FILE_POINT
       CHARACTER*180 BUFF180                                            
       INTEGER(4) LFNAM,NUM
       DIMENSION XYZ(4,5)

!C-------------------------------------------------------------
C....  READ COORDINATES OF STATION
        NUM = 0
C....//////////GNSS/A RESULT///////////////////////
       WRITE (FILE_POINT,'(256X)')
       READ (15,'(A)') FILE_POINT
       OPEN (11,FILE=FILE_POINT,STATUS = 'OLD',ERR=999)
       DO NUM = 1,4
       READ (11,'(A180)',END=199) BUFF180
         READ(BUFF180,'(12X,3(F15.3,2X),2F15.3)',ERR=199) XYZ(NUM,1:5)
       ENDDO
199    CONTINUE
       CLOSE(11)


999   CONTINUE
       RETURN
	END
