SUBROUTINE SORT3(N,RA,RB,RC)
      IMPLICIT REAL*8(A-H,O-Z)
! PURPOSE: SORT BY NA 
! CODE FROM THE COMMON NUMERICAL ALGORITHM SET OF VISUAL FORTRAN 2002 
! MODIFIED  2020
INTEGER(4) N
REAL*8 RA(N),RB(N),RC(N)
L=N/2+1
IR=N
DO
  IF(L>1) THEN
    L=L-1
    RRA=RA(L)
    RRB=RB(L)
    RRC=RC(L)
  ELSE
    RRA=RA(IR)
    RRB=RB(IR)
    RRC=RC(IR)
    RA(IR)=RA(1)
    RB(IR)=RB(1)
    RC(IR)=RC(1)
    IR=IR-1
    IF(IR==1) THEN
      RA(1)=RRA
      RB(1)=RRB
      RC(1)=RRC
      RETURN
    ENDIF
  ENDIF
  I=L
  J=L+L
  DO WHILE(J<=IR) 
    IF(J.LT.IR) THEN
      IF(RA(J)<RA(J+1)) J=J+1
    ENDIF
    IF(RRA<RA(J)) THEN
      RA(I)=RA(J)
      RB(I)=RB(J)
      RC(I)=RC(J)
      I=J
      J=J+J
    ELSE
      J=IR+1
    ENDIF
  END DO
  RA(I)=RRA
  RB(I)=RRB
  RC(I)=RRC
END DO
END SUBROUTINE SORT3
