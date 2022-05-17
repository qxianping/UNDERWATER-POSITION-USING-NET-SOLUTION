		SUBROUTINE MATRIX_MUL(X,Y,Z,L,M,N)
C...PURPOSE:  COMPUTE THE LXN MATRIX Z WHICH IS THE PRODUCT OF THE                                                                
C              THE LXM MATRIX X WITH THE MXN MATRIX Y                                                                 
C    TRANSPOSE OF                                                                                                                             
        IMPLICIT REAL*8(A-H,O-Z)                                        
        DIMENSION X(L,M),Y(M,N),Z(L,N)                                  
        DO 1 I=1,L                                                      
           DO 1 J=1,N                                                                
                  Z(I,J)=0.0D0                                              
              DO 1 K=1,M                                                
 1               Z(I,J)=Z(I,J)+X(I,K)*Y(K,J)                            
        RETURN                                                          
        END
		SUBROUTINE MATRIX_REV(X,Y,M,N)
C...PURPOSE:  COMPUTE THE LXN MATRIX Z WHICH IS THE PRODUCT OF THE                                                                
C              THE LXM MATRIX X WITH THE MXN MATRIX Y                                                                 
C    TRANSPOSE OF                                                                                                                             
        IMPLICIT REAL*8(A-H,O-Z)                                        
        DIMENSION X(M,N),Y(N,M)
        DO 1 I=1,M                                                      
           DO 1 J=1,N                                                                
 1               Y(J,I) = X(I,J)                            
        RETURN                                                          
        END

	 SUBROUTINE MATRIX_ADD(X,Y,Z,M,N)
C...PURPOSE:  COMPUTE THE MXN MATRIX Z WHICH IS THE PRODUCT OF THE                                                                
C             THE MATRIX X ADD THE MATRIX Y (Z=X+Y)                                                                 
        IMPLICIT REAL*8(A-H,O-Z)                                        
        DIMENSION X(M,N),Y(M,N),Z(M,N)
        DO 2 I=1,M                                                      
           DO 2 J=1,N                                                                
 2               Z(I,J) = 0.0D0
        DO 1 I=1,M                                                      
           DO 1 J=1,N                                                                
 1               Z(I,J) = X(I,J)+Y(I,J)                            
        RETURN                                                          
        END

	 SUBROUTINE MATRIX_MINUS(X,Y,Z,M,N)
C...PURPOSE:  COMPUTE THE MXN MATRIX Z WHICH IS THE PRODUCT OF THE                                                                
C             THE MATRIX X MINUS THE MATRIX Y (Z=X-Y)                                                                 
        IMPLICIT REAL*8(A-H,O-Z)                                        
        DIMENSION X(M,N),Y(M,N),Z(M,N)
        DO 2 I=1,M                                                      
           DO 2 J=1,N                                                                
 2               Z(I,J) = 0.0D0
        DO 1 I=1,M                                                      
           DO 1 J=1,N                                                                
 1               Z(I,J) = X(I,J)-Y(I,J)                            
        RETURN                                                          
        END


      SUBROUTINE M33T31 (A, B, C)
C
C...PURPOSE:  COMPUTE THE 3X1 MATRIX C WHICH IS THE PRODUCT OF THE
C             TRANSPOSE OF THE 3X3 MATRIX A WITH THE 3X1 MATRIX B
C
      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION A(9), B(3), C(3)
C
      C(1) = A(1)*B(1) + A(2)*B(2) + A(3)*B(3)
      C(2) = A(4)*B(1) + A(5)*B(2) + A(6)*B(3)
      C(3) = A(7)*B(1) + A(8)*B(2) + A(9)*B(3)
C
      RETURN
      END

