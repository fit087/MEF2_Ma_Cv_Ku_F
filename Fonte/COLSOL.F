      SUBROUTINE COLSOL (A,V,MAXA,NN,NWK,NNM,KKK,IOUT)                  COL00001
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . COL00002
C .                                                                   . COL00003
C .   P R O G R A M                                                   . COL00004
C .        TO SOLVE FINITE ELEMENT STATIC EQUILIBRIUM EQUATIONS IN    . COL00005
C .        CORE, USING COMPACTED STORAGE AND COLUMN REDUCTION SCHEME  . COL00006
C .                                                                   . COL00007
C .  - - INPUT VARIABLES - -                                          . COL00008
C .        A(NWK)    = STIFFNESS MATRIX STORED IN COMPACTED FORM      . COL00009
C .        V(NN)     = RIGHT-HAND-SIDE LOAD VECTOR                    . COL00010
C .        MAXA(NNM) = VECTOR CONTAINING ADDRESSES OF DIAGONAL        . COL00011
C .                    ELEMENTS OF STIFFNESS MATRIX IN A              . COL00012
C .        NN        = NUMBER OF EQUATIONS                            . COL00013
C .        NWK       = NUMBER OF ELEMENTS BELOW SKYLINE OF MATRIX     . COL00014
C .        NNM       = NN + 1                                         . COL00015
C .        KKK       = INPUT FLAG                                     . COL00016
C .            EQ. 1   TRIANGULARIZATION OF STIFFNESS MATRIX          . COL00017
C .            EQ. 2   REDUCTION AND BACK-SUBSTITUTION OF LOAD VECTOR . COL00018
C .        IOUT      = UNIT NUMBER USED FOR OUTPUT                    . COL00019
C .                                                                   . COL00020
C .  - - OUTPUT - -                                                   . COL00021
C .        A(NWK)    = D AND L - FACTORS OF STIFFNESS MATRIX          . COL00022
C .        V(NN)     = DISPLACEMENT VECTOR                            . COL00023
C .                                                                   . COL00024
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . COL00025
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               COL00026
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . COL00027
C .   THIS PROGRAM IS USED IN SINGLE PRECISION ARITHMETIC ON CRAY     . COL00028
C .   EQUIPMENT AND DOUBLE PRECISION ARITHMETIC ON IBM MACHINES,      . COL00029
C .   ENGINEERING WORKSTATIONS AND PCS. DEACTIVATE ABOVE LINE FOR     . COL00030
C .   SINGLE PRECISION ARITHMETIC.                                    . COL00031
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . COL00032
      DIMENSION A(NWK),V(0:NN),MAXA(NNM)                                COL00033
C                                                                       COL00034
C     PERFORM L*D*L(T) FACTORIZATION OF STIFFNESS MATRIX                COL00035
C                                                                       COL00036
      IF (KKK-2) 40,150,150                                             COL00037
   40 DO 140 N=1,NN                                                     COL00038
      KN=MAXA(N)                                                        COL00039
      KL=KN + 1                                                         COL00040
      KU=MAXA(N+1) - 1                                                  COL00041
      KH=KU - KL                                                        COL00042
      IF (KH) 110,90,50                                                 COL00043
   50 K=N - KH                                                          COL00044
      IC=0                                                              COL00045
      KLT=KU                                                            COL00046
      DO 80 J=1,KH                                                      COL00047
      IC=IC + 1                                                         COL00048
      KLT=KLT - 1                                                       COL00049
      KI=MAXA(K)                                                        COL00050
      ND=MAXA(K+1) - KI - 1                                             COL00051
      IF (ND) 80,80,60                                                  COL00052
   60 KK=MIN0(IC,ND)                                                    COL00053
      C=0.                                                              COL00054
      DO 70 L=1,KK                                                      COL00055
   70 C=C + A(KI+L)*A(KLT+L)                                            COL00056
      A(KLT)=A(KLT) - C                                                 COL00057
   80 K=K + 1                                                           COL00058
   90 K=N                                                               COL00059
      B=0.                                                              COL00060
      DO 100 KK=KL,KU                                                   COL00061
      K=K - 1                                                           COL00062
      KI=MAXA(K)                                                        COL00063
      C=A(KK)/A(KI)                                                     COL00064
      B=B + C*A(KK)                                                     COL00065
  100 A(KK)=C                                                           COL00066
      A(KN)=A(KN) - B                                                   COL00067
  110 IF (A(KN)) 120,120,140                                            COL00068
  120 WRITE (*,2000) N,A(KN) !kadu                                      COL00069
      GO TO 800                                                         COL00070
  140 CONTINUE                                                          COL00071
      GO TO 900                                                         COL00072
C                                                                       COL00073
C     REDUCE RIGHT-HAND-SIDE LOAD VECTOR                                COL00074
C                                                                       COL00075
  150 DO 180 N=1,NN                                                     COL00076
      KL=MAXA(N) + 1                                                    COL00077
      KU=MAXA(N+1) - 1                                                  COL00078
      IF (KU-KL) 180,160,160                                            COL00079
  160 K=N                                                               COL00080
      C=0.                                                              COL00081
      DO 170 KK=KL,KU                                                   COL00082
      K=K - 1                                                           COL00083
  170 C=C + A(KK)*V(K)                                                  COL00084
      V(N)=V(N) - C                                                     COL00085
  180 CONTINUE                                                          COL00086
C                                                                       COL00087
C     BACK-SUBSTITUTE                                                   COL00088
C                                                                       COL00089
      DO 200 N=1,NN                                                     COL00090
      K=MAXA(N)                                                         COL00091
  200 V(N)=V(N)/A(K)                                                    COL00092
      IF (NN.EQ.1) GO TO 900                                            COL00093
      N=NN                                                              COL00094
      DO 230 L=2,NN                                                     COL00095
      KL=MAXA(N) + 1                                                    COL00096
      KU=MAXA(N+1) - 1                                                  COL00097
      IF (KU-KL) 230,210,210                                            COL00098
  210 K=N                                                               COL00099
      DO 220 KK=KL,KU                                                   COL00100
      K=K - 1                                                           COL00101
  220 V(K)=V(K) - A(KK)*V(N)                                            COL00102
  230 N=N - 1                                                           COL00103
      GO TO 900                                                         COL00104
C                                                                       COL00105
  800 STOP                                                              COL00106
  900 RETURN                                                            COL00107
C                                                                       COL00108
 2000 FORMAT (//' STOP - STIFFNESS MATRIX NOT POSITIVE DEFINITE',//,    COL00109
     1          ' NONPOSITIVE PIVOT FOR EQUATION ',I8,//,               COL00110
     2          ' PIVOT = ',E20.12 )                                    COL00111
      END                                                               COL00112