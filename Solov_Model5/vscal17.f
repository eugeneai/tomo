C     calculation coefficients of decomposition over spherical 
C  haronics
C      September 6, 2003

c  schitaet coefficienti razlojeniya funcsii po spherical harmonics
c    continue program  vscal2.f     20 September, 2003

c  8 October 2003
c  vichel nulevuyu gsrmoniku iz modeli, vosstanovil koeff.,
c  prosummiroval i dobavil nulevuyu garmoniky   
c   proeksii vkluchayu v sleduyushei versii :  vscal8.f

c  10 October  L is changed from 0 to LL ( before it was from 1 to LL)
c-----------------------------------------------
c   20 October  Final version (previous version vscal11.f)
c  the line 1133  ubral ogranichenie (R .NE.0 >AND. .NE. 1) and 
c   in DATA (line 1457)  have chenged R form 0 to 1
c    discontinuety (razriv) in graf of projection is disapered.
c    (ischezlo)

C   Reconstruction of the vector field
C      27 Octorber 2003



      PROGRAMM COEFSPHER
C---------------------------------------------
      PARAMETER ( NALPHA   =6,
     #            NBETA    =6,
     #            NGAMMA   =6,
     #            NPP      =51,
     #            NKAPPA   =NPP, 
     #            LL       =5,
     #            NNMAX    =2,   !5,
     #            LN       =(LL+1)*(LL+1)*NNMAX,
     #            NITER    =10000,
     #            NRR      =NPP,
     #            NTETA    =41,
     #            NPHI     =41,
     #            RRE      =1.00, ! radius the region of reconstruction 
     #            RCC      =0.0,  ! radius central condactor 
     #            DD       =1.5000, ! distance from the top  
     #            NX     =NPP,  !41,
     #            NY     =NPP,  !41,
     #            NZ     =NPP,   !41,
     #            NU     =NPP,
     #            NV     =NPP,  ! NV=NPP
     #            NW     =51,
     #            NJ     =NALPHA*NBETA*NGAMMA,
     #            EOT    =0.0,
     #            INFNOT =1,
     #            IND    =1,
     #            KSEC   =NX/2+1 )    

  
      PARAMETER (PI = 3.1415926536)
      PARAMETER( EPS=1.0E-6)

      DIMENSION PAR(100),NPAR(100)

      DIMENSION IIND(LN)
      DIMENSION YM(LN)
      DIMENSION AM(LN,LN)
      DIMENSION XM(LN,3),XM2(LN,3),XM0(LN),XM1(LN)
      DIMENSION A1M(LN)



      DIMENSION PPM(NPP),UPPM(NPP)
      DIMENSION ALPHAM(NALPHA),GAMMAM(NGAMMA)

      DIMENSION BETAM(NBETA)

      DIMENSION BETKAPM(NKAPPA)

      DIMENSION XXM(NX),YYM(NY),ZZM(NZ)
      DIMENSION UM(NU),VM(NV),WM(NW)

      DIMENSION PROJ(NV*NJ)
      DIMENSION PROJ2(NV),PROJ3(NV)
      DIMENSION PROJJ(NV*NJ)
      DIMENSION BKAPPAM(NKAPPA)
      DIMENSION RRM(NRR),TETAM(NTETA),PHIM(NPHI)

      DIMENSION DR(NX*NY*NZ,3)

      DIMENSION DR11(NX*NY),DR21(NX*NY)
      DIMENSION DR12(NX*NY),DR22(NX*NY)
      DIMENSION DR13(NX*NY),DR23(NX*NY)



      DIMENSION GM(NRR*NTETA*NPHI,3)
      DIMENSION GGM(NRR*NTETA*NPHI,3)
      DIMENSION DDR(NX*NY*NZ,3)

      DIMENSION GM2(NRR)
      DIMENSION VLM(4)
      
c      DIMENSION AAAM(LLLN,LLLN)
c      DIMENSION DDDM(3*LLLN,3*LLLN)
c      DIMENSION XXXM(LLLN,3), YYYM(3*LLLN)
c      DIMENSION XXXM0(3*LLLN)

      DIMENSION YMVEC(LN,3),YMVEC3(3*LN)
      DIMENSION AMVEC(LN,LN),AMVEC3(3*LN,3*LN)
      DIMENSION XMVEC(3*LN),XMVEC0(3*LN),XMVEC1(3*LN)
      DIMENSION A1MVEC(3*LN)

      DIMENSION PSIM(NX*NZ)

      COMMON/PRIN/ IPRNT
      IPRNT=6

      OPEN(10,   FILE='MatrixRightPart')

      OPEN(16,   FILE='model.gnu')
      OPEN(14,   FILE='projec.gnu')
      OPEN(15,   FILE='proj_num.gnu')
      OPEN(17,   FILE='proj_www.gnu')

      OPEN(11,   FILE='proj_num.idl')
      OPEN(12,   FILE='proj_www.idl')


      OPEN(26,   FILE='model_X.gnu')
      OPEN(27,   FILE='model_Y.gnu')
      OPEN(28,   FILE='model_Z.gnu')

      OPEN(36,   FILE='model_X_rec.gnu')
      OPEN(37,   FILE='model_Y_rec.gnu')
      OPEN(38,   FILE='model_Z_rec.gnu')

c      OPEN(36,   FILE='model_rec.gnu')
      OPEN(18,   FILE='right_1.gnu')
      OPEN(19,   FILE='right_2.gnu')

      OPEN(20,   FILE='sol_exc1.gnu')
      OPEN(21,   FILE='sol_exc2.gnu')
      OPEN(22,   FILE='sol_exc3.gnu')

      OPEN(30,   FILE='cof_exc.gnu')
      OPEN(31,   FILE='cof_rec.gnu')

      OPEN(32,   FILE='cof_exc.idl')
      OPEN(33,   FILE='cof_rec.idl')

      OPEN(46,   FILE='vector_exc.gnu')
      OPEN(48,   FILE='vector_rec.gnu')

      OPEN(47,   FILE='vector_exc_y.gnu')
      OPEN(49,   FILE='vector_rec_y.gnu')

      OPEN(51,   FILE='vector_exc_x.gnu')
      OPEN(52,   FILE='vector_rec_x.gnu')


      OPEN(54,   FILE='phase1.gnu')
      OPEN(56,   FILE='phase2.gnu')
      OPEN(58,   FILE='phase3.gnu')

      OPEN(55,   FILE='phase11.gnu')
      OPEN(57,   FILE='phase22.gnu')
      OPEN(59,   FILE='phase33.gnu')
C---------------------------------------------
      NPAR(4) =NU
      NPAR(5) =NV
      NPAR(6) =NW

      NPAR(7)=NPHI
      NPAR(8)=NTETA

      NPAR(9)= NALPHA
      NPAR(10)=NBETA
      NPAR(18)=NGAMMA

      NPAR(11)=NJ
      NPAR(12)=NX
      NPAR(13)=NY
      NPAR(14)=NZ

      NPAR(19)=NPP
      NPAR(21)=NRR

      PAR(26)=EOT
      PAR(36)=DD
      PAR(37)=RRE
      PAR(42)=RCC


c      AA=0.
c      BB=-1.
c      CC=ATAN2(BB,AA)
c      write(*,*) 'CC=', CC*180/PI

C Proverka indeksov matrisi      
c      LL3=3

c      DO 4 I7=1,3
c         WRITE(*,*) 'I=',I7
c      DO 2 L=1,LL3
c      DO 2 M1=1,L+1
c         M=M1-1
c         IA=(I7-1)*(LL3**2+2*LL3)+ L**2+M
c         WRITE(*,*) 'IA=',IA
c 2    CONTINUE   
c      DO 3 L=1,LL3
c      DO 3 M1=1,L
c         IB=(I7-1)*(LL3**2+2*LL3)+ L**2+L+M1
c         WRITE(*,*) 'IB=',IB
c 3    CONTINUE    
c 4    CONTINUE 

c      return
c 222  continue
C

      CALL DATA(ALPHAM,BETAM,GAMMAM,XXM,YYM,ZZM,
     %          NKAPPA,BKAPPAM,TETAM,PHIM, 
     %          BETKAPM, UM,VM,WM,RRM,PPM,UPPM,NPAR,PAR)



      HPHI =PAR(12)
      HTETA=PAR(15)      


      HALPHA=PAR(38)
      HBETA =PAR(39)
      HGAMMA=PAR(40)
      HUPP  =PAR(41)

      LMAX=(LL+1)**2


      CALL XXMM(LL,NNMAX,XM)

      write(*,*) 'solution exact XM:'
      write(*,117) (XM(I,1), I=1, LN)
      write(*,117) (XM(I,2), I=1, LN)
      write(*,117) (XM(I,3), I=1, LN)



      CALL Vec_Line(XM,XMVEC,LN)


      CALL GNUFORM1 (3*LN,XMVEC,30)   !! plot cof_exc.gnu
      write(*,*) 'solution exact XMVEC'
      write(*,117) (XMVEC(I), I=1, 3*LN)

      CALL IDLFORM1 (3*LN,XMVEC,32)  !! plot cof_exc.idl

c   KK=1,2,3
      KK=1
      CALL GNUFORM1_1 (LN,KK,XM,20)   !! plot sol_exc1.gnu
      KK=2
      CALL GNUFORM1_1 (LN,KK,XM,21)
      KK=3
      CALL GNUFORM1_1 (LN,KK,XM,22)



cc      CALL Model_Infty(NNMAX,RRM,TETAM,PHIM,
cc     #                            NRR,NTETA,NPHI,GM)

c----------------------------
      IF(INFNOT) 2,4,6 
 2    CALL Summa_Harm_Vec(LL,NNMAX,RRM,TETAM,PHIM,
     #                            NRR,NTETA,NPHI,XM,GM)
      CALL POLAR_DECART(GM,XXM,YYM,ZZM,NX,NY,NZ,RRM,TETAM,PHIM,
     %                  NRR,NTETA,NPHI,DR)
      GOTO 8
 4    CALL Model_Infty2(NNMAX,RRM,TETAM,PHIM,
     #                            NRR,NTETA,NPHI,GM)
      CALL POLAR_DECART(GM,XXM,YYM,ZZM,NX,NY,NZ,RRM,TETAM,PHIM,
     %                  NRR,NTETA,NPHI,DR)
      GOTO 8
 6    CALL Model_Solov_2(ZZM,YYM,XXM,NZ,NY,NX,DR,PAR)

 8    CONTINUE

c-------------------

c      CALL PSIFUNM(ZZM,RRM,NZ,NRR,PSIM)

c      CALL GNUFORM_X_Y (ZZM,RRM,NZ,NX,PSIM,17)      


      KNZ=NZ/2
      KNY=NY/2
      KNX=NX/2

      CALL SPRINT_1(DR,DR11,XXM,YYM,ZZM,NX,NY,NZ,KNZ,1,26)
      CALL SPRINT_1(DR,DR12,XXM,YYM,ZZM,NX,NY,NZ,KNZ,2,27)
      CALL SPRINT_1(DR,DR13,XXM,YYM,ZZM,NX,NY,NZ,KNZ,3,28)

c      CALL GNUFORMVEC_Z (XXM,YYM,NX,NY,NZ,KNZ,DR,46,58)

c      CALL GNUFORMVEC_Y (XXM,ZZM,NX,NY,NZ,KNY,DR,47,56)

c      CALL GNUFORMVEC_X (YYM,ZZM,NX,NY,NZ,KNX,DR,51,54)


      KAN1=54
      KAN2=56
      KAN3=58

      CALL PHASE_PRINT(KSEC,IND,XXM,NX,NY,NZ,DR,KAN1,KAN2,KAN3)

c--------------------------------------------------
c      return
c 222  continue



      CALL ProjNum_Vec(DR,PROJ,PROJ2,UM,VM,
     %        WM, ALPHAM,BETAM,GAMMAM,BKAPPAM,  
     %         XXM,YYM,ZZM,NX,NY,NZ,NPAR,PAR,IERR)


cc      CALL ProjNum(DR,PROJ,PROJ2,UM,VM,
cc     %        WM, ALPHAM,BETAM,GAMMAM,BKAPPAM,  
cc     %         XXM,YYM,ZZM,NPAR,PAR,IERR)
      

      KJ=15
      DO 11 I4=1,NPP
      DO 11 II=KJ,KJ  ! 1,NJ
         CALL IALBETGAM(II,NALPHA,NBETA,I1,I2,I3)
         ALPHA=ALPHAM(I1)
         BETA =BETAM(I2)
         GAMMA=GAMMAM(I3)

         I321=(I3-1)*NBETA*NALPHA+(I2-1)*NALPHA+I1
         I4321=(I4-1)*NGAMMA*NBETA*NALPHA + I321
         PROJ2(I4)=PROJ(I4321)
 11   CONTINUE   

      write(*,*) 'alpha_n=',ALPHA*180/PI
      write(*,*) 'beta_n=',BETA*180/PI
      write(*,*) 'gamma_n=',GAMMA*180/PI
      CALL GNUFORM (NPP,1,1,PROJ2,15) 

      CALL IDLFORM1 (NPP,PROJ2,11)  !! plot proj_num.idl

      write(*,*) 'proj2',(PROJ2(I),I=1,NPP)


c      return
c 777  continue


c-----------------------------------------------


c      CALL Summa_W_Harm_Vec(LL,NJ,NNMAX,PPM,ALPHAM,BETAM,
c     #                   GAMMAM,NPP,NALPHA,NBETA,NGAMMA,XM,PROJ)

c---------------------------------------------------

C   noising of projections
      CALL Proj_Noise (PROJ, NPAR, PAR)
c--------------------------------------------------


      DO 14 I4=1,NPP
      DO 14 II=KJ,KJ  ! 1,NJ
         CALL IALBETGAM(II,NALPHA,NBETA,I1,I2,I3)
         I321=(I3-1)*NBETA*NALPHA+(I2-1)*NALPHA+I1
         I4321=(I4-1)*NGAMMA*NBETA*NALPHA + I321
         PROJ3(I4)=PROJ(I4321)
 14   CONTINUE   

      write(*,*) 'alpha_w=',ALPHA*180/PI
      write(*,*) 'beta_w=',BETA*180/PI
      write(*,*) 'gamma_w=',GAMMA*180/PI


      CALL GNUFORM (NPP,1,1,PROJ3,17) !!plot proj_www.gnu

      CALL IDLFORM1 (NPP,PROJ3,12)  !! plot proj_www.idl
c---------------------------------------------------


      CALL RightPart_WW_VEC(PROJ,LL,NJ,NNMAX,
     #           PPM,UPPM,ALPHAM,BETAM,GAMMAM,
     #           NPP,NALPHA,NBETA,NGAMMA,
     #           NPAR,PAR,YMVEC)


c      LLN=(LLL+1)*(LLL+1)*NNMAX
       write(18,*) 'right part:YMVEC(1)'
      write(18,118) (YMVEC(I,1), I=1, LN)
       write(18,*) 'right part:YMVEC(2)'
      write(18,118) (YMVEC(I,2), I=1, LN)
       write(18,*) 'right part:YMVEC(3)'
      write(18,118) (YMVEC(I,3), I=1, LN)

      CALL Vec_Line(YMVEC,YMVEC3,LN)

      write(10,*) 'right part:YMVEC3'
      write(10,117) (YMVEC3(I), I=1, 3*LN)
c-------------------------------------------------

      DO 81 INCP1=1,3
      DO 81 INCP2=1,INCP1     !1,3

         
      NCP1=INCP1
      NCP2=INCP2
      CALL Matrix_WW_VEC(LL,NNMAX,PPM,UPPM,ALPHAM,BETAM,GAMMAM,
     #         NPP,NALPHA,NBETA,NGAMMA,NPAR,PAR,AMVEC,NCP1,NCP2)

cc      write(*,*) 'Matrix AMVEC_
cc      DO 96 I=1,LN
cc      write(*,117) (AMVEC(I,J), J=1, LN)
cc 96   CONTINUE

cc      write(*,*) 'Matrix AMVEC_2'
cc      DO 93 I=1,LN
cc      write(*,117) (AMVEC(I,J,2), J=1, LN)
cc 93   CONTINUE

cc      write(*,*) 'Matrix AMVEC_3'
cc      DO 91 I=1,LN
cc      write(*,117) (AMVEC(I,J,3), J=1, LN)
cc 91   CONTINUE

      CALL Sort_Matr(AMVEC,LN,AMVEC3,NCP1,NCP2)
 81   CONTINUE

  
      write(10,*) 'Matrix AMVEC3'
      DO 97 I=1,3*LN
      write(10,117) (AMVEC3(I,J), J=1, 3*LN)
 97   CONTINUE
c-----------------------------------------------

c initial values for the solution
      DO 13 I=1,3*LN
         XMVEC0(I)=0.
 13   CONTINUE 

      CALL XYANG(AMVEC3,YMVEC3,XMVEC0,XMVEC1,A1MVEC,3*LN,3*LN,NITER)

      CALL GNUFORM1 (3*LN,XMVEC1,31)   !! cof_rec.gnu

      CALL IDLFORM1 (3*LN,XMVEC1,33)   !! cof_rec.idl



c      XMVEC1(18)=XMVEC1(17)
c      XMVEC1(44)=XMVEC1(42)

      WRITE(IPRNT,*) 'solution of system'       
      WRITE(IPRNT,117) (XMVEC1(I),I=1,3*LN) 

      CALL Vec_Line_Inv(XMVEC1,XM2,LN)
      write(*,*) 'solution exact XM2:'
      write(*,117) (XM2(I,1), I=1, LN)
      write(*,117) (XM2(I,2), I=1, LN)
      write(*,117) (XM2(I,3), I=1, LN)   

      CALL Summa_Harm_Vec(LL,NNMAX,RRM,TETAM,PHIM,
     #                            NRR,NTETA,NPHI,XM2,GGM)


      CALL POLAR_DECART(GGM,XXM,YYM,ZZM,NX,NY,NZ,RRM,TETAM,PHIM,
     %                  NRR,NTETA,NPHI,DDR)


      KAN1=55
      KAN2=57
      KAN3=59

      CALL PHASE_PRINT(KSEC,IND,XXM,NX,NY,NZ,DDR,KAN1,KAN2,KAN3)



      CALL CLEAN(ZZM,YYM,XXM,NZ,NY,NX,DDR)

      CALL NORMER3(DR,DDR,NX,NY,NZ,1,SVZKA1)
      write(*,*) 'SVAZKA1', SVZKA1

      CALL NORMER3(DR,DDR,NX,NY,NZ,2,SVZKA2)
      write(*,*) 'SVAZKA2', SVZKA2

      CALL NORMER3(DR,DDR,NX,NY,NZ,3,SVZKA3)
      write(*,*) 'SVAZKA3', SVZKA3

      CALL SPRINT_1(DDR,DR21,XXM,YYM,ZZM,NX,NY,NZ,KNX,1,36)
c      CALL NORMER3(DR11,DR21,NX,NY,NZ,SVZKA1)
      CALL NORMER2(DR11,DR21,NX,NY,SVZKA1)
      write(*,*) 'SVAZKA1', SVZKA1

      CALL SPRINT_1(DDR,DR22,XXM,YYM,ZZM,NX,NY,NZ,KNX,2,37)
c      CALL NORMER3(DR12,DR22,NX,NY,NZ,SVZKA2)
      CALL NORMER2(DR12,DR22,NX,NY,SVZKA2)
      write(*,*) 'SVAZKA2', SVZKA2

      CALL SPRINT_1(DDR,DR23,XXM,YYM,ZZM,NX,NY,NZ,KNX,3,38)
c      CALL NORMER3(DR13,DR23,NX,NY,NZ,SVZKA3)
      CALL NORMER2(DR13,DR23,NX,NY,SVZKA3)
      write(*,*) 'SVAZKA3', SVZKA3

c      CALL GNUFORMVEC_Z (XXM,YYM,NX,NY,NZ,KNZ,DDR,48,59)

c      CALL GNUFORMVEC_Y (XXM,ZZM,NX,NY,NZ,KNY,DDR,49,57)

c      CALL GNUFORMVEC_X (YYM,ZZM,NX,NY,NZ,KNX,DDR,52,55)


c      KAN1=55
c      KAN2=57
c      KAN3=59

c      CALL PHASE_PRINT(KSEC,IND,XXM,NX,NY,NZ,DDR,KAN1,KAN2,KAN3)

 117  FORMAT(12F7.2)
 118  FORMAT(5E20.3)
      STOP
      END


      SUBROUTINE PHASE_PRINT(KSEC,IND,XM,NX,NY,NZ,VELC,KAN1,KAN2,KAN3)
      REAL VELC(NX*NY*NZ,3)
      REAL XM(NX)
      REAL FASE1(NX*NZ)
      
c  IND=-1  -- section XZ
c  IND= 0  -- section XY
c  IND= 1  -- section YZ

      IF(IND) 2,4,6

 2    CONTINUE
      DO 10 I3=1,NZ
      DO 10 I2=KSEC,KSEC
      DO 10 I1=1,NX
         I13=(I3-1)*NX+I1
         I123=(I3-1)*NY*NX+(I2-1)*NX+I1
         RVE3= VELC(I123,3)
         RVE1= VELC(I123,1)
         IF(RVE1 .NE. 0.) THEN 
            FASE1(I13)=ATAN2(RVE3,RVE1)
         ELSE
            FASE1(I13)=0.  !PI/2.
         ENDIF   
 10   CONTINUE
      KZ=NZ/4+1
      CALL GNUFORM_1 (XM,NX,KZ,KZ,FASE1,KAN1)
      KZ=NZ/2+3 
      CALL GNUFORM_1 (XM,NX,KZ,KZ,FASE1,KAN2)
      KZ=2*NZ/3+2
      CALL GNUFORM_1 (XM,NX,KZ,KZ,FASE1,KAN3)
      GOTO 18

 4    CONTINUE
      DO 12 I3=KSEC,KSEC
      DO 12 I2=1,NY
      DO 12 I1=1,NX
         I12=(I2-1)*NX+I1
         I123=(I3-1)*NY*NX+(I2-1)*NX+I1
         RVE3= VELC(I123,2)
         RVE1= VELC(I123,1)
         IF(RVE1 .NE. 0.) THEN 
            FASE1(I12)=ATAN2(RVE3,RVE1)
         ELSE
            FASE1(I12)=0.  !PI/2.
         ENDIF   
 12   CONTINUE
      KY=NY/4+1
      CALL GNUFORM_1 (XM,NX,KY,KY,FASE1,KAN1)
      KY=NY/2+3 
      CALL GNUFORM_1 (XM,NX,KY,KY,FASE1,KAN2)
      KY=2*NY/3+2
      CALL GNUFORM_1 (XM,NX,KY,KY,FASE1,KAN3)
      GOTO 18

 6    CONTINUE
      DO 14 I3=1,NZ
      DO 14 I2=1,NY
      DO 14 I1=KSEC,KSEC
         I23=(I3-1)*NY+I2
         I123=(I3-1)*NY*NX+(I2-1)*NX+I1
         RVE3= VELC(I123,3)

         SIGNUM=SIGN(1.,RVE3)

         RVE1= VELC(I123,2)
         IF(RVE1 .NE. 0.) THEN 
            FASE1(I23)=ATAN2(RVE3,RVE1)
         ELSE
            FASE1(I23)=SIGNUM*PI/2.   !0.  !PI/2.
         ENDIF   
 14   CONTINUE
      KZ=NZ/4+1
      CALL GNUFORM_1 (XM,NX,KZ,KZ,FASE1,KAN1)
      KZ=NZ/2+3 
      CALL GNUFORM_1 (XM,NX,KZ,KZ,FASE1,KAN2)
      KZ=2*NZ/3+2
      CALL GNUFORM_1 (XM,NX,KZ,KZ,FASE1,KAN3)
      GOTO 18

 18   CONTINUE
      RETURN
      END

      SUBROUTINE GNUFORM_1 (XM,NX,NY,KK,GC,KAN)
      DIMENSION XM(*),GC(*)
      DO 10 J=KK,NY
      DO 10 I=1,NX   
         IJ=(J-1)*NX+I
      WRITE(KAN,100) XM(I),GC(IJ)
c      WRITE(KAN, *)
C      WRITE(KAN,'(A1)') char(13)//char(10)
10    CONTINUE 
100   FORMAT(2E12.3)
C110   FORMAT(A4)
      RETURN
      END

      SUBROUTINE GNUFORMVEC (KSEC,IND,XM,YM,ZM,NX,NY,NZ,GC,KAN1,KAN2)
      DIMENSION XM(NX),YM(NY),ZM(NZ)
      COMPLEX GC(NX*NY*NZ,3)
c  IND=-1  -- section XZ
c  IND= 0  -- section XY
c  IND= 1  -- section YZ

100   FORMAT(4E15.3)
      
      IF(IND) 2,4,6

 2    CONTINUE
      DO 10 I3=1,NZ
         ZZ=ZM(I3)
      DO 10 I2=KSEC,KSEC
         YY=YM(I2)
      DO 10 I1=1,NX
         XX=XM(I1)
      I123= (I3-1)*NY*NX+(I2-1)*NX+I1

      WRITE(KAN1,100) XX,ZZ,REAL(GC(I123,1)),REAL(GC(I123,3))
      WRITE(KAN1,100)
      WRITE(KAN2,100) XX,ZZ,AIMAG(GC(I123,1)),AIMAG(GC(I123,3))
      WRITE(KAN2,100) 
C      WRITE(KAN,'(A1)') char(13)//char(10)
10    CONTINUE 
      GOTO 18

 4    CONTINUE
      DO 12 I3=KSEC,KSEC
         ZZ=ZM(I3)
      DO 12 I2=1,NY
         YY=YM(I2)
      DO 12 I1=1,NX
         XX=XM(I1)
      I123= (I3-1)*NY*NX+(I2-1)*NX+I1

      WRITE(KAN1,100) XX,YY,REAL(GC(I123,1)),REAL(GC(I123,2))
      WRITE(KAN1,100)
      WRITE(KAN2,100) XX,YY,AIMAG(GC(I123,1)),AIMAG(GC(I123,2))
      WRITE(KAN2,100) 
C      WRITE(KAN,'(A1)') char(13)//char(10)
12    CONTINUE 
      GOTO 18

 6    CONTINUE
      DO 14 I3=1,NZ
         ZZ=ZM(I3)
      DO 14 I2=1,NY
         YY=YM(I2)
      DO 14 I1=KSEC,KSEC
         XX=XM(I1)
      I123= (I3-1)*NY*NX+(I2-1)*NX+I1

      WRITE(KAN1,100) YY,ZZ,REAL(GC(I123,2)),REAL(GC(I123,3))
      WRITE(KAN1,100)
      WRITE(KAN2,100) YY,ZZ,AIMAG(GC(I123,2)),AIMAG(GC(I123,3))
      WRITE(KAN2,100) 
C      WRITE(KAN,'(A1)') char(13)//char(10)
 14   CONTINUE 
 18   CONTINUE

      RETURN
      END


      SUBROUTINE CLEAN(ZZM,YYM,XXM,NZ,NY,NX,DDR)
      DIMENSION ZZM(NZ),YYM(NY),XXM(NX)
      DIMENSION DDR(NX*NY*NZ,3)

      DO 86 I3=1,NZ
         ZZ=ZZM(I3)
      DO 86 I2=1,NY
         YY=YYM(I2)
      DO 86 I1=1,NX
         XX=XXM(I1)
         I321=(I3-1)*NX*NY+(I2-1)*NX +I1
         CALL EMAX_VEC(NX*NY*NZ,DDR,EMAX,IERR)
 86   CONTINUE
      DO 88 J=1,3
      DO 88 I3=1,NZ
         ZZ=ZZM(I3)
      DO 88 I2=1,NY
         YY=YYM(I2)
      DO 88 I1=1,NX
         XX=XXM(I1)
         I321=(I3-1)*NX*NY+(I2-1)*NX +I1
         IF(ABS(DDR(I321,J)) .LT. 7*EMAX/100) DDR(I321,J)=0.
 88   CONTINUE 
      RETURN
      END

      SUBROUTINE EMAX_VEC(N,F,E,IERR)
      DIMENSION F(N,3)
c      K=1
      E=F(1,1)
      DO 2 J=1,3
      DO 2 I=1,N
      IF(F(I,J).LE.E) GO TO 2
      E=F(I,J)
c      K=I
   2  CONTINUE
      RETURN
      END

      SUBROUTINE EMAX(N,F,E,K,IERR)
      DIMENSION F(N)
      K=1
      E=F(1)
      DO 1 I=1,N
      IF(F(I).LE.E)GO TO 1
      E=F(I)
      K=I
   1  CONTINUE
      RETURN
      END
C
      SUBROUTINE EMIN(N,F,E,K,IERR)
      DIMENSION F(N)
      K=1
      E=F(1)
      DO 1 I=1,N
      IF(F(I).GE.E)GO TO 1
      E=F(I)
      K=I
   1  CONTINUE
      RETURN
      END

      SUBROUTINE PSIFUNM(ZZM,RRM,NZ,NRR,PSIM)
      DIMENSION PSIM(NZ*NRR)
      DIMENSION ZZM(NZ),RRM(NRR)
      PSI0=1.0
      R0=0.7
      RCC=0.
      AL=0.2

      DO 2 I=1,NRR
         RR=RRM(I)
      DO 2 J=1,NZ
         ZZ=ZZM(J)
         IJ=(I-1)*NZ+J
         PSIM(IJ)=PSIFUN(RR,ZZ,PSI0,RR0,RCC,AL)
 2    CONTINUE   

      RETURN
      END

      FUNCTION PSIFUN(RRR,ZZ,PSI0,R0,RCC,AL)
c вычисление PSI функции      
c      PSI0=1.0
c      R0=0.5
c      RCC=0.
cc      Z0=0.5
c      AL=0.2
      RR=RRR-RCC
      PSIFUN=PSI0*RR**2*(2*R0**2-RR**2-4*AL*ZZ**2)/R0**4
c      PSIFUN=PSI0*(RR-R0)**2*(2*R0**2-(RR-R0)**2-4*AL*ZZ**2)/R0**4
      RETURN
      END

      SUBROUTINE Model_Solov_2(ZZM,YYM,XXM,NZ,NY,NX,GM,PAR)
c
      REAL ZZM(NZ),YYM(NY),XXM(NX),PAR(*)
      REAL GM(NX*NY*NZ,3)

      PARAMETER (PI = 3.1415926536)

      RR1=PAR(37)   

      PSI0=1.0
      RR0=0.7  !0.5 
      RCC=PAR(42)  !0.05
      AL =0.2     ! 0.5

      DO 16 I3=1,NZ
         ZZ=ZZM(I3)    
 
      DO 16 I2=1,NY
         YY=YYM(I2)

      DO 16 I1=1,NX
         XX=XXM(I1)
         I321=(I3-1)*NX*NY+(I2-1)*NX+I1

         RR3=SQRT(XX**2+YY**2+ZZ**2)


         IF(RR3 .LE. RR1) THEN 

            RR2=SQRT(XX**2+YY**2)
            RR=RR2-RCC
            PHI=ATAN2(YY,XX)
            CPHI=COS(PHI)
            SPHI=SIN(PHI)
        

            PSIF= PSIFUN(RR2,ZZ,PSI0,RR0,RCC,AL)

            IF(PSIF .GE. 0. .AND. RR .GE. RCC) THEN    

               GR=1./(2*PI*RR)*8*PSI0*AL/RR0**4*RR**2*ZZ

               GZ=1./(2*PI*RR)*(-4*PSI0*RR**3/RR0**4+
     #              2*RR*(2*PSI0/RR0**2-4*AL*PSI0*ZZ**2/RR0**4))

               GM(I321,1)=GR*CPHI
               GM(I321,2)=GR*SPHI
               GM(I321,3)=GZ
            ELSE
               GM(I321,1)=0.
               GM(I321,2)=0.
               GM(I321,3)=0.
            ENDIF 
         ELSE
         GM(I321,1)=0.
         GM(I321,2)=0.
         GM(I321,3)=0.
      ENDIF   
c      write(*,*) 'I1=',I1,'I2=',I2,'I3=',I3
 16   CONTINUE
      RETURN
      END

      SUBROUTINE Model_Solov_2222(ZZM,YYM,XXM,NZ,NY,NX,GM)
c
      REAL ZZM(NZ),YYM(NY),XXM(NX)
      REAL GM(NX*NY*NZ,3)

      PARAMETER (PI = 3.1415926536)

      PSI0=1.0
      RR0=0.7
      RCC=0.
      AL =0.2
c      RCC=0.2

      DO 16 I3=1,NZ
         ZZ=ZZM(I3)    
 
      DO 16 I2=1,NY
         YY=YYM(I2)

      DO 16 I1=1,NX
         XX=XXM(I1)
         I321=(I3-1)*NX*NY+(I2-1)*NX+I1

         RR3=SQRT(XX**2+YY**2+ZZ**2)
         IF(RR3 .GT. 1.) GO TO 16 
         RR2=SQRT(XX**2+YY**2)
         RR=RR2-RCC
         PHI=ATAN2(YY,XX)
         CPHI=COS(PHI)
         SPHI=SIN(PHI)
        
c         CPHI=XX/RR
c         SPHI=YY/RR
c         RR=ABS(XX)

         PSIF= PSIFUN(RR2,ZZ,PSI0,RR0,RCC,AL)


c         IF(PSIF .GT. 0.) THEN    

         IF(PSIF .GE. 0. .AND. RR .GE. RCC) THEN    

c            GR=1./(2*PI*RR)*8*PSI0*AL/RR0**4*RR**2*ZZ

c            GZ=1./(2*PI*RR)*(-4*PSI0*RR**3/RR0**4+
c     #           2*RR*(2*PSI0/RR0**2-4*AL*PSI0*ZZ**2/RR0**4))

            GR=1./(2*PI*RR)*8*PSI0*AL/RR0**4*RR**2*ZZ

            GZ=1./(2*PI*RR)*(-4*PSI0*RR**3/RR0**4+
     #           2*RR*(2*PSI0/RR0**2-4*AL*PSI0*ZZ**2/RR0**4))

         GM(I321,1)=GR*CPHI
         GM(I321,2)=GR*SPHI
         GM(I321,3)=GZ
      ELSE
         GM(I321,1)=0.
         GM(I321,2)=0.
         GM(I321,3)=0.
      ENDIF 
c      IF(RR3 .GE. 0.9) THEN
c         GM(I321,1)=0.
c         GM(I321,2)=0.
c         GM(I321,3)=0
c      ENDIF   
c      write(*,*) 'I1=',I1,'I2=',I2,'I3=',I3
 16   CONTINUE
      RETURN
      END

      SUBROUTINE TORUS(XXM,YYM,ZZM,NX,NY,NZ,PHIM,NPHI,TORM,PAR)
      DIMENSION PHIM(NPHI)
      DIMENSION XXM(NX),YYM(NY),ZZM(NZ)
      DIMENSION PAR(100)
      DIMENSION TORM(NX*NY*NZ,3)
      AA0=0.3
      RR1=0.5
      RR0=0.
      NRR=21
      HRR=(RR1-RR0)/(NRR-1)
      
      DO 4 IR=1,NRR
         RR=RR0+(IR-1)*HRR
      DO 4 I=1,NPHI
         PHI=PHIM(I)
      DO 4 J=1,NPHI
         ALPHA=PHIM(J)

         XX=(AA0+RR*COS(ALPHA))*COS(PHI)
         YY=(AA0+RR*COS(ALPHA))*SIN(PHI)
         ZZ=RR*SIN(ALPHA)

         CALL INDNW3(XXM,YYM,ZZM,NX,NY,NZ,XX,YY,ZZ,PAR,
     %                  J1,J2,J3,J4)
C

      TORM(J4,1)=XX**2+YY**2+ZZ**2   
      TORM(J4,2)=XX**2+YY**2-ZZ**2
      TORM(J4,3)=XX**2-YY**2+ZZ**2

 4    CONTINUE
      RETURN
      END



      SUBROUTINE NORMER3(BT1,AR1,NX,NY,NZ,KK,SVZKA)
CL-  it is culculation of a relative error of 3D massive
CL- BT1 is exact massive
CL- AR1 is reconstructed massive
      DIMENSION BT1(NX*NY*NZ,3),AR1(NX*NY*NZ,3)
C
      CHARACTER*8 ITEX
      COMMON/PRIN/IPRNT
C
      ITEX='NORER3'
      IERR=0
      S1=0.
      S2=0.
C
      DO 2 I3=1,NZ
      DO 2 I2=1,NY
      DO 2 I1=1,NX
C
      J=(I3-1)*NY*NX+(I2-1)*NX+I1
      SR1=(BT1(J,KK)-AR1(J,KK))**2
      SR2=BT1(J,KK)**2
C
      S1=S1+SR1
      S2=S2+SR2
C
   2  CONTINUE
C
      IF(S2.NE.0.) GO TO 4
      IERR=1
      CALL ERRPRN(ITEX,IERR)
      WRITE(IPRNT,100) S2
      RETURN
C
   4  CONTINUE
C
      SVZKA=S1/S2*100.
      RETURN
 100  FORMAT(4X,' the norm of exact massive is equal 0',E10.3)
      END

      SUBROUTINE NORMER2(BT1,AR1,NX,NY,SVZKA)
CL-  it is culculation of a relative error of 2D massive
CL- BT1 is exact massive
CL- AR1 is reconstructed massive
      DIMENSION BT1(*),AR1(*)
C
      CHARACTER*8 ITEX
      COMMON/PRIN/IPRNT
C
      ITEX='NORER3'
      IERR=0
      S1=0.
      S2=0.
C
c      DO 2 I3=1,NZ
      DO 2 I2=1,NY
      DO 2 I1=1,NX
C
c      J=(I3-1)*NY*NX+(I2-1)*NX+I1
      J=(I2-1)*NX+I1
      SR1=(BT1(J)-AR1(J))**2
      SR2=BT1(J)**2
C
      S1=S1+SR1
      S2=S2+SR2
C
   2  CONTINUE
C
      IF(S2.NE.0.) GO TO 4
      IERR=1
      CALL ERRPRN(ITEX,IERR)
      WRITE(IPRNT,100) S2
      RETURN
C
   4  CONTINUE
C
      SVZKA=S1/S2*100.
      RETURN
 100  FORMAT(4X,' the norm of exact massive is equal 0',E10.3)
      END

      SUBROUTINE VLMM_RE (L,M,M1,ALPHA,BETA,GAMMA,VLMRE,VLMIM)

C          Computer the V harmonics
C          VLMRE - real part
C          VLMIM - image part
      PARAMETER (PI=3.1415926)
      INTEGER L,M
      REAL ALPHA,BETA,GAMMA
      REAL VLMRE,VLMIM
      
      DM=1.
      IF(M1 .EQ. 0) DM=0.5

c      KK=MOD((M-M1),2)
      KK=MOD(M1,2)
      SIGNUM=1.
      IF(KK .NE. 0) SIGNUM=-1.

c      KK=MOD((M+M1),2)
c      SIG1=1.
c      IF(KK .NE. 0) SIG1=-1.

c      KK=MOD((M1-M),2)
c      SIG2=1.
c      IF(KK .NE. 0) SIG2=-1.



c      CC1=COS(M*ALPHA + M1*(GAMMA+PI/2))
c      CC2=COS(M*ALPHA - M1*(GAMMA+PI/2))

c      SS1=SIN(M*ALPHA + M1*(GAMMA+PI/2))
c      SS2=SIN(M*ALPHA - M1*(GAMMA+PI/2))

      CC1=COS(M*ALPHA + M1*GAMMA)
      CC2=COS(M*ALPHA - M1*GAMMA)

      SS1= SIN(M*ALPHA + M1*GAMMA)
      SS2= SIN(M*ALPHA - M1*GAMMA)

      PLM1 = PLM1M2(L,M, M1,BETA)
      PLM2 = PLM1M2(L,M,-M1,BETA)

      VLMRE=DM*(CC1*PLM1 + SIGNUM*CC2*PLM2)
      VLMIM=DM*(SS1*PLM1 + SIGNUM*SS2*PLM2)

c      VLMRE=DM*(CC1*PLM1*(1+SIG2) + CC2*PLM2*(1+SIG1))
c      VLMIM=DM*(SS1*PLM1*(1+SIG2) + SS2*PLM2*(1+SIG1))

      RETURN
      END

      SUBROUTINE WWCS(N,L,M,ALPHA,BETA,GAMMA,PP,WWLMRE,WWLMIM)
c culculation  W(nlm;alha,beta,gamma) -functions
      S1=0.
      S2=0.
      DO 10 I=1,L+1
         M1=I-1
         CALL VLMM_RE (L,M,M1,ALPHA,BETA,GAMMA,VLMRE,VLMIM)
         GGLM= GNLMP(N,L,M1,PP)

         S1=S1+VLMRE*GGLM
         S2=S2+VLMIM*GGLM
 10   CONTINUE
      WWLMRE=S1
      WWLMIM=S2
      RETURN
      END


      SUBROUTINE Model_Infty3(NNMAX,XXM,YYM,ZZM,
     #                            NX,NY,NZ,GM)
c
      REAL XXM(NX),YYM(NY),ZZM(NZ)
      REAL GM(NX*NY*NZ,3)
c      REAL XM((LL+1)*(LL+1)*NNMAX,3)

c      COMMON/PRIN/IPRNT
      PARAMETER (PI = 3.1415926536)

      AA=0.4
      AA2=AA**2
      BB=0.3
      BB2=BB**2

      CC1=-20.
      CC2=10.
      CC3=10.
c      AMP=1./(2.*PI*SQRT(2.*PI)*AA*BB*CC)

c      DO 18 KK=1,3

      DO 10 I3=1,NZ
         ZZ=ZZM(I3)
c         RR2=RR**2
      DO 10 I2=1,NY
         YY=YYM(I2)
c         CTETA=COS(TETA)
c         STETA=SIN(TETA)
      DO 10 I1=1,NX
         XX=XXM(I1)
c         CPHI=COS(PHI)
c         SPHI=SIN(PHI)
         RR2=XX**2+YY**2+ZZ**2
         I21=(I2-1)*NX+I1
         I321=(I3-1)*NY*NX+I21
 
         GM(I321,1)=CC1*YY*ZZ*
     #               EXP(-RR2/(2*AA2)) 
         GM(I321,2)=CC2*XX*ZZ*
     #               EXP(-RR2/(2*AA2)) 
         GM(I321,3)=CC3*XX*YY*
     #               EXP(-RR2/(2*AA2))
 10   CONTINUE
      RETURN
      END

      SUBROUTINE Model_Infty2(NNMAX,RRM,TETAM,PHIM,
     #                            NRR,NTETA,NPHI,GM)
c
      REAL RRM(NRR),TETAM(NTETA),PHIM(NPHI)
      REAL GM(NRR*NTETA*NPHI,3)
c      REAL XM((LL+1)*(LL+1)*NNMAX,3)

c      COMMON/PRIN/IPRNT
      PARAMETER (PI = 3.1415926536)

      AA=0.3
      AA2=AA**2
      BB=0.3
      BB2=BB**2

      CC1=20.
      CC2=-10.
      CC3=-10.

c      AMP=1./(2.*PI*SQRT(2.*PI)*AA*BB*CC)

c      DO 18 KK=1,3

      DO 10 I3=1,NRR
         RR=RRM(I3)
         RR2=RR**2
      DO 10 I2=1,NTETA
         TETA=TETAM(I2)
         CTETA=COS(TETA)
         STETA=SIN(TETA)

      DO 10 I1=1,NPHI
         PHI=PHIM(I1)
         CPHI=COS(PHI)
         SPHI=SIN(PHI)
         
         I21=(I2-1)*NPHI+I1
         I321=(I3-1)*NTETA*NPHI+I21
 
         GM(I321,1)=CC1*RR2*STETA*CTETA*SPHI*
     #               EXP(-RR2/(2*AA2)) 
         GM(I321,2)=CC2*RR2*STETA*CTETA*CPHI*
     #               EXP(-RR2/(2*AA2)) 
         GM(I321,3)=CC3*RR2*STETA*STETA*CPHI*SPHI*
     #               EXP(-RR2/(2.*AA2))
 10   CONTINUE
      RETURN
      END

      SUBROUTINE Model_Infty(NNMAX,RRM,TETAM,PHIM,
     #                            NRR,NTETA,NPHI,GM)
c
      REAL RRM(NRR),TETAM(NTETA),PHIM(NPHI)
      REAL GM(NRR*NTETA*NPHI,3)
c      REAL XM((LL+1)*(LL+1)*NNMAX,3)

c      COMMON/PRIN/IPRNT
      PARAMETER (PI = 3.1415926536)

      AA=0.3
      AA2=AA**2
      BB=0.3
      BB2=BB**2
c      AMP=1./(2.*PI*SQRT(2.*PI)*AA*BB*CC)

c      DO 18 KK=1,3

      DO 10 I3=1,NRR
         RR=RRM(I3)
         RR2=RR**2
      DO 10 I2=1,NTETA
         TETA=TETAM(I2)
         CTETA=COS(TETA)
         STETA=SIN(TETA)

      DO 10 I1=1,NPHI
         PHI=PHIM(I1)
         CPHI=COS(PHI)
         SPHI=SIN(PHI)
         
         I21=(I2-1)*NPHI+I1
         I321=(I3-1)*NTETA*NPHI+I21
 
         GM(I321,1)=8.*RR2*SIN(2.*TETA)*SPHI*EXP(-RR2/(2*AA2)) 
         GM(I321,2)=-8.*RR2*SIN(2.*TETA)*CPHI*EXP(-RR2/(2*AA2)) 
c         GM(I321,3)=0.25*EXP(-(RR2*STETA**2/(2.*BB2)))
         GM(I321,3)=2.*RR2*STETA*STETA*SIN(2.*PHI)*
     #               EXP(-(RR2*STETA**2/(2.*BB2)))
 10   CONTINUE
c 18   CONTINUE
      RETURN
      END


      SUBROUTINE Proj_Noise (PROJ, NPAR, PAR)
      DIMENSION PROJ(*),NPAR(*),PAR(*)

c      NPAR(9)= NALPHA
c      NPAR(10)=NBETA
c      NPAR(18)=NGAMMA

      NJ = NPAR(11)
      NPP= NPAR(19)

      EOT= PAR(26)

      IDUM=-115
      NJPP= NJ*NPP

      DO 12 I4=1,NJPP
      PROJ(I4)=PROJ(I4)+ EOT*GASDEV(IDUM)       
  12  CONTINUE

      RETURN
      END



      SUBROUTINE GNUFORMVEC_X (YM,ZM,NX,NY,NZ,KNX,GC,KAN1,KAN2)
      DIMENSION YM(*),ZM(*)
      REAL GC(NX*NY*NZ,3)

      REAL PHAZE(NY)
      
      K1=2
      K2=3
      DO 10 I3=1,NZ
         ZZ=ZM(I3)
      DO 10 I2=1,NY
         YY=YM(I2)
      DO 10 I1=KNX,KNX    !1,NX
c         XX=XM(I1)
      I123= (I3-1)*NY*NX+(I2-1)*NX+I1

      WRITE(KAN1,100) YY,ZZ,GC(I123,K1),GC(I123,K2)
      WRITE(KAN1,100)
10    CONTINUE 

      DO 14 I3=1,NZ
         ZZ=ZM(I3)
      DO 12 I2=1,NY
         YY=YM(I2)

      DO 12 I1=KNX,KNX        !1,NX
         I123= (I3-1)*NY*NX+(I2-1)*NX+I1

      AA=GC(I123,K1)
      BB=GC(I123,K2)
c      PHAZE(I2)=ATAN2(BB,AA)
      IF(AA .NE. 0.0 .AND. BB .NE. 0.0) THEN  
         PHAZE(I2)=ATAN2(BB,AA)
      ELSE
         PHAZE(I2)=0.
      ENDIF
12    CONTINUE    
      WRITE(KAN2,110) (PHAZE(I),I=1,NY)
      WRITE(KAN2,*)
14    CONTINUE 
100   FORMAT(4E15.3)
110   FORMAT(E10.3)
      RETURN
      END

      SUBROUTINE GNUFORMVEC_Y (XM,ZM,NX,NY,NZ,KNY,GC,KAN1,KAN2)
      DIMENSION XM(*),ZM(*)
      REAL GC(NX*NY*NZ,3)

      REAL PHAZE(NX)

      K1=1
      K2=3
      DO 10 I3=1,NZ
         ZZ=ZM(I3)
      DO 10 I2=KNY,KNY   !1,NY
c         YY=YM(I2)
      DO 10 I1=1,NX
         XX=XM(I1)
      I123= (I3-1)*NY*NX+(I2-1)*NX+I1

      WRITE(KAN1,100) XX,ZZ,GC(I123,K1),GC(I123,K2)
      WRITE(KAN1,100)
10    CONTINUE 

      DO 14 I3=1,NZ
         ZZ=ZM(I3)
      DO 12 I1=1,NX
         XX=XM(I1)

      DO 12 I2=KNY,KNY        !1,NY
         I123= (I3-1)*NY*NX+(I2-1)*NX+I1
      AA=GC(I123,K1)
      BB=GC(I123,K2)
c      PHAZE(I1)=ATAN2(BB,AA)
      IF(AA .NE. 0.0 .AND. BB .NE. 0.0) THEN  
         PHAZE(I1)=ATAN2(BB,AA)
      ELSE
         PHAZE(I1)=0.
      ENDIF
12    CONTINUE    
      WRITE(KAN2,110) (PHAZE(I),I=1,NX)
      WRITE(KAN2,*)
14    CONTINUE 

100   FORMAT(4E15.3)
110   FORMAT(E10.3)
      RETURN
      END



      SUBROUTINE GNUFORMVEC_Z (XM,YM,NX,NY,NZ,KNZ,GC,KAN1,KAN2)
      DIMENSION XM(*),YM(*)
      REAL GC(NX*NY*NZ,3)

      REAL PHAZE(NX)

      K1=1
      K2=2
      DO 10 I3=KNZ,KNZ    !1,NZ
      DO 10 I2=1,NY
         YY=YM(I2)
      DO 10 I1=1,NX
         XX=XM(I1)
         I123= (I3-1)*NY*NX+(I2-1)*NX+I1

      WRITE(KAN1,100) XX,YY,GC(I123,K1),GC(I123,K2)
      WRITE(KAN1,100)
10    CONTINUE 

      DO 14 I2=1,NY
         YY=YM(I2)
      DO 12 I3=1,KNZ,KNZ

      DO 12 I1=1,NX
         XX=XM(I1)

         I123= (I3-1)*NY*NX+(I2-1)*NX+I1
      AA=GC(I123,K1)
      BB=GC(I123,K2)
      IF(AA .NE. 0.0 .AND. BB .NE. 0.0) THEN  
         PHAZE(I1)=ATAN2(BB,AA)
      ELSE
         PHAZE(I1)=0.
      ENDIF
12    CONTINUE    
      WRITE(KAN2,110) (PHAZE(I),I=1,NX)
      WRITE(KAN2,*)
14    CONTINUE 

100   FORMAT(4E15.3)
110   FORMAT(E10.3)
      RETURN
      END


      SUBROUTINE SPRINT_1(DR,DR1,XXM,YYM,ZZM,NX,NY,NZ,KNZ,KK,KAN)
c print of some section  
c KK - number of component of vector field       
      DIMENSION DR(NX*NY*NZ,3),DR1(NX*NY)
      DO 8 I3=KNZ,KNZ
      DO 8 I2=1,NY
      DO 8 I1=1,NX
         I21=(I2-1)*NX+I1
         I321=(I3-1)*NX*NY+(I2-1)*NX +I1
         DR1(I21)=DR(I321,KK)
 8    CONTINUE   

      CALL GNUFORM (NX,NY,1,DR1,KAN)
      RETURN
      END

      SUBROUTINE Sort_Matr(AAM,LN,DDM,NCP1,NCP2)
      DIMENSION AAM(LN,LN)
      DIMENSION DDM(3*LN,3*LN)

c      DO 2 I=1,3*LN
c      DO 2 J=1,3*LN
c         DDM(I,J)=0.
c 2    CONTINUE   

      DO 4 J=1,LN
      DO 4 I=1,LN
         DDM((NCP1-1)*LN+I,(NCP2-1)*LN+J)=AAM(I,J)
 4    CONTINUE

      RETURN
      END

      SUBROUTINE Vec_Line(XM,YM,LD)
      DIMENSION XM(LD,3), YM(3*LD)

      DO 4 K=1,3
      DO 4 I=1,LD
         YM((K-1)*LD+I)=XM(I,K)
 4    CONTINUE
      RETURN
      END

      SUBROUTINE Vec_Line_Inv(YM,XM,LD)
      DIMENSION XM(LD,3), YM(3*LD)

      DO 4 K=1,3
      DO 4 I=1,LD
         XM(I,K)=YM((K-1)*LD+I)
 4    CONTINUE
      RETURN
      END

      SUBROUTINE WW_WW_VEC(L1,M1,N1,L2,M2,N2,
     #        PPM,UPPM,NPP,HUPP,
     #        NBETA,NGAMMA,NALPHA,
     #        ALPHAM,BETAM,GAMMAM,NPAR,PAR,VLM,NCP1,NCP2)

      DIMENSION ALPHAM(*),BETAM(*),GAMMAM(*)
      DIMENSION PPM(*),UPPM(*)
      DIMENSION NPAR(*),PAR(*)
      DIMENSION VLMINT(4)      
      DIMENSION VLM(4)      

      VLM(1)=0.
      VLM(2)=0.
      VLM(3)=0.
      VLM(4)=0.

      DO 8 I1=1,L1+1
         K1=I1-1

         S1=0.
         S2=0.
         S3=0.
         S4=0.
      DO 6 I2=1,L2+1
         K2=I2-1

      CALL INT_GNLMP(PPM,UPPM,NPP,N1,L1,K1,N2,L2,K2,
     %               HUPP,GG12)


      CALL INT_VLMM_VEC(L1,M1,K1,L2,M2,K2,ALPHAM,BETAM,GAMMAM,
     %       NBETA,NGAMMA,NALPHA,NPAR,PAR,VLMINT,NCP1,NCP2)

      S1=S1+VLMINT(1)*GG12
      S2=S2+VLMINT(2)*GG12
      S3=S3+VLMINT(3)*GG12
      S4=S4+VLMINT(4)*GG12
 6    CONTINUE
      VLM(1)=VLM(1)+S1
      VLM(2)=VLM(2)+S2
      VLM(3)=VLM(3)+S3
      VLM(4)=VLM(4)+S4
 8    CONTINUE
      RETURN
      END

      SUBROUTINE INT_GNLMP(PPM,UPPM,NPP,N1,L1,K1,N2,L2,K2,
     %                     HUPP,GG12)
      DIMENSION PPM(NPP),UPPM(NPP)

      DIMENSION UG1M(NPP),UG2M(NPP)
      DIMENSION GG1M(NPP),GG2M(NPP)

      DO 4 I1=1,NPP
         PP=PPM(I1)
         GG1M(I1)=GNLMP(N1,L1,K1,PP)
         GG2M(I1)=GNLMP(N2,L2,K2,PP)
 4    CONTINUE

c      CALL GNUFORM (NPP,1,1,GG1M,20) 

      DO 6 I1=1,NPP
         UPP=UPPM(I1)
c         UG1M(I1)= FUNL1(UPPM,GG1M,1,NPP,UPP,IERR)
c         UG2M(I1)= FUNL1(UPPM,GG2M,1,NPP,UPP,IERR)
         UG1M(I1)= FUNL1(PPM,GG1M,1,NPP,UPP,IERR)
         UG2M(I1)= FUNL1(PPM,GG2M,1,NPP,UPP,IERR)
 6    CONTINUE   

c      CALL GNUFORM (NPP,1,1,UG1M,22)

      S=0.
      DO 8 I1=1,NPP
         WW1=1.0
         IF(I1 .EQ.1 .OR. I1 .EQ. NPP) WW1=0.5
         S=S+UG1M(I1)*UG2M(I1)*WW1
 8    CONTINUE   
      GG12= S*HUPP
      RETURN
      END

      SUBROUTINE INT_VLMM_VEC(L1,M1,K1,L2,M2,K2,ALPHAM,BETAM,GAMMAM,
     %  NBETA,NGAMMA,NALPHA,NPAR,PAR,VLMINT,NCP1,NCP2)

      DIMENSION ALPHAM(*),BETAM(*),GAMMAM(*)

      DIMENSION NPAR(*),PAR(*)
      DIMENSION BM1(NBETA*NGAMMA,4),VLMINT(4)
      DIMENSION BM2(NBETA,4)
      DIMENSION WUnit(3)

c      NALPHA=NPAR(9)
c      NBETA =NPAR(10)
c      NGAMMA=NPAR(18)

      HALPHA=PAR(38)
      HBETA =PAR(39)
      HGAMMA=PAR(40)
      
      DO 12 I2=1,NBETA
         BETA=BETAM(I2)

      DO 12 I3=1,NGAMMA

         I23=(I2-1)*NGAMMA+I3

         GAMMA=GAMMAM(I3)
         WW3=1.0
         IF(I3 .EQ.1 .OR. I3 .EQ. NGAMMA) WW3=0.5
         S1=0.
         S2=0.
         S3=0.
         S4=0.
      DO 10 I1=1,NALPHA
         ALPHA=ALPHAM(I1)

c         WUnit(1)=COS(ALPHA)*SIN(BETA)
c         WUnit(2)=SIN(ALPHA)*SIN(BETA)
c         WUnit(3)=COS(BETA)

         WUnit(1)=-SIN(BETA)*COS(GAMMA)
         WUnit(2)= SIN(BETA)*SIN(GAMMA)
         WUnit(3)= COS(BETA)

         WW1=1.0
         IF(I1 .EQ.1 .OR. I1 .EQ. NALPHA) WW1=0.5

         CALL VLMM_RE (L1,M1,K1,ALPHA,BETA,GAMMA,VLMRE1,VLMIM1)
         CALL VLMM_RE (L2,M2,K2,ALPHA,BETA,GAMMA,VLMRE2,VLMIM2)

         VLCC=VLMRE1*VLMRE2*WUnit(NCP1)*WUnit(NCP2)
         VLSC=VLMIM1*VLMRE2*WUnit(NCP1)*WUnit(NCP2)
         VLCS=VLMRE1*VLMIM2*WUnit(NCP1)*WUnit(NCP2)
         VLSS=VLMIM1*VLMIM2*WUnit(NCP1)*WUnit(NCP2)

         S1=S1+VLCC*WW1
         S2=S2+VLSC*WW1
         S3=S3+VLCS*WW1
         S4=S4+VLSS*WW1
 10   CONTINUE
         BM1(I23,1)=S1*HALPHA
         BM1(I23,2)=S2*HALPHA
         BM1(I23,3)=S3*HALPHA
         BM1(I23,4)=S4*HALPHA
 12   CONTINUE   

      DO 16 I2=1,NBETA
         BETA=BETAM(I2)
         S1=0.
         S2=0.
         S3=0.
         S4=0.
      DO 14 I3=1,NGAMMA
         I23=(I2-1)*NGAMMA+I3

         WW3=1.0
         IF(I3 .EQ.1 .OR. I3 .EQ. NGAMMA) WW3=0.5
         S1=S1+BM1(I23,1)*WW3
         S2=S2+BM1(I23,2)*WW3
         S3=S3+BM1(I23,3)*WW3
         S4=S4+BM1(I23,4)*WW3
 14   CONTINUE

         BM2(I2,1)=S1*HGAMMA
         BM2(I2,2)=S2*HGAMMA
         BM2(I2,3)=S3*HGAMMA
         BM2(I2,4)=S4*HGAMMA
 16   CONTINUE

         VLMINT(1)=0.
         VLMINT(2)=0.
         VLMINT(3)=0.
         VLMINT(4)=0.
      DO 18 I2=1,NBETA
         BETA=BETAM(I2)
         WW2=1.0
         IF(I2 .EQ.1 .OR. I2 .EQ. NBETA) WW2=0.5

         VLMINT(1)=VLMINT(1)+BM2(I2,1)*SIN(BETA)*HBETA*WW2
         VLMINT(2)=VLMINT(2)+BM2(I2,2)*SIN(BETA)*HBETA*WW2
         VLMINT(3)=VLMINT(3)+BM2(I2,3)*SIN(BETA)*HBETA*WW2
         VLMINT(4)=VLMINT(4)+BM2(I2,4)*SIN(BETA)*HBETA*WW2

 18   CONTINUE
      RETURN
      END


      SUBROUTINE Matrix_WW_VEC(LL,NNMAX,PPM,UPPM,ALPHAM,BETAM,GAMMAM,
     #                NPP,NALPHA,NBETA,NGAMMA,NPAR,PAR,AM,NCP1,NCP2)

      DIMENSION AM((LL+1)*(LL+1)*NNMAX,(LL+1)*(LL+1)*NNMAX)
      DIMENSION PPM(NPP),ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
      DIMENSION UPPM(NPP)
      DIMENSION NPAR(*),PAR(*)
      DIMENSION VLM(4)


      HALPHA=PAR(38)
      HBETA =PAR(39)
      HGAMMA=PAR(40)
      HUPP  =PAR(41)
    

c      write(*,*) 'Matrix_WW',HALPHA,HBETA,HGAMMA

      LMAX=(LL+1)*(LL+1)


C   cos * cos  (1,1)

      DO 30  N1=1,NNMAX
      DO 30  L1=1,LL+1
         L11=L1-1

      DO 30  M1=1,L11+1
         IIC1=(N1-1)*LMAX + (L11*L11+M1)

      DO 20  N=1,NNMAX
      DO 20  L=1,LL+1
         L22=L-1
      DO 20  M=1,L22+1
         IIC=(N-1)*LMAX + (L22*L22+M)

c        write(*,*) 'Matrix_W: 30-20', IIC,IIC1


      CALL WW_WW_VEC(L11,M1-1,N1,L22,M-1,N,
     #        PPM,UPPM,NPP,HUPP,
     #        NBETA,NGAMMA,NALPHA,
     #        ALPHAM,BETAM,GAMMAM,NPAR,PAR,VLM,NCP1,NCP2)

      AM(IIC,IIC1)=VLM(1)
 20   CONTINUE
 30   CONTINUE
c---------------------------------------------------

C   sin * sin  (2,2)

      DO 50  N1=1,NNMAX
      DO 50  L1=1,LL+1
         L11=L1-1
c        write(*,*) 'Matrix_W: 50 L1=  ',L1
      DO 50  M1=1,L11
         IIS1=(N1-1)*LMAX + (L11*L11+(M1+1)+L11)
      DO 40  N=1,NNMAX
      DO 40  L=1,LL+1
         L22=L-1
      DO 40  M=1,L22
         IIS=(N-1)*LMAX + (L22*L22+(M+1)+L22)

      CALL WW_WW_VEC(L11,M1,N1,L22,M,N,
     #        PPM,UPPM,NPP,HUPP,
     #        NBETA,NGAMMA,NALPHA,
     #        ALPHAM,BETAM,GAMMAM,NPAR,PAR,VLM,NCP1,NCP2)

        AM(IIS,IIS1)=VLM(4)
 40   CONTINUE
 50   CONTINUE
c-----------------------------------------------

C  cos * sin  (1,2)

      DO 70  N1=1,NNMAX
      DO 70  L1=1,LL+1
         L11=L1-1

      DO 70  M1=1,L11+1
         IIC1=(N1-1)*LMAX + (L11*L11+M1)
      DO 60  N=1,NNMAX
      DO 60  L=1,LL+1
         L22=L-1
      DO 60  M=1,L22
         IIS=(N-1)*LMAX + (L22*L22+(M+1)+L22)

      CALL WW_WW_VEC(L11,M1-1,N1,L22,M,N,
     #        PPM,UPPM,NPP,HUPP,
     #        NBETA,NGAMMA,NALPHA,
     #        ALPHAM,BETAM,GAMMAM,NPAR,PAR,VLM,NCP1,NCP2)

         AM(IIS,IIC1)=VLM(3)
 60   CONTINUE
 70   CONTINUE
c------------------------------------------

C  sin * cos   (2,1)

      DO 90  N1=1,NNMAX
      DO 90  L1=1,LL+1
         L11=L1-1

      DO 90  M1=1,L11
         IIS1=(N1-1)*LMAX + (L11*L11+(M1+1)+L11)

      DO 80  N=1,NNMAX
      DO 80  L=1,LL+1
         L22=L-1
      DO 80  M=1,L22+1
         IIC=(N-1)*LMAX + (L22*L22+M)

      CALL WW_WW_VEC(L11,M1,N1,L22,M-1,N,
     #        PPM,UPPM,NPP,HUPP,
     #        NBETA,NGAMMA,NALPHA,
     #        ALPHAM,BETAM,GAMMAM,NPAR,PAR,VLM,NCP1,NCP2)

         AM(IIC,IIS1)=VLM(2)
 80   CONTINUE
 90   CONTINUE
 112  CONTINUE
      RETURN
      END

      SUBROUTINE RightPart_WW_VEC(PROJ,LL,NJ,NNMAX,
     #           PPM,UPPM,ALPHAM,BETAM,GAMMAM,
     #           NPP,NALPHA,NBETA,NGAMMA,
     #           NPAR,PAR,YM)

      DIMENSION PROJ(NPP*NALPHA*NBETA*NGAMMA)

      DIMENSION YM((LL+1)*(LL+1)*NNMAX,3)
      DIMENSION PPM(NPP),UPPM(NPP)
      DIMENSION ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
      DIMENSION NPAR(*),PAR(*)

      DIMENSION BM1(NBETA*NGAMMA,3),BM2(NBETA,3)
      DIMENSION WUnit(3)



      HALPHA=PAR(38)
      HBETA =PAR(39)
      HGAMMA=PAR(40)
      UM1   =PAR(43)

      LMAX=(LL+1)*(LL+1)


c--------------------------------
c  cos

      DO 30  N=1,NNMAX
      DO 30  L=1,LL+1
         L22=L-1
      DO 30  M=1,L22+1
         IIC=(N-1)*LMAX + (L22*L22+M)


      DO 12 I2=1,NBETA
         BETA=BETAM(I2)

      DO 12 I3=1,NGAMMA
         GAMMA=GAMMAM(I3)

         I23=(I2-1)*NGAMMA+I3

         S1=0.
         S2=0.
         S3=0.

      DO 10 I1=1,NALPHA
         ALPHA=ALPHAM(I1)

c         WUnit(1)=COS(ALPHA)*SIN(BETA)
c         WUnit(2)=SIN(ALPHA)*SIN(BETA)
c         WUnit(3)=COS(BETA)

         WUnit(1)=-SIN(BETA)*COS(GAMMA)
         WUnit(2)= SIN(BETA)*SIN(GAMMA)
         WUnit(3)= COS(BETA)

         WW1=1.0
         IF(I1 .EQ.1 .OR. I1 .EQ. NALPHA) WW1=0.5

         CALL INT_Proj_GNLM_P(PROJ,PPM,UPPM,NPP,
     #        N,L22,M-1,
     #        ALPHA,BETA,GAMMA,    
     #        NJ,I1,I2,I3,NPAR,PAR,YYRE,YYIM)

         S1=S1+YYRE*WW1*WUnit(1)
         S2=S2+YYRE*WW1*WUnit(2)
         S3=S3+YYRE*WW1*WUnit(3)
 10   CONTINUE
         BM1(I23,1)=S1*HALPHA
         BM1(I23,2)=S2*HALPHA
         BM1(I23,3)=S3*HALPHA
 12   CONTINUE   

      DO 16 I2=1,NBETA
         S1=0.
         S2=0.
         S3=0.
      DO 14 I3=1,NGAMMA

         I23=(I2-1)*NGAMMA+I3

         WW3=1.0
         IF(I3 .EQ.1 .OR. I3 .EQ. NGAMMA) WW3=0.5
         S1=S1+BM1(I23,1)*WW3
         S2=S2+BM1(I23,2)*WW3
         S3=S3+BM1(I23,3)*WW3
 14   CONTINUE
         BM2(I2,1)=S1*HGAMMA
         BM2(I2,2)=S2*HGAMMA
         BM2(I2,3)=S3*HGAMMA
 16   CONTINUE

         S1=0. 
         S2=0. 
         S3=0. 
      DO 18 I2=1,NBETA
         BETA=BETAM(I2)
         WW2=1.0
         IF(I2 .EQ.1 .OR. I2 .EQ. NBETA) WW2=0.5
         S1=S1+BM2(I2,1)*SIN(BETA)*WW2
         S2=S2+BM2(I2,2)*SIN(BETA)*WW2
         S3=S3+BM2(I2,3)*SIN(BETA)*WW2

 18   CONTINUE
      YM(IIC,1)=S1*HBETA
      YM(IIC,2)=S2*HBETA
      YM(IIC,3)=S3*HBETA
 30   CONTINUE
c--------------------------------------
c     sin
      DO 50  N=1,NNMAX
      DO 50  L=1,LL+1
         L22=L-1
      DO 50  M=1,L22
         IIS=(N-1)*LMAX + (L22*L22+(M+1)+L22)

      DO 22 I2=1,NBETA
         BETA=BETAM(I2)
      DO 22 I3=1,NGAMMA
         GAMMA=GAMMAM(I3)

         I23=(I2-1)*NGAMMA+I3

c         WW3=1.0
c         IF(I3 .EQ.1 .OR. I3 .EQ. NGAMMA) WW3=0.5

         S1=0.
         S2=0.
         S3=0.
      DO 20 I1=1,NALPHA
         ALPHA=ALPHAM(I1)

c         WUnit(1)=COS(ALPHA)*SIN(BETA)
c         WUnit(2)=SIN(ALPHA)*SIN(BETA)
c         WUnit(3)=COS(BETA)

         WUnit(1)=-SIN(BETA)*COS(GAMMA)
         WUnit(2)= SIN(BETA)*SIN(GAMMA)
         WUnit(3)= COS(BETA)

         WW1=1.0
         IF(I1 .EQ.1 .OR. I1 .EQ. NALPHA) WW1=0.5

         CALL INT_Proj_GNLM_P(PROJ,PPM,UPPM,NPP,
     #        N,L22,M,
     #        ALPHA,BETA,GAMMA,    
     #        NJ,I1,I2,I3,NPAR,PAR,YYRE,YYIM)

         S1=S1+YYIM*WW1*WUnit(1)
         S2=S2+YYIM*WW1*WUnit(2)
         S3=S3+YYIM*WW1*WUnit(3)
 20   CONTINUE
         BM1(I23,1)=S1*HALPHA
         BM1(I23,2)=S2*HALPHA
         BM1(I23,3)=S3*HALPHA
 22   CONTINUE   

      DO 26 I2=1,NBETA
         BETA=BETAM(I2)

         S1=0.
         S2=0.
         S3=0.
      DO 24 I3=1,NGAMMA

         I23=(I2-1)*NGAMMA+I3

         WW3=1.0
         IF(I3 .EQ.1 .OR. I3 .EQ. NGAMMA) WW3=0.5
         S1=S1+BM1(I23,1)*WW3
         S2=S2+BM1(I23,2)*WW3
         S3=S3+BM1(I23,3)*WW3
 24   CONTINUE
         BM2(I2,1)=S1*HGAMMA
         BM2(I2,2)=S2*HGAMMA
         BM2(I2,3)=S3*HGAMMA
 26   CONTINUE

         S1=0. 
         S2=0. 
         S3=0. 
      DO 28 I2=1,NBETA
         BETA=BETAM(I2)
         WW2=1.0
         IF(I2 .EQ.1 .OR. I2 .EQ. NBETA) WW2=0.5
         S1=S1+BM2(I2,1)*SIN(BETA)*WW2
         S2=S2+BM2(I2,2)*SIN(BETA)*WW2
         S3=S3+BM2(I2,3)*SIN(BETA)*WW2

 28   CONTINUE
      YM(IIS,1)=S1*HBETA
      YM(IIS,2)=S2*HBETA
      YM(IIS,3)=S3*HBETA
 50   CONTINUE
      RETURN
      END

      SUBROUTINE POLAR_DECART_VEC(GM,XXM,YYM,ZZM,NX,NY,NZ,
     %       RRM,TETAM,PHIM,NRR,NTETA,NPHI,DR)

      DIMENSION GM(NRR*NTETA*NPHI,3)
      DIMENSION XXM(NX),YYM(NY),ZZM(NZ)
      DIMENSION RRM(NRR),TETAM(NTETA),PHIM(NPHI)
      DIMENSION DR(NX*NY*NZ,3)

      DIMENSION GMM1(NRR*NTETA*NPHI)
      DIMENSION GMM2(NRR*NTETA*NPHI)
      DIMENSION GMM3(NRR*NTETA*NPHI)



      DO 2 I=1,NRR*NTETA*NPHI
         GMM1(I)=0.
         GMM2(I)=0.
         GMM3(I)=0.
 2    CONTINUE   
         

    
      DO 4 I3=1,NRR
      DO 4 I2=1,NTETA
      DO 4 I1=1,NPHI

         I21=(I2-1)*NPHI+I1
         I321=(I3-1)*NTETA*NPHI+I21

         GMM1(I321)=GM(I321,1)
         GMM2(I321)=GM(I321,2)
         GMM3(I321)=GM(I321,3)
 4    CONTINUE   


      SS1=0.
      SS2=0.
      SS3=0.
      DO 6 I3=1,NZ
         ZZ=ZZM(I3)
      DO 6 I2=1,NY
         YY=YYM(I2)
      DO 6 I1=1,NX
         XX=XXM(I1)
         I321=(I3-1)*NX*NY+(I2-1)*NX +I1
         CALL POLAR(XX,YY,ZZ,TETA,PHI)
         RRR=SQRT(XX**2+YY**2+ZZ**2)
         CALL REGN(PHIM,TETAM,RRM,NPHI,NTETA,NRR,PHI,TETA,RRR,EI)

         ST=SIN(TETA)
         CT=COS(TETA)
         SP=SIN(PHI)
         CP=COS(PHI)


         IF(EI .GT. 0) THEN

         SS1=FUNL3(PHIM,TETAM,RRM,GMM1,NPHI,NTETA,NRR,
     #               PHI,TETA,RRR,IERR)
         SS2=FUNL3(PHIM,TETAM,RRM,GMM2,NPHI,NTETA,NRR,
     #               PHI,TETA,RRR,IERR)
         SS3=FUNL3(PHIM,TETAM,RRM,GMM3,NPHI,NTETA,NRR,
     #               PHI,TETA,RRR,IERR)
         ENDIF   
         DR(I321,1)=ST*CP*SS1+CT*CP*SS2-SP*SS3
         DR(I321,2)=ST*SP*SS1+CT*SP*SS2+CP*SS3
         DR(I321,3)=CT*SS1   -ST*SS2

 6    CONTINUE
      RETURN
      END

      SUBROUTINE POLAR_DECART(GM,XXM,YYM,ZZM,NX,NY,NZ,RRM,TETAM,PHIM,
     %                  NRR,NTETA,NPHI,DR)

      DIMENSION GM(NRR*NTETA*NPHI,3)
      DIMENSION XXM(NX),YYM(NY),ZZM(NZ)
      DIMENSION RRM(NRR),TETAM(NTETA),PHIM(NPHI)
      DIMENSION DR(NX*NY*NZ,3)

      DIMENSION GMM(NRR*NTETA*NPHI)

      DO 10 KK=1,3

      DO 2 I=1,NRR*NTETA*NPHI
         GMM(I)=0.
 2    CONTINUE   
         

      DO 4 I3=1,NRR
      DO 4 I2=1,NTETA
      DO 4 I1=1,NPHI

         I21=(I2-1)*NPHI+I1
         I321=(I3-1)*NTETA*NPHI+I21

         GMM(I321)=GM(I321,KK)
 4    CONTINUE   

      DO 6 I3=1,NZ
         ZZ=ZZM(I3)
      DO 6 I2=1,NY
         YY=YYM(I2)
      DO 6 I1=1,NX
         XX=XXM(I1)
         I321=(I3-1)*NX*NY+(I2-1)*NX +I1
         CALL POLAR(XX,YY,ZZ,TETA,PHI)
         RRR=SQRT(XX**2+YY**2+ZZ**2)
         CALL REGN(PHIM,TETAM,RRM,NPHI,NTETA,NRR,PHI,TETA,RRR,EI)

         IF(EI .GT. 0) THEN
            DR(I321,KK)=FUNL3(PHIM,TETAM,RRM,GMM,NPHI,NTETA,NRR,
     #               PHI,TETA,RRR,IERR)
         ENDIF   
 6    CONTINUE
 10   CONTINUE
      RETURN
      END

c      FUNCTION REGN_Hole(X1,Y1,RCC,EI)
c      RR
c      IF(X1.LT.XM(1).OR.X1.GT.XM(NX)) GO TO 10
c      IF(Y1.LT.YM(1).OR.Y1.GT.YM(NY)) GO TO 10
cc      IF(Z1.LT.ZM(1).OR.Z1.GT.ZM(NZ)) GO TO 10
c      EI=1.
c      RETURN
c  10  EI=-1.
c      RETURN
c      END

      SUBROUTINE ProjNum_Vec(DR,PROJ,PROJ2,UM,VM,
     %        WM, ALPHAM,BETAM,GAMMAM,BKAPPAM,  
     %         XM,YM,ZM,NX,NY,NZ,NPAR,PAR,IERR)
C
      DIMENSION  ALPHAM(*),BETAM(*),GAMMAM(*)
      DIMENSION  BKAPPAM(*)
      DIMENSION  UM(*),VM(*),WM(*),XM(*),YM(*),ZM(*)
      DIMENSION  PROJ(*),PROJ2(*)
      DIMENSION DR(NX*NY*NZ,3)

      DIMENSION DDR(NX*NY*NZ)
      DIMENSION WUnit(3)
      DIMENSION  NPAR(*),PAR(*)
      CHARACTER*8 ITEX
      PARAMETER (PI = 3.1415926536)

      NALPHA = NPAR(9)
      NBETA  = NPAR(10)
      NGAMMA = NPAR(18)
      NJ     = NPAR(11)

c      NX    = NPAR(12)
c      NY    = NPAR(13)
c      NZ    = NPAR(14)
      NU    = NPAR(4)
      NV    = NPAR(5)
      NW    = NPAR(6)

      DD  = PAR(36)
      RRE = PAR(37)
      RCC = PAR(42)
      UM1 = PAR(43)

      IERR=0
      ITEX='ProjNum'
      IF(NW.LT.100) GO TO 20
      IERR=1
      CALL ERRPRN(ITEX,IERR)
      RETURN
C
  20  CONTINUE
C
      BKAPMAX=ASIN(RRE/DD)
c      BKAPMAX=ATAN2(RRE,DD)

      DO 85 I4=1,NJ
         CALL IALBETGAM(I4,NALPHA,NBETA,I1,I2,I3)
         ALPHA=ALPHAM(I1)
         BETA =BETAM(I2)
         GAMMA=GAMMAM(I3)
         I321=(I3-1)*NBETA*NALPHA+(I2-1)*NALPHA+I1

         CALPHA=COS(ALPHA)
         SALPHA=SIN(ALPHA)
         CBETA=COS(BETA)
         SBETA=SIN(BETA)
         CGAMMA=COS(GAMMA)
         SGAMMA=SIN(GAMMA)

c         AA11= CALPHA*CBETA*CGAMMA - SALPHA*SGAMMA
c         AA12= SALPHA*CBETA*CGAMMA + CALPHA*SGAMMA
c         AA13=-SBETA*CGAMMA
c         AA21=-CALPHA*CBETA*SGAMMA - SALPHA*CGAMMA
c         AA22=-SALPHA*CBETA*SGAMMA + CALPHA*CGAMMA
c         AA23= SBETA*SGAMMA
c         AA31= CALPHA*SBETA
c         AA32= SALPHA*SBETA
c         AA33= CBETA

c         WUnit(1)=CALPHA*SBETA
c         WUnit(2)=SALPHA*SBETA
c         WUnit(3)=CBETA
         WUnit(1)=-SBETA*CGAMMA
         WUnit(2)= SBETA*SGAMMA
         WUnit(3)= CBETA

      DO 6 I3=1,NZ
      DO 6 I2=1,NY
      DO 6 I1=1,NX
         J321=(I3-1)*NX*NY+(I2-1)*NX +I1

         DDR(J321)= DR(J321,1)*WUnit(1)
     #            + DR(J321,2)*WUnit(2)
     #            + DR(J321,3)*WUnit(3)
 6    CONTINUE   

c      DO 65 J1=1,NV
c         VV=VM(J1)
c         QQQ=RRE**2-VV**2
c         IF(QQQ .LT. 0) QQQ=0.

      DO 65 J1=1,NU
         UU=UM(J1)
         QQQ=RRE**2-UU**2
         IF(QQQ .LT. 0) QQQ=0.

         WW0=SQRT(QQQ)
         HW =2.*WW0/(NW-1)

         S=0.
         SS=0.
      DO 55 J3=1,NW
         WW =-WW0+(J3-1)*HW
c         WWDD=1.+WW/DD

c--   T (al,BET,GAM)-------------------------------------------

      XX=
     % (CALPHA*CBETA*CGAMMA - SALPHA*SGAMMA)*UU + 
c     %    (-CALPHA*CBETA*SGAMMA - SALPHA*CGAMMA)*VV+
     %    CALPHA*SBETA*WW


      YY=
     % (SALPHA*CBETA*CGAMMA + CALPHA*SGAMMA)*UU + 
c     %    (-SALPHA*CBETA*SGAMMA + CALPHA*CGAMMA)*VV+
     %    SALPHA*SBETA*WW

      ZZ=
     % -SBETA*CGAMMA*UU + 
c     %     SBETA*SGAMMA*VV +
     %            CBETA*WW

c      XX=
c     % (CALPHA*CBETA*CGAMMA - SALPHA*SGAMMA)*UU + 
cc     %    (-CALPHA*CBETA*SGAMMA - SALPHA*CGAMMA)*VV+
c     %    -SBETA*CGAMMA*WW


c      YY=
c     % (-CALPHA*CBETA*SGAMMA - SALPHA*CGAMMA)*UU + 
cc     %    (-SALPHA*CBETA*SGAMMA + CALPHA*CGAMMA)*VV+
c     %    SBETA*SGAMMA*WW

c      ZZ=
c     % CALPHA*SBETA*UU + 
cc     %     SBETA*SGAMMA*VV +
c     %            CBETA*WW

      RRR=SQRT(XX**2+YY**2)

      CALL REGN(XM,YM,ZM,NX,NY,NZ,XX,YY,ZZ,EI)

      IF(EI .LT. 0.) GO TO 55

      S=FUNL3(XM,YM,ZM,DDR,NX,NY,NZ,XX,YY,ZZ,IERR)

      WW1=1.
      IF(J3.EQ.1.OR.J3.EQ.NW) WW1=0.5
      SS=SS+S*WW1

 55   CONTINUE
      
      IF(RRR .GT. RCC) THEN
         PROJ2(J1)=SS*HW 
      ELSE
         PROJ2(J1)=0. 
      ENDIF
   
 65   CONTINUE

c      DO 75 J1=1,NV
      DO 75 J1=1,NU
         I4321=(J1-1)*NGAMMA*NBETA*NALPHA + I321
         PROJ(I4321)=PROJ2(J1)       
  75  CONTINUE
  85  CONTINUE
      RETURN
      END


      SUBROUTINE Summa_W_Harm_Vec(LL,NJ,NNMAX,PPM,ALPHAM,BETAM,
     #                   GAMMAM,NPP,NALPHA,NBETA,NGAMMA,XM,GM)
c  summation of W  harmonics 
      REAL PPM(NPP),ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
      REAL GM(NPP*NALPHA*NBETA*NGAMMA)
      REAL XM((LL+1)*(LL+1)*NNMAX,3)
c      REAL WWLMRE(3),WWLMIM(3)
      REAL WUnit(3)

c      COMMON/PRIN/IPRNT
c      PARAMETER (PI = 3.1415926536)

      LMAX=(LL+1)*(LL+1)

c      DO 10 I4=1,NPP
c         PP=PPM(I4)
c      DO 10 I3=1,NGAMMA
c         GAMMA=GAMMAM(I3)
c      DO 10 I2=1,NBETA
c         BETA=BETAM(I2)
c      DO 10 I1=1,NALPHA
c         ALPHA=ALPHAM(I1)


      DO 10 II=1,NJ

         CALL IALBETGAM(II,NALPHA,NBETA,I1,I2,I3)

         ALPHA=ALPHAM(I1)
         BETA =BETAM(I2)
         GAMMA=GAMMAM(I3)
         I321=(I3-1)*NBETA*NALPHA+(I2-1)*NALPHA+I1

c         WUnit(1)=COS(ALPHA)*SIN(BETA)
c         WUnit(2)=SIN(ALPHA)*SIN(BETA)
c         WUnit(3)=COS(BETA)

         WUnit(1)=-SIN(BETA)*COS(GAMMA)
         WUnit(2)= SIN(BETA)*SIN(GAMMA)
         WUnit(3)= COS(BETA)


      DO 10 I4=1,NPP
         PP=PPM(I4)

         I4321=(I4-1)*NGAMMA*NBETA*NALPHA + I321

c         I4321=(II-1)*NPP + I4
 

        S1=0.
        DO 4  N=1,NNMAX
        DO 4  L=1,LL+1
           L22=L-1
        DO 4  M=1,L22+1
           IIC=(N-1)*LMAX + (L22*L22+M)

           CALL WWCS(N,L22,M-1,ALPHA,BETA,GAMMA,PP,WWLMRE,WWLMIM)

c           CALL YLM_RE (L,M-1,TETA,PHI,YLMRE,YLMIM)
           S1=S1+(XM(IIC,1)*WUnit(1)+
     #            XM(IIC,2)*WUnit(2)+
     #            XM(IIC,3)*WUnit(3))*WWLMRE
 4      CONTINUE

        DO 6 N=1,NNMAX
        DO 6 L=1,LL+1
           L22=L-1
        DO 6 M=1,L22
           IIS=(N-1)*LMAX + (L22*L22+(M+1)+L22)

           CALL WWCS(N,L22,M,ALPHA,BETA,GAMMA,PP,WWLMRE,WWLMIM)

c           CALL YLM_RE (L,M,TETA,PHI,YLMRE,YLMIM)
           S1=S1+(XM(IIS,1)*WUnit(1)+
     #            XM(IIS,2)*WUnit(2)+
     #            XM(IIS,3)*WUnit(3))*WWLMIM
 6      CONTINUE
        GM(I4321)=S1
 10   CONTINUE
      RETURN
      END

      SUBROUTINE WWCC(N,L,M,ALPHA,BETA,GAMMA,PP,WWLMRE)
c culculation  W(nlm;alha,beta,gamma) -functions

      DIMENSION WWLMRE(3)

      S1=0.
      DO 10 I=1,L+1
         M1=I-1
         CALL VLMM_RE (L,M,M1,ALPHA,BETA,GAMMA,VLMRE,VLMIM)
         GGLM= GNLMP(N,L,M1,PP)

         S1=S1+VLMRE*GGLM
 10   CONTINUE
      WWLMRE(1)=S1*COS(ALPHA)*SIN(BETA)
      WWLMRE(2)=S1*SIN(ALPHA)*SIN(BETA)
      WWLMRE(3)=S1*COS(BETA)
      RETURN
      END

      SUBROUTINE WWSS(N,L,M,ALPHA,BETA,GAMMA,PP,WWLMIM)
c culculation  W(nlm;alha,beta,gamma) -functions
      DIMENSION WWLMIM(3)

      S1=0.
      DO 10 I=1,L+1
         M1=I-1
         CALL VLMM_RE (L,M,M1,ALPHA,BETA,GAMMA,VLMRE,VLMIM)
         GGLM= GNLMP(N,L,M1,PP)

         S1=S1+VLMIM*GGLM
 10   CONTINUE
      WWLMIM(1)=S1*COS(ALPHA)*SIN(BETA)
      WWLMIM(2)=S1*SIN(ALPHA)*SIN(BETA)
      WWLMIM(3)=S1*COS(BETA)
      RETURN
      END

      FUNCTION GNLMP(N,L,M,PP)

      KK=MOD(L+M,2)

      IF(KK .EQ. 0) THEN 
         GGC= GCOF(N,L,M)   
         PPP=PP**M
         PP2=1-PP**2
         IF(PP2 .GE. 0.) THEN 
            PP32=PP2*SQRT(PP2)
         ELSE
            PP32=0.
         ENDIF   
         ZZ=1-2*PP**2
         LM2=(L-M)/2
         A=M
         B=3./2.
         ACOB=ACOBI(N+LM2,A,B,ZZ)

         GNLMP=GGC*PPP*PP32*ACOB
      ELSE
         GNLMP=0.
      ENDIF   
      RETURN
      END

      SUBROUTINE INT_Proj_GNLM_P(PROJ,PPM,UPPM,NPP,
     #        N1,L1,M1,
     #        ALPHA,BETA,GAMMA,
     #        NJ,I1,I2,I3,NPAR,PAR,YYRE,YYIM)
c   integration over p variable of 
c   g(p,al,bet,gam)*g(n,l,m)
c 

      DIMENSION PPM(NPP),UPPM(NPP)
      DIMENSION GG1M(NPP),GG2m(NPP)
      DIMENSION PROJ(NPP*NJ)

      DIMENSION PROJ2(NPP),PROJ3(NPP)
      DIMENSION NPAR(*),PAR(*)


      NALPHA=NPAR(9)
      NBETA =NPAR(10)
      NGAMMA=NPAR(18)
      HUPP  =PAR(41)

      DO 6 I4=1,NPP
         PP=PPM(I4)
         I321=(I3-1)*NBETA*NALPHA+(I2-1)*NALPHA+I1
         I4321=(I4-1)*NGAMMA*NBETA*NALPHA + I321

         PROJ2(I4)=PROJ(I4321)
 6    CONTINUE   

      DO 8 I4=1,NPP
         UPP=UPPM(I4)
c         PROJ3(I4)= FUNL1(UPPM,PROJ2,1,NPP,UPP,IERR)
         PROJ3(I4)= FUNL1(PPM,PROJ2,1,NPP,UPP,IERR)
 8    CONTINUE   

      S1=0.
      S2=0.
      DO 10 I4=1,NPP
         UPP=UPPM(I4)
         CALL WWCS(N1,L1,M1,ALPHA,BETA,GAMMA,UPP,WWLMRE,WWLMIM)

         WW1=1.0
         IF(I4 .EQ.1 .OR. I4 .EQ. NPP) WW1=0.5
         S1=S1 + WWLMRE * PROJ3(I4)*WW1
         S2=S2 + WWLMIM * PROJ3(I4)*WW1

 10   CONTINUE   
      YYRE= S1*HUPP
      YYIM= S2*HUPP
      RETURN
      END


      SUBROUTINE RightPart_WW(PROJ,LL,NJ,NNMAX,
     #           PPM,UPPM,ALPHAM,BETAM,GAMMAM,
     #           NPP,NALPHA,NBETA,NGAMMA,
     #           NPAR,PAR,YM)

      DIMENSION PROJ(NPP*NALPHA*NBETA*NGAMMA)

      DIMENSION YM((LL+1)*(LL+1)*NNMAX)
      DIMENSION PPM(NPP),UPPM(NPP)
      DIMENSION ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
      DIMENSION NPAR(*),PAR(*)
      DIMENSION BM1(NBETA*NGAMMA),BM2(NBETA)

      HALPHA=PAR(38)
      HBETA =PAR(39)
      HGAMMA=PAR(40)
c      HUPP  =PAR(41)

      LMAX=(LL+1)*(LL+1)


c--------------------------------
c  cos

      DO 30  N=1,NNMAX
      DO 30  L=1,LL+1
         L22=L-1
      DO 30  M=1,L22+1
         IIC=(N-1)*LMAX + (L22*L22+M)


      DO 12 I2=1,NBETA
         BETA=BETAM(I2)

      DO 12 I3=1,NGAMMA
         GAMMA=GAMMAM(I3)

         I23=(I2-1)*NGAMMA+I3

c         WW3=1.0
c         IF(I3 .EQ.1 .OR. I3 .EQ. NGAMMA) WW3=0.5

         S1=0.
      DO 10 I1=1,NALPHA
         ALPHA=ALPHAM(I1)
         WW1=1.0
         IF(I1 .EQ.1 .OR. I1 .EQ. NALPHA) WW1=0.5

         CALL INT_Proj_GNLM_P(PROJ,PPM,UPPM,NPP,
     #        N,L22,M-1,
     #        ALPHA,BETA,GAMMA,    
     #        NJ,I1,I2,I3,NPAR,PAR,YYRE,YYIM)

         S1=S1+YYRE*WW1
 10   CONTINUE
         BM1(I23)=S1*HALPHA
 12   CONTINUE   

      DO 16 I2=1,NBETA
         S1=0.
      DO 14 I3=1,NGAMMA

         I23=(I2-1)*NGAMMA+I3

         WW3=1.0
         IF(I3 .EQ.1 .OR. I3 .EQ. NGAMMA) WW3=0.5
         S1=S1+BM1(I23)*WW3
 14   CONTINUE
         BM2(I2)=S1*HGAMMA
 16   CONTINUE

         S1=0. 
      DO 18 I2=1,NBETA
         BETA=BETAM(I2)
         WW2=1.0
         IF(I2 .EQ.1 .OR. I2 .EQ. NBETA) WW2=0.5
         S1=S1+BM2(I2)*SIN(BETA)*WW2

 18   CONTINUE
      YM(IIC)=S1*HBETA
 30   CONTINUE
c--------------------------------------
c     sin
      DO 50  N=1,NNMAX
      DO 50  L=1,LL+1
         L22=L-1
      DO 50  M=1,L22
         IIS=(N-1)*LMAX + (L22*L22+(M+1)+L22)

      DO 22 I2=1,NBETA
         BETA=BETAM(I2)
      DO 22 I3=1,NGAMMA
         GAMMA=GAMMAM(I3)

         I23=(I2-1)*NGAMMA+I3

         WW3=1.0
         IF(I3 .EQ.1 .OR. I3 .EQ. NGAMMA) WW3=0.5
         S1=0.
      DO 20 I1=1,NALPHA
         ALPHA=ALPHAM(I1)
         WW1=1.0
         IF(I1 .EQ.1 .OR. I1 .EQ. NALPHA) WW1=0.5

         CALL INT_Proj_GNLM_P(PROJ,PPM,UPPM,NPP,
     #        N,L22,M,
     #        ALPHA,BETA,GAMMA,    
     #        NJ,I1,I2,I3,NPAR,PAR,YYRE,YYIM)

         S1=S1+YYIM*WW1
 20   CONTINUE
         BM1(I23)=S1*HALPHA
 22   CONTINUE   

      DO 26 I2=1,NBETA
         BETA=BETAM(I2)
         S1=0.
      DO 24 I3=1,NGAMMA

         I23=(I2-1)*NGAMMA+I3

         WW3=1.0
         IF(I3 .EQ.1 .OR. I3 .EQ. NGAMMA) WW3=0.5
         S1=S1+BM1(I23)*WW3
 24   CONTINUE
         BM2(I2)=S1*HGAMMA
 26   CONTINUE
         S1=0. 
      DO 28 I2=1,NBETA
         BETA=BETAM(I2)
         WW2=1.0
         IF(I2 .EQ.1 .OR. I2 .EQ. NBETA) WW2=0.5
         S1=S1+BM2(I2)*SIN(BETA)*WW2

 28   CONTINUE
      YM(IIS)=S1*HBETA
 50   CONTINUE
      RETURN
      END

      SUBROUTINE Matrix_WW(LL,NNMAX,PPM,UPPM,ALPHAM,BETAM,GAMMAM,
     #                    NPP,NALPHA,NBETA,NGAMMA,NPAR,PAR,AM)

      DIMENSION AM((LL+1)*(LL+1)*NNMAX,(LL+1)*(LL+1)*NNMAX)
      DIMENSION PPM(NPP),ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
      DIMENSION UPPM(NPP)
      DIMENSION NPAR(*),PAR(*)
      DIMENSION VLM(4)


      HALPHA=PAR(38)
      HBETA =PAR(39)
      HGAMMA=PAR(40)
      HUPP  =PAR(41)
    

c      write(*,*) 'Matrix_WW',HALPHA,HBETA,HGAMMA

      LMAX=(LL+1)*(LL+1)

C   cos * cos  (1,1)

      DO 30  N1=1,NNMAX
      DO 30  L1=1,LL+1
         L11=L1-1

      DO 30  M1=1,L11+1
         IIC1=(N1-1)*LMAX + (L11*L11+M1)

      DO 20  N=1,NNMAX
      DO 20  L=1,LL+1
         L22=L-1
      DO 20  M=1,L22+1
         IIC=(N-1)*LMAX + (L22*L22+M)

      CALL WW_WW(L11,M1-1,N1,L22,M-1,N,
     #        PPM,UPPM,NPP,HUPP,
     #        NBETA,NGAMMA,NALPHA,
     #        ALPHAM,BETAM,GAMMAM,NPAR,PAR,VLM)

      AM(IIC,IIC1)=VLM(1)
 20   CONTINUE
 30   CONTINUE
c---------------------------------------------------

C   sin * sin  (2,2)

      DO 50  N1=1,NNMAX
      DO 50  L1=1,LL+1
         L11=L1-1
c        write(*,*) 'Matrix_W: 50 L1=  ',L1
      DO 50  M1=1,L11
         IIS1=(N1-1)*LMAX + (L11*L11+(M1+1)+L11)
      DO 40  N=1,NNMAX
      DO 40  L=1,LL+1
         L22=L-1
      DO 40  M=1,L22
         IIS=(N-1)*LMAX + (L22*L22+(M+1)+L22)

      CALL WW_WW(L11,M1,N1,L22,M,N,
     #        PPM,UPPM,NPP,HUPP,
     #        NBETA,NGAMMA,NALPHA,
     #        ALPHAM,BETAM,GAMMAM,NPAR,PAR,VLM)

        AM(IIS,IIS1)=VLM(4)
 40   CONTINUE
 50   CONTINUE
c-----------------------------------------------

C  cos * sin  (1,2)

      DO 70  N1=1,NNMAX
      DO 70  L1=1,LL+1
         L11=L1-1

      DO 70  M1=1,L11+1
         IIC1=(N1-1)*LMAX + (L11*L11+M1)
      DO 60  N=1,NNMAX
      DO 60  L=1,LL+1
         L22=L-1
      DO 60  M=1,L22
         IIS=(N-1)*LMAX + (L22*L22+(M+1)+L22)

      CALL WW_WW(L11,M1-1,N1,L22,M,N,
     #        PPM,UPPM,NPP,HUPP,
     #        NBETA,NGAMMA,NALPHA,
     #        ALPHAM,BETAM,GAMMAM,NPAR,PAR,VLM)

         AM(IIS,IIC1)=VLM(3)
 60   CONTINUE
 70   CONTINUE
c------------------------------------------

C  sin * cos   (2,1)

      DO 90  N1=1,NNMAX
      DO 90  L1=1,LL+1
         L11=L1-1

      DO 90  M1=1,L11
         IIS1=(N1-1)*LMAX + (L11*L11+(M1+1)+L11)

      DO 80  N=1,NNMAX
      DO 80  L=1,LL+1
         L22=L-1
      DO 80  M=1,L22+1
         IIC=(N-1)*LMAX + (L22*L22+M)

      CALL WW_WW(L11,M1,N1,L22,M-1,N,
     #        PPM,UPPM,NPP,HUPP,
     #        NBETA,NGAMMA,NALPHA,
     #        ALPHAM,BETAM,GAMMAM,NPAR,PAR,VLM)

         AM(IIC,IIS1)=VLM(2)
 80   CONTINUE
 90   CONTINUE
      RETURN
      END


      SUBROUTINE WW_WW(L1,M1,N1,L2,M2,N2,
     #        PPM,UPPM,NPP,HUPP,
     #        NBETA,NGAMMA,NALPHA,
     #        ALPHAM,BETAM,GAMMAM,NPAR,PAR,VLM)

      DIMENSION ALPHAM(*),BETAM(*),GAMMAM(*)
      DIMENSION PPM(*),UPPM(*)
      DIMENSION NPAR(*),PAR(*)
      DIMENSION VLMINT(4)      
      DIMENSION VLM(4)      


      VLM(1)=0.
      VLM(2)=0.
      VLM(3)=0.
      VLM(4)=0.

      DO 8 I1=1,L1+1
         K1=I1-1

         S1=0.
         S2=0.
         S3=0.
         S4=0.
      DO 6 I2=1,L2+1
         K2=I2-1

      CALL INT_GNLMP(PPM,UPPM,NPP,N1,L1,K1,N2,L2,K2,
     %               HUPP,GG12)


      CALL INT_VLMM(L1,M1,K1,L2,M2,K2,ALPHAM,BETAM,GAMMAM,
     %       NBETA,NGAMMA,NALPHA,NPAR,PAR,VLMINT)

      S1=S1+VLMINT(1)*GG12
      S2=S2+VLMINT(2)*GG12
      S3=S3+VLMINT(3)*GG12
      S4=S4+VLMINT(4)*GG12
 6    CONTINUE
      VLM(1)=VLM(1)+S1
      VLM(2)=VLM(2)+S2
      VLM(3)=VLM(3)+S3
      VLM(4)=VLM(4)+S4
 8    CONTINUE
      RETURN
      END

      SUBROUTINE Summa_Harm_Vec(LL,NNMAX,RRM,TETAM,PHIM,
     #                            NRR,NTETA,NPHI,XM,GM)
c  summation of scalar spherical harmonics 
      REAL RRM(NRR),TETAM(NTETA),PHIM(NPHI)
      REAL GM(NRR*NTETA*NPHI,3)
      REAL XM((LL+1)*(LL+1)*NNMAX,3)

c      COMMON/PRIN/IPRNT
c      PARAMETER (PI = 3.1415926536)

      LMAX=(LL+1)*(LL+1)

      DO 18 KK=1,3

      DO 10 I3=1,NRR
         RR=RRM(I3)
      DO 10 I2=1,NTETA
         TETA=TETAM(I2)
      DO 10 I1=1,NPHI
         PHI=PHIM(I1)

         I21=(I2-1)*NPHI+I1
         I321=(I3-1)*NTETA*NPHI+I21
 
        S1=0.
        DO 4  N=1,NNMAX
        DO 4  L=1,LL+1
           L22=L-1
        DO 4  M=1,L22+1
           IIC=(N-1)*LMAX + (L22*L22+M)
           CALL YLM_RE (L22,M-1,TETA,PHI,YLMRE,YLMIM)
           S1=S1+XM(IIC,KK)*YLMRE* HNLR(L22,N,RR)
 4      CONTINUE
        
        DO 6 N=1,NNMAX
        DO 6 L=1,LL+1
           L22=L-1
        DO 6 M=1,L22
           IIS=(N-1)*LMAX + (L22*L22+(M+1)+L22)
           CALL YLM_RE (L22,M,TETA,PHI,YLMRE,YLMIM)
           S1=S1+XM(IIS,KK)*YLMIM* HNLR(L22,N,RR)  
 6      CONTINUE
        GM(I321,KK)=S1
 10   CONTINUE
 18   CONTINUE
      RETURN
      END

      SUBROUTINE Summa_Scal_Harm(LL,NNMAX,RRM,TETAM,PHIM,
     #                            NRR,NTETA,NPHI,XM,GM)
c  summation of scalar spherical harmonics 
      REAL RRM(NRR),TETAM(NTETA),PHIM(NPHI)
      REAL GM(NRR*NTETA*NPHI)
      REAL XM((LL+1)*(LL+1)*NNMAX)

c      COMMON/PRIN/IPRNT
c      PARAMETER (PI = 3.1415926536)

      LMAX=(LL+1)*(LL+1)

      DO 10 I3=1,NRR
         RR=RRM(I3)
      DO 10 I2=1,NTETA
         TETA=TETAM(I2)
      DO 10 I1=1,NPHI
         PHI=PHIM(I1)

         I21=(I2-1)*NPHI+I1
         I321=(I3-1)*NTETA*NPHI+I21
 
        S1=0.
        DO 4  N=1,NNMAX
        DO 4  L=1,LL+1
           L22=L-1
        DO 4  M=1,L22+1
           IIC=(N-1)*LMAX + (L22*L22+M)
           CALL YLM_RE (L22,M-1,TETA,PHI,YLMRE,YLMIM)
           S1=S1+XM(IIC)*YLMRE* HNLR(L22,N,RR)   
 4      CONTINUE
        
        DO 6 N=1,NNMAX
        DO 6 L=1,LL+1
           L22=L-1
        DO 6 M=1,L22
           IIS=(N-1)*LMAX + (L22*L22+(M+1)+L22)
           CALL YLM_RE (L22,M,TETA,PHI,YLMRE,YLMIM)
           S1=S1+XM(IIS)*YLMIM* HNLR(L22,N,RR)  
 6      CONTINUE
        GM(I321)=S1
 10   CONTINUE
      RETURN
      END

      SUBROUTINE INT_VLMM(L1,M1,K1,L2,M2,K2,ALPHAM,BETAM,GAMMAM,
     %  NBETA,NGAMMA,NALPHA,NPAR,PAR,VLMINT)

      DIMENSION ALPHAM(*),BETAM(*),GAMMAM(*)

      DIMENSION NPAR(*),PAR(*)
      DIMENSION BM1(NBETA*NGAMMA,4),VLMINT(4)
      DIMENSION BM2(NBETA,4)

c      NALPHA=NPAR(9)
c      NBETA =NPAR(10)
c      NGAMMA=NPAR(18)

      HALPHA=PAR(38)
      HBETA =PAR(39)
      HGAMMA=PAR(40)
      
      DO 12 I2=1,NBETA
         BETA=BETAM(I2)

      DO 12 I3=1,NGAMMA

         I23=(I2-1)*NGAMMA+I3

         GAMMA=GAMMAM(I3)
         WW3=1.0
         IF(I3 .EQ.1 .OR. I3 .EQ. NGAMMA) WW3=0.5
         S1=0.
         S2=0.
         S3=0.
         S4=0.
      DO 10 I1=1,NALPHA
         ALPHA=ALPHAM(I1)
         WW1=1.0
         IF(I1 .EQ.1 .OR. I1 .EQ. NALPHA) WW1=0.5

         CALL VLMM_RE (L1,M1,K1,ALPHA,BETA,GAMMA,VLMRE1,VLMIM1)
         CALL VLMM_RE (L2,M2,K2,ALPHA,BETA,GAMMA,VLMRE2,VLMIM2)

         VLCC=VLMRE1*VLMRE2
         VLSC=VLMIM1*VLMRE2
         VLCS=VLMRE1*VLMIM2
         VLSS=VLMIM1*VLMIM2

         S1=S1+VLCC*WW1
         S2=S2+VLSC*WW1
         S3=S3+VLCS*WW1
         S4=S4+VLSS*WW1
 10   CONTINUE
         BM1(I23,1)=S1*HALPHA
         BM1(I23,2)=S2*HALPHA
         BM1(I23,3)=S3*HALPHA
         BM1(I23,4)=S4*HALPHA
 12   CONTINUE   

      DO 16 I2=1,NBETA
         BETA=BETAM(I2)
         S1=0.
         S2=0.
         S3=0.
         S4=0.
      DO 14 I3=1,NGAMMA
         I23=(I2-1)*NGAMMA+I3

         WW3=1.0
         IF(I3 .EQ.1 .OR. I3 .EQ. NGAMMA) WW3=0.5
         S1=S1+BM1(I23,1)*WW3
         S2=S2+BM1(I23,2)*WW3
         S3=S3+BM1(I23,3)*WW3
         S4=S4+BM1(I23,4)*WW3
 14   CONTINUE

         BM2(I2,1)=S1*HGAMMA
         BM2(I2,2)=S2*HGAMMA
         BM2(I2,3)=S3*HGAMMA
         BM2(I2,4)=S4*HGAMMA
 16   CONTINUE

         VLMINT(1)=0.
         VLMINT(2)=0.
         VLMINT(3)=0.
         VLMINT(4)=0.
      DO 18 I2=1,NBETA
         BETA=BETAM(I2)
         WW2=1.0
         IF(I2 .EQ.1 .OR. I2 .EQ. NBETA) WW2=0.5
         VLMINT(1)=VLMINT(1)+BM2(I2,1)*SIN(BETA)*HBETA*WW2
         VLMINT(2)=VLMINT(2)+BM2(I2,2)*SIN(BETA)*HBETA*WW2
         VLMINT(3)=VLMINT(3)+BM2(I2,3)*SIN(BETA)*HBETA*WW2
         VLMINT(4)=VLMINT(4)+BM2(I2,4)*SIN(BETA)*HBETA*WW2

 18   CONTINUE
      RETURN
      END


      FUNCTION FUNL1(XM,GM,KK,N1,X1,IERR)
C-         1D line interpolation
      DIMENSION XM(N1)
      DIMENSION GM(N1)
      CHARACTER*8 ITEX
      COMMON /PRIN/IPRNT
C
      ITEX='FUNL1'
      A17=9.E17
C
      IERR=0
      R1=(XM(N1)-XM(1))*1.E-5
      IF(X1-XM(1)+R1)  2,10,10
   2  IERR=1
      CALL ERRPRN(ITEX,IERR)

      write(*,*) 'FUNL1', X1,XM(1),R1
      RETURN
  10  IF(XM(N1)+R1-X1) 4,17,17
   4  IERR=2
      CALL ERRPRN(ITEX,IERR)
      RETURN
  17  CONTINUE
C
      K1=N1-1
      DO 21 I1=1,K1
      I11=I1
      IF(X1-XM(I1+1)) 27,21,21
  21  CONTINUE
C
  27  CONTINUE
C
      IF(GM(I11)+A17) 28,29,29
  28  IERR=5
      CALL ERRPRN(ITEX,IERR)
C
  29  IF(GM(I11+1)+A17) 30,42,42
  30  IERR=6
      CALL ERRPRN(ITEX,IERR)
  42  CONTINUE
C
      J1=(KK-1)*N1+I11
      J2=J1+1
      HX1=X1-XM(I11)
C     HX2=XM(I11+1)-X1
      HX=XM(I11+1)-XM(I11)
C
      FY1=GM(J1)+(GM(J2)-GM(J1))/HX*HX1
C
      FUNL1=FY1
      RETURN
      END


      SUBROUTINE Summa_YLM_Harm_W(LL,NNMAX,RRM,TETAM,PHIM,
     #                            NRR,NTETA,NPHI,XM,GM)
c  summation of scalar spherical harmonics 
      REAL RRM(NRR),TETAM(NTETA),PHIM(NPHI)
      REAL GM(NRR*NTETA*NPHI)
      REAL XM((LL+1)*(LL+1)*NNMAX)

c      COMMON/PRIN/IPRNT
      PARAMETER (PI = 3.1415926536)

c      CALL IntOverSher(GM,GM2,RRM,TETAM,PHIM,NRR,NTETA,NPHI,
c     #      HPHI,HTETA)


      LMAX=(LL+1)*(LL+1)

      DO 10 I3=1,NRR
         RR=RRM(I3)
      DO 10 I2=1,NTETA
         TETA=TETAM(I2)
      DO 10 I1=1,NPHI
         PHI=PHIM(I1)

         I21=(I2-1)*NPHI+I1
         I321=(I3-1)*NTETA*NPHI+I21
 
        S1=0.
        DO 4  N=1,NNMAX
        DO 4  L=1,LL+1
           L1=L-1
        DO 4  M=1,L1+1
           IIC=(N-1)*LMAX + (L1*L1+M)
           CALL YLM_RE (L1,M-1,TETA,PHI,YLMRE,YLMIM)
           S1=S1+XM(IIC)*YLMRE * HNLR(L1,N,RR)   
 4      CONTINUE
        
        DO 6 N=1,NNMAX
        DO 6 L=1,LL+1
           L1=L-1
        DO 6 M=1,L1
           IIS=(N-1)*LMAX + (L1*L1+(M+1)+L1)
           CALL YLM_RE (L1,M,TETA,PHI,YLMRE,YLMIM)
           S1=S1+XM(IIS)*YLMIM * HNLR(L1,N,RR)  
 6      CONTINUE
        GM(I321)=S1
 10   CONTINUE
      RETURN
      END

      SUBROUTINE XXMM(LL,NNMAX,XM)
c  coefficients of decomposition are written as 
c  vector XM 

      DIMENSION XM((LL+1)*(LL+1)*NNMAX,3)

      LMAX=(LL+1)*(LL+1)

      DO 10  N=1,NNMAX
      DO 10  L=1,LL+1
         L22=L-1
      DO 10  M=1,L22+1
         IIC=(N-1)*LMAX + (L22*L22+M)
c         write(*,*) 'XXMM-cos',L,M,IIC
         AAC1=AACOS1(L22,M-1,N)
         AAC2=AACOS2(L22,M-1,N)
         AAC3=AACOS3(L22,M-1,N)
         XM(IIC,1)=AAC1
         XM(IIC,2)=AAC2
         XM(IIC,3)=AAC3

 10   CONTINUE

      DO 20  N=1,NNMAX
      DO 20  L=1,LL+1
         L22=L-1
      DO 20  M=1,L22
         IIS=(N-1)*LMAX + (L22*L22+(M+1)+L22)
c         write(*,*) 'XXMM-sin',L,M,IIS
         AAS1=AASIN1(L22,M,N)
         AAS2=AASIN2(L22,M,N)
         AAS3=AASIN3(L22,M,N)
         XM(IIS,1)=AAS1
         XM(IIS,2)=AAS2
         XM(IIS,3)=AAS3
 20   CONTINUE
      RETURN
      END

      FUNCTION AACOS1(L,M,N)
      AACOS1=1.
      IF(L .EQ.0 .AND. M .EQ.0) AACOS1=0.5
      RETURN
      END

      FUNCTION AASIN1(L,M,N)
      AASIN1=-1.
      RETURN
      END

      FUNCTION AACOS2(L,M,N)
      AACOS2=1.   
c      AACOS2=0.
      IF(L .EQ.0 .AND. M .EQ.0) AACOS2=0.     !  0.1
      RETURN
      END

      FUNCTION AASIN2(L,M,N)
      AASIN2=-0.1 
c      AASIN2=0.
      RETURN
      END

      FUNCTION AACOS3(L,M,N)
      AACOS3=0.5  
c      AACOS3=0.
      IF(L .EQ.0 .AND. M .EQ.0) AACOS3=0.2  
      RETURN
      END

      FUNCTION AASIN3(L,M,N)
      AASIN3=-0.5 
c      AASIN3=0.
      RETURN
      END

      SUBROUTINE Summa_YLM_Harm(LL,NNMAX,RRM,TETAM,PHIM,
     #                            NRR,NTETA,NPHI,XM,GM)
c  summation of scalar spherical harmonics 
      REAL RRM(NRR),TETAM(NTETA),PHIM(NPHI)
      REAL GM(NRR*NTETA*NPHI)
      REAL XM(LL*(LL+2)*NNMAX)

c      COMMON/PRIN/IPRNT
      PARAMETER (PI = 3.1415926536)

      LMAX=(LL+1)*(LL+1)

      DO 10 I3=1,NRR
         RR=RRM(I3)
      DO 10 I2=1,NTETA
         TETA=TETAM(I2)
      DO 10 I1=1,NPHI
         PHI=PHIM(I1)

         I21=(I2-1)*NPHI+I1
         I321=(I3-1)*NTETA*NPHI+I21
 
        S1=0.
        DO 4  N=1,NNMAX
        DO 4  L=1,LL+1
           L1=L-1
        DO 4  M=1,L1+1
           IIC=(N-1)*LMAX + (L1*L1+M)
           CALL YLM_RE (L1,M-1,TETA,PHI,YLMRE,YLMIM)
           S1=S1+XM(IIC)*YLMRE * HNLR(L1,N,RR)   
 4      CONTINUE
        
        DO 6 N=1,NNMAX
        DO 6 L=1,LL+1
           L1=L-1
        DO 6 M=1,L1
           IIS=(N-1)*LMAX + (L1*L1+(M+1)+L1)
           CALL YLM_RE (L1,M,TETA,PHI,YLMRE,YLMIM)
           S1=S1+XM(IIS)*YLMIM * HNLR(L1,N,RR)  
 6      CONTINUE
        GM(I321)=S1
 10   CONTINUE
      RETURN
      END

      FUNCTION HNLR(L,N,RR)
      COMMON/PRIN/ IPRNT

c      IF(RR .NE. 0 .AND. RR .NE. 1) GOTO 12
c         WRITE(IPRNT,*) 'Procedure HNLR ERROR: RR=0 or 1 '
c         RETURN
c 12   CONTINUE   
c  AL=2.5
      Z=1-2*RR**2
      A=L+1./2.
      B=1.
      ACOB=ACOBI(N,A,B,Z)
      HNLR=RR**L*(1.0-RR**2)*ACOB

c  AL=3.5
c      Z=1-2*RR**2
c      A=L+1./2.
c      B=2.
c      ACOB=ACOBI(N,A,B,Z)
c      HNLR=RR**L*(1.0-RR**2)**2*ACOB

c  AL=4.5
c      Z=1-2*RR**2
c      A=L+1./2.
c      B=3.
c      ACOB=ACOBI(N,A,B,Z)
c      HNLR=RR**L*(1.0-RR**2)**3*ACOB
c------------------------------------
c      A=0.
c      B=L+1./2.
c      Z=2*RR**2-1.
c      ACOB=ACOBI(N,A,B,Z)
c      CC=SQRT(4.0*N+2.0*L+3.0)
c      HNLR=CC*RR**L*ACOB
c-------------------------------------
C   Laguerra
c      Z=RR**2
c      LAG=1+L+0.5-Z
c     HNLR=RR**L*EXP(-Z)*LAG


      RETURN
      END

      SUBROUTINE POLAR(X,Y,Z,TET,FI)
C procedure perfom a passage to spherical coordinates TET, FI
C  by use decart X,Y,Z coordinate
C SIGNY,SIGNZ- signs of  X,Y,Z
C TET,FI- spherical coordinates
      PARAMETER (PI = 3.1415926536)

      SIGNY=SIGN(1.,Y)
      SIGNZ=SIGN(1.,Z)
      IF(X.NE.0)GO TO 30
      FI=PI*0.5*SIGNY
      GO TO 35
   30 FI=ATAN2(Y,X)
   35 CONTINUE
      IF(Z.NE.0)GO TO 40
      TET=PI*0.5*SIGNZ
      GO TO 45
   40 TET=ATAN2(SQRT(X*X+Y*Y),Z)
   45 CONTINUE
      RETURN
      END

      SUBROUTINE REGN2(XM,YM,NX,NY,X1,Y1,EI)
      DIMENSION XM(*),YM(*)
      IF(X1.LT.XM(1).OR.X1.GT.XM(NX)) GO TO 10
      IF(Y1.LT.YM(1).OR.Y1.GT.YM(NY)) GO TO 10
c      IF(Z1.LT.ZM(1).OR.Z1.GT.ZM(NZ)) GO TO 10
      EI=1.
      RETURN
  10  EI=-1.
      RETURN
      END

      SUBROUTINE MODEL1(GM,XM,YM,ZM,NX,NY,NZ)
      DIMENSION GM(NX*NY*NZ)
      DIMENSION XM(NX),YM(NY),ZM(NZ)


      AA=0.3        !0.2
      BB=0.6        !0.3
      CC=0.3        !0.4
      DO 4 I3=1,NZ
         ZZ=ZM(I3)
      DO 4 I2=1,NY
         YY=YM(I2)
      DO 4 I1=1,NX
         XX=XM(I1)
         I321=(I3-1)*NX*NY+(I2-1)*NX +I1
         GM(I321)=EXP(-((XX/AA)**2+(YY/BB)**2+(ZZ/CC)**2))
 4    CONTINUE   
      RETURN
      END

      FUNCTION GCOF(N,L,M)
      PARAMETER (PI = 3.1415926536)

      SIGNUM=1
      KK=MOD(L,2)
      IF(KK .NE. 0) SIGNUM=-1.

      N11=N+(L-M)/2
      A11=FACTRL(N11)
c      write(*,*) 'A11=',A11
c---------------------------
      A22=FACTRL(L-M)
      A33=FACTRL(L+M)

c      write(*,*) 'A22=',A22
c      write(*,*) 'A33=',A33
c-------------------------
      N22=(L-M)/2+N+2
      S=1.
      DO 4 I=1, 2*N22-1,2
         S=S*I
 4    CONTINUE   
c      write(*,*) 'S=',S
      GN22=SQRT(PI)*S/(2**N22)
c      write(*,*) 'GN22=',GN22

c--------------------------------
      BB=(N+1)*A11*SQRT((2*L+1)*A22*A33)
c            write(*,*) 'BB=',BB
      A44=FACTRL((L-M)/2)
      A55=FACTRL((L+M)/2)
      A66=2**(L+1)
c---------------------------------
c            write(*,*) 'A44=',A44
c            write(*,*) 'A55=',A55
c            write(*,*) 'A66=',A66

      GCOF=SIGNUM*BB/A44/A55/A66/GN22

      RETURN
      END

      SUBROUTINE Summa_W_Harm(LL,NJ,NNMAX,PPM,ALPHAM,BETAM,
     #                   GAMMAM,NPP,NALPHA,NBETA,NGAMMA,XM,GM)
c  summation of W  harmonics 
      REAL PPM(NPP),ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
      REAL GM(NPP*NALPHA*NBETA*NGAMMA)
      REAL XM(LL*(LL+2)*NNMAX)

c      COMMON/PRIN/IPRNT
c      PARAMETER (PI = 3.1415926536)

      LMAX=LL*(LL+2)

c      DO 10 I4=1,NPP
c         PP=PPM(I4)
c      DO 10 I3=1,NGAMMA
c         GAMMA=GAMMAM(I3)
c      DO 10 I2=1,NBETA
c         BETA=BETAM(I2)
c      DO 10 I1=1,NALPHA
c         ALPHA=ALPHAM(I1)

      DO 10 I4=1,NPP
         PP=PPM(I4)
      DO 10 II=1,NJ
         CALL IALBETGAM(II,NALPHA,NBETA,I1,I2,I3)

         I321=(I3-1)*NBETA*NALPHA+(I2-1)*NALPHA+I1
         I4321=(I4-1)*NGAMMA*NBETA*NALPHA + I321
 
         ALPHA=ALPHAM(I1)
         BETA =BETAM(I2)
         GAMMA=GAMMAM(I3)

        S1=0.
        DO 4  N=1,NNMAX
        DO 4  L=1,LL
        DO 4  M=1,L+1
           IIC=(N-1)*LMAX + (L*L+M-1)
           CALL WWCS(N,L,M-1,ALPHA,BETA,GAMMA,PP,WWLMRE,WWLMIM)

c           CALL YLM_RE (L,M-1,TETA,PHI,YLMRE,YLMIM)
           S1=S1+XM(IIC)*WWLMRE
 4      CONTINUE
        
        DO 6 N=1,NNMAX
        DO 6 L=1,LL
        DO 6 M=1,L
           IIS=(N-1)*LMAX + (L*L+M+L)
           CALL WWCS(N,L,M,ALPHA,BETA,GAMMA,PP,WWLMRE,WWLMIM)

c           CALL YLM_RE (L,M,TETA,PHI,YLMRE,YLMIM)
           S1=S1+XM(IIS)*WWLMIM
 6      CONTINUE
        GM(I4321)=S1
 10   CONTINUE
      RETURN
      END


      FUNCTION PLM1M2(L,M1,M2,BET)
C   Varshalovich p.77 eq.(5)
      REAL BET
      REAL Y1,Y2,Y3,Y4
      PARAMETER (PI=3.1415926)

      CBET=COS(BET)
      KK1=MAX(0, M2-M1)
      KK2=MIN(L-M1, L+M2)
      Y1=SIGNUM(M1-M2)*SQRT(FACTRL(L+M1)*FACTRL(L-M1)) *
     #   SQRT(FACTRL(L+M2)*FACTRL(L-M2))
      KKK1=MAX(KK1,KK2)
      KKK2=MIN(KK1,KK2)

      SS=0.
      DO 4 K=KKK2, KKK1
      Y2=SIN(BET/2.)**(M1-M2+2*K)*COS(BET/2.)**(2*L+M2-M1-2*K)
      Y3=SIGNUM(K)  
      Y4=FACTRL(L+M2-K)*FACTRL(M1-M2+K)*FACTRL(L-M1-K)*FACTRL(K)
      SS=SS+Y2*Y3/Y4
 4    CONTINUE
      PLM1M2=SS*Y1
      RETURN
      END


      SUBROUTINE DATA(ALPHAM,BETAM,GAMMAM,XM,YM,ZM,
     #               NKAPPA,BKAPPAM,TETAM,PHIM,
     #               BETKAPM, UM,VM,WM,RRM,PPM,UPPM,NPAR,PAR)

      DIMENSION ALPHAM(*),BETAM(*)
      DIMENSION GAMMAM(*)
      DIMENSION PPM(*),UPPM(*)
      DIMENSION TETAM(*),PHIM(*),RRM(*)




      DIMENSION XM(*),YM(*),ZM(*)
      DIMENSION UM(*),VM(*),WM(*)

      DIMENSION PAR(*),NPAR(*)

      DIMENSION BETKAPM(NKAPPA)
      DIMENSION BKAPPAM(NKAPPA)

      PARAMETER ( PI = 3.141592653)




      DD = PAR(36)


      RRE= PAR(37)
      RCC= PAR(42)
      UM1= PAR(43)

      NU    =NPAR(4)
      NV    =NPAR(5)
      NW    =NPAR(6)

      NPHI  =NPAR(7)
      NTETA  =NPAR(8)

      NALPHA=NPAR(9)
      NBETA =NPAR(10)
      NGAMMA=NPAR(18)
c      NKAPPA=NPAR(20)

      NJ   =NPAR(11)
      NPP   =NPAR(19)
      NX    =NPAR(12)
      NY    =NPAR(13)
      NZ    =NPAR(14)
      NRR   =NPAR(21)


      BETA0=0.
      BETA1=PI
      HBETA=(BETA1-BETA0)/(NBETA-1)
c      HBETA=(BETA1-BETA0)/NBETA

c      DO 10 I=1,NJ
cc      DO 10 J=1,NKAPPA
c         BETAM(I)=BETA0+(I-1)*HBETA
c 10   CONTINUE

      ALPHA0=-PI
      ALPHA1= PI
      HALPHA=(ALPHA1-ALPHA0)/(NALPHA-1)
c      HALPHA=(ALPHA1-ALPHA0)/NALPHA

      GAMMA0=0.
      GAMMA1=2.*PI
      HGAMMA=(GAMMA1-GAMMA0)/(NGAMMA-1)
c      HGAMMA=(GAMMA1-GAMMA0)/NGAMMA

      X0=-1.
      X1= 1.
      HX=(X1-X0)/(NX-1)

      Y0=-1.
      Y1= 1.
      HY=(Y1-Y0)/(NY-1)

      Z0=-1.
      Z1= 1.
      HZ=(Z1-Z0)/(NZ-1)

      U0=RCC  !  -1.
      U1= 1.
      HU=(U1-U0)/(NU-1)
      
      V0=-1.
      V1= 1.
      HV=(V1-V0)/(NV-1)
 
      W0=-1.
      W1= 1.
      HW=(W1-W0)/(NW-1)

      RR0=0.   !0.01
      RR1=1.   !0.99
      HRR=(RR1-RR0)/(NRR-1)
c      HRR=(RR1-RR0)/NRR

      DO 24 I=1,NRR
         RRM(I)=RR0+(I-1)*HRR
 24   CONTINUE


      PP0=RR0      !0.01
      PP1=RR1         !0.99
      HPP=(PP1-PP0)/(NPP-1)
c      HPP=(PP1-PP0)/NPP

c      BKAPPA0 = ATAN2(RRE,DD)
      BKAPPA0 = ASIN(RRE/DD)
      BKAPPA1 = ASIN(RCC/DD)  ! угол затенения
      UM1 = DD*SIN(BKAPPA1)
c----------------------
      HKAPPA=2.*BKAPPA0/(NKAPPA-1)
c      HKAPPA=(BKAPPA0-BKAPPA1)/(NKAPPA-1)

        DO 3 I=1,NKAPPA
         BKAPPAM(I)=-BKAPPA0+(I-1)*HKAPPA
c         BKAPPAM(I)=(I-1)*HKAPPA
 3    CONTINUE
c---------------------------
      DO 25 I=1,NPP
         BKAPPA=BKAPPAM(I)
c         PPM(I)=PP0+(I-1)*HPP
         PPM(I)=DD*SIN(BKAPPA)
         
cc         VM(I)=PPM(I)   !  VM nije zacomentirival
         UM(I)=PPM(I)   !  UM nije zacomentirival

 25   CONTINUE
  

      UPP0=PPM(1)
      UPPN=PPM(NPP)
      HUPP=(UPPN-UPP0)/(NPP-1)
      DO 26 I=1,NPP
         UPPM(I)=UPP0+(I-1)*HUPP
 26   CONTINUE



      PI100=PI/100.
      TETA0=0.    !PI100
      TETA1=PI    !PI-PI100
      HTETA=(TETA1-TETA0)/(NTETA-1)
      
      PHI0=-PI
      PHI1=PI
      HPHI=(PHI1-PHI0)/(NPHI-1)
c      HPHI=(PHI1-PHI0)/NPHI



      DO 4 I=1,NTETA
         TETAM(I)=TETA0+(I-1)*HTETA
 4    CONTINUE

      DO 6 I=1,NPHI
         PHIM(I)=PHI0+(I-1)*HPHI
 6    CONTINUE

      DO 8 I=1,NALPHA
         ALPHAM(I)=ALPHA0+(I-1)*HALPHA
 8    CONTINUE

      DO 10 I=1,NBETA
         BETAM(I)=BETA0+(I-1)*HBETA
 10   CONTINUE

      DO 11 I=1,NGAMMA
         GAMMAM(I)=GAMMA0+(I-1)*HGAMMA
 11   CONTINUE

      DO 12 I=1,NX
         XM(I)=X0+(I-1)*HX
 12   CONTINUE

      DO 14 I=1,NY
         YM(I)=Y0+(I-1)*HY
 14   CONTINUE

      DO 16 I=1,NZ
         ZM(I)=Z0+(I-1)*HZ
 16   CONTINUE

c      DO 18 I=1,NU
c         UM(I)=U0+(I-1)*HU
c 18   CONTINUE

      DO 20 I=1,NV
         VM(I)=V0+(I-1)*HV
 20   CONTINUE

      DO 22 I=1,NW
         WM(I)=W0+(I-1)*HW
 22   CONTINUE

c      NPAR(4) =NU
c      NPAR(5) =NV
c      NPAR(6) =NW

c      NPAR(7) =NPHI
c      NPAR(8) =NTET
c      NPAR(9) =NALPHA
c      NPAR(10)=NBETA

c      NPAR(12)=NX
c      NPAR(13)=NY
c      NPAR(14)=NZ

c------------------------
      PAR(1)=U0
      PAR(2)=U1
      PAR(3)=HU
      PAR(4)=V0
      PAR(5)=V1
      PAR(6)=HV
      PAR(7)=W0
      PAR(8)=W1
      PAR(9)=HW

c      PAR(10)=PHI0
c      PAR(11)=PHI1
      PAR(12)=HPHI

c      PAR(13)=TET0
c      PAR(14)=TET1
      PAR(15)=HTETA

      PAR(16)=X0
      PAR(17)=X1
      PAR(18)=HX

      PAR(19)=Y0
      PAR(20)=Y1
      PAR(21)=HY

      PAR(22)=Z0
      PAR(23)=Z1
      PAR(24)=HZ

      PAR(33)=RR0
      PAR(34)=RR1
      PAR(35)=HRR
      PAR(38)=HALPHA
      PAR(39)=HBETA
      PAR(40)=HGAMMA
      PAR(41)=HUPP

      RETURN
      END

      SUBROUTINE YLM_RE (L,M,TETA,PHI,YLMRE,YLMIM)
C          Computer the spherical harmonics Y(L,M)
C          YLMRE - real part
C          YLMIM - image part
      PARAMETER (PI=3.1415926)
      INTEGER L,M
      REAL TETA,PHI
      REAL YLMRE,YLMIM

      MM=ABS(M)
      KK=MOD(MM,2)
      SIGNUM=1.
      IF(M .LT. 0 .AND. KK .EQ. 1) SIGNUM=-1.
      IF(M .LT. 0 .AND. KK .EQ. 0) SIGNUM=1.
      X=COS(TETA) 

      PNN=PNORM(L,MM)
      PLAG=PLGNDR(L,MM,X)
      YLMRE= SIGNUM*PNN*PLAG*COS(MM*PHI)
      YLMIM= SIGNUM*PNN*PLAG*SIN(MM*PHI)

c      CYLM_1=CMPLX(YLMRE,YLMIM)
      RETURN
      END


      FUNCTION PNORM(L,MM)
C           Computers the norm of associated Legendre polynomial 
      INTEGER L,M
      INTEGER K
      REAL S
      PARAMETER (PI=3.1415926)
      M=ABS(MM)   !!!!!!!! vnes ispravlenie 11 August 2004
      K=L-M
      S=1. 
      DO 4 I=1, 2*M
      J=K+I 
      S=S/FLOAT(J)
 4    CONTINUE
      PNORM=SQRT((2*L+1)*S/(4.*PI))
      RETURN
      END
      FUNCTION PLGNDR(L,M,X)
C                                                     m
C       computers the associated Legendre polinomial P (x). 
C                                                     l
C       Here m and l are integer satisfying  0 <= m <= l, while x 
C       -1 <= x <=1

      INTEGER L,M
      REAL PLGNDR,X
      INTEGER I,LL
      REAL FACT, PLL,PMM,PMM1,SOMX2
      IF(M .LT. 0 .OR. ABS(X) .GT. 1.) PAUSE 'BAD ARGUMENTS IN PLGNDR'
      PMM=1.                    ! Computer Pmm      
      IF(M .GT. 0) THEN
            SOMX2=SQRT((1.-X)*(1.+X))
        FACT=1.
        DO 11 I=1,M
            PMM=-PMM*FACT*SOMX2
            FACT=FACT+2
 11     CONTINUE
      ENDIF
      IF(L .EQ. M) THEN
            PLGNDR=PMM
      ELSE
            PMMP1=X*(2*M+1)*PMM      !Compute Pm+1,m     
      IF(L .EQ. M+1) THEN
            PLGNDR=PMMP1
      ELSE
       DO 12 LL=M+2,L
            PLL=(X*(2*LL-1)*PMMP1-(LL+M-1)*PMM)/(LL-M)
            PMM=PMMP1
            PMMP1=PLL
 12    CONTINUE
      PLGNDR=PLL              !*PNORM(L,M)
      ENDIF
      ENDIF
      RETURN
      END  

      SUBROUTINE XYANG(am,gm,fm1,fm2,a1m,n1,n2,niter)
cl... am - 
cl....gm - 
cl ...fm1 -
cl ...a1m -
      DIMENSION am(*),gm(*),fm1(*),fm2(*),a1m(*)
      CHARACTER*8 itex
      INTEGER IPRNT

      COMMON/PRIN/ IPRNT
     
      CC=1.
      itex = 'xyang'
      WRITE(IPRNT,*) 'number iterations:',niter
      ierr = 0
      DO 40 i = 1,niter
      WRITE(IPRNT,*) 'iteration I=:',i
      DO 30 i1 = 1,n1
      DO 10 i2 = 1,n2
      j = (i2-1)*n1+i1
      a1m(i2) = am(j)
10    CONTINUE
      CALL scal(fm1,a1m,n2,q1)
      CALL scal(a1m,a1m,n2,q2)
      IF(q2 .NE. 0) GO TO 12
      WRITE(IPRNT,*) 'xyang: error ',i1
      RETURN
12    CONTINUE
      DO 14 i2 = 1,n2
      fm2 (i2) = fm1 (i2) - CC*(q1-gm(i1))/q2*a1m(i2)
C      IF(fm2(i2) .LT. 0.) fm2(i2) = 0.
14    fm1(i2) = fm2(i2)
30    CONTINUE
40    CONTINUE
      RETURN
      END   

      SUBROUTINE scal(a,b,n,q)
c...     scale production 
       DIMENSION a(*),b(*)
      s = 0.
      DO 2 i = 1 , n
2     s = s + a(i) * b(i)
      q = s
      RETURN
      END

      FUNCTION SIGNUM (LA)
      NN=ABS(MOD (LA,2))
        SIGNUM=1.
        IF(NN .EQ. 1) SIGNUM=-1.
      RETURN
      END

      FUNCTION gammln(xx)
      REAL gammln,xx
C   return the value lnG(xx)
      INTEGER j
      DOUBLE PRECISION ser, stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp /76.18009172947146d0, -86.50532032941677d0,
     # 24.01409824083091d0, -1.231739572450155d0, .1208650973866179d-2,
     # -.5395239384953d-5, 2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      DO 11 j=1,6
        y=y+1.d0
        ser=ser+cof(6)/y
 11   CONTINUE
      gammln=tmp+log(stp+ser/x)
      RETURN
      END

      FUNCTION FACTRL(n)
      INTEGER n
      REAL factrl
C return the value n! as a floating-point number
      INTEGER j,ntop
      REAL a(33),gammln
      SAVE ntop,a
      DATA ntop, a(1) /0,1./   !table initialized with 0! only
      if(n .lt. 0) then
          pause 'negativefactorial in "factrl" '
      else if (n .le. ntop) then
              factrl=a(n+1)
      else if (n .le. 32) then
        DO 11 j=ntop+1,n
         a(j+1)=j*a(j)
 11     CONTINUE
      ntop=n
         FACTRL=a(n+1)
      else
         FACTRL= exp(gammln(n+1.))
      endif
      RETURN
      END 

      FUNCTION ACOBI(N,A,B,Z)
C procedure calculate Yacobi polinome at point  Z
C by summation Hiper-Geometrical series
C look Beitmen p.172
C Z=COS(BET),BET- an angle of rotation over new X axis
C look Nikiforov p.122.
c      INTEGER A,B
      REAL A,B
      INTEGER RN,RI
      RN=N
      P=1.
      DO 1I=1,N
      RI=I
    1 P=P*(RN+A+1.-RI)/RI
      S1=1.
      S=1.
      X=(1.-Z)/2.
      DO 2 I=1,N
      RI=I
      S1=S1*(-RN+RI-1.)*(RN+A+B+RI)/(A+RI)*X/RI
      S=S+S1
    2 CONTINUE
      ACOBI=P*S
      RETURN
      END

      COMPLEX FUNCTION CYLM (L,M,TETA,PHI)
C          Computation the spherical harmonics L, M
C          YLMRE - real part
C          YLMIM - image part
      PARAMETER (PI=3.1415926)
      INTEGER L,M
      REAL TETA,PHI
      REAL YLMRE,YLMIM
      COMPLEX CCEXP
      CCEXP(X)=CMPLX(COS(X),SIN(X))

      MM=ABS(M)
      KK=MOD(MM,2)
      X=COS(TETA) 

      PNN=PNORM(L,MM)
      PLAG=PLGNDR(L,MM,X)

      IF(M .LT. 0) THEN 
           IF(KK .EQ. 1) THEN
              SIGNUM=-1.
           ELSE
              SIGNUM=1.
           ENDIF  
        CYLM= SIGNUM*PNN*PLAG*CONJG(CCEXP(MM*PHI))
      ELSE
        CYLM= PNN*PLAG*CCEXP(MM*PHI)  
      ENDIF

      RETURN
      END

      SUBROUTINE GNUFORM3 (NX,NY,KK,GC,KAN)
      DIMENSION GC(*)
      DO 10 I3=KK,KK
      DO 10 I2=1,NY
      WRITE(KAN,100) (GC((I3-1)*NX*NY+(I2-1)*NX),I1=1,NX)
      WRITE(KAN, *)
C      WRITE(KAN,'(A1)') char(13)//char(10)
10    CONTINUE 
100   FORMAT(E10.3)
C110   FORMAT(A4)
      RETURN
      END

      SUBROUTINE GNUFORM (NX,NY,KK,GC,KAN)
      DIMENSION GC(*)
      DO 10 J=KK,NY
      WRITE(KAN,100) (GC((J-1)*NX+I),I=1,NX)
      WRITE(KAN, *)
C      WRITE(KAN,'(A1)') char(13)//char(10)
10    CONTINUE 
100   FORMAT(E10.3)
C110   FORMAT(A4)
      RETURN
      END

      SUBROUTINE GNUFORM1 (NX,GC,KAN)
      DIMENSION GC(*)
      DO 10 I=1,NX
      WRITE(KAN,100) I,GC(I)
c      WRITE(KAN, *)
C      WRITE(KAN,'(A1)') char(13)//char(10)
10    CONTINUE 
100   FORMAT(I3,E12.3)
C110   FORMAT(A4)
      RETURN
      END

      SUBROUTINE IDLFORM (NX,GC,KAN)
      DIMENSION GC(*)
      DO 10 I=1,NX
      WRITE(KAN,100) GC(I)
10    CONTINUE 
100   FORMAT(5(1x,E10.3))
      RETURN
      END

      SUBROUTINE IDLFORM1 (NX,GC,KAN)
      DIMENSION GC(*)
      DO 10 I=1,NX
      WRITE(KAN,100) GC(I)
10    CONTINUE 
100   FORMAT(E12.3)
      RETURN
      END

      SUBROUTINE GNUFORM1_1 (NX,KK,GC,KAN)
      DIMENSION GC(NX,3)

c      write(*,*) 'GNUFORM1_1',(GC(I,3),I=1,NX)
      DO 10 I=1,NX
      WRITE(KAN,100) I,GC(I,KK)
c      WRITE(KAN, *)
C      WRITE(KAN,'(A1)') char(13)//char(10)
10    CONTINUE 
100   FORMAT(I3,E12.3)
C110   FORMAT(A4)
      RETURN
      END

      SUBROUTINE GNUFORM_X_Y (XM,YM,NX,NY,GC,KAN1)
      DIMENSION XM(*),YM(*)
      REAL GC(NX*NY)

      DO 12 I2=1,NY
         YY=YM(I2)
      DO 10 I1=1,NX
         XX=XM(I1)
         I12= (I2-1)*NX+I1

      WRITE(KAN1,100) XX,YY,GC(I12)

c      WRITE(KAN1,100)
10    CONTINUE 
      WRITE(17,'(A1)') char(13)//char(10)
 12   CONTINUE
100   FORMAT(4E15.3)
      RETURN
      END


      SUBROUTINE IALBETGAM(JJ,N1,N2,I1,I2,I3)
C calculation I1, I2, I3--  indexes by the number of  JJ,
C where JJ is calculated as  JJ=(I3-1)*N1*N2 + (I2-1)*N1 + I1
      
      RJJ=JJ
      RN12=FLOAT(N1*N2)
      RMODUL=AMOD(RJJ,RN12)
      IC=1
      IF(RMODUL.EQ.0) IC=0
      I3=INT(RJJ/RN12)+IC
      I12=JJ-(I3-1)*N1*N2
     
      RJJ=I12
      RN1=FLOAT(N1)
      RMODUL=AMOD(RJJ,RN1)
      IC=1
      IF(RMODUL.EQ.0) IC=0
      I2=INT(RJJ/N1)+IC
      I1=I12-(I2-1)*N1
      RETURN
      END

      SUBROUTINE REGN(XM,YM,ZM,NX,NY,NZ,X1,Y1,Z1,EI)
      DIMENSION XM(1),YM(1),ZM(1)
      IF(X1.LT.XM(1).OR.X1.GT.XM(NX)) GO TO 10
      IF(Y1.LT.YM(1).OR.Y1.GT.YM(NY)) GO TO 10
      IF(Z1.LT.ZM(1).OR.Z1.GT.ZM(NZ)) GO TO 10
      EI=1.
      RETURN
  10  EI=-1.
      RETURN
      END

      FUNCTION FUNL3(XM,YM,ZM,GM,N1,N2,N3,X1,Y1,Z1,IERR)
CL- 3D interpolation 

      DIMENSION XM(*),YM(*),ZM(*)
      DIMENSION GM(*)
C
      COMMON/PRIN/IPRNT
      CHARACTER*8 ITEX
C
      ITEX='FUNL3'
      IERR=0
C
      A17=9.E17
C
      R1=(XM(N1)-XM(1))*1.E-5
      R2=(YM(N2)-YM(1))*1.E-5
      R3=(ZM(N3)-ZM(1))*1.E-5
C
      IF(X1-XM(1)+R1)  2,10,10
   2  IERR=1
      CALL ERRPRN(ITEX,IERR)
      RETURN
C
  10  IF(XM(N1)+R1-X1) 4,11,11
   4  IERR=2
      CALL ERRPRN(ITEX,IERR)
      RETURN
C
  11  IF(Y1-YM(1)+R2)  6,12,12
   6  IERR=3
      CALL ERRPRN(ITEX,IERR)
      RETURN
C
  12  IF(YM(N2)+R2-Y1) 8,13,13
   8  IERR=4
      CALL ERRPRN(ITEX,IERR)
      RETURN
C
  13  IF(Z1-ZM(1)+R3)  14,15,15
  14  IERR=5
      CALL ERRPRN(ITEX,IERR)
      RETURN
C
  15  IF(ZM(N3)+R3-Z1) 17,18,18
  17  IERR=6
      CALL ERRPRN(ITEX,IERR)
      RETURN
C
  18  CONTINUE
C
      K1=N1-1
      DO 21 I1=1,K1
      I11=I1
      IF(X1-XM(I1+1)) 22,21,21
  21  CONTINUE
C
  22  K2=N2-1
      DO 23 I2=1,K2
      I22=I2
      IF(Y1-YM(I2+1)) 24,23,23
  23  CONTINUE
C
  24  K3=N3-1
      DO 25 I3=1,K3
      I33=I3
      IF(Z1-ZM(I3+1)) 27,25,25
  25  CONTINUE
C
  27  CONTINUE
      J1=(I33-1)*N1*N2+(I22-1)*N1+I11
      J2=(I33-1)*N1*N2+(I22-1)*N1+I11+1
      J3=(I33-1)*N1*N2+(I22+1-1)*N1+I11+1
      J4=(I33-1)*N1*N2+(I22+1-1)*N1+I11
      J5=(I33+1-1)*N1*N2+(I22-1)*N1+I11
      J6=(I33+1-1)*N1*N2+(I22-1)*N1+I11+1
      J7=(I33+1-1)*N1*N2+(I22+1-1)*N1+I11+1
      J8=(I33+1-1)*N1*N2+(I22+1-1)*N1+I11
C
      IF(GM(J1)+A17) 29,39,39
  29  GM(J1)=-A17
C
  39  IF(GM(J4)+A17) 30,40,40
  30  GM(J4)=-A17
C
  40  IF(GM(J8)+A17) 32,42,42
  32  GM(J8)=-A17
C
  42  IF(GM(J5)+A17) 34,44,44
  34  GM(J5)=-A17
C
  44  IF(GM(J2)+A17) 36,46,46
  36  GM(J2)=-A17
C
  46  IF(GM(J3)+A17) 38,48,48
  38  GM(J3)=-A17
C
  48  IF(GM(J7)+A17) 50,52,52
  50  GM(J7)=-A17
C
  52  IF(GM(J6)+A17) 54,56,56
  54  GM(J6)=-A17
C
  56  CONTINUE
C
      HX1=X1-XM(I11)
      HX2=XM(I11+1)-X1
      HX=XM(I11+1)-XM(I11)
C
      HY1=Y1-YM(I22)
      HY2=YM(I22+1)-Y1
      HY=YM(I22+1)-YM(I22)
C
      HZ1=Z1-ZM(I33)
      HZ2=ZM(I33+1)-Z1
      HZ=ZM(I33+1)-ZM(I33)
C
      FY1=GM(J1)+(GM(J4)-GM(J1))/HY*HY1
      FY2=GM(J2)+(GM(J3)-GM(J2))/HY*HY1
      FXY1=FY1+(FY2-FY1)/HX*HX1
C
      FY1=GM(J5)+(GM(J8)-GM(J5))/HY*HY1
      FY2=GM(J6)+(GM(J7)-GM(J6))/HY*HY1
      FXY2=FY1+(FY2-FY1)/HX*HX1
C
      FUN=FXY1+(FXY2-FXY1)/HZ*HZ1
C
      FUNL3=FUN
      RETURN
      END

      SUBROUTINE ERRPRN(ITEX,IERR)
      CHARACTER*8 ITEX
      COMMON/PRIN/IPRNT
      WRITE(IPRNT,100) ITEX,IERR
 100  FORMAT('MODUL<< ',A8,' >>-code error-IERR= ',I6)
      RETURN
      END


      FUNCTION FUNL2(XM,YM,GM,N1,N2,X1,Y1,IERR)
C
      DIMENSION XM(*),YM(*)
      DIMENSION GM(*)
      CHARACTER*8 ITEX
      COMMON /PRIN/IPRNT
C
      ITEX='FUNL2'
      A17=9.E17
C
      IERR=0
      R1=(XM(N1)-XM(1))*1.E-5
      R2=(YM(N2)-YM(1))*1.E-5
      IF(X1-XM(1)+R1)  2,10,10
   2  IERR=1
      CALL ERRPRN(ITEX,IERR)
      RETURN
  10  IF(XM(N1)+R1-X1) 4,11,11
   4  IERR=2
      CALL ERRPRN(ITEX,IERR)
      RETURN
  11  IF(Y1-YM(1)+R2)  6,12,12
   6  IERR=3
      CALL ERRPRN(ITEX,IERR)
      RETURN
  12  IF(YM(N2)+R2-Y1) 8,17,17
   8  IERR=4
      CALL ERRPRN(ITEX,IERR)
      RETURN
  17  CONTINUE
C
      K1=N1-1
      DO 21 I1=1,K1
      I11=I1
      IF(X1-XM(I1+1)) 22,21,21
  21  CONTINUE
C
  22  K2=N2-1
      DO 23 I2=1,K2
      I22=I2
      IF(Y1-YM(I2+1)) 27,23,23
  23  CONTINUE
C
  27  CONTINUE
      J1=(I22-1)*N1+I11
      J2=(I22-1)*N1+I11+1
      J3=(I22+1-1)*N1+I11+1
      J4=(I22+1-1)*N1+I11
C
      IF(GM(J1)+A17) 28,29,29
  28  IERR=5
      CALL ERRPRN(ITEX,IERR)
C
  29  IF(GM(J4)+A17) 30,34,34
  30  IERR=6
      CALL ERRPRN(ITEX,IERR)
C
  34  IF(GM(J2)+A17) 32,36,36
  32  IERR=7
      CALL ERRPRN(ITEX,IERR)
C
  36  IF(GM(J3)+A17) 40,42,42
  40  IERR=8
      CALL ERRPRN(ITEX,IERR)
C
  42  CONTINUE
C
      HX1=X1-XM(I11)
      HX2=XM(I11+1)-X1
      HX=XM(I11+1)-XM(I11)
C
      HY1=Y1-YM(I22)
      HY2=YM(I22+1)-Y1
      HY=YM(I22+1)-YM(I22)
C
      FY1=GM(J1)+(GM(J4)-GM(J1))/HY*HY1
      FY2=GM(J2)+(GM(J3)-GM(J2))/HY*HY1
      FXY1=FY1+(FY2-FY1)/HX*HX1
C
      FUN=FXY1
      FUNL2=FUN
      RETURN
      END


      SUBROUTINE ProjNum(DR,PROJ,PROJ2,UM,VM,
     %        WM, ALPHAM,BETAM,GAMMAM,BKAPPAM,  
     %         XM,YM,ZM,NPAR,PAR,IERR)
C
      DIMENSION  ALPHAM(*),BETAM(*),GAMMAM(*)
      DIMENSION  BKAPPAM(*)
      DIMENSION  UM(*),VM(*),WM(*),XM(*),YM(*),ZM(*)
      DIMENSION  DR(*),PROJ(*),PROJ2(*)
c      DIMENSION DR(NX*NY*NZ,3)

      DIMENSION  NPAR(*),PAR(*)
      CHARACTER*8 ITEX
      PARAMETER (PI = 3.1415926536)

      NALPHA = NPAR(9)
      NBETA  = NPAR(10)
      NGAMMA = NPAR(18)
      NJ     = NPAR(11)

      NX    = NPAR(12)
      NY    = NPAR(13)
      NZ    = NPAR(14)
      NU    = NPAR(4)
      NV    = NPAR(5)
      NW    = NPAR(6)

      DD  = PAR(36)
      RRE = PAR(37)

      IERR=0
      ITEX='ProjNum'
      IF(NW.LT.100) GO TO 20
      IERR=1
      CALL ERRPRN(ITEX,IERR)
      RETURN
C
  20  CONTINUE
C
      BKAPMAX=ASIN(RRE/DD)
c      BKAPMAX=ATAN2(RRE,DD)

      DO 85 I4=1,NJ
         CALL IALBETGAM(I4,NALPHA,NBETA,I1,I2,I3)
         ALPHA=ALPHAM(I1)
         BETA =BETAM(I2)
         GAMMA=GAMMAM(I3)
         I321=(I3-1)*NBETA*NALPHA+(I2-1)*NALPHA+I1

         CALPHA=COS(ALPHA)
         SALPHA=SIN(ALPHA)
         CBETA=COS(BETA)
         SBETA=SIN(BETA)
         CGAMMA=COS(GAMMA)
         SGAMMA=SIN(GAMMA)

      DO 65 J1=1,NV
         VV=VM(J1)
         QQQ=RRE**2-VV**2
         IF(QQQ .LT. 0) QQQ=0.
         WW0=SQRT(QQQ)
         HW =2.*WW0/(NW-1)

         S=0.
         SS=0.
      DO 55 J3=1,NW
         WW =-WW0+(J3-1)*HW
c         WWDD=1.+WW/DD

c--   T (al,BET,GAM)-------------------------------------------
      XX=
c     % (CALPHA*CBETA*CGAMMA - SALPHA*SGAMMA)*UU + 
     %    (-CALPHA*CBETA*SGAMMA - SALPHA*CGAMMA)*VV+
     %    CALPHA*SBETA*WW

      YY=
c     % (SALPHA*CBETA*CGAMMA + CALPHA*SGAMMA)*UU + 
     %    (-SALPHA*CBETA*SGAMMA + CALPHA*CGAMMA)*VV+
     %    SALPHA*SBETA*WW

      ZZ=
c     % -SBETA*CGAMMA*UU + 
     %     SBETA*SGAMMA*VV +
     %            CBETA*WW


      CALL REGN(XM,YM,ZM,NX,NY,NZ,XX,YY,ZZ,EI)

      IF(EI .LT. 0.) GO TO 55

      S=FUNL3(XM,YM,ZM,DR,NX,NY,NZ,XX,YY,ZZ,IERR)

      WW1=1.
      IF(J3.EQ.1.OR.J3.EQ.NW) WW1=0.5
      SS=SS+S*WW1

  55  CONTINUE

      PROJ2(J1)=SS*HW !*WWDD
     
  65  CONTINUE

      DO 75 J1=1,NV
         I4321=(J1-1)*NGAMMA*NBETA*NALPHA + I321
         PROJ(I4321)=PROJ2(J1)       
  75  CONTINUE
  85  CONTINUE
      RETURN
      END


      SUBROUTINE Summa_W_Harm33(LL,NJ,NNMAX,PPM,ALPHAM,BETAM,
     #                   GAMMAM,NPP,NALPHA,NBETA,NGAMMA,XM,GM)
c  summation of W  harmonics 
      REAL PPM(NPP),ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
      REAL GM(NPP*NALPHA*NBETA*NGAMMA)
      REAL XM((LL+1)*(LL+1)*NNMAX,3)
      REAL WWLMRE(3),WWLMIM(3)

c      COMMON/PRIN/IPRNT
c      PARAMETER (PI = 3.1415926536)

      LMAX=(LL+1)*(LL+1)

c      DO 10 I4=1,NPP
c         PP=PPM(I4)
c      DO 10 I3=1,NGAMMA
c         GAMMA=GAMMAM(I3)
c      DO 10 I2=1,NBETA
c         BETA=BETAM(I2)
c      DO 10 I1=1,NALPHA
c         ALPHA=ALPHAM(I1)


      DO 10 II=1,NJ

         CALL IALBETGAM(II,NALPHA,NBETA,I1,I2,I3)

         ALPHA=ALPHAM(I1)
         BETA =BETAM(I2)
         GAMMA=GAMMAM(I3)
         I321=(I3-1)*NBETA*NALPHA+(I2-1)*NALPHA+I1

      DO 10 I4=1,NPP
         PP=PPM(I4)

         I4321=(I4-1)*NGAMMA*NBETA*NALPHA + I321

c         I4321=(II-1)*NPP + I4
 

        S1=0.
        DO 4  N=1,NNMAX
        DO 4  L=1,LL+1
           L22=L-1
        DO 4  M=1,L22+1
           IIC=(N-1)*LMAX + (L22*L22+M)
           CALL WWCC(N,L22,M-1,ALPHA,BETA,GAMMA,PP,WWLMRE)

c           CALL YLM_RE (L,M-1,TETA,PHI,YLMRE,YLMIM)
           S1=S1+XM(IIC,1)*WWLMRE(1)+XM(IIC,2)*WWLMRE(2)+
     #           XM(IIC,3)*WWLMRE(3)
 4      CONTINUE
        
        DO 6 N=1,NNMAX
        DO 6 L=1,LL+1
           L22=L-1
        DO 6 M=1,L22
           IIS=(N-1)*LMAX + (L22*L22+(M+1)+L22)
           CALL WWSS(N,L22,M,ALPHA,BETA,GAMMA,PP,WWLMIM)

c           CALL YLM_RE (L,M,TETA,PHI,YLMRE,YLMIM)
           S1=S1+XM(IIS,1)*WWLMIM(1)+XM(IIS,2)*WWLMIM(2)+
     #           XM(IIS,3)*WWLMIM(3)
 6      CONTINUE
        GM(I4321)=S1
 10   CONTINUE
      RETURN
      END

      FUNCTION GASDEV(IDUM)
      INTEGER IDUM
      REAL GASDEV
      INTEGER ISET
      REAL FAC, GSET, RSQ, V1,V2,RAN1
      SAVE ISET,GSET
      DATA ISET/0/
      IF (ISET .EQ. 0) THEN
 1       V1=2.*ran1(IDUM)-1.
         V2=2.*ran1(IDUM)-1.
         RSQ=V1**2+V2**2
         IF(RSQ .GE. 1. .OR. RSQ .EQ. 0.) GOTO 1
         FAC=SQRT(-2.*ALOG(RSQ)/RSQ)
         GSET=V1*FAC
         GASDEV=V2*FAC
         ISET=1
       ELSE
          GASDEV=GSET
          ISET=0
       ENDIF
       RETURN
       END

      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END

      SUBROUTINE INDNW3(XM,YM,ZM,NX,NY,NZ,XNEW,YNEW,ZNEW,PAR,
     %                  J1,J2,J3,J4)
C
      DIMENSION XM(*),YM(*),ZM(*),PAR(*)
C
      HX=PAR(18)
      HY=PAR(21)      
      HZ=PAR(24)


      KX=NX-1
      DO 46 I1=1,KX
      I11=I1
      IF(XNEW-XM(I1+1)) 48,46,46
  46  CONTINUE
C
  48  KY=NY-1
      DO 53 I2=1,KY
      I22=I2
      IF(YNEW-YM(I2+1)) 54,53,53
  53  CONTINUE
C
  54  KZ=NZ-1
      DO 55 I3=1,KZ
      I33=I3
      IF(ZNEW-ZM(I3+1)) 57,55,55
  55  CONTINUE
C
  57  CONTINUE
C
      LEFTX=1
      LEFTY=1
      LEFTZ=1
C
      IF(ABS(XNEW-XM(I11)).LE.0.5*HX) LEFTX=0
      IF(ABS(YNEW-YM(I22)).LE.0.5*HY) LEFTY=0
      IF(ABS(ZNEW-ZM(I33)).LE.0.5*HZ) LEFTZ=0
C
      I111=I11+LEFTX
      I222=I22+LEFTY
      I333=I33+LEFTZ
C
      J1=I111
      J2=I222
      J3=I333
      J4=(I333-1)*NX*NY+(I222-1)*NX+I111
C
      RETURN
      END
