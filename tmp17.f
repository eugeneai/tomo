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
     #            LL       =0,
     #            NNMAX    =1,   !5,
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
     #            INFNOT =-1,
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
      

c      DIMENSION DDDM(3*LLLN,3*LLLN)
c      DIMENSION XXXM(LLLN,3), YYYM(3*LLLN)
c      DIMENSION XXXM0(3*LLLN)

      DIMENSION YMVEC(LN,3),YMVEC3(3*LN)
      DIMENSION AMVEC(LN,LN),AMVEC3(3*LN,3*LN)
      DIMENSION XMVEC(3*LN),XMVEC0(3*LN),XMVEC1(3*LN)
      DIMENSION A1MVEC(3*LN)

      DIMENSION PSIM(NX*NZ)

      DIMENSION AAM(9),BM(9),CM(81)


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
!--------------------------------------------------
      IA=3
      JA=1
      IB=1
      JB=3

      ALPHA=PI/6.
      BETA=PI/4.
      GAMMA=PI/6.

      AAM(1)= COS(ALPHA)*SIN(BETA)  
      AAM(2)= SIN(ALPHA)*SIN(BETA) 
      AAM(3)= COS(BETA)        
      BM(1)= COS(ALPHA)*SIN(BETA)
      BM(2)= SIN(ALPHA)*SIN(BETA)
      BM(3)= COS(BETA)

!         AAM(1)=-SIN(BETA)*COS(GAMMA)
!         AAM(2)= SIN(BETA)*SIN(GAMMA)
!         AAM(3)= COS(BETA)
         BM(1)=-SIN(BETA)*COS(GAMMA)
         BM(2)= SIN(BETA)*SIN(GAMMA)
         BM(3)= COS(BETA)


      CALL Tensor_Product_Vec(AAM,IA,JA,BM,IB,JB,CM)

      GOTO 777
!--------------------------------------------------------

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

 777  CONTINUE
 117  FORMAT(12F7.2)
 118  FORMAT(5E20.3)
      STOP
      END
!------------------------------------------------------------
      SUBROUTINE Tensor_Product_Vec(AM,IA,JA,BM,IB,JB,CM)
! tensor product of two matrix 

      INTEGER IA,JA,IB,JB
      REAL AM(IA*JA), BM(IB*JB), CM(IA*IB*JA*JB)
!local
      INTEGER I1,J1,I2,J2,II,JJ      

      DO 8 J1=1,JA 
      DO 8 I1=1,IA
         IJ1=(J1-1)*IA+I1
         DO 4 J2=1,JB
         DO 4 I2=1,IB
            IJ2=(J2-1)*IB+I2

            JJ=(J1-1)*JB+J2
            II=(I1-1)*IB+I2

!            IJ=(JJ-1)*IA*IB+II   ! вывод по столцам
            IJ=(II-1)*JA*JB+JJ    ! вывод по строчкам
!            CM(II,JJ)=AM(IJ1)*BM(IJ2)   
            CM(IJ)=AM(IJ1)*BM(IJ2)   
 4          CONTINUE   
 8       CONTINUE   

      WRITE(*,*) 'Tensor product'
      DO 10 I=1,IA*IB
         WRITE(*,*) (CM((I-1)*JA*JB+J),J=1,JA*JB)
         WRITE(*,*) "-----------------------------------"
 10   CONTINUE
! 100  FORMAT(F8.1)
      RETURN
      END   


