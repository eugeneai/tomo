
!------------------------------------------------------------
      MODULE CHANALS
        IMPLICIT NONE
        INTEGER :: KAN(100)
      END MODULE CHANALS
!------------------------------------------------------------
      MODULE PARAM_INPUT
        IMPLICIT NONE
        REAL(4), allocatable :: PARAM(:)
      END MODULE PARAM_INPUT
!------------------------------------------------------------
      MODULE FREQUENCY
        IMPLICIT NONE
        INTEGER, PARAMETER  :: dp = kind(1.0d0) !working precision=double 
        REAL(dp),DIMENSION(6) :: FREQ  
      END MODULE FREQUENCY

!------------------------------------------------------------
      MODULE MODEL_GAUSS
        IMPLICIT NONE
        REAL(4), allocatable :: GA(:),GAH(:),GB(:),GBH(:),GC(:),GCH(:),GAMPL(:),& 
                                GTET1(:),GPHI1(:) 
      END MODULE MODEL_GAUSS
!-------------------------------------------------------------
      MODULE PARAMTR
        IMPLICIT NONE
        INTEGER, PARAMETER  :: dp = kind(1.0d0) !working precision=double 
!        INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single 
        REAL(dp),  PARAMETER :: AX=20.0  !12.0
        REAL(dp), allocatable :: PAR(:)
        INTEGER, PARAMETER :: QMAX=30   ! the number spherical Bessel functions
        INTEGER, PARAMETER :: SMAXIMUM=20, LMAXIMUM=11 
        REAL(dp), allocatable :: ROOTB(:,:)
      END MODULE PARAMTR  
!------------------------------------------------------------
      MODULE MODEL_PARALL
        IMPLICIT NONE
        REAL(4), allocatable :: PA(:),PAH(:),PB(:),PBH(:),PC(:),PCH(:),PAMPL(:),&
                                PTET1(:),PPHI1(:)  
      END MODULE MODEL_PARALL
!-------------------------------------------------------------
      MODULE MODEL_CYLINDER
        IMPLICIT NONE
        REAL(4), allocatable :: CA(:),CAH(:),CB(:),CBH(:),CC(:),CCH(:),CAMPL(:),&
                                CTET1(:),CPHI1(:)  
      END MODULE MODEL_CYLINDER
!-------------------------------------------------------------
      MODULE COEFFICIENTS
        IMPLICIT NONE
        REAL(8), allocatable :: COEF(:)
      END MODULE COEFFICIENTS
!-------------------------------------------------------------
      MODULE MY_INTERFACES
        IMPLICIT NONE
        INTERFACE
           SUBROUTINE ProjNum_Vec_Decart(DR,PROJ,UM,VM,WM,ALPHA,BETA,GAMMAM, &
                                 XM,YM,ZM,NX,NY,NZ,NU,NV,NW,NGAMMA,IERR)
             INTEGER, PARAMETER  :: dp = kind(1.0d0) !working precision=double 
             INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
             INTEGER, INTENT(IN) :: NX,NY,NZ,NU,NV,NW,NGAMMA,IERR
             REAL(dp),INTENT(IN) :: DR(NX*NY*NZ,3)
             REAL(dp),INTENT(IN) :: XM(NX),YM(NY),ZM(NZ)
             REAL(dp),INTENT(OUT) :: PROJ(NU*NV)
             REAL(dp),INTENT(IN) :: UM(NU),VM(NV),WM(NW)  
             REAL(dp),INTENT(IN) :: GAMMAM(NGAMMA)
             REAL(dp),INTENT(IN) :: ALPHA,BETA
           END SUBROUTINE ProjNum_Vec_Decart
        END INTERFACE

!        INTERFACE
!           FUNCTION PSI_FUN(L1,M1,QMAX,RR,TETA,PHI)
!             INTEGER, PARAMETER  :: dp = kind(1.0d0) !working precision=double 
!             INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
!             INTEGER, INTENT(IN) :: L1,M1,QMAX
!             REAL(dp),INTENT(IN) :: RR,TETA,PHI
!!             COMPLEX(dp),INTENT(OUT) :: PSI_FUN
!           END FUNCTION PSI_FUN
!        END INTERFACE
      END MODULE MY_INTERFACES
!-------------------------------------------------------------

      MODULE MAIN_ARRAYS
        IMPLICIT NONE
        INTEGER, PARAMETER  :: dp = kind(1.0d0) !working precision=double 
!        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 

        REAL(dp), allocatable::PPM(:),TETM(:),PHIM(:),RRM(:)
        REAL(dp), allocatable::XM(:),YM(:)
        REAL(dp), allocatable::ZM(:)
        REAL(dp), allocatable::UM(:),VM(:),WM(:)
        REAL(dp), allocatable::XNU(:),YNU(:),ZNU(:),RNU(:)
        REAL(dp), allocatable::ALPHAM(:),BETAM(:),GAMMAM(:)
!        REAL(dp), allocatable::DR(:),DI(:)
!        REAL(dp), allocatable::FDR(:),FDI(:)
!        REAL(dp), allocatable::PRJ3DRE(:),PRJ3DIM(:)
!        REAL(dp), allocatable::PROJRE2D(:),PROJIM2D(:)  !!!
        REAL(dp), allocatable:: REM(:)
        REAL(dp), allocatable::FURDR(:)
!        REAL(dp), allocatable::DRRE(:),DRIM(:)
!        REAL(dp), allocatable::FURPRJ3DRE(:),FURPRJ3DIM(:)
      END MODULE MAIN_ARRAYS
!-------------------------------------------------------------
!      MODULE LAMBDA
!        IMPLICIT NONE
!        REAL(4) :: SLAMBDA
!      END MODULE LAMBDA
!-------------------------------------------------------------


      PROGRAM COEFSPHER

        USE MAIN_ARRAYS

        USE FREQUENCY

        USE PARAM_INPUT

        USE CHANALS
        
!        USE SHTOOLS        

        USE PARAMTR

        USE MY_INTERFACES

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

        INTEGER :: LN
        INTEGER :: MX,MY,MZ
        INTEGER :: NX,NY,NZ,NPP
        INTEGER :: NU,NV,NW    
        INTEGER :: NTETA,NPHI,NRR
        INTEGER :: NNMAX,NTEN
!        INTEGER :: INTRP,LZER1
        INTEGER :: IMOD,KKK
        INTEGER :: NJ,NITER,KNZ,KNY,KNX
        INTEGER :: NII,NJJ

        INTEGER :: I,II,ISIGNUM   
        INTEGER :: I1,I2,I3,I4,I21,I321,I4321
        INTEGER :: J,J1,J2,J21,J3
        INTEGER :: IERR,DUMMY,ERR_ALLOC
        INTEGER :: NDAT,NDAT1,NDAT2,NDAT3
        INTEGER :: KSEC,IND,JP
!        INTEGER*4 MM   ! число членов в ряде для рассеянного поля 

        REAL(wp) :: EOT,RRE
        REAL(wp) :: HX,HY,HZ,HPHI,HTET,HPP
        REAL(wp) :: HNUX,HNUY,HNUZ 
        REAL(wp) :: PHI,TET
        REAL(wp) :: SKZER,SLZER,SDELTA
        REAL(wp) :: UNU,VNU,WNU
        REAL(wp) :: SVZKA1,SVZKA2

        REAL(wp) :: ALPHA,BETA,GAMMA,LAMBDA
        REAL(wp) :: HALPHA,HBETA,HGAMMA
        INTEGER :: NALPHA,NBETA,NGAMMA
        INTEGER :: JJ,KJ,NN1,NN2

 !       REAL(4) ZR1,ZI1,Inv_Ref_Ind
        REAL(wp) :: SM1,SM2,SM3
        REAL(wp) :: S_Bessel_J,INT_BESSEL_FUNC
!        REAL(4) :: JJI
! ВРЕМЕННО: ------------------------------
        REAL(wp) ARE,AIM,CEILING
        INTEGER :: NN,MM,II1,NIT,LL,MMM ,ERR
        INTEGER :: LMAX,SMAX,NB,NB1,NZZ
        INTEGER :: L1,M1,L2,M2,LM2,JJJ,JPP,LMS,K1,K2,L3
        COMPLEX(wp) :: CKLM
        REAL(wp) :: AA(2),BB(3),QQ,AAM(6),BBM(6), AU(2),AV(2)
        REAL(wp) :: CCM(24)
        REAL(wp) :: X,XX,YY,ZZ,XXX
        REAL(wp) :: ALL
        REAL(wp) :: XXNU,YYNU,ZZNU
        REAL(wp) :: DRR,DAX,GN
        REAL(wp) :: SBES,BBBB,CCCC
        COMPLEX(wp) :: uuee,F_L2M2_LM,PSI_2D,CYY,CYY1
        COMPLEX(wp) :: CYYM(3),CYYM1(3),SS
        REAL(wp) :: R0,RR,KK,GNN,DER
        REAL(wp) :: rUnit(3), tetaUnit(3), fiUnit(3)
        COMPLEX(wp),PARAMETER :: ci = (0.,1.)
        COMPLEX(wp) :: CCEXP
!        REAL(wp) :: ALPHA,BETA
        REAL(wp) :: RL1,RM1,RM2
        INTEGER :: NPROJ    ! the projection number
!functions 
        REAL(wp) :: djmn
        COMPLEX(wp) :: FCoefSWF,CYLM,ylm,D_WIGNER
        COMPLEX(wp) :: PSI_FUN
        REAL(wp) :: PLM1M2,DerSpherBessel
        REAL(wp) :: BessJ,PLGNDR_NORM,PLGNDR
!===============================================================
!	integer, parameter ::	maxdeg = 100, maxgrid = 400
!	character*80 ::		infile
!	real*8 ::		cilm(2, maxdeg+1, maxdeg+1), header(8), grid(maxgrid, maxgrid), &
!				griddh(maxgrid, maxgrid), cilm2(2, maxdeg+1, maxdeg+1), interval, &
!				error, maxerror
!	integer ::		n, nlat, nlong, lmax2,l,m

!==============================================================	
        REAL(wp), allocatable:: XXM(:,:),XXM2(:,:),XXM0(:),XXM1(:)
        REAL(wp), allocatable:: XMVEC(:),XMVEC0(:),XMVEC1(:)
        REAL(wp), allocatable:: YMVEC(:,:),YMVEC3(:) 
        REAL(wp), allocatable:: DR(:),DI(:)
        REAL(wp), allocatable:: DGMRE(:,:),DGMIM(:,:)
        REAL(wp), allocatable:: GMRE(:,:),GMIM(:,:)
        REAL(wp), allocatable:: AMVEC(:,:), AMVEC3(:,:)
!        REAL(wp), allocatable:: ALPHAM(:),BETAM(:),GAMMAM(:)
        REAL(wp), allocatable:: DR1(:),DI1(:)
        REAL(wp), allocatable:: ACOEF(:)
        REAL(wp), allocatable:: RPART2(:)
!        COMPLEX(wp), allocatable:: RPART(:)
        REAL(wp), allocatable:: ACOF2(:)
        REAL(wp), allocatable:: ACOFRE(:),ACOFIM(:)
        REAL(wp), allocatable::PRJ3DRE(:),PRJ3DIM(:)
        REAL(wp), allocatable:: PROJRE1D(:),PROJIM1D(:)
        REAL(wp), allocatable::DREX(:),DRIX(:)
        REAL(wp), allocatable::FDR(:),FDI(:),FDR2(:)
!        REAL(wp), allocatable::SVZKARE(:),SVZKAIM(:)
        REAL(wp), allocatable::PRJ3DRE2(:),PRJ3DIM2(:)
        REAL(wp), allocatable::DRRE(:),DRIM(:)
        REAL(wp), allocatable::FM(:)
        REAL(wp), allocatable::PRFURRE(:),PRFURIM(:)
        REAL(wp), allocatable::FURPRJ3DRE(:),FURPRJ3DIM(:)
        REAL(wp), allocatable:: PROJRE2D(:) ,PROJIM2D(:)
        REAL(wp), allocatable::SMATRE2(:)  !,AM(:)
        REAL(wp), allocatable::SMATRE(:),SMATIM(:)
        REAL(wp), allocatable::SMA1(:)
        REAL(wp), allocatable::SMA2(:)
        REAL(wp), allocatable::SMAT(:)
        REAL(wp), allocatable::BES(:)
        REAL(wp), allocatable:: BMODRE(:),BMODIM(:),BMODRE1(:),BMODIM1(:)
        REAL(wp), allocatable:: SMODRE(:),SMODIM(:)
        REAL(wp), allocatable:: BMOD(:),BMOD1(:)
        REAL(wp), allocatable:: PMATRE(:),PMATIM(:)
        INTEGER, allocatable:: AM(:,:)


        COMPLEX(wp), allocatable::PMATC(:), SMATC(:)
        COMPLEX(wp), allocatable::PRFURC(:)
        COMPLEX(wp), allocatable::MLMM(:),MLMM1(:)
        COMPLEX(wp), allocatable:: BMODC(:,:)  ,BMODC1(:)
        COMPLEX(wp), allocatable:: PROJC(:)


        OPEN(11,FILE ='Data_Input.dat', ERR=300)


        OPEN(12,FILE ='Data_Gauss.dat', ERR=400)
        OPEN(14,FILE ='Data_Parall.dat',ERR=410)
        OPEN(15,FILE ='Data_Cylinder.dat',ERR=420)

        OPEN(19,FILE='modelYZ.idl')
        OPEN(22,FILE='modelXY.idl')
        OPEN(25,FILE='modelXZ.idl')

        OPEN(21,FILE='modelYZ_1D.gnu')
        OPEN(24,FILE='modelXY_1D.gnu')
        OPEN(27,FILE='modelXZ_1D.gnu')

        OPEN(20,FILE='modelYZ_EXC.gnu')
        OPEN(23,FILE='modelXY_EXC.gnu')
        OPEN(26,FILE='modelXZ_EXC.gnu')

        OPEN(29,FILE='fourierYZ.idl')
        OPEN(30,FILE='cof_exc.gnu')
!        OPEN(30,FILE='fourierYZ.gnu')
        OPEN(31,FILE='cof_rec.gnu')

        OPEN(32,FILE='fourierXY.idl')
        OPEN(33,FILE='fourierXY.gnu')
        OPEN(34,FILE='fourierXY_1D.gnu')

        OPEN(35,FILE='fourierXZ.idl')
        OPEN(36,FILE='fourierXZ.gnu')
        OPEN(37,FILE='fourierXZ_1D.gnu')

        OPEN(39,FILE='proj_furRE1.gnu')
        OPEN(40,FILE='difrag_projRE.gnu')
!        OPEN(41,FILE='proj_furRE3.gnu')
        OPEN(42,FILE='proj_furIM1.gnu')
        OPEN(43,FILE='difrag_projIM.gnu')

        OPEN(44,FILE='proj_EXAC_RE.gnu')
        OPEN(45,FILE='proj_EXAC_IM.gnu')

        OPEN(46,FILE='vector_exc_xy.gnu')
        OPEN(47,FILE='vector_exc_xz.gnu')
        OPEN(51,FILE='vector_exc_yz.gnu')

        OPEN(54,FILE='phaze_exc_zy.gnu')
        OPEN(56,FILE='phaze_exc_zx.gnu')
        OPEN(58,FILE='phaze_exc_yx.gnu')

        OPEN(74,FILE='modelYZ_REC.gnu')
        OPEN(75,FILE='modelXY_REC.gnu')
        OPEN(76,FILE='modelXZ_REC.gnu')

        OPEN(77,FILE='modelYZ_REC_1D.gnu')
        OPEN(78,FILE='modelXY_REC_1D.gnu')
        OPEN(79,FILE='modelXZ_REC_1D.gnu')

        OPEN(80,FILE='Model.gnu')
        OPEN(81,FILE='Projection.gnu')
        OPEN(82,FILE='Projection2.gnu')

!        OPEN(93,FILE='model_cyl_exc3D_6.tor',form='binary',STATUS='unknown')
!        OPEN(94,FILE='model_cyl_rec3D_6.tor',form='binary',STATUS='unknown')

!        OPEN(83,FILE='PRJ3DRE_3D_6.tor',form='binary',STATUS='unknown')
!        OPEN(84,FILE='PRJ3DIM_3D_6.tor',form='binary',STATUS='unknown')
!        OPEN(85,FILE='INTENSITY_6.tor',form='binary',STATUS='unknown')
        OPEN(86,FILE='INTENSITY_2.gnu')
        OPEN(87,FILE='87.gnu')

        OPEN(48,FILE='FurProjec_RE.dat')
        OPEN(49,FILE='FurProjec_IM.dat')

        OPEN(50,FILE='Projec_RE.dat')
!        OPEN(51,FILE='Projec_IM.dat')
!        OPEN(58,FILE='phaze_exc_z.gnu')
        OPEN(59,FILE='phaze_rec_z.gnu')

!----------------------------------------------------------

        KAN(9)  =9
        KAN(11) =11
        KAN(12) =12
        KAN(14) =14
        KAN(15) =15
        KAN(19) =19
        KAN(20) =20
        KAN(21) =21
        KAN(22) =22
        KAN(23) =23
        KAN(24) =24
        KAN(25) =25
        KAN(26) =26
        KAN(27) =27
        KAN(29) =29
        KAN(30) =30
        KAN(31) =31
        KAN(32) =32
        KAN(33) =33
        KAN(34) =34
        KAN(35) =35
        KAN(36) =36
        KAN(37) =37
        KAN(39) =39
        KAN(40) =40
        KAN(42) =42 
        KAN(43) =43
        KAN(44) =44
        KAN(45) =45

        KAN(46) =46
        KAN(47) =47
        KAN(51) =51

        KAN(54) =54
        KAN(56) =56      
        KAN(58) =58

        KAN(74) =74
        KAN(75) =75
        KAN(76) =76
        KAN(77) =77
        KAN(78) =78
        KAN(79) =79
        KAN(80) =80
        KAN(81) =81
        KAN(82) =82
        KAN(83) =83
        KAN(84) =84
        KAN(85) =85
        KAN(86) =86
        KAN(87) =87

        KAN(48) =48
        KAN(49) =49
        KAN(50) =50
!        KAN(51) =51

!-----------------------------------
!        EOTT=SIGNUM(0.)
!        write(*,*) 'EOTT=',EOTT
!-----------------------------------



        CALL SIZE_PARAM(NDAT)  !_sp

        ALLOCATE (PARAM(NDAT),STAT=ERR_ALLOC)  
        IF(ERR_ALLOC .NE. 0) THEN
           IERR=1
           GOTO 999
        ENDIF  

!        write(*,*) 'NDAT=',NDAT

        CALL INPUT_PARAM(PARAM,NDAT)           !_sp
	
        write(*,*) 'PARAM_INPUT'
        write(*,*) (PARAM(I),I=1,NDAT)

        EOT   = PARAM(1)
        NX    = PARAM(2)
        NY    = PARAM(3)
        NZ    = PARAM(4)
        RRE   = PARAM(5)
        NPHI  = PARAM(6)
        NTETA = PARAM(7)
        NRR   = PARAM(8)
        NTEN  = PARAM(9)
        NNMAX = PARAM(10)
        LL    = PARAM(11)
        NPP   = PARAM(12)
        NITER = PARAM(13)
        NALPHA = PARAM(14)
        NBETA  = PARAM(15)
        NGAMMA = PARAM(16)
        NU     = PARAM(17)
        NV     = PARAM(18)
        NW     = PARAM(19)

        LN=(LL+1)*(LL+1)*NNMAX

!------------------------------------
!        MVV = MUU
!        MWW = MUU
!        NUU = 2**MUU+1
!        NVV = 2**MVV+1
!        NWW = 2**MWW+1
!        NU  = NUU
!        NV  = NVV
!        NW  = NWW

!        MX  = MUU
!        MY  = MVV
!        MZ  = MWW
        
!        NX  = 2**MX+1
!        NY  = 2**MY+1
!        NZ  = 2**MZ+1 
!        NPP = NX

!        NALPHA=1
!        NBETA =1
!        NGAMMA=1
        NJ=NALPHA*NBETA*NGAMMA

!        NJ  = NTET*NPHI
        JJ  = 16   !  JJ=(L+1)**2,  L -is power of spherical harmonics  
!        QMAX=1
!        LMAX=QMAX-1   ! LMAX < or = QMAX
        NB  = QMAX


!        LMAXIMUM=11 
!        SMAXIMUM=20 
        LMAX=LMAXIMUM - 5  
        SMAX=5 ! not larger 20

        LM2=(LMAX+1)**2
        LMS=(LMAX+1)**2*SMAX
        L3=2  ! для проверки матрицы  
!        YM((LL+1)*(LL+1)*NNMAX,NTEN*NTEN)
        ALLOCATE (PAR(50), PPM(NPP), RRM(NRR),  &       
                  TETM(NTETA),PHIM(NPHI),       &
                  XM(NX),YM(NY),ZM(NZ),        &
                  UM(NU),VM(NV),WM(NW),        &
                  XXM(LN,NTEN*NTEN),XXM2(LN,NTEN*NTEN),XXM0(LN),XXM1(LN), &     
                  YMVEC(LN,NTEN*NTEN), YMVEC3(LN*NTEN*NTEN),   &
!                  XMVEC(NTEN*NTEN*LN),XMVEC0(NTEN*LN),XMVEC1(NTEN*LN), &
                  AMVEC(LN,LN), AMVEC3(NTEN*NTEN*LN,NTEN*NTEN*LN),     &
                  XMVEC(NTEN*NTEN*LN),XMVEC0(NTEN*NTEN*LN),   &
                  XMVEC1(NTEN*NTEN*LN),&

                  GMRE(NPP*NTETA*NPHI,NTEN*NTEN),            &
                  DGMRE(NX*NY*NZ,NTEN*NTEN),      &

                  DR(NX*NY*NZ),DI(NX*NY*NZ),   &
!                  ROOTB(SMAXIMUM,LMAXIMUM),    &
!               DGMIM(NX*NY*NZ,NTEN),           &
!               GMIM(NX*NY*NZ,NTEN),             &
!                  DR1(NX*NY*NZ),DI1(NX*NY*NZ),   &
!                   ACOEF((LMAX+1)*SMAX),          &
!                  RPART2(2*(LMAX+1)**2*SMAX),          &
!                  RPARTRE((LMAX+1)**2*SMAX),          &
!                  RPARTIM((LMAX+1)**2*SMAX),          &
!                  RPART((LMAX+1)**2*SMAX),          &
!                  ACOFRE((LMAX+1)**2*SMAX),          &
!                  ACOFIM((LMAX+1)**2*SMAX),          &
!                  ACOF2(2*(LMAX+1)**2*SMAX),          &
                  ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA), &
!                  DR2(NX*NY*NZ),DI2(NX*NY*NZ),   &
!                  DREX((LMAX+1)**2*SMAX**2),   &
!                  DRIX((LMAX+1)*SMAX),   &
                  PRJ3DRE(NU*NGAMMA*NALPHA*NBETA),   &
                  PRJ3DIM(NU*NGAMMA*NALPHA*NBETA),   &
                  PROJRE1D(NU),PROJIM1D(NU),           &    !
!                  PROJC(NUU*NVV),           &    
!                  PROJIM2D(NUU*NVV),           &    !
!                  PROJIM2D(NUU*NGAMMA),           &    !
!                  FDR(NX*NY*NZ),               &
!                  FDI(NX*NY*NZ),               &
!                  FDR2(2*NX*NY*NZ),            &
!                  DRRE(NUU*NVV*NWW),           &
!                  DRIM(NUU*NVV*NWW),           &
!                  FM(NX*NB),                   &
!                  BES(NB),                     &
!                  SVZKARE(NITER),              & 
!                  SVZKAIM(NITER),              & 
!                  PRFURRE(NUU*NVV*(LMAX+1)**2),   &
!                  PRFURIM(NUU*NVV*(LMAX+1)**2),   &
!                  PRFURC(NUU*NVV*(LMAX+1)**2),   &

!                  FURPRJ3DRE(NUU*NVV*NJ),   &
!                  FURPRJ3DIM(NUU*NVV*NJ),   &
!                  SMATRE2(4*(LMAX+1)**2*SMAX*(LMAX+1)**2*SMAX),       &
!                  SMATRE((LMAX+1)**2*(LMAX+1)**2),       &
!                  SMATIM((LMAX+1)**2*(LMAX+1)**2),       &
!                  SMATC((LMAX+1)**2*SMAX*(LMAX+1)**2*SMAX),       &  !!!!!!!!!!

!                  SMAT((LMAX+1)**2),       &
!                  SMA1(LMAX+1),       &
!                  SMA2(2*(LMAX+1)**2),       &
!                  PMATC((LMAX+1)**2*NPP),       &
!                  PMATRE((LMAX+1)**2*NPP),       &
!                  PMATIM((LMAX+1)**2*NPP),       &
!                  AM(4*L3*(L3+2),4*L3*(L3+2)),                    &
!                  BMODRE(NPP*NTETA*NPHI),            &
!                  BMODIM(NPP*NTETA*NPHI),            &     
!                  SMODRE(NTETA*NPHI),                &
!                  SMODIM(NTETA*NPHI),                 &
!                  BMODC1(NPP*NPHI),            &
!                  BMODRE1(NPP*NPHI),                 &
!                  BMODIM1(NPP*NPHI),                 &
!                  BMODC(NX*NY*NZ,3),                 &
!                  BMOD(NPP),BMOD1(NPP), &
!                  MLMM(3),MLMM1(3),  &
                  STAT=ERR_ALLOC)
        IF(ERR_ALLOC .NE. 0) THEN
           IERR=2
           GOTO 888
        ENDIF   
!-------------------------------------------------------
        XX=-1.0
        YY=-1.0
        ZZ= 1.0
!        CALL POLAR_2D(XX,YY,ARE,AIM)
!        CALL POLAR(XX,YY,ZZ,TET,PHI)
!        write(*,*) 'qqqqqqq',PHI*180./PI,TET*180/PI 

!------------------------------------------------------

! complex exp definition
!      CCEXP(X)=CMPLX(COS(X),SIN(X))

      CALL DATA(RRE,PPM,RRM,XM,YM,ZM,NPP,NRR,NX,NY,NZ,PHIM,TETM,  & 
              NPHI,NTETA,ALPHAM,BETAM,GAMMAM,NALPHA,NBETA,NGAMMA, &
              HX,HY,HZ,HPHI,HTET,UM,VM,WM,NU,NV,NW)



!=================================================================

      CALL XXMM(LL,NNMAX,XXM,NTEN)  
!      write(*,*) 'solution exact XXM:'
!      write(*,117) (XXM(I,1), I=1, LN)
!      write(*,117) (XXM(I,2), I=1, LN)
!      write(*,117) (XXM(I,3), I=1, LN)
!      write(*,117) (XXM(I,4), I=1, LN)
!      write(*,117) (XXM(I,5), I=1, LN)
!      write(*,117) (XXM(I,6), I=1, LN)
!      write(*,117) (XXM(I,7), I=1, LN)
!      write(*,117) (XXM(I,8), I=1, LN)
!      write(*,117) (XXM(I,9), I=1, LN)

      CALL Vec_Line(XXM,XMVEC,LN,NTEN)

      CALL GNUFORM_0 (NTEN**2*LN,XMVEC,KAN(30))   !! plot cof_exc.gnu
 
     write(*,*) 'solution exact XMVEC',LL,LN, NTEN
      write(*,117) (XMVEC(I), I=1, NTEN**2*LN)


!не забыть поменял RRM на UM
      CALL Summa_Harm_Vec(LL,NNMAX,RRM,TETM,PHIM, &
                         NRR,NTETA,NPHI,NTEN,XXM,GMRE)
!!не забыть поменял RRM на UM
      CALL POLAR_DECART(GMRE,XM,YM,ZM,NX,NY,NZ,RRM,TETM,PHIM, &
                              NRR,NTETA,NPHI,NTEN,DGMRE)
      KNZ=NZ/2
      KNY=NY/2
      KNX=NX/2
      KKK=1
!  IND=-1  -- section XZ
!  IND= 0  -- section XY
!  IND= 1  -- section YZ
      IND=0
      CALL SPRINT_3D(DGMRE,XM,YM,ZM,NX,NY,NZ,IND,KNY,NTEN,1,KAN(26))
      IND=0
      CALL SPRINT_3D(DGMRE,XM,YM,ZM,NX,NY,NZ,IND,KNZ,NTEN,2,KAN(23))
      IND=0
      CALL SPRINT_3D(DGMRE,XM,YM,ZM,NX,NY,NZ,IND,KNX,NTEN,3,KAN(20))

      CALL GNUFORMVEC_Z(XM,YM,NX,NY,NZ,KNZ,DGMRE,NTEN,KAN(46),KAN(58))

      CALL GNUFORMVEC_Y(XM,ZM,NX,NY,NZ,KNY,DGMRE,NTEN,KAN(47),KAN(56))

      CALL GNUFORMVEC_X(YM,ZM,NX,NY,NZ,KNX,DGMRE,NTEN,KAN(51),KAN(54))
!----------------


      CALL ProjNum_TEN(DGMRE,PRJ3DIM,PPM,  &
              ALPHAM,BETAM,GAMMAM,NALPHA,NBETA,NGAMMA, &
              XM,YM,ZM,NX,NY,NZ,NPP,NW,NTEN)
       
      CALL Proj_Summa_W_Harm(LL,NJ,NNMAX,PPM,ALPHAM,BETAM, &
              GAMMAM,NPP,NALPHA,NBETA,NGAMMA,XXM,NTEN,PRJ3DRE)

       KJ=35

       DO 11 I4=1,NU
       DO 11 II=KJ,KJ  ! 1,NJ
!       DO 11 I4=1,NU
          CALL IALBETGAM(II,NALPHA,NBETA,I1,I2,I3)
          ALPHA=ALPHAM(I1)
          BETA =BETAM(I2)
          GAMMA=GAMMAM(I3)

          I321=(I3-1)*NBETA*NALPHA+(I2-1)*NALPHA+I1
          I4321=(I4-1)*NGAMMA*NBETA*NALPHA + I321
          PROJIM1D(I4)=PRJ3DIM(I4321)
          PROJRE1D(I4)=PRJ3DRE(I4321)
11     CONTINUE   

!         write(*,*) 'ALPHA',PRJ3DIM

!         write(*,*) 'BETA',PRJ3DIM

!          write(*,*) 'GAMMA',GAMMAM

!          write(*,*) 'I1,I2,I3',I1,I2,I3
          write(*,*) 'alpha_n=',ALPHA*180/PI, 'I1=',I1
          write(*,*) 'beta_n=',BETA*180/PI, 'I2=', I2
          write(*,*) 'gamma_n=',GAMMA*180/PI, 'I3=',I3
          CALL GNUFORM1(UM,NU,PROJIM1D,KAN(45)) 
          CALL GNUFORM1(UM,NU,PROJRE1D,KAN(44)) 


!          GOTO 777
!-------------------------

!          CALL RightPart_WW_TEN(PRJ3DIM,LL,NJ,NNMAX,  &
!                PPM,PPM,ALPHAM,BETAM,GAMMAM,        &
!                NPP,NALPHA,NBETA,NGAMMA,             &
!                NTEN,PAR,YMVEC)

!          CALL Vec_Line(YMVEC,YMVEC3,LN,NTEN)


          CALL RightPart_Full(LL,NNMAX,NTEN,     &
                    PRJ3DIM,PPM,ALPHAM,BETAM,GAMMAM,  &
                    NPP,NALPHA,NBETA,NGAMMA,PAR,YMVEC3)

!       write(*,*) 'right part:YMVEC(1)'
!       write(*,117) (YMVEC(I,1), I=1, LN)
!       write(*,*) 'right part:YMVEC(2)'
!       write(*,117) (YMVEC(I,2), I=1, LN)
!       write(*,*) 'right part:YMVEC(3)'
!       write(*,117) (YMVEC(I,3), I=1, LN)

       write(*,*) 'right part:YMVEC3'
       write(*,117) (YMVEC3(I), I=1, NTEN*NTEN*LN)

!       GOTO 777


       CALL Matrix_WW_Full(LL,NNMAX,PPM,ALPHAM,BETAM,GAMMAM, &
                     NPP,NALPHA,NBETA,NGAMMA,NTEN,PAR,AMVEC3)

! initial values for the solution
       DO 13 I=1,NTEN**2*LN
          XMVEC0(I)=0.
13     CONTINUE 

!          NN1=NTEN**2*LN
!          NN2=NTEN**2*LN
          CALL XYANG(AMVEC3,YMVEC3,XMVEC0,NTEN**2*LN,NTEN**2*LN,NITER)
!         write(*,*) '111111111111111111111'
!         write(*,*) YMVEC3

!       GOTO 777

      CALL Vec_Line_Inv(XMVEC0,XXM2,LN,NTEN)

      write(*,*) 'solution reconstr. XMVEC0'
      write(*,117) (XMVEC0(I), I=1, NTEN**2*LN)


777   CONTINUE
      STOP
117  FORMAT(9F7.2)
118  FORMAT(5E20.3)
100   FORMAT(E10.3)
130   FORMAT(4E15.3)
110   FORMAT(1X,8F10.3)   
300   PRINT *, 'Can not open Data_Input.dat'
400   PRINT *, 'Can not open Data_Models.dat'
410   PRINT *, 'Can not open Data_Parall.dat'
420   PRINT *, 'Can not open Data_Cylinder.dat'
999   PRINT *, 'Can not ALLOCATE in a Main program IERR=',IERR
888   PRINT *, 'Can not ALLOCATE in a Main program IERR=',IERR
      STOP
      END PROGRAM COEFSPHER
!====================================================


      SUBROUTINE Matrix_WW_Full(LL,NNMAX,PPM,ALPHAM,BETAM,GAMMAM, &
                     NPP,NALPHA,NBETA,NGAMMA,NTEN,PAR,AM)


        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

        INTEGER :: LL,NNMAX,NPP,NALPHA,NBETA,NGAMMA,NTEN,NI,NJ
        REAL(wp) :: PAR(50)   !!!!!!!
        REAL(wp) :: PPM(NPP),ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
!        REAL(wp) :: RRM(NPP)   !!! !NRR=NPP
        REAL(wp) :: AM((LL+1)**2*NNMAX*NTEN**2,(LL+1)**2*NNMAX*NTEN**2)

!local
        INTEGER :: LMAX,N1,N2,L1,L2,M1,M2,L11,L22
        INTEGER :: I1,J1,IJ1,I2,J2,IJ2
        INTEGER :: MU_W,NU_W,MU_V,NU_V
        INTEGER :: IIC1,IIC2,IIS1,IIS2
        REAL(wp) :: VLM
        REAL(wp) :: HALPHA,HBETA,HGAMMA,HPP
!functions
        REAL(wp) :: Matrix_Term_WW, Matrix_Term_VV, Matrix_Term_WV

!      CALL IALBET(NI,NTEN,NI1,NI2)
!      CALL IALBET(NJ,NTEN,NJ1,NJ2)

      HALPHA=PAR(38)
      HBETA =PAR(39)
      HGAMMA=PAR(40)
      HPP  =PAR(35)
    
      LMAX=(LL+1)*(LL+1)

!---------------------------------------
!   cos * cos  (1,1)

      DO 32 J1=1,NTEN
      DO 32 I1=1,NTEN
         IJ1=(J1-1)*NTEN+I1
      DO 32 J2=1,NTEN
      DO 32 I2=1,NTEN
         IJ2=(J2-1)*NTEN+I2


      DO 30  N1=1,NNMAX
      DO 30  L1=1,LL+1
         L11=L1-1
      DO 30  M1=1,L11+1
         IIC1=(N1-1)*LMAX + (L11*L11+M1)
         MU_W=LMAX*(IJ1-1)*NNMAX+IIC1  

!      DO 28 J2=1,NTEN
!      DO 28 I2=1,NTEN
!         IJ2=(J2-1)*NTEN+I2
      DO 30  N2=1,NNMAX
      DO 30  L2=1,LL+1
         L22=L2-1
      DO 30  M2=1,L22+1
         IIC2=(N2-1)*LMAX + (L22*L22+M2)

         NU_W=LMAX*(IJ2-1)*NNMAX+IIC2  

         VLM=Matrix_Term_WW(I1,J1,I2,J2,LL,N1,L11,M1-1,N2,L22,M2-1,NNMAX,  &
             PPM,ALPHAM,BETAM,GAMMAM,NPP,NALPHA,NBETA,NGAMMA,NTEN,PAR)

!      CALL WW_WW_TEN(L11,M1-1,N1,L22,M-1,N,       &
!              PPM,RRM,NPP,HPP,                    &
!              NBETA,NGAMMA,NALPHA,                  &
!              ALPHAM,BETAM,GAMMAM,PAR,VLM,NI,NJ,NTEN)

         AM(MU_W,NU_W)=VLM

!      write(*,*) '1111111111111111111111111'
!      write(29,*) MU_W,NU_W
 30   CONTINUE
 32   CONTINUE
 
!---------------------------------------------------
!  cos * sin  (1,2)

      DO 50 J1=1,NTEN
      DO 50 I1=1,NTEN
         IJ1=(J1-1)*NTEN+I1

      DO 50 J2=1,NTEN
      DO 50 I2=1,NTEN
         IJ2=(J2-1)*NTEN+I2

      DO 48  N1=1,NNMAX
      DO 48  L1=1,LL+1
         L11=L1-1
      DO 48  M1=1,L11+1
         IIC1=(N1-1)*LMAX + (L11*L11+M1)

         MU_W=LMAX*(IJ1-1)*NNMAX+IIC1  

!      DO 48 J2=1,NTEN
!      DO 48 I2=1,NTEN
!         IJ2=(J2-1)*NTEN+I2
      DO 48  N2=1,NNMAX
      DO 48  L2=1,LL+1
         L22=L2-1
      DO 48  M2=1,L22
         IIS2=(N2-1)*LMAX + (L22*L22+(M2+1)+L22)

         NU_V=LMAX*(IJ2-1)*NNMAX+IIS2  

         VLM=Matrix_Term_WV(I1,J1,I2,J2,LL,N1,L11,M1-1,N2,L22,M2,NNMAX,  &
             PPM,ALPHAM,BETAM,GAMMAM,NPP,NALPHA,NBETA,NGAMMA,NTEN,PAR)

!      CALL WW_WW_TEN(L11,M1-1,N1,L22,M,N,       &
!              PPM,RRM,NPP,HPP,                    &
!              NBETA,NGAMMA,NALPHA,                  &
!              ALPHAM,BETAM,GAMMAM,PAR,VLM,NI,NJ,NTEN)

         AM(MU_W,NU_V)=VLM

!      write(*,*) '1111111111111111111111111'
!      write(29,*) MU_W,NU_V
 48   CONTINUE
 50   CONTINUE

!------------------------------------------

!   sin * sin  (2,2)

      DO 40 J1=1,NTEN
      DO 40 I1=1,NTEN
         IJ1=(J1-1)*NTEN+I1
      DO 40 J2=1,NTEN
      DO 40 I2=1,NTEN
         IJ2=(J2-1)*NTEN+I2

      DO 38  N1=1,NNMAX
      DO 38  L1=1,LL+1
         L11=L1-1
      DO 38  M1=1,L11
         IIS1=(N1-1)*LMAX + (L11*L11+(M1+1)+L11)
!         IIS1=(N1-1)*LMAX + L11*L11+M1+L11

         MU_V=LMAX*(IJ1-1)*NNMAX+IIS1  

!      DO 38 J2=1,NTEN
!      DO 38 I2=1,NTEN
!         IJ2=(J2-1)*NTEN+I2
      DO 38  N2=1,NNMAX
      DO 38  L2=1,LL+1
         L22=L2-1
      DO 38  M2=1,L22
         IIS2=(N2-1)*LMAX + (L22*L22+(M2+1)+L22)
!         IIS2=(N2-1)*LMAX + L22*L22+M2+L22

         NU_V=LMAX*(IJ2-1)*NNMAX+IIS2  

         VLM=Matrix_Term_VV(I1,J1,I2,J2,LL,N1,L11,M1,N2,L22,M2,NNMAX,  &
             PPM,ALPHAM,BETAM,GAMMAM,NPP,NALPHA,NBETA,NGAMMA,NTEN,PAR)

!      CALL WW_WW_TEN(L11,M1,N1,L22,M,N,       &
!              PPM,RRM,NPP,HPP,                    &
!              NBETA,NGAMMA,NALPHA,                  &
!              ALPHAM,BETAM,GAMMAM,PAR,VLM,NI,NJ,NTEN)

      AM(MU_V,NU_V)=VLM

!      write(*,*) '1111111111111111111111111'
!      write(*,*) MU_V,NU_V
38    CONTINUE
40    CONTINUE
!-----------------------------------------------

!  sin * cos   (2,1)
      DO 60 J1=1,NTEN
      DO 60 I1=1,NTEN
         IJ1=(J1-1)*NTEN+I1
      DO 60 J2=1,NTEN
      DO 60 I2=1,NTEN
         IJ2=(J2-1)*NTEN+I2

      DO 58  N1=1,NNMAX
      DO 58  L1=1,LL+1
         L11=L1-1
      DO 58  M1=1,L11
         IIS1=(N1-1)*LMAX + (L11*L11+(M1+1)+L11)
!         IIS1=(N1-1)*LMAX + L11*L11+M1+L11

         MU_V=LMAX*(IJ1-1)*NNMAX+IIS1  

!      DO 58 J2=1,NTEN
!      DO 58 I2=1,NTEN
!         IJ2=(J2-1)*NTEN+I2
      DO 58  N2=1,NNMAX
      DO 58  L2=1,LL+1
         L22=L2-1
      DO 58  M2=1,L22+1
         IIC2=(N2-1)*LMAX + (L22*L22+M2)

         NU_W=LMAX*(IJ2-1)*NNMAX+IIC2  

         VLM=Matrix_Term_WV(I1,J1,I2,J2,LL,N1,L11,M1,N2,L22,M2-1,NNMAX,  &
             PPM,ALPHAM,BETAM,GAMMAM,NPP,NALPHA,NBETA,NGAMMA,NTEN,PAR)

!      CALL WW_WW_TEN(L11,M1,N1,L22,M-1,N,       &
!              PPM,RRM,NPP,HPP,                    &
!              NBETA,NGAMMA,NALPHA,                  &
!              ALPHAM,BETAM,GAMMAM,PAR,VLM,NI,NJ,NTEN)

         AM(MU_V,NU_W)=VLM
!      write(*,*) '1111111111111111111111111'
!      write(29,*) MU_V,NU_W
 58   CONTINUE
 60   CONTINUE
777   CONTINUE
      RETURN
      END SUBROUTINE Matrix_WW_Full
!---------------------------------------------------

      FUNCTION Matrix_Term_WW(IC1,JC1,IC2,JC2,LL,N1,L1,M1,N2,L2,M2,NNMAX,  &
           PPM,ALPHAM,BETAM,GAMMAM,NPP,NALPHA,NBETA,NGAMMA,NTEN,PAR)

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 
       
        INTEGER :: IC1,JC1,IC2,JC2
        INTEGER :: LL,N1,L1,M1,N2,L2,M2
        INTEGER :: NNMAX,NPP,NTEN !NJ,NPAR
        INTEGER :: NALPHA,NBETA,NGAMMA
!        REAL(wp) :: PROJ(NPP*NALPHA*NBETA*NGAMMA)
 !       REAL(wp) :: YM((LL+1)*(LL+1)*NNMAX,NTEN)
!        REAL(wp) :: YM(2,NTEN*NTEN)
!        REAL(wp) :: YM(2)
        REAL(wp) :: PPM(NPP),UPPM(NPP)
        REAL(wp) :: ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
        REAL(wp) :: PAR(50)

!local 
!        REAL(wp) :: BM1(NBETA*NGAMMA,NTEN*NTEN),BM2(NBETA*NGAMMA,NTEN*NTEN)
        REAL(wp) :: BM1(NBETA*NGAMMA)   !,BM2(NBETA*NGAMMA,NTEN,NTEN)
!        REAL(wp) :: CM1(NBETA,NTEN*NTEN), CM2(NBETA,NTEN*NTEN)
        REAL(wp) :: CM1(NBETA)
!        REAL(wp) :: WUnit(3)
        INTEGER :: LMAX,L22,IIC,IIS
        INTEGER :: I,J,IJ,I1,I2,I3,I23
!        REAL(wp) :: S1(NTEN,NTEN),S2(NTEN,NTEN),WW1,WW2,WW3
        REAL(wp) :: S1,WW1,WW2,WW3
        REAL(wp) :: ALPHA,BETA,GAMMA
        REAL(wp) :: HALPHA,HBETA,HGAMMA
        REAL(wp) :: YYRE
        REAL(wp) :: YM
!functions
        REAL(wp) :: Matrix_Term_WW, INT_Matrix_P_WW



!      DIMENSION BM1(NBETA*NGAMMA,3),BM2(NBETA,3)
!      DIMENSION WUnit(3)

!        write(*,*) '444444444',PROJ

      HALPHA=PAR(38)
      HBETA =PAR(39)
      HGAMMA=PAR(40)
!      UM1   =PAR(43)

      LMAX=(LL+1)*(LL+1)

 !--------------------------------
!  cos*cos

      DO 12 I2=1,NBETA
         BETA=BETAM(I2)

      DO 12 I3=1,NGAMMA
         GAMMA=GAMMAM(I3)
         I23=(I2-1)*NGAMMA+I3

          S1=0.
          DO 10 I1=1,NALPHA
          ALPHA=ALPHAM(I1)

          WW1=1.0
          IF(I1 .EQ.1 .OR. I1 .EQ. NALPHA) WW1=0.5

          YYRE=INT_Matrix_P_WW(PPM,NPP,N1,L1,M1,N2,L2,M2,IC1,JC1,IC2,JC2, &
                            ALPHA,BETA,GAMMA)

          S1=S1+YYRE*WW1
10       CONTINUE
         
          BM1(I23)=S1*HALPHA  
12    CONTINUE   

          DO 16 I2=1,NBETA

            S1=0.
            DO 14 I3=1,NGAMMA
               I23=(I2-1)*NGAMMA+I3
               
               WW3=1.0
               IF(I3 .EQ.1 .OR. I3 .EQ. NGAMMA) WW3=0.5
                  S1=S1+BM1(I23)*WW3   
14          CONTINUE
               CM1(I2)=S1*HGAMMA
16    CONTINUE

      S1=0.
      DO 18 I2=1,NBETA
         BETA=BETAM(I2)
         WW2=1.0
         IF(I2 .EQ.1 .OR. I2 .EQ. NBETA) WW2=0.5

            S1=S1+CM1(I2)*SIN(BETA)*WW2

18    CONTINUE

      YM=S1*HBETA
      Matrix_Term_WW=S1*HBETA
      RETURN
      END FUNCTION Matrix_Term_WW
!--------------------------------------
      FUNCTION Matrix_Term_WV(IC1,JC1,IC2,JC2,LL,N1,L1,M1,N2,L2,M2,NNMAX,  &
                PPM,ALPHAM,BETAM,GAMMAM,NPP,NALPHA,NBETA,NGAMMA,NTEN,PAR)

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 
       
        INTEGER :: IC1,JC1,IC2,JC2
        INTEGER :: LL,N1,L1,M1,N2,L2,M2
        INTEGER :: NNMAX,NJ,NPP,NTEN !,NPAR
        INTEGER :: NALPHA,NBETA,NGAMMA
!        REAL(wp) :: PROJ(NPP*NALPHA*NBETA*NGAMMA)
 !       REAL(wp) :: YM((LL+1)*(LL+1)*NNMAX,NTEN)
!        REAL(wp) :: YM(2,NTEN*NTEN)
!        REAL(wp) :: YM(2)
        REAL(wp) :: PPM(NPP),UPPM(NPP)
        REAL(wp) :: ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
        REAL(wp) :: PAR(50)

!local 
!        REAL(wp) :: BM1(NBETA*NGAMMA,NTEN*NTEN),BM2(NBETA*NGAMMA,NTEN*NTEN)
        REAL(wp) :: BM1(NBETA*NGAMMA)   !,BM2(NBETA*NGAMMA,NTEN,NTEN)
!        REAL(wp) :: CM1(NBETA,NTEN*NTEN), CM2(NBETA,NTEN*NTEN)
        REAL(wp) :: CM1(NBETA)
!        REAL(wp) :: WUnit(3)
        INTEGER :: LMAX,L22,IIC,IIS
        INTEGER :: I,J,IJ,I1,I2,I3,I23
!        REAL(wp) :: S1(NTEN,NTEN),S2(NTEN,NTEN),WW1,WW2,WW3
        REAL(wp) :: S1,WW1,WW2,WW3
        REAL(wp) :: ALPHA,BETA,GAMMA
        REAL(wp) :: HALPHA,HBETA,HGAMMA
        REAL(wp) :: YYRE
        REAL(wp) :: YM
!functions
        REAL(wp) :: Matrix_Term_WV, INT_Matrix_P_WV



!      DIMENSION BM1(NBETA*NGAMMA,3),BM2(NBETA,3)
!      DIMENSION WUnit(3)

!        write(*,*) '444444444',PROJ

      HALPHA=PAR(38)
      HBETA =PAR(39)
      HGAMMA=PAR(40)
!      UM1   =PAR(43)

      LMAX=(LL+1)*(LL+1)

 !--------------------------------
!  cos*sin

      DO 12 I2=1,NBETA
         BETA=BETAM(I2)

      DO 12 I3=1,NGAMMA
         GAMMA=GAMMAM(I3)
         I23=(I2-1)*NGAMMA+I3

          S1=0.
          DO 10 I1=1,NALPHA
          ALPHA=ALPHAM(I1)

          WW1=1.0
          IF(I1 .EQ.1 .OR. I1 .EQ. NALPHA) WW1=0.5

          YYRE=INT_Matrix_P_WV(PPM,NPP,N1,L1,M1,N2,L2,M2,IC1,JC1,IC2,JC2, &
                            ALPHA,BETA,GAMMA)

          S1=S1+YYRE*WW1
10       CONTINUE
         
          BM1(I23)=S1*HALPHA  
12    CONTINUE   

          DO 16 I2=1,NBETA

            S1=0.
            DO 14 I3=1,NGAMMA
               I23=(I2-1)*NGAMMA+I3
               
               WW3=1.0
               IF(I3 .EQ.1 .OR. I3 .EQ. NGAMMA) WW3=0.5
                  S1=S1+BM1(I23)*WW3   
14          CONTINUE
               CM1(I2)=S1*HGAMMA
16    CONTINUE

      S1=0.
      DO 18 I2=1,NBETA
         BETA=BETAM(I2)
         WW2=1.0
         IF(I2 .EQ.1 .OR. I2 .EQ. NBETA) WW2=0.5

            S1=S1+CM1(I2)*SIN(BETA)*WW2

18    CONTINUE

      YM=S1*HBETA
      Matrix_Term_WV=S1*HBETA
      RETURN
      END FUNCTION Matrix_Term_WV
!--------------------------------------
      FUNCTION Matrix_Term_VV(IC1,JC1,IC2,JC2,LL,N1,L1,M1,N2,L2,M2,NNMAX, &
                PPM,ALPHAM,BETAM,GAMMAM,NPP,NALPHA,NBETA,NGAMMA,NTEN,PAR)

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 
       
        INTEGER :: IC1,JC1,IC2,JC2
        INTEGER :: LL,N1,L1,M1,N2,L2,M2
        INTEGER :: NNMAX,NJ,NPP,NTEN !,NPAR
        INTEGER :: NALPHA,NBETA,NGAMMA
!        REAL(wp) :: PROJ(NPP*NALPHA*NBETA*NGAMMA)
 !       REAL(wp) :: YM((LL+1)*(LL+1)*NNMAX,NTEN)
!        REAL(wp) :: YM(2,NTEN*NTEN)
!        REAL(wp) :: YM(2)
        REAL(wp) :: PPM(NPP),UPPM(NPP)
        REAL(wp) :: ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
        REAL(wp) :: PAR(50)

!local 
!        REAL(wp) :: BM1(NBETA*NGAMMA,NTEN*NTEN),BM2(NBETA*NGAMMA,NTEN*NTEN)
        REAL(wp) :: BM1(NBETA*NGAMMA)   !,BM2(NBETA*NGAMMA,NTEN,NTEN)
!        REAL(wp) :: CM1(NBETA,NTEN*NTEN), CM2(NBETA,NTEN*NTEN)
        REAL(wp) :: CM1(NBETA)
!        REAL(wp) :: WUnit(3)
        INTEGER :: LMAX,L22,IIC,IIS
        INTEGER :: I,J,IJ,I1,I2,I3,I23
!        REAL(wp) :: S1(NTEN,NTEN),S2(NTEN,NTEN),WW1,WW2,WW3
        REAL(wp) :: S1,WW1,WW2,WW3
        REAL(wp) :: ALPHA,BETA,GAMMA
        REAL(wp) :: HALPHA,HBETA,HGAMMA
        REAL(wp) :: YYRE
        REAL(wp) :: YM
!functions
        REAL(wp) :: Matrix_Term_VV, INT_Matrix_P_VV



!      DIMENSION BM1(NBETA*NGAMMA,3),BM2(NBETA,3)
!      DIMENSION WUnit(3)

!        write(*,*) '444444444',PROJ

      HALPHA=PAR(38)
      HBETA =PAR(39)
      HGAMMA=PAR(40)
!      UM1   =PAR(43)

      LMAX=(LL+1)*(LL+1)

 !--------------------------------
!  cos*sin

      DO 12 I2=1,NBETA
         BETA=BETAM(I2)

      DO 12 I3=1,NGAMMA
         GAMMA=GAMMAM(I3)
         I23=(I2-1)*NGAMMA+I3

          S1=0.
          DO 10 I1=1,NALPHA
          ALPHA=ALPHAM(I1)

          WW1=1.0
          IF(I1 .EQ.1 .OR. I1 .EQ. NALPHA) WW1=0.5

          YYRE=INT_Matrix_P_VV(PPM,NPP,N1,L1,M1,N2,L2,M2,IC1,JC1,IC2,JC2, &
                            ALPHA,BETA,GAMMA)

          S1=S1+YYRE*WW1
10       CONTINUE
         
          BM1(I23)=S1*HALPHA  
12    CONTINUE   

          DO 16 I2=1,NBETA

            S1=0.
            DO 14 I3=1,NGAMMA
               I23=(I2-1)*NGAMMA+I3
               
               WW3=1.0
               IF(I3 .EQ.1 .OR. I3 .EQ. NGAMMA) WW3=0.5
                  S1=S1+BM1(I23)*WW3   
14          CONTINUE
               CM1(I2)=S1*HGAMMA
16    CONTINUE

      S1=0.
      DO 18 I2=1,NBETA
         BETA=BETAM(I2)
         WW2=1.0
         IF(I2 .EQ.1 .OR. I2 .EQ. NBETA) WW2=0.5

            S1=S1+CM1(I2)*SIN(BETA)*WW2

18    CONTINUE

      YM=S1*HBETA
      Matrix_Term_VV=S1*HBETA
      RETURN
      END FUNCTION Matrix_Term_VV
!--------------------------------------


      FUNCTION INT_Matrix_P_WW(PPM,NPP,N1,L1,M1,N2,L2,M2,IC1,JC1,IC2,JC2,     &
                             ALPHA,BETA,GAMMA)
!   integration over p variable W*W
!   
!        USE PARAM_INPUT
        USE PARAMTR
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 
       
        INTEGER :: NPP,N1,L1,M1,N2,L2,M2,IC1,JC1,IC2,JC2
!        INTEGER :: I1,I2,I3,NJ
        REAL(wp) :: ALPHA,BETA,GAMMA
!        REAL(wp) :: YYRE
        REAL(wp) :: PPM(NPP)
!        REAL(wp) :: PAR(50)
!local
!        REAL(wp) :: PROJ3(NPP) !,PROJ2(NPP),
        REAL(wp) :: PP
        REAL(wp) :: WW1,WWLM1,WWLM2
        REAL(wp) :: S1,HPP
!        INTEGER :: NALPHA,NBETA,NGAMMA
        INTEGER :: I4  !,I321,I4321
        REAL(wp) :: WUnit(3)
!functions
!        REAL(wp) :: FUNL1
        REAL(wp) :: WWC, INT_Matrix_P_WW
 
!      NALPHA=PARAM(14)
!      NBETA =PARAM(15)
!      NGAMMA=PARAM(16)
      HPP  =PAR(35)

!      DO 6 I4=1,NPP
!         I321=(I3-1)*NBETA*NALPHA+(I2-1)*NALPHA+I1
!         I4321=(I4-1)*NGAMMA*NBETA*NALPHA + I321
!         PROJ3(I4)=PROJ(I4321)
! 6    CONTINUE   

!      DO 8 I4=1,NPP
!         UPP=UPPM(I4)
!         PROJ3(I4)= FUNL1(PPM,PROJ2,1,NPP,UPP,IERR)
! 8    CONTINUE   

         WUnit(1)=-SIN(BETA)*COS(GAMMA)
         WUnit(2)= SIN(BETA)*SIN(GAMMA)
         WUnit(3)= COS(BETA)


      S1=0.
      DO 10 I4=1,NPP
         PP=PPM(I4)

         WWLM1= WWC(WUnit,IC1,JC1,N1,L1,M1,ALPHA,BETA,GAMMA,PP)
         WWLM2= WWC(WUnit,IC2,JC2,N2,L2,M2,ALPHA,BETA,GAMMA,PP)

         WW1=1.0
         IF(I4 .EQ.1 .OR. I4 .EQ. NPP) WW1=0.5
         S1=S1 + WWLM1 * WWLM2 * WW1

 10   CONTINUE   
         INT_Matrix_P_WW= S1*HPP
      RETURN
      END FUNCTION INT_Matrix_P_WW
!------------------------------------------------------------
      FUNCTION INT_Matrix_P_WV(PPM,NPP,N1,L1,M1,N2,L2,M2,IC1,JC1,IC2,JC2,     &
                             ALPHA,BETA,GAMMA)
!   integration over p variable W*V
!   
!        USE PARAM_INPUT
        USE PARAMTR
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 
       
        INTEGER :: NPP,N1,L1,M1,N2,L2,M2,IC1,JC1,IC2,JC2
        REAL(wp) :: ALPHA,BETA,GAMMA
        REAL(wp) :: PPM(NPP)
!local
        REAL(wp) :: PP
        REAL(wp) :: WW1,WWLM1,WWLM2
        REAL(wp) :: S1,HPP
        INTEGER :: I4 
        REAL(wp) :: WUnit(3) 
!functions
        REAL(wp) :: WWC,VVS, INT_Matrix_P_WV
 
        HPP  =PAR(35)

!      DO 8 I4=1,NPP
!         UPP=UPPM(I4)
!         PROJ3(I4)= FUNL1(PPM,PROJ2,1,NPP,UPP,IERR)
! 8    CONTINUE   

         WUnit(1)=-SIN(BETA)*COS(GAMMA)
         WUnit(2)= SIN(BETA)*SIN(GAMMA)
         WUnit(3)= COS(BETA)

      S1=0.
      DO 10 I4=1,NPP
         PP=PPM(I4)

         WWLM1= WWC(WUnit,IC1,JC1,N1,L1,M1,ALPHA,BETA,GAMMA,PP)
         WWLM2= VVS(WUnit,IC2,JC2,N2,L2,M2,ALPHA,BETA,GAMMA,PP)

         WW1=1.0
         IF(I4 .EQ.1 .OR. I4 .EQ. NPP) WW1=0.5
         S1=S1 + WWLM1 * WWLM2 * WW1

 10   CONTINUE   
         INT_Matrix_P_WV= S1*HPP
      RETURN
      END FUNCTION INT_Matrix_P_WV
!------------------------------------------------------------
      FUNCTION INT_Matrix_P_VV(PPM,NPP,N1,L1,M1,N2,L2,M2,IC1,JC1,IC2,JC2,     &
                             ALPHA,BETA,GAMMA)
!   integration over p variable V*V
!   
!        USE PARAM_INPUT
        USE PARAMTR
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 
       
        INTEGER :: NPP,N1,L1,M1,N2,L2,M2,IC1,JC1,IC2,JC2
        REAL(wp) :: ALPHA,BETA,GAMMA
        REAL(wp) :: PPM(NPP)
!local
        REAL(wp) :: PP
        REAL(wp) :: WW1,WWLM1,WWLM2
        REAL(wp) :: S1,HPP
        INTEGER :: I4  
        REAL(wp) :: WUnit(3)
!functions
        REAL(wp) :: VVS, INT_Matrix_P_VV
 
        HPP  =PAR(35)

!      DO 8 I4=1,NPP
!         UPP=UPPM(I4)
!         PROJ3(I4)= FUNL1(PPM,PROJ2,1,NPP,UPP,IERR)
! 8    CONTINUE   

         WUnit(1)=-SIN(BETA)*COS(GAMMA)
         WUnit(2)= SIN(BETA)*SIN(GAMMA)
         WUnit(3)= COS(BETA)

      S1=0.
      DO 10 I4=1,NPP
         PP=PPM(I4)

         WWLM1= VVS(WUnit,IC1,JC1,N1,L1,M1,ALPHA,BETA,GAMMA,PP)
         WWLM2= VVS(WUnit,IC2,JC2,N2,L2,M2,ALPHA,BETA,GAMMA,PP)

         WW1=1.0
         IF(I4 .EQ.1 .OR. I4 .EQ. NPP) WW1=0.5
         S1=S1 + WWLM1 * WWLM2 * WW1

 10   CONTINUE   
         INT_Matrix_P_VV= S1*HPP
      RETURN
      END FUNCTION INT_Matrix_P_VV
!------------------------------------------------------------

      SUBROUTINE Proj_Summa_W_Harm(LL,NJ,NNMAX,PPM,ALPHAM,BETAM, &
                        GAMMAM,NPP,NALPHA,NBETA,NGAMMA,XXM,NTEN,GM)
! вычисление проекций по формуле (15)
!  summation of W  harmonics 
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 

        INTEGER :: LL,NJ,NNMAX,NPP,NALPHA,NBETA,NGAMMA
        INTEGER :: NTEN
        REAL(wp) :: ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
        REAL(wp) :: PPM(NPP)
        REAL(wp) :: XXM((LL+1)*(LL+1)*NNMAX,NTEN*NTEN)
        REAL(wp) :: GM(NPP*NALPHA*NBETA*NGAMMA)
!local 
        INTEGER :: II,I1,I2,I3,I321
        INTEGER :: I4,I4321,I,J,IJ
        INTEGER :: N,L,L22,M,IIC,IIS,LMAX
        REAL(wp) :: ALPHA,BETA,GAMMA
        REAL(wp) :: WWLMRE,WWLMIM
        REAL(wp) :: S(NTEN),PP,S1
        REAL(wp) :: WUnit(3)
!functions 
        REAL(wp) :: WWC,VVS

      LMAX=(LL+1)*(LL+1)

      DO 10 I4=1,NPP
         PP=PPM(I4)
      DO 10 II=1,NJ
         CALL IALBETGAM(II,NALPHA,NBETA,I1,I2,I3)

         ALPHA=ALPHAM(I1)
         BETA =BETAM(I2)
         GAMMA=GAMMAM(I3)
         I321=(I3-1)*NBETA*NALPHA+(I2-1)*NALPHA+I1
         I4321=(I4-1)*NGAMMA*NBETA*NALPHA + I321

         WUnit(1)=-SIN(BETA)*COS(GAMMA)
         WUnit(2)= SIN(BETA)*SIN(GAMMA)
         WUnit(3)= COS(BETA)


!      DO 10 I4=1,NPP
!            PP=PPM(I4)
!         I4321=(II-1)*NPP + I4
!         write(*,*) '11111111111',ALPHA,BETA,GAMMA

!         WUnit(1)=COS(ALPHA)*SIN(BETA)
!         WUnit(2)=SIN(ALPHA)*SIN(BETA)
!         WUnit(3)=COS(BETA)

        S1=0.
        DO 4  N=1,NNMAX
        DO 4  L=1,LL+1
           L22=L-1
        DO 4  M=1,L22+1
           IIC=(N-1)*LMAX + (L22*L22+M)

!         WWLMRE= WWC(IC,JC,N,L22,M-1,ALPHA,BETA,GAMMA,PP)

           DO 5 J=1,NTEN
           DO 5 I=1,NTEN   
              IJ=(J-1)*NTEN+I
!              S1=S1+XXM(IIC,IJ)*WWLMRE*WUnit(I)*WUnit(J)
              WWLMRE=WWC(WUnit,I,J,N,L22,M-1,ALPHA,BETA,GAMMA,PP)
              S1=S1+XXM(IIC,IJ)*WWLMRE
5          CONTINUE
4       CONTINUE

        DO 6 N=1,NNMAX
        DO 6 L=1,LL+1
           L22=L-1
        DO 6 M=1,L22
           IIS=(N-1)*LMAX + L22*L22+L22+M+1

!           WWLMIM= VVS(N,L22,M,ALPHA,BETA,GAMMA,PP)
!           CALL WWCS(N,L22,M,ALPHA,BETA,GAMMA,PP,WWLMRE,WWLMIM)
!           write(*,*) 'Proj_Summa_W_Harm:IMAGE',WWLMIM

           DO 9 J=1,NTEN
           DO 9 I=1,NTEN   
              IJ=(J-1)*NTEN+I
              WWLMIM= VVS(WUnit,I,J,N,L22,M,ALPHA,BETA,GAMMA,PP)
!              S1=S1+XXM(IIS,IJ)*WWLMIM*WUnit(I)*WUnit(J)
              S1=S1+XXM(IIS,IJ)*WWLMIM

9          CONTINUE
6       CONTINUE
        GM(I4321)=S1
10    CONTINUE
      RETURN
      END SUBROUTINE Proj_Summa_W_Harm
!------------------------------------------------------

      SUBROUTINE RightPart_Full(LL,NNMAX,NTEN,     &
                    PROJ,PPM,ALPHAM,BETAM,GAMMAM,  &
                    NPP,NALPHA,NBETA,NGAMMA,PAR,YMM)

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

        INTEGER :: LL,NNMAX,NTEN
        REAL(wp) :: YMM((LL+1)**2*NNMAX*NTEN*NTEN)
        REAL(wp) :: PROJ(NPP*NALPHA*NBETA*NGAMMA)
        REAL(wp) :: PPM(NPP),ALPHAM(NALPHA)
        REAL(wp) ::BETAM(NBETA),GAMMAM(NGAMMA)
        INTEGER :: NPP,NALPHA,NBETA,NGAMMA
        REAL(wp) :: PAR(50) !,WUnit(3)
!local
        INTEGER :: IC,JC,N,L,M,MU_W,MU_V,I,J,IJ
        INTEGER :: LMAX,NJ,L22,IIC,IIS
        REAL(wp) :: YM
!functions
        REAL(wp) :: RightPart_Term_RE,RightPart_Term_IM

        LMAX=(LL+1)*(LL+1)
!        NJ=NALPHA*NBETA*NGAMMA

!--------------------------------
!cos

        DO 60 J=1,NTEN
        DO 60 I=1,NTEN
           IJ=(J-1)*NTEN+I

        DO 20  N=1,NNMAX
        DO 20  L=1,LL+1
           L22=L-1
        DO 20  M=1,L22+1
           IIC=(N-1)*LMAX + (L22*L22+M)
           MU_W=LMAX*(IJ-1)*NNMAX+IIC  
!           YM= RightPart_Term_RE(PROJ,IC,JC,LL,N,L22,M-1,NNMAX,NJ,  &
!                PPM,ALPHAM,BETAM,GAMMAM,        &
!                NPP,NALPHA,NBETA,NGAMMA,             &
!                NTEN,PAR)
!           DO 20 J=1,NTEN
!           DO 20 I=1,NTEN
!           IJ=(J-1)*NTEN+I
!           MU_W=LMAX*(IJ-1)*NNMAX+IIC  

           YM= RightPart_Term_RE(PROJ,I,J,LL,N,L22,M-1,NNMAX,  &
                PPM,ALPHAM,BETAM,GAMMAM,        &
                NPP,NALPHA,NBETA,NGAMMA,NTEN,PAR)

           YMM(MU_W)=YM
20         CONTINUE
!30      CONTINUE
!------------------------------
!  sin
!        DO 60 J=1,NTEN
!        DO 60 I=1,NTEN
!           IJ=(J-1)*NTEN+I
        DO 50  N=1,NNMAX
        DO 50  L=1,LL+1
           L22=L-1
        DO 50  M=1,L22
           IIS=(N-1)*LMAX + L22*L22+L22+M+1
           MU_V=LMAX*(IJ-1)*NNMAX+IIS
!           write(*,*) '111111', MU_V
!           YM=RightPart_Term_IM(PROJ,IC,JC,LL,N,L22,M,NNMAX,NJ,  &
!                PPM,ALPHAM,BETAM,GAMMAM,        &
!                NPP,NALPHA,NBETA,NGAMMA,             &
!                NTEN,PAR)

!           DO 40 J=1,NTEN
!           DO 40 I=1,NTEN
!           IJ=(J-1)*NTEN+I
!           MU_V=LMAX*(IJ-1)*NNMAX+IIS

           YM=RightPart_Term_IM(PROJ,I,J,LL,N,L22,M,NNMAX,  &
                PPM,ALPHAM,BETAM,GAMMAM,        &
                NPP,NALPHA,NBETA,NGAMMA,NTEN,PAR)

              YMM(MU_V)=YM
50         CONTINUE
!           write(*,*) 'RightPart_Full:M=', IJ
60      CONTINUE   
      RETURN
      END SUBROUTINE RightPart_Full
!--------------------------------------------------

     FUNCTION RightPart_Term_RE(PROJ,IC,JC,LL,N,L,M,NNMAX,  &
                PPM,ALPHAM,BETAM,GAMMAM,        &
                NPP,NALPHA,NBETA,NGAMMA,         &
                NTEN,PAR)

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 
       
        INTEGER :: IC,JC,LL,N,L,M,NNMAX,NJ,NPP,NTEN !,NPAR
        INTEGER :: NALPHA,NBETA,NGAMMA
        REAL(wp) :: PROJ(NPP*NALPHA*NBETA*NGAMMA)
        REAL(wp) :: PPM(NPP),UPPM(NPP)
        REAL(wp) :: ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
        REAL(wp) :: PAR(50)  !,WUnit(3)
!local 
        REAL(wp) :: BM1(NBETA*NGAMMA)
        REAL(wp) :: CM1(NBETA)
        INTEGER :: LMAX,L22,IIC,IIS
        INTEGER :: I,J,IJ,I1,I2,I3,I23
        REAL(wp) :: S1,WW1,WW2,WW3
        REAL(wp) :: ALPHA,BETA,GAMMA
        REAL(wp) :: HALPHA,HBETA,HGAMMA
        REAL(wp) :: YYRE
        REAL(wp) :: YM
!functions
        REAL(wp) :: INT_Proj_GNLM_RE,INT_Proj_GNLM_IM,RightPart_Term_RE

      HALPHA=PAR(38)
      HBETA =PAR(39)
      HGAMMA=PAR(40)
!      UM1   =PAR(43)

      LMAX=(LL+1)*(LL+1)

 !--------------------------------
!  cos
          DO 12 I2=1,NBETA
             BETA=BETAM(I2)
          DO 12 I3=1,NGAMMA
             GAMMA=GAMMAM(I3)
             I23=(I2-1)*NGAMMA+I3
!----
          S1=0.
          DO 10 I1=1,NALPHA
             ALPHA=ALPHAM(I1)

             WW1=1.0
             IF(I1 .EQ.1 .OR. I1 .EQ. NALPHA) WW1=0.5

             YYRE=INT_Proj_GNLM_RE(PROJ,PPM,NPP,N,L,M,IC,JC, &
                      I1,I2,I3,NJ,ALPHA,BETA,GAMMA,PAR)

!            YYIM=INT_Proj_GNLM_IM(PROJ,PPM,NPP,N,L,M,IC,JC, &
!                      I1,I2,I3,NJ,ALPHA,BETA,GAMMA,PAR)

!            CALL INT_Proj_GNLM_P(PROJ,PPM,NPP,N,L,M,I1,I2,I3,     &
!                      NJ,ALPHA,BETA,GAMMA,PAR,YYRE,YYIM)

!            CALL INT_Proj_GNLM_P(PROJ,PPM,UPPM,NPP,         &
!                             N,L22,M-1,                          & 
!                             ALPHA,BETA,GAMMA,                   &
!                             NJ,I1,I2,I3,PAR,YYRE,YYIM)

             S1=S1+YYRE*WW1
10        CONTINUE
!            BM1(I23,IC,JC)=S1*HALPHA
             BM1(I23)=S1*HALPHA
12        CONTINUE   
!----
          DO 16 I2=1,NBETA
          S1=0.
          DO 14 I3=1,NGAMMA
             I23=(I2-1)*NGAMMA+I3
             WW3=1.0
             IF(I3 .EQ.1 .OR. I3 .EQ. NGAMMA) WW3=0.5
!                  S1=S1+BM1(I23,IC,JC)*WW3
             S1=S1+BM1(I23)*WW3
14        CONTINUE
             CM1(I2)=S1*HGAMMA
16        CONTINUE
!----
          S1=0.
          DO 18 I2=1,NBETA
             BETA=BETAM(I2)
             WW2=1.0
             IF(I2 .EQ.1 .OR. I2 .EQ. NBETA) WW2=0.5
             S1=S1+CM1(I2)*SIN(BETA)*WW2
18        CONTINUE
!            YM=S1*HBETA
      RightPart_Term_RE=S1*HBETA
      RETURN
      END FUNCTION RightPart_Term_RE
!--------------------------------------
     FUNCTION RightPart_Term_IM(PROJ,IC,JC,LL,N,L,M,NNMAX,  &
                PPM,ALPHAM,BETAM,GAMMAM,        &
                NPP,NALPHA,NBETA,NGAMMA,         &
                NTEN,PAR)

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 
       
        INTEGER :: IC,JC,LL,N,L,M,NNMAX,NJ,NPP,NTEN !,NPAR
        INTEGER :: NALPHA,NBETA,NGAMMA
        REAL(wp) :: PROJ(NPP*NALPHA*NBETA*NGAMMA)
 !       REAL(wp) :: YM((LL+1)*(LL+1)*NNMAX,NTEN)
!        REAL(wp) :: YM(2,NTEN*NTEN)
!        REAL(wp) :: YM(2)
        REAL(wp) :: PPM(NPP),UPPM(NPP)
        REAL(wp) :: ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
        REAL(wp) :: PAR(50) !,WUnit(3)

!local 
!        REAL(wp) :: BM1(NBETA*NGAMMA,NTEN*NTEN),BM2(NBETA*NGAMMA,NTEN*NTEN)
        REAL(wp) :: BM1(NBETA*NGAMMA)  !,BM2(NBETA*NGAMMA,NTEN,NTEN)
!        REAL(wp) :: CM1(NBETA,NTEN*NTEN), CM2(NBETA,NTEN*NTEN)
        REAL(wp) :: CM1(NBETA), CM2(NBETA)
!        REAL(wp) :: WUnit(3)
        INTEGER :: LMAX,L22,IIC,IIS
        INTEGER :: I,J,IJ,I1,I2,I3,I23
!        REAL(wp) :: S1(NTEN,NTEN),S2(NTEN,NTEN),WW1,WW2,WW3
        REAL(wp) :: S1,WW1,WW2,WW3
        REAL(wp) :: ALPHA,BETA,GAMMA
        REAL(wp) :: HALPHA,HBETA,HGAMMA
        REAL(wp) :: YYIM
!        REAL(wp) :: YM
!functions
        REAL(wp) :: INT_Proj_GNLM_RE,INT_Proj_GNLM_IM,RightPart_Term_IM

      HALPHA=PAR(38)
      HBETA =PAR(39)
      HGAMMA=PAR(40)
!      UM1   =PAR(43)

      LMAX=(LL+1)*(LL+1)

 !--------------------------------
!  cos
          DO 12 I2=1,NBETA
             BETA=BETAM(I2)

          DO 12 I3=1,NGAMMA
             GAMMA=GAMMAM(I3)
             I23=(I2-1)*NGAMMA+I3

          S1=0.
          DO 10 I1=1,NALPHA
             ALPHA=ALPHAM(I1)
             WW1=1.0
             IF(I1 .EQ.1 .OR. I1 .EQ. NALPHA) WW1=0.5

             YYIM=INT_Proj_GNLM_IM(PROJ,PPM,NPP,N,L,M,IC,JC, &
                      I1,I2,I3,NJ,ALPHA,BETA,GAMMA,PAR)

!            CALL INT_Proj_GNLM_P(PROJ,PPM,NPP,N,L,M,I1,I2,I3,     &
!                      NJ,ALPHA,BETA,GAMMA,PAR,YYRE,YYIM)

!            CALL INT_Proj_GNLM_P(PROJ,PPM,UPPM,NPP,         &
!                             N,L22,M-1,                          & 
!                             ALPHA,BETA,GAMMA,                   &
!                             NJ,I1,I2,I3,PAR,YYRE,YYIM)
            S1=S1+YYIM*WW1
10       CONTINUE
         
!            BM1(I23,IC,JC)=S1*HALPHA
            BM1(I23)=S1*HALPHA
12    CONTINUE   

      DO 16 I2=1,NBETA

      S1=0.
      DO 14 I3=1,NGAMMA
         I23=(I2-1)*NGAMMA+I3
         WW3=1.0
         IF(I3 .EQ.1 .OR. I3 .EQ. NGAMMA) WW3=0.5
!                  S1=S1+BM1(I23,IC,JC)*WW3
         S1=S1+BM1(I23)*WW3
         
14    CONTINUE
         CM1(I2)=S1*HGAMMA
16    CONTINUE

      S1=0.
      DO 18 I2=1,NBETA
         BETA=BETAM(I2)
         WW2=1.0
         IF(I2 .EQ.1 .OR. I2 .EQ. NBETA) WW2=0.5
            S1=S1+CM1(I2)*SIN(BETA)*WW2
18    CONTINUE
!            YM=S1*HBETA
      RightPart_Term_IM=S1*HBETA
      RETURN
      END FUNCTION RightPart_Term_IM
!--------------------------------------

      FUNCTION INT_Proj_GNLM_RE(PROJ,PPM,NPP,N1,L1,M1,IC,JC,     &
                             I1,I2,I3,NJ,ALPHA,BETA,GAMMA,PAR)
!   integration over p variable of REAL part
!   g(p,al,bet,gam)*WWC(n,l,m)

        USE PARAM_INPUT

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

       
        INTEGER :: NPP,N1,L1,M1,IC,JC
        INTEGER :: I1,I2,I3,NJ
        REAL(wp) :: ALPHA,BETA,GAMMA
!        REAL(wp) :: YYRE
        REAL(wp) :: PPM(NPP)
        REAL(wp) :: PROJ(NPP*NJ)
        REAL(wp) :: PAR(50)

!local
        REAL(wp) :: PROJ3(NPP) !,PROJ2(NPP),
        REAL(wp) :: PP
        REAL(wp) :: WW1,WWLMRE,WWLMIM
        REAL(wp) :: S1,S2,HUPP
        INTEGER :: NALPHA,NBETA,NGAMMA
        INTEGER :: I4,I321,I4321
        REAL(wp) :: WUnit(3)
!functions
!        REAL(wp) :: FUNL1
        REAL(wp) :: WWC,INT_Proj_GNLM_RE

      NALPHA=PARAM(14)
      NBETA =PARAM(15)
      NGAMMA=PARAM(16)
      HUPP  =PAR(35)

      DO 6 I4=1,NPP
         I321=(I3-1)*NBETA*NALPHA+(I2-1)*NALPHA+I1
         I4321=(I4-1)*NGAMMA*NBETA*NALPHA + I321
         PROJ3(I4)=PROJ(I4321)
 6    CONTINUE   

!      I321=(I3-1)*NBETA*NALPHA+(I2-1)*NALPHA+I1
!      DO 6 I4=1,NPP
!         I4321=(I321-1)*NPP + I4
!         PROJ3(I4)=PROJ(I4321)
! 6    CONTINUE  

        WUnit(1)=-SIN(BETA)*COS(GAMMA)
        WUnit(2)= SIN(BETA)*SIN(GAMMA)
        WUnit(3)= COS(BETA)

      S1=0.
      S2=0.
      DO 10 I4=1,NPP
         PP=PPM(I4)

         WWLMRE= WWC(WUnit,IC,JC,N1,L1,M1,ALPHA,BETA,GAMMA,PP)
!         WWLMIM= VVS(N1,L1,M1,ALPHA,BETA,GAMMA,PP)
!         CALL WWCS(N1,L1,M1,ALPHA,BETA,GAMMA,UPP,WWLMRE,WWLMIM)

         WW1=1.0
         IF(I4 .EQ.1 .OR. I4 .EQ. NPP) WW1=0.5
         S1=S1 + WWLMRE * PROJ3(I4)*WW1
!         S2=S2 + WWLMIM * PROJ3(I4)*WW1

 10   CONTINUE   
         INT_Proj_GNLM_RE= S1*HUPP
!      YYIM= S2*HUPP

!      write(*,*) '2222222222',YYRE,YYIM
      RETURN
      END FUNCTION INT_Proj_GNLM_RE
!------------------------------------------------------------
       FUNCTION INT_Proj_GNLM_IM(PROJ,PPM,NPP,N1,L1,M1,IC,JC,     &
                        I1,I2,I3,NJ,ALPHA,BETA,GAMMA,PAR)
!   integration over p variable of IMARGE part
!   g(p,al,bet,gam)*g(n,l,m)

        USE PARAM_INPUT

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

       
        INTEGER :: NPP,N1,L1,M1,IC,JC
        INTEGER :: I1,I2,I3,NJ
        REAL(wp) :: ALPHA,BETA,GAMMA
!        REAL(wp) :: YYRE
        REAL(wp) :: PPM(NPP)
        REAL(wp) :: PROJ(NPP*NJ)
        REAL(wp) :: PAR(50)
!local
        REAL(wp) :: PROJ3(NPP) !,PROJ2(NPP),
        REAL(wp) :: PP
        REAL(wp) :: WW1,WWLMRE,WWLMIM
        REAL(wp) :: S1,S2,HUPP
        INTEGER :: NALPHA,NBETA,NGAMMA
        INTEGER :: I4,I321,I4321
        REAL(wp) :: WUnit(3)
!functions
!        REAL(wp) :: FUNL1
        REAL(wp) :: VVS,INT_Proj_GNLM_IM

      NALPHA=PARAM(14)
      NBETA =PARAM(15)
      NGAMMA=PARAM(16)
      HUPP  =PAR(35)

      DO 6 I4=1,NPP
         I321=(I3-1)*NBETA*NALPHA+(I2-1)*NALPHA+I1
         I4321=(I4-1)*NGAMMA*NBETA*NALPHA + I321
         PROJ3(I4)=PROJ(I4321)
 6    CONTINUE   

!      I321=(I3-1)*NBETA*NALPHA+(I2-1)*NALPHA+I1
!      DO 6 I4=1,NPP
!         I4321=(I321-1)*NPP + I4
!         PROJ3(I4)=PROJ(I4321)
! 6    CONTINUE  

        WUnit(1)=-SIN(BETA)*COS(GAMMA)
        WUnit(2)= SIN(BETA)*SIN(GAMMA)
        WUnit(3)= COS(BETA)

      S1=0.
      S2=0.
      DO 10 I4=1,NPP
         PP=PPM(I4)

!         WWLMRE= WWC(N1,L1,M1,ALPHA,BETA,GAMMA,PP)
         WWLMIM= VVS(WUnit,IC,JC,N1,L1,M1,ALPHA,BETA,GAMMA,PP)
!         CALL WWCS(N1,L1,M1,ALPHA,BETA,GAMMA,UPP,WWLMRE,WWLMIM)

         WW1=1.0
         IF(I4 .EQ.1 .OR. I4 .EQ. NPP) WW1=0.5
         S1=S1 + WWLMIM * PROJ3(I4)*WW1
!         S2=S2 + WWLMIM * PROJ3(I4)*WW1

 10   CONTINUE   
      INT_Proj_GNLM_IM= S1*HUPP
!      YYIM= S2*HUPP

!      write(*,*) '2222222222',YYRE,YYIM
      RETURN
      END FUNCTION INT_Proj_GNLM_IM
!------------------------------------------------------------

      FUNCTION WWC(WUnit,IC,JC,N,L,M,ALPHA,BETA,GAMMA,PP)
! culculation  W(nlm;alha,beta,gamma)*U(nlm) -functions
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 

        REAL(wp) :: WUnit(3)
        INTEGER :: IC,JC,N,L,M
        REAL(wp) :: ALPHA,BETA,GAMMA,PP
! local
        INTEGER ::I,M1
        REAL(wp) :: S1,GGLM
        REAL(wp) :: VLMRE
! functions
        REAL(wp) :: GNLMP,WLMM_RE
        REAL(wp) :: WWC

!        WUnit(1)=-SIN(BETA)*COS(GAMMA)
!        WUnit(2)= SIN(BETA)*SIN(GAMMA)
!        WUnit(3)= COS(BETA)

        S1=0.
        DO 10 I=1,L+1
           M1=I-1

           VLMRE= WLMM_RE(L,M,M1,ALPHA,BETA,GAMMA)
           GGLM= GNLMP(N,L,M1,PP)

           S1=S1+GGLM*VLMRE
10      CONTINUE
      WWC=S1*WUnit(IC)*WUnit(JC)

      RETURN
      END FUNCTION WWC
!------------------------------------------------
      FUNCTION VVS(WUnit,IC,JC,N,L,M,ALPHA,BETA,GAMMA,PP)
! culculation  V(nlm;alha,beta,gamma)*U(lmn) -functions
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 

        REAL(wp) :: WUnit(3)
        INTEGER :: IC,JC,N,L,M
        REAL(wp) :: ALPHA,BETA,GAMMA,PP
! local
        INTEGER ::I,M1
        REAL(wp) :: S1,GGLM,VLMIM
! functions
        REAL(wp) :: GNLMP,VLMM_IM
        REAL(wp) :: VVS

!        WUnit(1)=-SIN(BETA)*COS(GAMMA)
!        WUnit(2)= SIN(BETA)*SIN(GAMMA)
!        WUnit(3)= COS(BETA)

        S1=0.
        DO 10 I=1,L+1
           M1=I-1

           VLMIM=VLMM_IM(L,M,M1,ALPHA,BETA,GAMMA)
           GGLM= GNLMP(N,L,M1,PP)

           S1=S1+GGLM*VLMIM

10      CONTINUE
      VVS=S1*WUnit(IC)*WUnit(JC)
      RETURN
      END FUNCTION VVS
!------------------------------------------------

      FUNCTION WLMM_RE(L,M,M1,ALPHA,BETA,GAMMA)

!          Computer the W harmonics

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

        INTEGER :: L,M,M1
        REAL(wp) :: ALPHA,BETA,GAMMA
!        REAL(wp) :: VLMRE
! local
        INTEGER ::KK
        REAL(wp) :: DM,SIGNUM,CC1,CC2
        REAL(wp) :: PLM1,PLM2
! functions
        REAL(wp) :: PLM1M2, WLMM_RE

        DM=1.
        IF(M1 .EQ. 0) DM=0.5

        KK=MOD(M1,2)
        SIGNUM=1.
        IF(KK .NE. 0) SIGNUM=-1.

        CC1=COS(M*ALPHA + M1*GAMMA)
        CC2=COS(M*ALPHA - M1*GAMMA)

        PLM1 = PLM1M2(L,M, M1,BETA)
        PLM2 = PLM1M2(L,M,-M1,BETA)

        WLMM_RE=DM*(CC1*PLM1 + SIGNUM*CC2*PLM2)

      RETURN
      END FUNCTION WLMM_RE
!----------------------------------------------
      FUNCTION VLMM_IM(L,M,M1,ALPHA,BETA,GAMMA)

!          Computer the V harmonics

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

        INTEGER :: L,M,M1
        REAL(wp) :: ALPHA,BETA,GAMMA
! local
        INTEGER ::KK
        REAL(wp) :: DM,SIGNUM,SS1,SS2
        REAL(wp) :: PLM1,PLM2
! functions
        REAL(wp) :: PLM1M2, VLMM_IM

        DM=1.
        IF(M1 .EQ. 0) DM=0.5

        KK=MOD(M1,2)
        SIGNUM=1.
        IF(KK .NE. 0) SIGNUM=-1.

        SS1= SIN(M*ALPHA + M1*GAMMA)
        SS2= SIN(M*ALPHA - M1*GAMMA)

        PLM1 = PLM1M2(L,M, M1,BETA)
        PLM2 = PLM1M2(L,M,-M1,BETA)

        VLMM_IM=DM*(SS1*PLM1 + SIGNUM*SS2*PLM2)

      RETURN
      END FUNCTION VLMM_IM
!----------------------------------------------

      SUBROUTINE Sort_Matr(AAM,LL,LN,DDM,NTEN,NNMAX,NII,NJJ)
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 

        INTEGER :: LL,NTEN,NII,NJJ,NNMAX
        REAL(wp) :: AAM(LN,LN),DDM(NTEN*LN,NTEN*LN)
!local 
        INTEGER :: I,J,II,JJ,III,JJJ,LN

!      DIMENSION AAM(LN,LN)
!      DIMENSION DDM(3*LN,3*LN)
     
!      LN=(LL+1)*(LL+1)*NNMAX

        DO 4 J=1,LN
           JJ=(NJJ-1)*NTEN*LN+J
        DO 4 I=1,LN
           II=(NII-1)*NTEN*LN+I
           DDM(II,JJ)=AAM(I,J)
!         write(*,*) 'Sort_Matr',II,JJ
4       CONTINUE
6     CONTINUE
      RETURN
      END SUBROUTINE Sort_Matr
!--------------------------------------------------------

      FUNCTION INDEX_COS(LL,L,M,N,NNMAX,NI,NJ)

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

        INTEGER :: LL,L,N,NNMAX,NI,NJ
!local
        INTEGER :: LMAX 
        INTEGER :: M
        INTEGER :: IIC
!functions 
        REAL(wp) :: INDEX_COS

        LMAX=(LL+1)*(LL+1)

!   cos 

        IIC=(N-1)*LMAX + L*L+M+1

        INDEX_COS=IIC
      RETURN
      END FUNCTION INDEX_COS
       
!---------------------------------------------------

      FUNCTION INDEX_SIN(LL,L,M,N,NNMAX,NI,NJ)

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

        INTEGER :: LL,L,NNMAX,NI,NJ

!local
        INTEGER :: LMAX 
        INTEGER :: N,M
        INTEGER :: IIS
!functions 
        REAL(wp) :: INDEX_SIN

        LMAX=(LL+1)*(LL+1)

!   sin
        IF(M .NE. 0) THEN  
           IIS=(N-1)*LMAX + L*L+L+M+1
        ENDIF
   
        INDEX_SIN=IIS
      RETURN
      END FUNCTION INDEX_SIN
!------------------------------------------------- 



      FUNCTION VLMM_IM_TEN(L,M,N,ALPHA,BETA,GAMMA,PP)

!          Computer the V harmonics
!           image part
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

        INTEGER :: L,M,N !,NJ1,NJ2
        REAL(wp) :: ALPHA,BETA,GAMMA,PP
! local
        INTEGER ::K,KK,I
        REAL(wp) :: DM,SIGNUM,SS1,SS2
        REAL(wp) :: PLM1,PLM2
        REAL(wp) :: VLMIM,S,GN
! functions
        REAL(wp) :: PLM1M2,GNLMP,VLMM_IM_TEN 

        DM=0.5
        S=0.
        DO 4 I=1,L
           K=I-1
        IF(K .NE. 0) DM=1.

        KK=MOD(K,2)
        SIGNUM=1.
        IF(KK .NE. 0) SIGNUM=-1.

        SS1= SIN(M*ALPHA + K*GAMMA)
        SS2= SIN(M*ALPHA - K*GAMMA)

        PLM1 = PLM1M2(L,M, K,BETA)
        PLM2 = PLM1M2(L,M,-K,BETA)

        VLMIM=DM*(SS1*PLM1 + SIGNUM*SS2*PLM2)
        GN=GNLMP(N,L,K,PP)

        S=S+VLMIM*GN
4       CONTINUE
        VLMM_IM_TEN=S 
      RETURN
      END FUNCTION VLMM_IM_TEN
!----------------------------------------------


      SUBROUTINE INT_VLMM_VEC(L1,M1,K1,L2,M2,K2,ALPHAM,BETAM,GAMMAM, &
        NBETA,NGAMMA,NALPHA,NPAR,PAR,VLMINT,NCP1,NCP2)

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 

        INTEGER :: L1,M1,K1,L2,M2,K2
        INTEGER :: LL1,MM1,KK1,LL2,MM2,KK2
        INTEGER :: NALPHA,NBETA,NGAMMA
        INTEGER :: NPP,NCP1,NCP2
        REAL(wp) :: ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
        REAL(wp) :: VLMINT(4)
        INTEGER :: NPAR(100)  !!!!!
        REAL(wp) :: PAR(50)   !!!!!!!
!local
        INTEGER :: I1,I2,I3,I23
        REAL(wp) :: S1,S2,S3,S4
        REAL(wp) :: ALPHA,BETA,GAMMA,HALPHA,HBETA,HGAMMA
        REAL(wp) :: WW1,WW2,WW3
        REAL(wp) :: VLCC,VLSC,VLCS,VLSS
        REAL(wp) :: VLMRE1,VLMIM1,VLMRE2,VLMIM2
        REAL(wp) :: WLMRE1,WLMIM1,WLMRE2,WLMIM2
        REAL(wp) :: WUnit(3)
        REAL(wp) :: BM1(NBETA*NGAMMA,4),BM2(NBETA,4)   !!!!!!!!
!functions
        REAL(wp) :: WLMM_RE,WLMM_IM


!      DIMENSION ALPHAM(*),BETAM(*),GAMMAM(*)
!      DIMENSION NPAR(*),PAR(*)
!      DIMENSION BM1(NBETA*NGAMMA,4),VLMINT(4)
!      DIMENSION BM2(NBETA,4)
!      DIMENSION WUnit(3)

!      NALPHA=NPAR(9)
!      NBETA =NPAR(10)
!      NGAMMA=NPAR(18)

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

         WUnit(1)=COS(ALPHA)*SIN(BETA)
         WUnit(2)=SIN(ALPHA)*SIN(BETA)
         WUnit(3)=COS(BETA)

!         WUnit(1)=-SIN(BETA)*COS(GAMMA)
!         WUnit(2)= SIN(BETA)*SIN(GAMMA)
!         WUnit(3)= COS(BETA)

         WW1=1.0
         IF(I1 .EQ.1 .OR. I1 .EQ. NALPHA) WW1=0.5
! correspond to i
!         VLMRE1=WLMM_RE(L1,M1,K1,ALPHA,BETA,GAMMA,NCP1)
! correspond to j
!         VLMRE2=WLMM_RE(L2,M2,K2,ALPHA,BETA,GAMMA,NCP2)

!         VLMIM1=WLMM_IM(L1,M1,K1,ALPHA,BETA,GAMMA,NCP1)
!         VLMIM2=WLMM_IM(L2,M2,K2,ALPHA,BETA,GAMMA,NCP2)

!         CALL VLMM_RE (L1,M1,K1,ALPHA,BETA,GAMMA,VLMRE1,VLMIM1)
!         CALL VLMM_RE (L2,M2,K2,ALPHA,BETA,GAMMA,VLMRE2,VLMIM2)

         VLCC=VLMRE1*VLMRE2 !*WUnit(NCP1)*WUnit(NCP2)
         VLSC=VLMIM1*VLMRE2 !*WUnit(NCP1)*WUnit(NCP2)
         VLCS=VLMRE1*VLMIM2 !*WUnit(NCP1)*WUnit(NCP2)
         VLSS=VLMIM1*VLMIM2 !*WUnit(NCP1)*WUnit(NCP2)

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
      END SUBROUTINE INT_VLMM_VEC
!--------------------------------------------------------

      SUBROUTINE RightPart_WW_TEN(PROJ,LL,NJ,NNMAX,  &
                PPM,UPPM,ALPHAM,BETAM,GAMMAM,        &
                NPP,NALPHA,NBETA,NGAMMA,             &
                NTEN,PAR,YM)

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 
       
        INTEGER :: LL,NJ,NNMAX,NPP,NTEN,NPAR
        INTEGER :: NALPHA,NBETA,NGAMMA
        REAL(wp) :: PROJ(NPP*NALPHA*NBETA*NGAMMA)
        REAL(wp) :: YM((LL+1)*(LL+1)*NNMAX,NTEN)
        REAL(wp) :: PPM(NPP),UPPM(NPP)
        REAL(wp) :: ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
        REAL(wp) :: PAR(50)

!      DIMENSION PROJ(NPP*NALPHA*NBETA*NGAMMA)

!      DIMENSION YM((LL+1)*(LL+1)*NNMAX,3)
!      DIMENSION PPM(NPP),UPPM(NPP)
!      DIMENSION ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
!      DIMENSION NPAR(*),PAR(*)

!local 
        REAL(wp) :: BM1(NBETA*NGAMMA,NTEN),BM2(NBETA,NTEN)
        REAL(wp) :: WUnit(3)
        INTEGER :: LMAX,N,L,L22,M,IIC,IIS
        INTEGER :: I,J,IJ,I1,I2,I3,I23
        REAL(wp) :: S1(NTEN/3),S2(NTEN),S3,WW1,WW2,WW3
        REAL(wp) :: ALPHA,BETA,GAMMA
        REAL(wp) :: HALPHA,HBETA,HGAMMA,PP
        REAL(wp) :: YYRE,YYIM
!        REAL(wp) :: WWLMRE,WWLMIM

!      DIMENSION BM1(NBETA*NGAMMA,3),BM2(NBETA,3)
!      DIMENSION WUnit(3)

!        write(*,*) '444444444',PROJ

      HALPHA=PAR(38)
      HBETA =PAR(39)
      HGAMMA=PAR(40)
!      UM1   =PAR(43)



      LMAX=(LL+1)*(LL+1)

!--------------------------------
!  cos

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

            DO 3 J=1,NTEN/3
               S1(J)=0.
            DO 3 I=1,NTEN/3
               IJ=(J-1)*NTEN/3+I
               S2(IJ)=0.
3           CONTINUE   

         DO 10 I1=1,NALPHA
            ALPHA=ALPHAM(I1)

 !           WUnit(1)=COS(ALPHA)*SIN(BETA)
 !           WUnit(2)=SIN(ALPHA)*SIN(BETA)
 !           WUnit(3)=COS(BETA)

            WUnit(1)=-SIN(BETA)*COS(GAMMA)
            WUnit(2)= SIN(BETA)*SIN(GAMMA)
            WUnit(3)= COS(BETA)

!            write(*,*) 'RightPart_WW_VEC_11',L22,M-1

            WW1=1.0
            IF(I1 .EQ.1 .OR. I1 .EQ. NALPHA) WW1=0.5

!            CALL WWCS(N,L22,M-1,ALPHA,BETA,GAMMA,PP,WWLMRE,WWLMIM)


            CALL INT_Proj_GNLM_P(PROJ,PPM,UPPM,NPP,         &
                             N,L22,M-1,                          & 
                             ALPHA,BETA,GAMMA,                   &
                             NJ,I1,I2,I3,PAR,YYRE,YYIM)

            DO 4 J=1,NTEN/3
               S1(J)=0.
            DO 2 I=1,NTEN/3   
               IJ=(J-1)*NTEN/3+I
               S1(J)=S1(J)+YYRE*WUnit(I)*WW1
!               S1(J)=WW1
2           CONTINUE
               S2(IJ)=S2(IJ)+S1(J)*WUnit(J)
!               S2(IJ)=S2(IJ)+S1(J)
4           CONTINUE

10       CONTINUE
         
         DO 5 J=1,NTEN/3
         DO 5 I=1,NTEN/3   
            IJ=(J-1)*NTEN/3+I
            BM1(I23,IJ)=S2(IJ)*HALPHA
5        CONTINUE

12    CONTINUE   

!            write(*,*) '1111111111',BM1
      DO 16 I2=1,NBETA

         DO 7 J=1,NTEN/3
         DO 7 I=1,NTEN/3   
            IJ=(J-1)*NTEN/3+I
            S2(IJ)=0.
7        CONTINUE

            DO 14 I3=1,NGAMMA

               I23=(I2-1)*NGAMMA+I3
               
               WW3=1.0
               IF(I3 .EQ.1 .OR. I3 .EQ. NGAMMA) WW3=0.5
               DO 9 J=1,NTEN/3
               DO 9 I=1,NTEN/3   
                  IJ=(J-1)*NTEN/3+I
                  S2(IJ)=S2(IJ)+BM1(I23,IJ)*WW3
9              CONTINUE
14          CONTINUE
            DO 11 J=1,NTEN/3
            DO 11 I=1,NTEN/3   
               IJ=(J-1)*NTEN/3+I
               BM2(I2,IJ)=S2(IJ)*HGAMMA
11          CONTINUE
16    CONTINUE


      DO 13 J=1,NTEN/3
      DO 13 I=1,NTEN/3   
         IJ=(J-1)*NTEN/3+I
         S2(IJ)=0.
13    CONTINUE

      DO 18 I2=1,NBETA
         BETA=BETAM(I2)
         WW2=1.0
         IF(I2 .EQ.1 .OR. I2 .EQ. NBETA) WW2=0.5
         DO 17 J=1,NTEN/3
         DO 17 I=1,NTEN/3   
            IJ=(J-1)*NTEN/3+I
            S2(IJ)=S2(IJ)+BM2(I2,IJ)*SIN(BETA)*WW2
17       CONTINUE
18    CONTINUE
         DO 19 J=1,NTEN/3
         DO 19 I=1,NTEN/3   
            IJ=(J-1)*NTEN/3+I
            YM(IIC,IJ)=S2(IJ)*HBETA
19       CONTINUE

30    CONTINUE
!--------------------------------------
!    sin

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

         DO 31 J=1,NTEN/3
            S1(J)=0.
         DO 31 I=1,NTEN/3
            IJ=(J-1)*NTEN/3+I
            S2(IJ)=0.
31       CONTINUE   

      DO 20 I1=1,NALPHA
         ALPHA=ALPHAM(I1)
!         WUnit(1)=COS(ALPHA)*SIN(BETA)
!         WUnit(2)=SIN(ALPHA)*SIN(BETA)
!         WUnit(3)=COS(BETA)

         WUnit(1)=-SIN(BETA)*COS(GAMMA)
         WUnit(2)= SIN(BETA)*SIN(GAMMA)
         WUnit(3)= COS(BETA)
!         write(*,*) 'RightPart_WW_VEC_22',L22,M

         WW1=1.0
         IF(I1 .EQ.1 .OR. I1 .EQ. NALPHA) WW1=0.5

!         CALL WWCS(N,L22,M,ALPHA,BETA,GAMMA,PP,WWLMRE,WWLMIM)

         CALL INT_Proj_GNLM_P(PROJ,PPM,UPPM,NPP,         &
                     N,L22,M,                          & 
                     ALPHA,BETA,GAMMA,                   &
                     NJ,I1,I2,I3,PAR,YYRE,YYIM)
         
         DO 41 J=1,NTEN/3
               S1(J)=0.
         DO 21 I=1,NTEN/3   
            IJ=(J-1)*NTEN/3+I
            S1(J)=S1(J)+YYIM*WUnit(I)*WW1
!            S1(J)=WW1
21       CONTINUE
            S2(IJ)=S2(IJ)+S1(J)*WUnit(J)
!            S2(IJ)=S2(IJ)+S1(J)
41       CONTINUE

20    CONTINUE
         
         DO 51 J=1,NTEN/3
         DO 51 I=1,NTEN/3   
            IJ=(J-1)*NTEN/3+I
            BM1(I23,IJ)=S2(IJ)*HALPHA
51       CONTINUE

22    CONTINUE   

      DO 26 I2=1,NBETA
         BETA=BETAM(I2)

         DO 71 J=1,NTEN/3
         DO 71 I=1,NTEN/3   
            IJ=(J-1)*NTEN/3+I
            S2(IJ)=0.
71       CONTINUE

      DO 24 I3=1,NGAMMA
         I23=(I2-1)*NGAMMA+I3

         WW3=1.0
         IF(I3 .EQ.1 .OR. I3 .EQ. NGAMMA) WW3=0.5
         DO 91 J=1,NTEN/3
         DO 91 I=1,NTEN/3   
            IJ=(J-1)*NTEN/3+I
            S2(IJ)=S2(IJ)+BM1(I23,IJ)*WW3
91       CONTINUE
24    CONTINUE

         DO 111 J=1,NTEN/3
         DO 111 I=1,NTEN/3   
            IJ=(J-1)*NTEN/3+I
            BM2(I2,IJ)=S2(IJ)*HGAMMA
111      CONTINUE
26    CONTINUE


         DO 131 J=1,NTEN/3
         DO 131 I=1,NTEN/3   
            IJ=(J-1)*NTEN/3+I
            S2(IJ)=0.
131       CONTINUE

      DO 28 I2=1,NBETA
         BETA=BETAM(I2)
         WW2=1.0
         IF(I2 .EQ.1 .OR. I2 .EQ. NBETA) WW2=0.5
         DO 171 J=1,NTEN/3
         DO 171 I=1,NTEN/3   
            IJ=(J-1)*NTEN/3+I
            S2(IJ)=S2(IJ)+BM2(I2,IJ)*SIN(BETA)*WW2
171       CONTINUE

28    CONTINUE
         DO 191 J=1,NTEN/3
         DO 191 I=1,NTEN/3   
            IJ=(J-1)*NTEN/3+I
            YM(IIS,IJ)=S2(IJ)*HBETA
191       CONTINUE
50    CONTINUE
      RETURN
    END SUBROUTINE RightPart_WW_TEN
!---------------------------------------------------------------------------

      SUBROUTINE RightPart_WW_VEC(PROJ,LL,NJ,NNMAX,  &
                PPM,UPPM,ALPHAM,BETAM,GAMMAM,        &
                NPP,NALPHA,NBETA,NGAMMA,             &
                NTEN,PAR,YM)

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 
       
        INTEGER :: LL,NJ,NNMAX,NPP,NTEN,NPAR
        INTEGER :: NALPHA,NBETA,NGAMMA
        REAL(wp) :: PROJ(NPP*NALPHA*NBETA*NGAMMA)
        REAL(wp) :: YM((LL+1)*(LL+1)*NNMAX,NTEN)
        REAL(wp) :: PPM(NPP),UPPM(NPP)
        REAL(wp) :: ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
        REAL(wp) :: PAR(50)

!      DIMENSION PROJ(NPP*NALPHA*NBETA*NGAMMA)

!      DIMENSION YM((LL+1)*(LL+1)*NNMAX,3)
!      DIMENSION PPM(NPP),UPPM(NPP)
!      DIMENSION ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
!      DIMENSION NPAR(*),PAR(*)

!local 
        REAL(wp) :: BM1(NBETA*NGAMMA,NTEN),BM2(NBETA,NTEN)
        REAL(wp) :: WUnit(3)
        INTEGER :: LMAX,N,L,L22,M,IIC,IIS
        INTEGER :: I,J,IJ,I1,I2,I3,I23
        REAL(wp) :: S1(NTEN),S2(NTEN*NTEN),S3,WW1,WW2,WW3
        REAL(wp) :: ALPHA,BETA,GAMMA
        REAL(wp) :: HALPHA,HBETA,HGAMMA
        REAL(wp) :: YYRE,YYIM
!      DIMENSION BM1(NBETA*NGAMMA,3),BM2(NBETA,3)
!      DIMENSION WUnit(3)

!        write(*,*) '444444444',PROJ

      HALPHA=PAR(38)
      HBETA =PAR(39)
      HGAMMA=PAR(40)
!      UM1   =PAR(43)



      LMAX=(LL+1)*(LL+1)

!--------------------------------
!  cos

      DO 30  N=1,NNMAX
      DO 30  L=1,LL+1
         L22=L-1
      DO 30  M=1,L22+1
         IIC=(N-1)*LMAX + (L22*L22+M)

         DO 12 I3=1,NGAMMA
            GAMMA=GAMMAM(I3)

         DO 12 I2=1,NBETA
            BETA=BETAM(I2)

!      DO 12 I3=1,NGAMMA
!         GAMMA=GAMMAM(I3)
!         I23=(I2-1)*NGAMMA+I3
            I23=(I3-1)*NBETA+I2

            DO 3 J=1,NTEN/3
               S1(J)=0.
            DO 3 I=1,NTEN/3
               IJ=(J-1)*NTEN/3+I
               S2(IJ)=0.
3           CONTINUE   

         DO 10 I1=1,NALPHA
            ALPHA=ALPHAM(I1)

!            WUnit(1)=COS(ALPHA)*SIN(BETA)
!            WUnit(2)=SIN(ALPHA)*SIN(BETA)
!            WUnit(3)=COS(BETA)

         WUnit(1)=-SIN(BETA)*COS(GAMMA)
         WUnit(2)= SIN(BETA)*SIN(GAMMA)
         WUnit(3)= COS(BETA)

!            write(*,*) 'RightPart_WW_VEC_11',L22,M-1

            WW1=1.0
            IF(I1 .EQ.1 .OR. I1 .EQ. NALPHA) WW1=0.5

            CALL INT_Proj_GNLM_P(PROJ,PPM,UPPM,NPP,         &
                             N,L22,M-1,                          & 
                             ALPHA,BETA,GAMMA,                   &
                             NJ,I1,I2,I3,PAR,YYRE,YYIM)

            DO 4 J=1,NTEN/3
            DO 2 I=1,NTEN/3   
               IJ=(J-1)*NTEN/3+I
               S1(J)=S1(J)+YYRE*WW1*WUnit(I)
!               S1(J)=1.
2           CONTINUE
               S2(IJ)=S2(IJ)+S1(J)*WUnit(J)
!               S2(IJ)=S2(IJ)+S1(J)
4           CONTINUE

10       CONTINUE
         
         DO 5 J=1,NTEN/3
         DO 5 I=1,NTEN/3   
            IJ=(J-1)*NTEN/3+I
            BM1(I23,IJ)=S2(IJ)*HALPHA
5        CONTINUE

12    CONTINUE   

!            write(*,*) '1111111111',BM1
      DO 16 I3=1,NGAMMA

         DO 7 J=1,NTEN/3
         DO 7 I=1,NTEN/3   
            IJ=(J-1)*NTEN/3+I
            S2(IJ)=0.
7        CONTINUE

            DO 14 I2=1,NBETA
               BETA=BETAM(I2)
               I23=(I3-1)*NBETA+I2
               
               WW2=1.0
               IF(I2 .EQ.1 .OR. I2 .EQ. NBETA) WW2=0.5
               DO 9 J=1,NTEN/3
               DO 9 I=1,NTEN/3   
                  IJ=(J-1)*NTEN/3+I
                  S2(IJ)=S2(IJ)+BM1(I23,IJ)*SIN(BETA)*WW2
9              CONTINUE
14          CONTINUE
            DO 11 J=1,NTEN/3
            DO 11 I=1,NTEN/3   
               IJ=(J-1)*NTEN/3+I
               BM2(I3,IJ)=S2(IJ)*HBETA
11          CONTINUE
16    CONTINUE


      DO 13 J=1,NTEN/3
      DO 13 I=1,NTEN/3   
         IJ=(J-1)*NTEN/3+I
         S2(IJ)=0.
13    CONTINUE

      DO 18 I3=1,NGAMMA
         WW3=1.0
         IF(I3 .EQ.1 .OR. I3 .EQ. NGAMMA) WW3=0.5
         DO 17 J=1,NTEN/3
         DO 17 I=1,NTEN/3   
            IJ=(J-1)*NTEN/3+I
            S2(IJ)=S2(IJ)+BM2(I3,IJ)*WW3
17       CONTINUE
18    CONTINUE
         DO 19 J=1,NTEN/3
         DO 19 I=1,NTEN/3   
            IJ=(J-1)*NTEN/3+I
            YM(IIC,IJ)=S2(IJ)*HGAMMA
19       CONTINUE

30    CONTINUE
!--------------------------------------
!    sin
      DO 50  N=1,NNMAX
      DO 50  L=1,LL+1
         L22=L-1
      DO 50  M=1,L22 
         IIS=(N-1)*LMAX + (L22*L22+(M+1)+L22)
 
      DO 22 I3=1,NGAMMA
         GAMMA=GAMMAM(I3)

      DO 22 I2=1,NBETA
         BETA=BETAM(I2)

!      DO 12 I3=1,NGAMMA
!         GAMMA=GAMMAM(I3)
!         I23=(I2-1)*NGAMMA+I3
         I23=(I3-1)*NBETA+I2

         DO 31 J=1,NTEN/3
            S1(J)=0.
         DO 31 I=1,NTEN/3
            IJ=(J-1)*NTEN/3+I
            S2(IJ)=0.
31       CONTINUE   

      DO 20 I1=1,NALPHA
         ALPHA=ALPHAM(I1)

         WUnit(1)=COS(ALPHA)*SIN(BETA)
         WUnit(2)=SIN(ALPHA)*SIN(BETA)
         WUnit(3)=COS(BETA)

!         WUnit(1)=-SIN(BETA)*COS(GAMMA)
!         WUnit(2)= SIN(BETA)*SIN(GAMMA)
!         WUnit(3)= COS(BETA)
!         write(*,*) 'RightPart_WW_VEC_22',L22,M


         WW1=1.0
         IF(I1 .EQ.1 .OR. I1 .EQ. NALPHA) WW1=0.5

         CALL INT_Proj_GNLM_P(PROJ,PPM,UPPM,NPP,         &
                     N,L22,M,                          & 
                     ALPHA,BETA,GAMMA,                   &
                     NJ,I1,I2,I3,PAR,YYRE,YYIM)
         
         DO 41 J=1,NTEN/3
         DO 21 I=1,NTEN/3   
            IJ=(J-1)*NTEN/3+I
            S1(J)=S1(J)+YYIM*WW1*WUnit(I)
!            S1(J)=1.
21       CONTINUE
            S2(IJ)=S2(IJ)+S1(J)*WUnit(J)
!            S2(IJ)=S2(IJ)+S1(J)
41       CONTINUE

20    CONTINUE
         
         DO 51 J=1,NTEN/3
         DO 51 I=1,NTEN/3   
            IJ=(J-1)*NTEN/3+I
            BM1(I23,IJ)=S2(IJ)*HALPHA
51       CONTINUE

22    CONTINUE   

      DO 26 I3=1,NGAMMA

         DO 71 J=1,NTEN/3
         DO 71 I=1,NTEN/3   
            IJ=(J-1)*NTEN/3+I
            S2(IJ)=0.
71       CONTINUE

      DO 24 I2=1,NBETA
         BETA=BETAM(I2)
         I23=(I3-1)*NBETA+I2

         WW2=1.0
         IF(I2 .EQ.1 .OR. I2 .EQ. NBETA) WW2=0.5
         DO 91 J=1,NTEN/3
         DO 91 I=1,NTEN/3   
            IJ=(J-1)*NTEN/3+I
            S2(IJ)=S2(IJ)+BM1(I23,IJ)*SIN(BETA)*WW2
91       CONTINUE
24    CONTINUE

         DO 111 J=1,NTEN/3
         DO 111 I=1,NTEN/3   
            IJ=(J-1)*NTEN/3+I
            BM2(I3,IJ)=S2(IJ)*HBETA
111      CONTINUE
26    CONTINUE


         DO 131 J=1,NTEN/3
         DO 131 I=1,NTEN/3   
            IJ=(J-1)*NTEN/3+I
            S2(IJ)=0.
131       CONTINUE

      DO 28 I3=1,NGAMMA
         WW3=1.0
         IF(I3 .EQ.1 .OR. I3 .EQ. NGAMMA) WW3=0.5
         DO 171 J=1,NTEN/3
         DO 171 I=1,NTEN/3   
            IJ=(J-1)*NTEN/3+I
            S2(IJ)=S2(IJ)+BM2(I3,IJ)*WW3
171       CONTINUE

28    CONTINUE
         DO 191 J=1,NTEN/3
         DO 191 I=1,NTEN/3   
            IJ=(J-1)*NTEN/3+I
            YM(IIS,IJ)=S2(IJ)*HGAMMA
191       CONTINUE
50    CONTINUE
      RETURN
    END SUBROUTINE RightPart_WW_VEC
!--------------------------------------

      FUNCTION GNLMP(N,L,M,PP)
!  Вычисление функции U(lmn)  формула (20)
! статья Мат. Мод. 2005
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        INTEGER :: N,L,M
        REAL(wp) :: PP
!local 
        INTEGER :: KK,LM2 
        REAL(wp) :: GGC,PPP,PP2,PP32,ZZ,A,B,ACOB
!functions
        REAL(wp) :: GNLMP,GCOF,ACOBI
        REAL(sp) :: ACOBI_SP

!         write(*,*) 'GNLMP', N,L,M,PP

        KK=MOD(L+M,2)

      IF(KK .EQ. 0) THEN 
         GGC=GCOF(N,L,M)   
         PPP=PP**M
         PP2=1-PP**2
         IF(PP2 .GT. 0.) THEN 
            PP32=PP2*SQRT(PP2)
         ELSE
            PP32=0.
         ENDIF   
         ZZ=1.-2.*PP**2
         LM2=(L-M)/2.
         A=M
         B=3./2.
         ACOB=ACOBI(N+LM2,A,B,ZZ)

         GNLMP=GGC*PPP*PP32*ACOB

      ELSE
         GNLMP=0.
      ENDIF   
      RETURN
      END FUNCTION GNLMP
!-----------------------------------------------------

      SUBROUTINE VLMM_RE_IM(L,M,M1,ALPHA,BETA,GAMMA,VLMRE,VLMIM)

!          Computer the V harmonics
!          VLMRE - real part
!          VLMIM - image part
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

        INTEGER :: L,M,M1
        REAL(wp) :: ALPHA,BETA,GAMMA
        REAL(wp) :: VLMRE,VLMIM   
! local
        INTEGER ::KK
        REAL(wp) :: DM,SIGNUM,CC1,CC2,SS1,SS2
        REAL(wp) :: PLM1,PLM2
! functions
        REAL(wp) :: PLM1M2 

        DM=1.
        IF(M1 .EQ. 0) DM=0.5

!      KK=MOD((M-M1),2)
        KK=MOD(M1,2)
        SIGNUM=1.
        IF(KK .NE. 0) SIGNUM=-1.

!      KK=MOD((M+M1),2)
!      SIG1=1.
!      IF(KK .NE. 0) SIG1=-1.

!      KK=MOD((M1-M),2)
!      SIG2=1.
!      IF(KK .NE. 0) SIG2=-1.

!      CC1=COS(M*ALPHA + M1*(GAMMA+PI/2))
!      CC2=COS(M*ALPHA - M1*(GAMMA+PI/2))

!      SS1=SIN(M*ALPHA + M1*(GAMMA+PI/2))
!      SS2=SIN(M*ALPHA - M1*(GAMMA+PI/2))

        CC1=COS(M*ALPHA + M1*GAMMA)
        CC2=COS(M*ALPHA - M1*GAMMA)

        SS1= SIN(M*ALPHA + M1*GAMMA)
        SS2= SIN(M*ALPHA - M1*GAMMA)

!        write(*,*) 'VLMM_RE',L,M,M1

        PLM1 = PLM1M2(L,M, M1,BETA)
        PLM2 = PLM1M2(L,M,-M1,BETA)

!       write(*,*) 'VLMM_RE',PLM1,PLM2

        VLMRE=DM*(CC1*PLM1 + SIGNUM*CC2*PLM2)
        VLMIM=DM*(SS1*PLM1 + SIGNUM*SS2*PLM2)

!      VLMRE=DM*(CC1*PLM1*(1+SIG2) + CC2*PLM2*(1+SIG1))
!      VLMIM=DM*(SS1*PLM1*(1+SIG2) + SS2*PLM2*(1+SIG1))

      RETURN
      END SUBROUTINE VLMM_RE_IM
!----------------------------------------------

      SUBROUTINE ProjNum_TEN(DR,PROJ,PPM,  &
              ALPHAM,BETAM,GAMMAM,NALPHA,NBETA,NGAMMA, &
              XM,YM,ZM,NX,NY,NZ,NPP,NW,NTEN)

!        USE PARAM_INPUT

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 
      
        INTEGER :: NX,NY,NZ,NPP,NW,NTEN
        INTEGER :: NALPHA,NBETA,NGAMMA
        REAL(wp) :: DR(NX*NY*NZ,NTEN*NTEN)
        REAL(wp) :: ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
!        REAL(wp) ::  UM(NU),WM(NW),PPM(NPP)
        REAL(wp) ::  PPM(NPP)
        REAL(wp) :: XM(NX),YM(NY),ZM(NZ)
        REAL(wp) :: PROJ(NPP*NALPHA*NBETA*NGAMMA)

!local 
        REAL(wp) :: WUnit(3)
        REAL(wp) :: PP,ALPHA,BETA,GAMMA 
        REAL(wp) :: RRE
        INTEGER :: NJ  
        INTEGER :: I4,J4,I1,I2,I3,I321,I4321
        INTEGER :: IERR  
        CHARACTER*8 ITEX
!functions
        REAL(wp) :: FUNL3,ProjNum_FUN

      NJ = NALPHA*NBETA*NGAMMA   
     
      IERR=0
      ITEX='ProjNum_TEN'
      IF(NPP.LT.100) GOTO 2
      IERR=1
      CALL ERRPRN(ITEX,IERR)
      RETURN
2     CONTINUE

      DO 85 J4=1,NPP
         PP=PPM(J4)
      DO 85 I4=1,NJ 
         CALL IALBETGAM(I4,NALPHA,NBETA,I1,I2,I3)
         ALPHA=ALPHAM(I1)
         BETA =BETAM(I2)
         GAMMA=GAMMAM(I3)
         I321=(I3-1)*NBETA*NALPHA+(I2-1)*NALPHA+I1
         I4321=(J4-1)*NJ+I321

         WUnit(1)=-SIN(BETA)*COS(GAMMA)
         WUnit(2)= SIN(BETA)*SIN(GAMMA)
         WUnit(3)= COS(BETA)

!      DO 85 J4=1,NPP
!         PP=PPM(J4) 
!         I4321=(I4-1)*NPP+J4

         PROJ(I4321)= ProjNum_FUN(WUnit,DR,PP,ALPHA,BETA,GAMMA,   & 
                       XM,YM,ZM,NX,NY,NZ,NPP,NW,NTEN)
85    CONTINUE
      RETURN
      END SUBROUTINE ProjNum_TEN
!------------------------------------------

      FUNCTION ProjNum_FUN(WUnit,DR,PP,ALPHA,BETA,GAMMA,   & 
                       XM,YM,ZM,NX,NY,NZ,NPP,NW,NTEN)
        USE PARAM_INPUT

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 
      
        INTEGER :: NX,NY,NZ,NPP,NW,NTEN
        REAL(wp) :: DR(NX*NY*NZ,NTEN*NTEN)
        REAL(wp) :: XM(NX),YM(NY),ZM(NZ)
        REAL(wp) :: PP,ALPHA,BETA,GAMMA 
        REAL(wp) :: WUnit(3)
!local 
        REAL(wp) :: DDR(NX*NY*NZ)
        REAL(wp) :: CALPHA,SALPHA,CBETA,SBETA,SGAMMA,CGAMMA
        REAL(wp) :: AA11,AA12,AA13,AA21,AA22,AA23,AA31,AA32,AA33
        REAL(wp) :: UU,VV,WW,HW,QQQ,WW0,S,SS
        REAL(wp) :: XX,YY,ZZ
        REAL(wp) :: RRR,RRE,WW1,EI

        INTEGER :: I,I1,I2,I3,I321,J,IJ,J3
        INTEGER :: IERR  
        CHARACTER*8 ITEX
!functions
        REAL(wp) :: FUNL3,ProjNum_FUN

        RRE= PARAM(5) 
     
        IERR=0
        ITEX='ProjNum_TEN'
        IF(NPP.LT.100) GOTO 2
        IERR=1
        CALL ERRPRN(ITEX,IERR)
        RETURN
2       CONTINUE

        CALPHA=COS(ALPHA)
        SALPHA=SIN(ALPHA)
        CBETA=COS(BETA)
        SBETA=SIN(BETA)
        CGAMMA=COS(GAMMA)
        SGAMMA=SIN(GAMMA)

! direct matrix of the active rotation
        AA11= CALPHA*CBETA*CGAMMA - SALPHA*SGAMMA
        AA21= SALPHA*CBETA*CGAMMA + CALPHA*SGAMMA
        AA31=-SBETA*CGAMMA
        AA12=-CALPHA*CBETA*SGAMMA - SALPHA*CGAMMA
        AA22=-SALPHA*CBETA*SGAMMA + CALPHA*CGAMMA
        AA32= SBETA*SGAMMA
        AA13= CALPHA*SBETA
        AA23= SALPHA*SBETA
        AA33= CBETA

!  unit vector along the Z' axis
! integration is always performed along Z'  

!         WUnit(1)= CALPHA*SBETA  
!         WUnit(2)= SALPHA*SBETA 
!         WUnit(3)= CBETA        

!        WUnit(1)=-SBETA*CGAMMA
!        WUnit(2)= SBETA*SGAMMA
!        WUnit(3)= CBETA

        DO 4 I3=1,NZ
        DO 4 I2=1,NY
        DO 4 I1=1,NX
           I321=(I3-1)*NX*NY+(I2-1)*NX +I1
           DDR(I321)=0.
4       CONTINUE

        DO 30 I3=1,NZ
        DO 30 I2=1,NY
        DO 30 I1=1,NX
           I321=(I3-1)*NX*NY+(I2-1)*NX +I1

           DO 24 J=1,NTEN   
           DO 24 I=1,NTEN
              IJ=(J-1)*NTEN+I
              DDR(I321)= DDR(I321)+DR(I321,IJ)*WUnit(I)*WUnit(J)
24         CONTINUE
30      CONTINUE

        VV=0.                       !!!!!!!!!!!!!!!!!!!
        QQQ=RRE**2-PP**2
        IF(QQQ .LT. 0) QQQ=0.
        WW0=SQRT(QQQ)
        HW =2.*WW0/(NW-1)

        S=0.
        SS=0.
        DO 55 J3=1,NW
           WW =-WW0+(J3-1)*HW

!      -1
!--   T  (AL,BET,GAM)-------------------------------------------

!         XX=  AA11*UU + AA21*VV + AA31*WW
!         YY=  AA12*UU + AA22*VV + AA32*WW
!         ZZ=  AA13*UU + AA23*VV + AA33*WW

           XX=  AA11*PP + AA12*VV + AA13*WW
           YY=  AA21*PP + AA22*VV + AA23*WW
           ZZ=  AA31*PP + AA32*VV + AA33*WW

           CALL REGN(XM,YM,ZM,NX,NY,NZ,XX,YY,ZZ,EI)

           IF(EI .LT. 0.) GO TO 55

           S=FUNL3(XM,YM,ZM,DDR,NX,NY,NZ,XX,YY,ZZ,IERR)

           WW1=1.
           IF(J3.EQ.1.OR.J3.EQ.NW) WW1=0.5
           SS=SS+S*WW1
 55     CONTINUE
      
      ProjNum_FUN=SS*HW 
      END FUNCTION ProjNum_FUN
!------------------------------------------

      FUNCTION GCOF(N,L,M)
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

        INTEGER :: L,N,M
!local 
        INTEGER ::KK,N11,N22,I
        REAL(wp) :: SIGNUM,A11,A22,A33,S,GN22
        REAL(wp) :: BB,A44,A55,A66
        REAL(wp) :: RL,RN,RM
!functions
        REAL(wp) :: FACTRL,GCOF

        RL=L
        RN=N
        RM=M
        
      SIGNUM=1.
      KK=MOD(L,2)
      IF(KK .NE. 0) SIGNUM=-1.

!      N11=N+(L-M)/2
      N11=RN+(RL-RM)/2.
      A11=FACTRL(N11)
!      write(*,*) 'A11=',A11
!---------------------------
      A22=FACTRL(L-M)
      A33=FACTRL(L+M)

!      write(*,*) 'A22=',A22
!      write(*,*) 'A33=',A33
!-------------------------
!      N22=(L-M)/2+N+2
      N22=(RL-RM)/2.0+RN+2.0
      S=1.
      DO 4 I=1, 2*N22-1,2
         S=S*I
4     CONTINUE   
!      write(*,*) 'S=',S
      GN22=SQRT(PI)*S/(2**N22)
!      write(*,*) 'GN22=',GN22

!--------------------------------
      BB=(N+1)*A11*SQRT((2*L+1)*A22*A33)
!            write(*,*) 'BB=',BB
!      A44=FACTRL((L-M)/2)
      A44=FACTRL(INT((RL-RM)/2.))
!      A55=FACTRL((L+M)/2)
      A55=FACTRL(INT((RL+RM)/2.))
      A66=2**(L+1)
!---------------------------------
!            write(*,*) 'A44=',A44
!            write(*,*) 'A55=',A55
!            write(*,*) 'A66=',A66

      GCOF=SIGNUM*BB/A44/A55/A66/GN22

!      write(*,*) 'GCOF',GCOF

      RETURN
      END FUNCTION GCOF
!-----------------------------------------------
      FUNCTION HNLR(L,N,RR)
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 

        INTEGER :: L,N
        REAL(wp) :: RR
!local 
         REAL(wp) :: Z,A,B,ACOB 
!functions
        REAL(wp) :: HNLR,ACOBI

!      IF(RR .NE. 0 .AND. RR .NE. 1) GOTO 12
!         WRITE(IPRNT,*) 'Procedure HNLR ERROR: RR=0 or 1 '
!         RETURN
! 12   CONTINUE   

!  AL=2.5
      Z=1.-2.*RR**2
      A=L+1./2.
      B=1.
      ACOB=ACOBI(N,A,B,Z)
      HNLR=RR**L*(1.0-RR**2)*ACOB

!  AL=3.5
!      Z=1-2*RR**2
!      A=L+1./2.
!      B=2.
!      ACOB=ACOBI(N,A,B,Z)
!      HNLR=RR**L*(1.0-RR**2)**2*ACOB

!  AL=4.5
!      Z=1-2*RR**2
!      A=L+1./2.
!      B=3.
!      ACOB=ACOBI(N,A,B,Z)
!      HNLR=RR**L*(1.0-RR**2)**3*ACOB
!------------------------------------
!      A=0.
!      B=L+1./2.
!      Z=2*RR**2-1.
!      ACOB=ACOBI(N,A,B,Z)
!      CC=SQRT(4.0*N+2.0*L+3.0)
!      HNLR=CC*RR**L*ACOB
!-------------------------------------
!   Laguerra
!      Z=RR**2
!      LAG=1+L+0.5-Z
!     HNLR=RR**L*EXP(-Z)*LAG
      RETURN
      END FUNCTION HNLR
!-------------------------------------------------

      SUBROUTINE Summa_Harm_Vec(LL,NNMAX,RRM,TETAM,PHIM, &
                                 NRR,NTETA,NPHI,NTEN,XXM,GM)
!  summation of scalar spherical harmonics 
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 

        INTEGER :: LL,NNMAX,NRR,NTETA,NPHI,NTEN
        REAL(wp) :: RRM(NRR),TETAM(NTETA),PHIM(NPHI)
        REAL(wp) :: GM(NRR*NTETA*NPHI,NTEN*NTEN)
        REAL(wp) :: XXM((LL+1)*(LL+1)*NNMAX,NTEN*NTEN)
!local
        INTEGER :: LMAX,K1,K2,KK,I3,I2,I1,I21,I321
        INTEGER :: N,L,L22,M,IIC,IIS
        REAL(wp) :: RR,TETA,PHI,S1
        REAL(wp) :: YLMRE,YLMIM
!functions
        REAL(wp) :: HNLR

        LMAX=(LL+1)*(LL+1)

        DO 18 K2=1,NTEN  
        DO 18 K1=1,NTEN
           KK=(K2-1)*NTEN+K1
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
              CALL YLM_RE(L22,M-1,TETA,PHI,YLMRE,YLMIM)
              S1=S1+XXM(IIC,KK)*YLMRE* HNLR(L22,N,RR)
4          CONTINUE
        
           DO 6 N=1,NNMAX
           DO 6 L=1,LL+1
              L22=L-1
           DO 6 M=1,L22
              IIS=(N-1)*LMAX + (L22*L22+(M+1)+L22)
              CALL YLM_RE(L22,M,TETA,PHI,YLMRE,YLMIM)
              S1=S1+XXM(IIS,KK)*YLMIM* HNLR(L22,N,RR)  
6          CONTINUE
              GM(I321,KK)=S1
10         CONTINUE
18      CONTINUE
      RETURN
      END  SUBROUTINE Summa_Harm_Vec
!---------------------------------------------------
      SUBROUTINE YLM_RE (L,M,TETA,PHI,YLMRE,YLMIM)
!          Computer the spherical harmonics Y(L,M)
!          YLMRE - real part
!          YLMIM - image part
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

        INTEGER :: L,M
        REAL(wp) :: TETA,PHI
        REAL(wp) :: YLMRE,YLMIM
!local
        INTEGER :: MM,KK
        REAL(wp) :: SIGNUM,X,PNN,PLAG
!functions
        REAL(wp) :: PNORM,PLGNDR


      MM=ABS(M)
      KK=MOD(MM,2)
      SIGNUM=1.
      IF(M .LT. 0 .AND. KK .EQ. 1) SIGNUM=-1.
      IF(M .LT. 0 .AND. KK .EQ. 0) SIGNUM=1.
      X=COS(TETA) 

      PNN=PNORM(L,MM)/SQRT(2.*PI)
      PLAG=PLGNDR(L,MM,X)
      YLMRE= SIGNUM*PNN*PLAG*COS(MM*PHI)
      YLMIM= SIGNUM*PNN*PLAG*SIN(MM*PHI)
!      CYLM_1=CMPLX(YLMRE,YLMIM)
      RETURN
      END SUBROUTINE YLM_RE
!---------------------------------------------------
      SUBROUTINE XXMM(LL,NNMAX,XM,NTEN)
!  coefficients of decomposition are written as 
!  vector XM 
! В этой процедуре устанавливается тип тензора 
! (симметричный или антисимметричный)  

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        INTEGER :: LL,NNMAX,NTEN
        REAL(wp) :: XM((LL+1)*(LL+1)*NNMAX,NTEN)

!local 
        INTEGER ::LMAX,N,L,L22,M,IIC,IIS
        INTEGER :: I,J,IJ
        REAL(wp) :: AAC(NTEN),AAS(NTEN)
! functions
!        REAL(wp) :: AACOS1,AACOS2,AACOS3,AACOS4,AACOS5
!        REAL(wp) :: AASIN1,AASIN2,AASIN3,AASIN4,AASIN5

!        IF(NTEN .GT. 3) THEN  
!           WRITE(*,*) 'Procedure XXMM: NTEN > 3'
!           STOP
!        ENDIF

        LMAX=(LL+1)*(LL+1)
!        DO 30 I=1,NTEN
        DO 10  N=1,NNMAX
        DO 10  L=1,LL+1
           L22=L-1
        DO 10  M=1,L22+1
           IIC=(N-1)*LMAX + (L22*L22+M)
!         write(*,*) 'XXMM-cos',L,M,IIC
           CALL ACOS(L22,M-1,N,AAC)

!           DO 2 J=1,NTEN
           DO 2 I=1,NTEN
!              IJ=(J-1)*NTEN+I
              XM(IIC,I)=AAC(I)
2          CONTINUE
!           write(*,*) 'XXMM:IIC'
!           write(*,*) IIC,L22,M
10      CONTINUE

        DO 20  N=1,NNMAX
        DO 20  L=1,LL+1
           L22=L-1
        DO 20  M=1,L22
           IIS=(N-1)*LMAX + (L22*L22+(M+1)+L22)
!         write(*,*) 'XXMM-sin',L,M,IIS
           CALL ASIN(L22,M,N,AAS)

!           DO 4 J=1,NTEN
           DO 4 I=1,NTEN
!              IJ=(J-1)*NTEN+I
              XM(IIS,I)=AAS(I)
4          CONTINUE
!           write(*,*) 'XXMM:IIS'
!           write(*,*) IIS,L22,M
20      CONTINUE
!30    CONTINUE
      RETURN
      END SUBROUTINE XXMM
!-------------------------------------------

      SUBROUTINE  ACOS(L,M,N,AACOS)
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        INTEGER :: L,M,N
! local
        REAL(wp) :: AACOS(9)

        AACOS(1)=1.
        IF(L .EQ.0 .AND. M .EQ.0) AACOS(1)=0.5

        AACOS(2)=1.   
        IF(L .EQ.0 .AND. M .EQ.0) AACOS(2)=0.     

        AACOS(3)=0.5   
        IF(L .EQ.0 .AND. M .EQ.0) AACOS(3)=0.2 

        AACOS(4)=0.0   
 !       IF(L .EQ.0 .AND. M .EQ.0) AACOS(4)=0.5 

        AACOS(5)=0.0   
 !       IF(L .EQ.0 .AND. M .EQ.0) AACOS(5)=0.5 

        AACOS(6)=0.0   
!        IF(L .EQ.0 .AND. M .EQ.0) AACOS(6)=1. 

        AACOS(7)=0.   
        IF(L .EQ.0 .AND. M .EQ.0) AACOS(7)=0. 

        AACOS(8)=0.   
        IF(L .EQ.0 .AND. M .EQ.0) AACOS(8)=0. 

        AACOS(9)=0.   
        IF(L .EQ.0 .AND. M .EQ.0) AACOS(9)=0. 

      RETURN
      END SUBROUTINE ACOS
!----------------------------------------------------------
      SUBROUTINE  ASIN(L,M,N,AASIN)
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        INTEGER :: L,M,N
! local
        REAL(wp) :: AASIN(9)

        AASIN(1)=-1.

        AASIN(2)= 0.   

        AASIN(3)=-0.5   

        AASIN(4)= 0.
 
        AASIN(5)=0.   

        AASIN(6)=0.   

        AASIN(7)=0.   

        AASIN(8)=0.   

        AASIN(9)=0.   

      RETURN  
      END SUBROUTINE ASIN
!----------------------------------------------------------

      SUBROUTINE Vec_Line(XM,YM,LN,NTEN)
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        INTEGER:: LN,NTEN
        REAL(wp):: XM(LN,NTEN), YM(NTEN*LN)
! local 
        INTEGER:: K,I,II

        DO 4 K=1,NTEN**2
        DO 4 I=1,LN
           II=(K-1)*LN+I
           YM(II)=XM(I,K)      
!           IF(II .GT. NTEN2*LN*0.6) THEN   !можно будет убрать
!              YM(II)=0.0                 !можно будет убрать  
!           ENDIF                          !можно будет убрать  
4       CONTINUE
      RETURN
      END SUBROUTINE Vec_Line
!---------------------------------------------
      SUBROUTINE Vec_Line_Inv(YM,XM,LN,NTEN)
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 

        INTEGER:: LN,NTEN
        REAL(wp):: XM(LN,NTEN), YM(NTEN*LN)
! local 
        INTEGER:: K,I,II

        DO 4 K=1,NTEN**2
        DO 4 I=1,LN
           II=(K-1)*LN+I
           XM(I,K)=YM(II)
4       CONTINUE
      RETURN
      END SUBROUTINE Vec_Line_Inv
!----------------------------------------------

      SUBROUTINE GNUFORM_0(NX,GC,KAN)
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
      
        INTEGER:: NX,KAN
        REAL(wp) :: GC(NX)
!local  
        INTEGER:: I

        DO 10 I=1,NX
           WRITE(KAN,100) I,GC(I)
!            WRITE(KAN, *)
!            WRITE(KAN,'(A1)') char(13)//char(10)
10     CONTINUE 
100   FORMAT(I3,E12.3)
!110   FORMAT(A4)
      RETURN
      END SUBROUTINE GNUFORM_0
!------------------------------------------
      SUBROUTINE GNUFORM1(XM,NX,GC,KAN)
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
      
        INTEGER:: NX,KAN
        REAL(wp) :: XM(NX),GC(NX)
!local  
        INTEGER:: I

        DO 10 I=1,NX
           WRITE(KAN,100) XM(I),GC(I)
!            WRITE(KAN, *)
!            WRITE(KAN,'(A1)') char(13)//char(10)
10     CONTINUE 
100   FORMAT(2E12.3)
!110   FORMAT(A4)
      RETURN
      END SUBROUTINE GNUFORM1
!------------------------------------------

      SUBROUTINE IALPHBETGAM(JJ,N1,N2,I1,I2,I3)
! calculation I1, I2, I3--  indexes by the number of  JJ,
! where JJ is calculated as  JJ=(I3-1)*N1*N2 + (I2-1)*N1 + I1
      
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=double 

        INTEGER :: JJ,N1,N2,I1,I2,I3

!  local
        REAL(wp) :: RJJ,RN12,RMODUL,RN1
        INTEGER :: IC,I12
      RJJ=JJ
      RN12=FLOAT(N1*N2)
!      RMODUL=AMOD(RJJ,RN12)
      RMODUL=MOD(RJJ,RN12)
      IC=1
      IF(RMODUL.EQ.0.) IC=0
      I3=INT(RJJ/RN12)+IC
      I12=JJ-(I3-1)*N1*N2
     
      RJJ=I12
      RN1=FLOAT(N1)
!      RMODUL=AMOD(RJJ,RN1)
      RMODUL=MOD(RJJ,RN1)
      IC=1
      IF(RMODUL.EQ.0.) IC=0
      I2=INT(RJJ/N1)+IC
      I1=I12-(I2-1)*N1
      RETURN
      END SUBROUTINE IALPHBETGAM
!----------------------------------------------------------
          SUBROUTINE POLAR2D(X,Y,FI)
!    procedure perfom a passage to polar coordinates  FI
!    by use decart X,Y coordinate
!    SIGNY- signs of  Y
!    FI- polar angle
          IMPLICIT NONE
          INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
          INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=double 
          REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

          REAL(wp) :: X,Y,XY,FI
          REAL(wp) :: SIGNY


          SIGNY=SIGN(1.0_wp,Y)
!          SIGNZ=SIGN(1.0_wp,Z)

          IF(X.NE.0.0_wp)GO TO 30
          IF(Y.NE.0.0_wp) THEN
!             FI=0.5*(1.-SIGNY)*PI+PI*0.5
             FI=PI/2.0*(2.0-SIGNY)
          ELSE
             FI=0.0_wp
          ENDIF   
          GO TO 45
30        FI=ATAN2(Y,X) 
          IF(FI .LT. 0.0) FI=FI+ PI* 2.0    ! for the case when  FI \in (0,2*PI)
 
45        CONTINUE
        RETURN
        END SUBROUTINE POLAR2D
!--------------------------------------------

      FUNCTION PLM1M2(L,M1,M2,BET)
!   Varshalovich (english) p.77 eq.(5)

      IMPLICIT NONE
      INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
      INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single 
      REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

      INTEGER :: L,M1,M2
      REAL(wp) :: BET
! functions
      REAL(wp) :: PLM1M2,SIGNUM,FACTRL

! local
      REAL(wp) :: Y1,Y2,Y3,Y4,CBET
      REAL(wp) :: SS
      INTEGER :: K,KK1,KK2,KKK1,KKK2 
!      INTEGER :: M22

!      write(*,*) 'PLM1M2',L,M1,M2


      CBET=COS(BET)
      KK1=MAX(0, M2-M1)
      KK2=MIN(L-M1, L+M2)

!      write(*,*) 'PLM1M2:KK1,KK2',KK1,KK2,L,M1,M2

      Y1=SIGNUM(M1-M2)*SQRT(FACTRL(L+M1)*FACTRL(L-M1)) *  &
         SQRT(FACTRL(L+M2)*FACTRL(L-M2))
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
      END FUNCTION PLM1M2
!-------------------------------------------------------------
      FUNCTION FACTRL(n)

      IMPLICIT NONE
      INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
      INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single 
!      REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 
      INTEGER :: n
! function
      REAL(wp) :: FACTRL
!  local
!      REAL(wp) :: factrl
! return the value n! as a floating-point number
      INTEGER :: j,ntop
      REAL(wp) :: a(33),gammln
      SAVE ntop,a
      DATA ntop, a(1) /0,1./   !table initialized with 0! only
      if(n .lt. 0) then
          pause 'negativefactorial in "factrl" '
      else if (n .le. ntop) then
              FACTRL=a(n+1)
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
      END FUNCTION FACTRL
!------------------------------------------------
      FUNCTION gammln(xx)
        IMPLICIT NONE
      INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
      INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single 

        REAL(wp) :: gammln,xx
!   return the value lnG(xx)

! local
        INTEGER :: j
!        REAL(wp) :: ser, stp,tmp,x,y,cof(6)
        DOUBLE PRECISION ser, stp,tmp,x,y,cof(6)
        SAVE cof,stp
        DATA cof,stp /76.18009172947146d0, -86.50532032941677d0, &
             24.01409824083091d0, -1.231739572450155d0, .1208650973866179d-2, &
             -.5395239384953d-5, 2.5066282746310005d0/
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
      END FUNCTION gammln
!------------------------------------------------------
      FUNCTION SIGNUM (LA)
        IMPLICIT NONE
      INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
      INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single 
        INTEGER :: LA
! function
        REAL(wp) :: SIGNUM

! local 
        INTEGER :: NN
        NN=ABS(MOD (LA,2))
        SIGNUM=1.
        IF(NN .EQ. 1) SIGNUM=-1.
      RETURN
      END FUNCTION SIGNUM
!----------------------------------------



      SUBROUTINE SIZE_MODELS(N1,N2,N3)
!     Read input file to determine data needs
        USE CHANALS
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single

        INTEGER :: N1,N2,N3
! local
        REAL(sp) :: DUMMY
!        OPEN(12,FILE ='Data_Gauss.dat', ERR=400)
!        OPEN(14,FILE ='Data_Parall.dat', ERR=410)

        READ(KAN(12),*)
        N1=0
10      READ(KAN(12),110,END=20) DUMMY
        N1=N1+1
        GO TO 10
20      CONTINUE
!    Having gone through Unit 12, I need to reposition it at the beginning
        REWIND(KAN(12))

!--------------------------------
        READ(KAN(14),*)
        N2=0
12      READ(KAN(14),110,END=30) DUMMY
        N2=N2+1
        GO TO 12
30      CONTINUE
!    Having gone through Unit 14, I need to reposition it at the beginning
        REWIND(KAN(14))

!--------------------------------
        READ(KAN(15),*)
        N3=0
14      READ(KAN(15),110,END=40) DUMMY
        N3=N3+1
        GO TO 14
40      CONTINUE
!    Having gone through Unit 15, I need to reposition it at the beginning
        REWIND(KAN(15))

110     FORMAT(1X,9F10.4)
        RETURN
!400     PRINT *, 'Can not open Data_Models.dat'
!410     PRINT *, 'Can not open Data_Parall.dat'
!      STOP
      END SUBROUTINE SIZE_MODELS
!---------------------------------------------------

      SUBROUTINE SIZE_PARAM(N)

      USE CHANALS
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single

!     Read input file to determine the data needs
        INTEGER :: N
        REAL(sp) :: DUMMY
!        OPEN(KAN(11),FILE ='Data_Input.dat', ERR=400)
        READ(KAN(11),*)
        N=0
10      READ(KAN(11),100,END=20) DUMMY
        N=N+1
        GO TO 10
20      CONTINUE
!    Having gone through Unit 11, I need to reposition it at the beginning
        REWIND(KAN(11))
100     FORMAT(1x,F10.4)   
        RETURN
400     PRINT *, 'Can not open Data_Input.dat'
      STOP
      END SUBROUTINE SIZE_PARAM
!---------------------------------------------------
      SUBROUTINE INPUT_PARAM(A,N)
!      Read array
      USE CHANALS
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single

        INTEGER :: N
        REAL(sp) :: A(N) 

        READ(KAN(11),*)
        DO 30 I=1,N
           READ(KAN(11),100) A(I)
30      CONTINUE
100   FORMAT((1X,F10.4))   
      RETURN
      END SUBROUTINE INPUT_PARAM
!-----------------------------------------------------

      SUBROUTINE XYANG(am,gm,fm1,n1,n2,niter)
!... am(n1,n2) - matrix 
!....gm(n1) - right part 
! ...fm1(n2) - initial estimation (input) and solution (output)
!-   niter -the number of iterations
      IMPLICIT NONE
      INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
      INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single 

      REAL(wp) :: am(n1*n2),gm(n1),fm1(n2)
      INTEGER :: n1,n2,niter

      REAL(wp) :: CC,q1,q2,a1m(n2),fm2(n2)
      INTEGER :: i,i1,i2,j

      CC=1.0
!      itex = 'xyang'
!      WRITE(*,*) 'number iterations:',niter
!      ierr = 0
      DO 40 i = 1,niter
!         WRITE(*,*) 'iteration I=:',i
      DO 30 i1 = 1,n1
      DO 10 i2 = 1,n2
      j = (i2-1)*n1+i1
      a1m(i2) = am(j)
!      a1m(i2) = am(i1,i2)
10    CONTINUE
      CALL scal(fm1,a1m,n2,q1)
      CALL scal(a1m,a1m,n2,q2)

!      write(*,*) 'XYANG:q2=',q2

      IF(q2 .NE. 0.0) GO TO 12
      WRITE(*,*) 'xyang: error ','i1=',i1,'i2=',i2,'j=',j
      RETURN
12    CONTINUE
      DO 14 i2 = 1,n2
      fm2 (i2) = fm1 (i2) - CC*(q1-gm(i1))/q2*a1m(i2)
!      IF(fm2(i2) .LT. 0.) fm2(i2) = 0.
14    fm1(i2) = fm2(i2)
30    CONTINUE
40    CONTINUE
      RETURN
      END SUBROUTINE XYANG  
!-----------------------------------------
      SUBROUTINE scal(a,b,n,q)
!  scale production 
      IMPLICIT NONE
      INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
      INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single 
      REAL(wp) :: a(n),b(n),q
      INTEGER :: n

      INTEGER :: i
      REAL(wp) :: s

      s = 0.0_wp
      DO 2 i = 1 , n
2     s = s + a(i) * b(i)
      q = s
      RETURN
      END SUBROUTINE scal
!-------------------------------------------------------
 
!---------------------------------------------------------------
      SUBROUTINE IALBET(JJ,N1,I1,I2)
! calculation I1, I2, --  indexes by the number of  JJ,
! where JJ is calculated as  JJ=(I2-1)*N1 + I1
        IMPLICIT NONE      
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=double 

        INTEGER :: JJ,N1,I1,I2,I32
!  local
        REAL(wp) :: RJJ, RN12,RMODUL
        INTEGER :: IC

      RJJ=JJ
      RN12=FLOAT(N1)
      RMODUL=MOD(RJJ,RN12)
      IC=1
      IF(RMODUL .EQ. 0.) IC=0
      I2=INT(RJJ/RN12)+IC
      I1=JJ-(I2-1)*N1
      RETURN
      END SUBROUTINE IALBET
!----------------------------------------------------------

      SUBROUTINE IALBETGAM(JJ,N1,N2,I1,I2,I3)
! calculation I1, I2, I3--  indexes by the number of  JJ,
! where JJ is calculated as  JJ=(I3-1)*N1*N2 + (I2-1)*N1 + I1
        IMPLICIT NONE      
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=double 

        INTEGER :: JJ,N1,N2,I1,I2,I3,I32
!  local
        REAL(wp) :: RJJ, RN12,RMODUL,RN1
        INTEGER :: IC,I12

      RJJ=JJ
      RN12=FLOAT(N1*N2)
      RMODUL=MOD(RJJ,RN12)
      IC=1
      IF(RMODUL .EQ. 0.) IC=0
      I3=INT(RJJ/RN12)+IC
      I12=JJ-(I3-1)*N1*N2
     
      RJJ=I12
      RN1=FLOAT(N1)
      RMODUL=MOD(RJJ,RN1)
      IC=1
      IF(RMODUL .EQ. 0.) IC=0
      I2=INT(RJJ/N1)+IC
      I1=I12-(I2-1)*N1
      RETURN
      END SUBROUTINE IALBETGAM
!----------------------------------------------------------

          SUBROUTINE POLAR2(X,Y,Z,TET,FI)
!    procedure perfom a passage to spherical coordinates TET, FI
!    by use decart X,Y,Z coordinate
!    SIGNY,SIGNZ- signs of  X,Y,Z
!    TET,FI- spherical coordinates
          IMPLICIT NONE
          INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
          INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
          REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

          REAL(wp) :: X,Y,Z,XY,TET,FI
          REAL(wp) :: SIGNY,SIGNZ
!          INTEGER :: SIGNY,SIGNZ

          SIGNY=SIGN(1.0,Y)
          SIGNZ=SIGN(1.0,Z)

          IF(X.NE.0.0)GO TO 30
          FI=PI*0.5*SIGNY
          GO TO 35
30        FI=ATAN2(Y,X)
35        CONTINUE
          IF(Z.NE.0)GO TO 40
          TET=PI*0.5*SIGNZ
          GO TO 45
40        TET=ATAN2(SQRT(X*X+Y*Y),Z)
45        CONTINUE

        RETURN
        END SUBROUTINE POLAR2
!--------------------------------------------


          SUBROUTINE POLAR(X,Y,Z,TET,FI)
!    procedure perfom a passage to spherical coordinates TET, FI
!    by use decart X,Y,Z coordinate
!    SIGNY,SIGNZ- signs of  X,Y,Z
!    TET,FI- spherical coordinates
          IMPLICIT NONE
          INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
          INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
          REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

          REAL(wp) :: X,Y,Z,XY,TET,FI
          REAL(wp) :: SIGNY,SIGNZ
!          INTEGER :: SIGNY,SIGNZ

          SIGNY=SIGN(1.0,Y)
          SIGNZ=SIGN(1.0,Z)

          IF(X.NE.0.0)GO TO 30
          IF(Y.NE.0.0) THEN
!             FI=0.5*(1.-SIGNY)*PI+PI*0.5
             FI=PI/2.0*(2.0-SIGNY)
          ELSE
             FI=0.0_wp
          ENDIF   
          GO TO 35
30        FI=ATAN2(Y,X) 
          IF(FI .LT. 0.0) FI=FI+ PI* 2.0    ! for the case when  FI \in (0,2*PI)
35        CONTINUE
          IF(Z.NE.0.0)GO TO 40
          TET=PI*0.5_wp  !*SIGNZ
          GO TO 45
40        XY=X*X+Y*Y
          IF(XY .NE. 0.0) THEN 
             TET=ATAN2(SQRT(XY),Z)  
!             TET=ATAN(SQRT(X*X+Y*Y)/Z) +PI/2.0 
          ELSE
             TET=ABS(SIGNZ-1.0_wp)*PI/2.0_wp
          ENDIF   
45        CONTINUE
        RETURN
        END SUBROUTINE POLAR
!--------------------------------------------
          SUBROUTINE POLAR_1(X,Y,Z,TET,FI)
!    procedure perfom a passage to spherical coordinates TET, FI
!    by use decart X,Y,Z coordinate
!    SIGNY,SIGNZ- signs of  X,Y,Z
!    TET,FI- spherical coordinates
          IMPLICIT NONE
          INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
          INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
          REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

          REAL(wp) :: X,Y,Z,XY,TET,FI
          INTEGER :: SIGNX,SIGNY,SIGNZ

          SIGNY=SIGN(1.0,X)
          SIGNY=SIGN(1.0,Y)
          SIGNZ=SIGN(1.0,Z)
          XY=X*X+Y*Y

          IF(X .NE. 0.0) THEN
             IF(Y .EQ. 0.0) THEN 
                FI=(1-SIGNX)*PI/2.0
             ELSE
                FI=ATAN2(Y,X) 
             ENDIF   
          ENDIF
          IF(X .EQ. 0.0) THEN
             IF(Y .NE. 0.0) THEN 
                FI=(2-SIGNY)*PI/2.0
             ELSE
                FI=0.0
             ENDIF   
          ENDIF

          IF(Z.NE.0.0) THEN
             IF(XY .NE. 0.0) THEN 
                TET=ATAN2(SQRT(XY),Z)
             ELSE
                TET=ABS(SIGNZ-1)*PI/2.
             ENDIF
          ENDIF
          IF(Z.EQ.0.0) THEN
              TET=PI/2.0
          ENDIF    
        RETURN
        END SUBROUTINE POLAR_1
!--------------------------------------------

        SUBROUTINE SPRINT_1(DR,NX,NY,NZ,KNZ,KK,NTEN,KAN)
! print of some section  
! KK - component number of the vector field   
          USE MAIN_ARRAYS
          IMPLICIT NONE
!          INTEGER, PARAMETER  :: dp = kind(1.0d0) !working precision=double 
          INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
          INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=double    

          INTEGER :: NX,NY,NZ,KNZ,KK,NTEN,KAN 
          REAL(wp) :: DR(NX*NY*NZ,NTEN*NTEN)
! local
          REAL(wp) :: DR1(NX*NY)
          INTEGER :: I1,I2,I3,I21,I321

          DO 8 I3=KNZ,KNZ
          DO 8 I2=1,NY
          DO 8 I1=1,NX
             I21=(I2-1)*NX+I1
             I321=(I3-1)*NX*NY+(I2-1)*NX +I1
             DR1(I21)=DR(I321,KK)
 8        CONTINUE   

!          DO 10 I2=1,NY
!             WRITE(KAN,100) (DR1((I2-1)*NX+I1),I1=1,NX)
!             WRITE(KAN, *)
!10       CONTINUE 
          CALL GNUFORM_M(XM,YM,NX,NY,DR1,KAN)
      RETURN
100   FORMAT(E10.3)
      END SUBROUTINE SPRINT_1
!------------------------------------------
      SUBROUTINE GNUFORM(NX,NY,GC,KAN)
!      SUBROUTINE GNUFORM (NX,NY,KK,GC,KAN)
!
          IMPLICIT NONE
          INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
          INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single  
          
          REAL(wp) :: GC(NX*NY)
          INTEGER :: NX,NY,KAN   !,KK
! local 
          INTEGER :: J,I

!          DO 10 J=KK,NY
          DO 10 J=1,NY
             WRITE(KAN,100) (GC((J-1)*NX+I),I=1,NX)
             WRITE(KAN, *)
!               WRITE(KAN,'(A1)') char(13)//char(10)
10       CONTINUE 
!             write(*,*) 'GNUFORM',GC
100    FORMAT(E10.3)
!110   FORMAT(A4)
      RETURN
      END SUBROUTINE GNUFORM
!----------------------------------------------------------------------
      SUBROUTINE GNUFORM_M (XM,YM,NX,NY,GC,KAN)
          IMPLICIT NONE
          INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
          INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single  
          
          REAL(wp) :: XM(NX),YM(NY),GC(NX*NY)
          INTEGER :: NX,NY,KAN 
! local 
          INTEGER :: J,I
          DO 12 J=1,NY
          DO 10 I=1,NX
             WRITE(KAN,100) XM(I),YM(J),GC((J-1)*NX+I)
10       CONTINUE 
             WRITE(KAN, *)
12       CONTINUE 
100   FORMAT(3E15.2)
      RETURN
      END SUBROUTINE GNUFORM_M
!----------------------------------------------------------------------


      SUBROUTINE DATA(RR,PPM,RRM,XM,YM,ZM,NPP,NRR,NX,NY,NZ,PHIM,TETM,  & 
                 NPHI,NTET,ALPHAM,BETAM,GAMMAM,NALPHA,NBETA,NGAMMA,    &
                 HX,HY,HZ,HPHI,HTET,UM,VM,WM,NU,NV,NW)
!CL.... -
!CL.... -
        USE PARAMTR
!        USE MAIN_ARRAYS
        IMPLICIT NONE
        
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp  

        INTEGER :: NPP,NRR,NX,NY,NZ,NPHI,NTET,NALPHA,NBETA,NGAMMA
        INTEGER :: NU,NV,NW
        INTEGER :: I
       
        REAL(wp) :: XM(NX),YM(NY),ZM(NZ)
        REAL(wp) :: UM(NU),VM(NV),WM(NW)
        REAL(wp) :: PHIM(NPHI),TETM(NTET)
        REAL(wp) :: PPM(NPP),RRM(NRR)
        REAL(wp) :: ALPHAM(NALPHA)
        REAL(wp) :: BETAM(NBETA)
        REAL(wp) :: GAMMAM(NGAMMA)

        REAL(wp) :: HX,HY,HZ,HPHI,HTET,HRE
        REAL(wp) :: HU,HV   !!!!!!!!!!!!!

        REAL(wp) ::  RR
        REAL(wp) ::  U0,U1,V0,V1,W0,W1,X0,X1,Y0,Y1,Z0,Z1
        REAL(wp) ::  TET0,TET1,PHI0,PHI1,RE1,RE2
        REAL(wp) ::  PP0,PP1,HPP
        REAL(wp) ::  ALPHA0,ALPHA1,BETA0,BETA1,HBETA,HALPHA
        REAL(wp) ::  GAMMA0,GAMMA1,HGAMMA

        U0  = -RR
        U1  =  RR
        V0  = -RR
        V1  =  RR
        W0  = -RR
        W1  =  RR
        X0  = -RR 
        X1  =  RR
        Y0  = -RR
        Y1  =  RR
        Z0  = -RR
        Z1  =  RR
        TET0=  0.    !PI/2.   
        TET1=  PI
        PHI0=  0.
        PHI1= 2.*PI  
!        RE1 =  RR
!        RE2 =  2.*RR 
!        NU=NX
!        NV=NY
!        NW=NZ

        PP0   = 0.0   !-RR 
        PP1   = RR  !*10.
!        HPP=(PP1-PP0)/FLOAT(NX)
        HPP=(PP1-PP0)/FLOAT(NPP-1)
        DO 10 I=1,NPP
           PPM(I)=PP0+(I-1)*HPP
           RRM(I)=PPM(I)
10      CONTINUE

!          HX=(X1-X0)/FLOAT(NX-2)
           HX=(X1-X0)/FLOAT(NX-1)
           HU=HX
        DO 20 I=1,NX
           XM(I)=X0+(I-1)*HX
           UM(I)=XM(I)
20      CONTINUE

!------------------
!           U0=-10*RR
!           U1= 10*RR
!           HU= (U1-U0)/(NX-2)
!        DO 21 I=1,NX
!           UM(I)=U0+(I-1)*HU
!21      CONTINUE
!----------------------
!           HY=(Y1-Y0)/FLOAT(NY-2)
           HY=(Y1-Y0)/FLOAT(NY-1)
           HV=HY
        DO 24 I=1,NY
           YM(I)=Y0+(I-1)*HY
           VM(I)=YM(I)
 24     CONTINUE

!           HZ=(Z1-Z0)/FLOAT(NZ-2)
           HZ=(Z1-Z0)/FLOAT(NZ-1)
        DO 26 I=1,NZ
           ZM(I)=Z0+(I-1)*HZ
           WM(I)=ZM(I)
26      CONTINUE

           HPHI=(PHI1-PHI0)/FLOAT(NPHI-1)
        DO 28 I=1,NPHI
28         PHIM(I)=PHI0+(I-1)*HPHI

           HTET=(TET1-TET0)/FLOAT(NTET-1)
        DO 30 I=1,NTET
30         TETM(I)=TET0+(I-1)*HTET


           ALPHA0= 0.0  !-PI
           ALPHA1= 2.0*PI
!           HALPHA=(ALPHA1-ALPHA0)/FLOAT(NALPHA-1)
           HALPHA=(ALPHA1-ALPHA0)/FLOAT(NALPHA)
        DO 32 I=1,NALPHA
32         ALPHAM(I)=ALPHA0+(I-1)*HALPHA

           BETA0=0.0
           BETA1=PI
!           HBETA=(BETA1-BETA0)/FLOAT(NBETA-1)
           HBETA=(BETA1-BETA0)/FLOAT(NBETA)
        DO 34 I=1,NBETA
34         BETAM(I)=BETA0+(I-1)*HBETA

           GAMMA0=0.0
           GAMMA1=2.0*PI
!           HGAMMA=(GAMMA1-GAMMA0)/FLOAT(NGAMMA-1)
           HGAMMA=(GAMMA1-GAMMA0)/FLOAT(NGAMMA)
        DO 36 I=1,NGAMMA
36         GAMMAM(I)=GAMMA0+(I-1)*HGAMMA

!      PAR(1)=U0
!      PAR(2)=U1
!      PAR(3)=HU
!      PAR(4)=V0
!      PAR(5)=V1
!      PAR(6)=HV
!      PAR(7)=W0
!      PAR(8)=W1
!!      PAR(9)=HW

!!c      PAR(10)=PHI0
!!c      PAR(11)=PHI1
!      PAR(12)=HPHI
!      PAR(15)=HTET

!!c      PAR(13)=TET0
!!c      PAR(14)=TET1

!      PAR(16)=X0
!      PAR(17)=X1
!      PAR(18)=HX

!      PAR(19)=Y0
!      PAR(20)=Y1
!      PAR(21)=HY

!      PAR(22)=Z0
!      PAR(23)=Z1
!      PAR(24)=HZ

!      PAR(33)=RR0
!      PAR(34)=RR1
      PAR(35)=HPP
      PAR(38)=HALPHA
      PAR(39)=HBETA
      PAR(40)=HGAMMA
!      PAR(41)=HUPP

!      write(*,*) 'DATA', HPP

      RETURN
      END SUBROUTINE DATA
!----------------------------------------------------------


      SUBROUTINE ERRPRN(ITEX,IERR)
      CHARACTER*8 ITEX
!      COMMON/PRIN/IPRNT
      WRITE(*,100) ITEX,IERR
 100  FORMAT(' MODUL''<< ',A8,' >>-IERR= ',I6)
      RETURN
      END SUBROUTINE ERRPRN
!------------------------------------------------



      SUBROUTINE SPRINT_3D_3(IND,KSEC1,KSEC2,DR,XM,YM,ZM,NX,NY,NZ, &
                            KAN1,KAN2,KAN3)
! print of some section  
      IMPLICIT NONE
      INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
      INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single 

      INTEGER :: NX,NY,NZ
      INTEGER :: IND,KSEC1,KSEC2,KAN1,KAN2,KAN3

      REAL(wp) DR(NX*NY*NZ)  
      REAL(wp) XM(NX),YM(NY),ZM(NZ)

! local
      REAL(wp) DR1(NX*NY*3)   
      INTEGER*4 I1,I2,I3,I21,I321     

!  IND=-1  -- section XZ
!  IND= 0  -- section XY
!  IND= 1  -- section YZ

      IF(IND) 2,4,6

 2    CONTINUE
      DO 10 I3=1,NZ
      DO 10 I2=KSEC1,KSEC1
      DO 10 I1=1,NX
         I21=(I3-1)*NX+I1
         I321=(I3-1)*NY*NX+(I2-1)*NX+I1
         DR1(I21)=DR(I321)
 
 10   CONTINUE
      CALL IDLFORM_2D (XM,ZM,NX,NZ,DR1,KAN1)
      CALL GNUFORM_2D (NX,NZ,DR1,KAN2)
      CALL GNUFORM_1D (KSEC2,NX,NZ,DR1,KAN3)
      GOTO 30

 4    CONTINUE
      DO 12 I3=KSEC1,KSEC1
      DO 12 I2=1,NY
      DO 12 I1=1,NX
         I21=(I2-1)*NX+I1
         I321=(I3-1)*NY*NX+(I2-1)*NX+I1
         DR1(I21)=DR(I321)
 12   CONTINUE
      CALL IDLFORM_2D (XM,YM,NX,NY,DR1,KAN1)
      CALL GNUFORM_2D (NX,NY,DR1,KAN2)
      CALL GNUFORM_1D (KSEC2,NX,NY,DR1,KAN3)
      GOTO 30

 6    CONTINUE
      DO 14 I3=1,NZ
      DO 14 I2=1,NY
      DO 14 I1=KSEC1,KSEC1
         I21=(I3-1)*NY+I2
         I321=(I3-1)*NY*NX+(I2-1)*NX+I1
         DR1(I21)=DR(I321)

 14   CONTINUE

      CALL IDLFORM_2D (YM,ZM,NY,NZ,DR1,KAN1)
      CALL GNUFORM_2D (NY,NZ,DR1,KAN2)
      CALL GNUFORM_1D (KSEC2,NY,NZ,DR1,KAN3)

 30   CONTINUE
      RETURN
      END SUBROUTINE SPRINT_3D_3
!-----------------------------------------------

      SUBROUTINE GNUFORM_2D (NX,NY,GC,KAN)

      IMPLICIT NONE
      INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
      INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single 

      INTEGER :: NX,NY,KAN
      INTEGER :: I,J

      REAL(wp) :: GC(NX*NY)

      DO 10 J=1,NY
      WRITE(KAN,100) (GC((J-1)*NX+I),I=1,NX)
      WRITE(KAN, *)
!      WRITE(KAN,'(A1)') char(13)//char(10)
10    CONTINUE 
100   FORMAT(E10.3)
!110   FORMAT(A4)
      CLOSE(KAN)
      RETURN
      END  SUBROUTINE GNUFORM_2D
!--------------------------------------------------

      SUBROUTINE IDLFORM_3D (XM,YM,ZM,NX,NY,NZ,GC,KAN)

      IMPLICIT NONE
      INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
      INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single 

      INTEGER*4 NX,NY,NZ,KAN
      INTEGER*4 I1,I2,I3,I321

      REAL(wp) GC(NX*NY*NZ)
      REAL(sp) XM(NX),YM(NY),ZM(NZ)
      
      DO 4 I3=1,NZ
      DO 4 I2=1,NY
      DO 4 I1=1,NX
         I321=(I3-1)*NX*NY+(I2-1)*NX+I1

      WRITE(KAN,100) GC(I321),XM(I1),YM(I2),ZM(I3)
!      WRITE(KAN) GC(I321)

4     CONTINUE

100   FORMAT(4(1X,E10.3))
      RETURN
      END SUBROUTINE IDLFORM_3D
!--------------------------------------------------

      SUBROUTINE IDLFORM_2D (XM,YM,NX,NY,GC,KAN)

      IMPLICIT NONE
      INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
      INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single 

      INTEGER*4 NX,NY,KAN
      INTEGER*4 I,J

      REAL(wp) GC(NX,NY)
      REAL(sp) XM(NX),YM(NY)
   
      WRITE(KAN,100) ((GC(I,J),I=1,NX),J=1,NY),&
                    (XM(I),I=1,NX),(YM(I),I=1,NY)
100   FORMAT(3(1X,E10.3))
      RETURN
      END SUBROUTINE IDLFORM_2D
!-------------------------------------------------

      SUBROUTINE GNUFORM_1D (KSEC,NX,NY,GC,KAN)

      IMPLICIT NONE
      INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
      INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single 

      INTEGER*4 KSEC,NX,NY,KAN
      INTEGER*4 I,J

      REAL(wp) GC(NX*NY)

      DO 10 J=KSEC,KSEC
      WRITE(KAN,100) (GC((J-1)*NX+I),I=1,NX)      
      WRITE(KAN, *)
!      WRITE(KAN,'(A1)') char(13)//char(10)
10    CONTINUE 
100   FORMAT(E10.3)
!110   FORMAT(A4)
      RETURN
      END SUBROUTINE GNUFORM_1D

!--------------------------------------------------
      SUBROUTINE GNUFORM_11D (NX,NY,K,GC,KAN)

!        USE CHANALS

      IMPLICIT NONE
      INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
      INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=double 
      INTEGER :: NX,NY,KAN
      INTEGER :: I,K
      REAL(wp) :: GC(NX*NY)
      REAL(sp) :: GCC(NX*NY)
      
!      DO 2 I=1,NX*NY
!         GCC(I)=SNGL(GC(I))
!2     CONTINUE

!      DO 10 J=1,NY
      WRITE(KAN,100) (GC((K-1)*NX+I),I=1,NX)
      WRITE(KAN, *)
!      WRITE(KAN,'(A1)') char(13)//char(10)
!10    CONTINUE 
100   FORMAT(E10.3)
!110   FORMAT(A4)
      CLOSE(KAN)
      RETURN
      END SUBROUTINE GNUFORM_11D
!--------------------------------------------------

      SUBROUTINE GNUFORM_12D (NX,NY,GC,KAN)

      IMPLICIT NONE
      INTEGER*4 KSEC,NX,NY,KAN
      INTEGER*4 I,J

      REAL(4) GC(NX*NY)

      DO 10 J=1,NX*NY
      WRITE(KAN,100) GC(J),J      
!      WRITE(KAN, *)
!      WRITE(KAN,'(A1)') char(13)//char(10)
10    CONTINUE 
100   FORMAT(E10.3,I3)
!110   FORMAT(A4)
      RETURN
      END SUBROUTINE GNUFORM_12D
!--------------------------------------------------------


      SUBROUTINE REGN(XM,YM,ZM,NX,NY,NZ,X1,Y1,Z1,EI)

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=double 
        REAL(wp) :: XM(NX),YM(NY),ZM(NZ)
!        REAL(wp) :: X11,Y11,Z11
        INTEGER :: NX,NY,NZ
        INTEGER :: EI
!local
        REAL(wp) :: X1,Y1,Z1

!        X1=SNGL(X11)
!        Y1=SNGL(Y11)
!        Z1=SNGL(Z11)

!        IF(X1.LT.XM(1).OR.X1.GT.XM(NX)) GO TO 10
!        IF(Y1.LT.YM(1).OR.Y1.GT.YM(NY)) GO TO 10
!        IF(Z1.LT.ZM(1).OR.Z1.GT.ZM(NZ)) GO TO 10
!        EI=1
!        RETURN
!10      EI=-1
        IF(X1 .GE. XM(1) .AND. X1 .LE. XM(NX) &
           .AND. Y1 .GE. YM(1) .AND. Y1 .LE. YM(NY) &
           .AND. Z1 .GE. ZM(1) .AND. Z1 .LE. ZM(NZ)) THEN
           EI=1
        ELSE
           EI=-1
        ENDIF
       RETURN
      END SUBROUTINE REGN
!----------------------------------------------    
      SUBROUTINE INDNW3(XM,YM,ZM,NX,NY,NZ,      &
                        HX,HY,HZ,XNEW,YNEW,ZNEW,&
                        J1,J2,J3,J4)
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 

        REAL(wp) :: XM(NX),YM(NY),ZM(NZ)
        REAL(wp) :: HX,HY,HZ
        REAL(wp) :: XNEW,YNEW,ZNEW
        INTEGER :: NX,NY,NZ
        INTEGER :: J1,J2,J3,J4

        INTEGER :: I1,I2,I3,I11,I22,I33
        INTEGER :: I111,I222,I333,KX,KY,KZ
        INTEGER :: LEFTX,LEFTY,LEFTZ

        KX=NX-1
        DO 46 I1=1,KX
           I11=I1
           IF(XNEW-XM(I1+1)) 48,46,46
46         CONTINUE

48         KY=NY-1
           DO 53 I2=1,KY
              I22=I2
              IF(YNEW-YM(I2+1)) 54,53,53
53         CONTINUE

54            KZ=NZ-1
           DO 55 I3=1,KZ
              I33=I3
              IF(ZNEW-ZM(I3+1)) 57,55,55
55         CONTINUE

57         CONTINUE

              LEFTX=1
              LEFTY=1
              LEFTZ=1
              
              IF(ABS(XNEW-XM(I11)).LE.0.5*HX) LEFTX=0
              IF(ABS(YNEW-YM(I22)).LE.0.5*HY) LEFTY=0
              IF(ABS(ZNEW-ZM(I33)).LE.0.5*HZ) LEFTZ=0

              I111=I11+LEFTX
              I222=I22+LEFTY
              I333=I33+LEFTZ

              J1=I111
              J2=I222
              J3=I333
              J4=(I333-1)*NX*NY+(I222-1)*NX+I111

      RETURN
      END SUBROUTINE INDNW3
!--------------------------------------------------


      SUBROUTINE INDNW2(XM,YM,NX,NY,HX,HY,XNEW,YNEW,J1,J2,J4)
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=singl 

        INTEGER :: NX,NY
        REAL(wp) :: XM(NX),YM(NY)
        REAL(wp) :: HX,HY,XNEW,YNEW
        INTEGER  :: J1,J2,J4
! local
        INTEGER :: I1,I2,I3,I11,I22
        INTEGER :: I111,I222,KX,KY
        INTEGER :: LEFTX,LEFTY

        KX=NX-1
        DO 46 I1=1,KX
           I11=I1
           IF(XNEW-XM(I1+1)) 48,46,46
46         CONTINUE

48         KY=NY-1
           DO 53 I2=1,KY
              I22=I2
              IF(YNEW-YM(I2+1)) 57,53,53
53         CONTINUE

57         CONTINUE

              LEFTX=1
              LEFTY=1
              
              IF(ABS(XNEW-XM(I11)).LE.0.5*HX) LEFTX=0
              IF(ABS(YNEW-YM(I22)).LE.0.5*HY) LEFTY=0

              I111=I11+LEFTX
              I222=I22+LEFTY

              J1=I111
              J2=I222
              J4=(I222-1)*NX+I111

      RETURN
      END SUBROUTINE INDNW2
!------------------------------------------------------------

      SUBROUTINE INDNEIGHBOR(XM,YM,ZM,NX,NY,NZ, &
                        XNEW,YNEW,ZNEW,J1,J2,J3,J4)

        USE PARAM_INPUT
        IMPLICIT NONE

        REAL(4) XM(NX),YM(NY),ZM(NZ)
        REAL(4) XNEW,YNEW,ZNEW
!        REAL(4) HX,HY,HZ
        INTEGER*4 NX,NY,NZ
        INTEGER*4 J1,J2,J3,J4

        INTEGER*4 I1,I2,I3

        REAL(4) RR,RRR
        REAL(4) XX,YY,ZZ
        REAL(4) SS1,SS2

        RR=PARAM(5)
        
        RRR=RR/5.
        

        J1=1
        J2=1
        J3=1

        SS1=(XNEW-XM(1))**2+(YNEW-YM(1))**2+(ZNEW-ZM(1))**2
        
        
        DO 40 I3=2,NZ
           ZZ=ZM(I3)

        DO 40 I2=2,NY
           YY=YM(I2)
           
        DO 40 I1=2,NX
           XX=XM(I1)

        SS2=(XNEW-XX)**2+(YNEW-YY)**2+(ZNEW-ZZ)**2

        IF(SS2 .LT. RRR) THEN

           IF(SS2 .GE. SS1) GOTO 40
           SS1=SS2
           J1=I1
           J2=I2
           J3=I3

        ENDIF
40      CONTINUE

        J4=(J3-1)*NX*NY+(J2-1)*NX+J1
   
      RETURN
      END SUBROUTINE INDNEIGHBOR
!------------------------------------------------------
      SUBROUTINE EMIN(N,F,E,K,IERR)

        IMPLICIT NONE

        REAL(4) F(N)
        REAL(4) E
        INTEGER*4 N,K,IERR 

        INTEGER*4 I 

        K=1
        E=F(1)
        DO 1 I=1,N
           IF(F(I).GE.E)GO TO 1
           E=F(I)
           K=I
1       CONTINUE
      RETURN
      END SUBROUTINE EMIN

!-------------------------------------------------------

      FUNCTION REGION(HW,XX)
        IMPLICIT NONE

        REAL(4) REGION
        REAL(4) HW,XX

        IF(XX .GE.-HW .AND. XX .LE. HW) THEN 
           REGION=1.
        ELSE
           REGION=-1.
        ENDIF

      RETURN
      END FUNCTION REGION
!----------------------------------------------------

      FUNCTION FUNL3(XM,YM,ZM,GM,N1,N2,N3,X1,Y1,Z1,IERR)
!       - 3D interpolation 

        IMPLICIT NONE 
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single

        REAL(wp) :: XM(N1),YM(N2),ZM(N3)
        REAL(wp) :: GM(N1*N2*N3)
!        REAL(wp) :: X11,Y11,Z11
        REAL(wp) :: FUNL3,FUN
        INTEGER :: N1,N2,N3,IERR

        REAL(wp) :: X1,Y1,Z1
        REAL(wp) :: A17,R1,R2,R3
        REAL(wp) :: HX1,HX2,HX,HY1,HY2,HY,HZ1,HZ2,HZ
        REAL(wp) :: FY1,FY2,FXY1,FXY2

        INTEGER :: I1,I2,I3,I11,I22,I33,K1,K2,K3
        INTEGER :: J1,J2,J3,J4,J5,J6,J7,J8
        CHARACTER*8 ITEX
!XM(N1)+R1-X1

!        X1=SNGL(X11)
!        Y1=SNGL(Y11)
!        Z1=SNGL(Z11)

        ITEX='FUNL3'
        IERR=0

        A17=9.E32_wp   !9.E17

        R1=(XM(N1)-XM(1))*1.E-10_wp          !*1.E-5
        R2=(YM(N2)-YM(1))*1.E-10_wp          ! *1.E-5
        R3=(ZM(N3)-ZM(1))*1.E-10_wp           !*1.E-5

        IF(X1-XM(1)+R1)  2,10,10
2       IERR=1
        CALL ERRPRN(ITEX,IERR)
        RETURN

10      IF(XM(N1)+R1-X1) 4,11,11
4       IERR=2
        CALL ERRPRN(ITEX,IERR)
        RETURN

11      IF(Y1-YM(1)+R2)  6,12,12
6       IERR=3
        CALL ERRPRN(ITEX,IERR)
        RETURN

12      IF(YM(N2)+R2-Y1) 8,13,13
8       IERR=4
        CALL ERRPRN(ITEX,IERR)
        RETURN

13      IF(Z1-ZM(1)+R3)  14,15,15
14      IERR=5
        CALL ERRPRN(ITEX,IERR)
        RETURN

15      IF(ZM(N3)+R3-Z1) 17,18,18
17      IERR=6
        CALL ERRPRN(ITEX,IERR)
        RETURN

18      CONTINUE
        K1=N1-1
        DO 21 I1=1,K1
           I11=I1
           IF(X1-XM(I1+1)) 22,21,21
21         CONTINUE

22         K2=N2-1
        DO 23 I2=1,K2
           I22=I2
           IF(Y1-YM(I2+1)) 24,23,23
23         CONTINUE

24         K3=N3-1
        DO 25 I3=1,K3
           I33=I3
           IF(Z1-ZM(I3+1)) 27,25,25
25         CONTINUE

27         CONTINUE
           J1=(I33-1)*N1*N2+(I22-1)*N1+I11
           J2=(I33-1)*N1*N2+(I22-1)*N1+I11+1
           J3=(I33-1)*N1*N2+(I22+1-1)*N1+I11+1
           J4=(I33-1)*N1*N2+(I22+1-1)*N1+I11
           J5=(I33+1-1)*N1*N2+(I22-1)*N1+I11
           J6=(I33+1-1)*N1*N2+(I22-1)*N1+I11+1
           J7=(I33+1-1)*N1*N2+(I22+1-1)*N1+I11+1
           J8=(I33+1-1)*N1*N2+(I22+1-1)*N1+I11

           IF(GM(J1)+A17) 29,39,39
29         GM(J1)=-A17

39         IF(GM(J4)+A17) 30,40,40
30         GM(J4)=-A17

40         IF(GM(J8)+A17) 32,42,42
32         GM(J8)=-A17

42         IF(GM(J5)+A17) 34,44,44
34         GM(J5)=-A17

44         IF(GM(J2)+A17) 36,46,46
36         GM(J2)=-A17

46         IF(GM(J3)+A17) 38,48,48
38         GM(J3)=-A17

48         IF(GM(J7)+A17) 50,52,52
50         GM(J7)=-A17

52         IF(GM(J6)+A17) 54,56,56
54         GM(J6)=-A17

56         CONTINUE

           HX1=X1-XM(I11)
           HX2=XM(I11+1)-X1
           HX=XM(I11+1)-XM(I11)

           HY1=Y1-YM(I22)
           HY2=YM(I22+1)-Y1
           HY=YM(I22+1)-YM(I22)

           HZ1=Z1-ZM(I33)
           HZ2=ZM(I33+1)-Z1
           HZ=ZM(I33+1)-ZM(I33)

           FY1=GM(J1)+(GM(J4)-GM(J1))/HY*HY1
           FY2=GM(J2)+(GM(J3)-GM(J2))/HY*HY1
           FXY1=FY1+(FY2-FY1)/HX*HX1

           FY1=GM(J5)+(GM(J8)-GM(J5))/HY*HY1
           FY2=GM(J6)+(GM(J7)-GM(J6))/HY*HY1
           FXY2=FY1+(FY2-FY1)/HX*HX1
           
           FUN=FXY1+(FXY2-FXY1)/HZ*HZ1

           FUNL3=FUN
      RETURN
      END FUNCTION FUNL3
!-----------------------------------------------------



      FUNCTION FACT(NN,HX)

      IMPLICIT NONE
!  output -1*3*5...*(2*NN-1)*HX**(-(NN+1))

      INTEGER*4 NN
      REAL(8) HX,FACT

      INTEGER*4 I
      REAL(8) SS

      SS=1.
      DO 4 I=1,2*NN-1,2 
         SS=SS*I
4     CONTINUE   

      FACT=-SS/HX**(NN+1)

      RETURN
      END FUNCTION FACT
!----------------------------------------------------

      FUNCTION PLGNDR(L,MM,X)
!                                                     m
!       computers the associated Legendre polinomial P (x). 
!                                                     l
!       Here m and l are integer satisfying  -l <= m <= l, while x 
!       -1 <= x <=1
        IMPLICIT NONE 
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single

      INTEGER :: L,MM
      REAL(wp) :: PLGNDR,X
! local
      INTEGER :: M
      INTEGER :: I,J,LL,K
      REAL(wp) :: FACT, PLL,PMM,PMMP1,SOMX2
      REAL(wp) :: NOR

      M=ABS(MM)
!--------------------------
! Norm for negative MM
!      K=L-M
!      NOR=1. 
!      DO 4 I=1, 2*M
!      J=K+I 
!      NOR=NOR/FLOAT(J)*(-1)**M
! 4    CONTINUE
!---------------------------
!      IF(M .LT. 0 .OR. ABS(X) .GT. 1.) THEN 
      IF(ABS(X) .GT. 1.) THEN 
         WRITE(*,*) 'BAD ARGUMENTS IN PLGNDR:  X', X
      STOP
      RETURN
      ENDIF
         PMM=1.                    ! Computer Pmm      
         IF(M .GT. 0) THEN
            SOMX2=SQRT((1.-X)*(1.+X))
            FACT=1.
            DO 11 I=1,M
               PMM=-PMM*FACT*SOMX2
               FACT=FACT+2
11          CONTINUE
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
12          CONTINUE
               PLGNDR=PLL              !*PNORM(L,M)
            ENDIF
         ENDIF
      IF(MM .LT. 0) THEN
         PLGNDR=PLGNDR*(-1)**M
      ENDIF   
      RETURN
      END FUNCTION PLGNDR 
!-----------------------------------------------------------

      FUNCTION PLGNDR_NORM(L,MM,X)
!                                                              m
!       computers the associated norming Legendre polinomial  P (x). 
!                                                              l
!       Here m and l are integer satisfying  -l <= m <= l, while x 
!       -1 <= x <=1
        IMPLICIT NONE 
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single

      INTEGER :: L,MM
      REAL(wp) :: PLGNDR_NORM,X
      INTEGER :: I,LL,M
      REAL(wp) :: FACT,PLL,PMM,PMMP1,SOMX2
! functions
      REAL(wp) :: PNORM

      M=ABS(MM)

      IF(ABS(X) .GT. 1.) THEN 
         WRITE(*,*) 'BAD ARGUMENTS IN PLGNDR_NORM: M or X',X
      STOP
      RETURN
      ENDIF
   
         PMM=1.                    ! Computer Pmm      
         IF(M .GT. 0) THEN
            SOMX2=SQRT((1.-X)*(1.+X))
            FACT=1.
            DO 11 I=1,M
               PMM=-PMM*FACT*SOMX2
               FACT=FACT+2
11          CONTINUE
         ENDIF
         IF(L .EQ. M) THEN
            PLGNDR_NORM=PMM *PNORM(L,M)
         ELSE
            PMMP1=X*(2*M+1)*PMM      !Compute Pm+1,m     
            IF(L .EQ. M+1) THEN
               PLGNDR_NORM=PMMP1*PNORM(L,M)
            ELSE
            DO 12 LL=M+2,L
               PLL=(X*(2*LL-1)*PMMP1-(LL+M-1)*PMM)/(LL-M)
               PMM=PMMP1
               PMMP1=PLL
12          CONTINUE
               PLGNDR_NORM=PLL*PNORM(L,M)
            ENDIF
         ENDIF
      IF(MM .LT. 0) THEN
         PLGNDR_NORM=PLGNDR_NORM*(-1)**M
      ENDIF 
      RETURN
      END FUNCTION PLGNDR_NORM 
!--------------------------------------------------------
      FUNCTION PNORM(L,M)
!           Computers the norm of associated Legendre polynomial 
        IMPLICIT NONE 
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 
      INTEGER :: L,M
      INTEGER :: K,I,J
      REAL(wp) :: S
! functions
      REAL(wp) :: PNORM

      K=L-M
      S=1. 
      DO 4 I=1, 2*M
      J=K+I 
      S=S/FLOAT(J)
 4    CONTINUE
!      PNORM=SQRT((2*L+1)*S/(4.*PI))
      PNORM=SQRT((2*L+1)*S/2.0)
      RETURN
      END FUNCTION PNORM
!-------------------------------------------------------------------

      FUNCTION CYLM (L,M,TETA,PHI)
!          Computation the spherical harmonics L, M
!          YLMRE - real part
!          YLMIM - image part
        IMPLICIT NONE 
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

      INTEGER :: L,M
      REAL(wp) :: TETA,PHI
!      REAL(wp) :: YLMRE,YLMIM
      COMPLEX(wp) :: CYLM
! functions
      REAL(wp) :: PNORM,PLGNDR

!local 
      INTEGER :: MM,KK
      REAL(wp) :: X,PNN,PLAG,SIGNUM
      COMPLEX(wp) :: CCEXP

      CCEXP(X)=CMPLX(COS(X),SIN(X))

      MM=ABS(M)
      KK=MOD(MM,2)
      X=COS(TETA) 

!      PNN=PNORM(L,MM)
      PNN=PNORM(L,MM)/SQRT(2.*PI)
      PLAG=PLGNDR(L,MM,X)

      CYLM= PNN*PLAG*CCEXP(MM*PHI)  

      IF(M .LT. 0) THEN 
           IF(KK .EQ. 1) THEN
              SIGNUM=-1.
           ELSE
              SIGNUM=1.
           ENDIF  
        CYLM= SIGNUM*PNN*PLAG*CCEXP(-MM*PHI)
!      ELSE
!        CYLM= PNN*PLAG*CCEXP(MM*PHI)  
      ENDIF

      RETURN
      END FUNCTION CYLM 
!-----------------------------------------------------

      FUNCTION FUNL2(GM,XM,YM,N1,N2,X1,Y1,IERR)
!       - 2D interpolation 
        
        IMPLICIT NONE 
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single

        INTEGER :: N1,N2,IERR
        REAL(wp) :: XM(N1),YM(N2)
        REAL(wp) :: GM(N1*N2)
! functions
        REAL(wp) :: FUNL2
! local
        REAL(wp) :: X1,Y1,FUN
        REAL(wp) :: A17,R1,R2
        REAL(wp) :: HX1,HX2,HX,HY1,HY2,HY !,HZ1,HZ2,HZ
        REAL(wp) :: FY1,FY2,FXY1,FXY2

        INTEGER :: I1,I2,I3,I11,I22,K1,K2
        INTEGER :: J1,J2,J3,J4    !,J5,J6,J7,J8
        CHARACTER*8 ITEX
!
      ITEX='FUNL2'
      A17=9.E17
!
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
!
      K1=N1-1
      DO 21 I1=1,K1
      I11=I1
      IF(X1-XM(I1+1)) 22,21,21
  21  CONTINUE
!
  22  K2=N2-1
      DO 23 I2=1,K2
      I22=I2
      IF(Y1-YM(I2+1)) 27,23,23
  23  CONTINUE
!
  27  CONTINUE
      J1=(I22-1)*N1+I11
      J2=(I22-1)*N1+I11+1
      J3=(I22+1-1)*N1+I11+1
      J4=(I22+1-1)*N1+I11
!
      IF(GM(J1)+A17) 28,29,29
  28  IERR=5
      CALL ERRPRN(ITEX,IERR)
!
  29  IF(GM(J4)+A17) 30,34,34
  30  IERR=6
      CALL ERRPRN(ITEX,IERR)
!
  34  IF(GM(J2)+A17) 32,36,36
  32  IERR=7
      CALL ERRPRN(ITEX,IERR)
!
  36  IF(GM(J3)+A17) 40,42,42
  40  IERR=8
      CALL ERRPRN(ITEX,IERR)
!
  42  CONTINUE
!
      HX1=X1-XM(I11)
      HX2=XM(I11+1)-X1
      HX=XM(I11+1)-XM(I11)
!
      HY1=Y1-YM(I22)
      HY2=YM(I22+1)-Y1
      HY=YM(I22+1)-YM(I22)
!
      FY1=GM(J1)+(GM(J4)-GM(J1))/HY*HY1
      FY2=GM(J2)+(GM(J3)-GM(J2))/HY*HY1
      FXY1=FY1+(FY2-FY1)/HX*HX1
!
      FUN=FXY1
      FUNL2=FUN
      RETURN
      END FUNCTION FUNL2
!---------------------------------------------------------
          SUBROUTINE POLAR_2D(X,Y,RR,FI)
!    procedure perfom a passage to spherical coordinates RR, FI
!    by the use of decart X,Y coordinate
!    SIGNY,SIGNZ- signs of  X,Y
!    RR- radial- and FI- angle ccordinates
          IMPLICIT NONE
          INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
          INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=double 
          REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

          REAL(wp) :: X,Y,RR,FI
          REAL(wp) :: SIGNY  !,SIGNZ

          SIGNY=SIGN(1.0_wp,Y)
          IF(X.NE.0.0_wp)GO TO 30
          IF(Y.NE.0.0_wp) THEN
!             FI=0.5*(1.-SIGNY)*PI+PI*0.5
             FI=PI/2.0*(2.0-SIGNY)
          ELSE
             FI=0.0_wp
          ENDIF   
          GO TO 35
30        FI=ATAN2(Y,X) 
          IF(FI .LT. 0.0) FI=FI+ PI* 2.0    ! for the case when  FI \in (0,2*PI)
35        CONTINUE
          RR=SQRT(X**2+Y**2)
        RETURN
        END SUBROUTINE POLAR_2D
!--------------------------------------------
      SUBROUTINE REGN2(XM,YM,NX,NY,X1,Y1,EI)
          IMPLICIT NONE
          INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
          INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
          REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

          INTEGER :: NX,NY,EI
          REAL(wp) :: XM(NX),YM(NY),X1,Y1

          IF(X1.LT.XM(1).OR.X1.GT.XM(NX)) GO TO 10
          IF(Y1.LT.YM(1).OR.Y1.GT.YM(NY)) GO TO 10
          EI=1
          RETURN
10        EI=-1
          RETURN
      END SUBROUTINE REGN2
!-----------------------------------------------------------
      SUBROUTINE POLAR_DECART2D_2(GM,XXM,YYM,NX,NY,RRM,PHIM,NRR,NPHI,DR)

        USE PARAMTR
          IMPLICIT NONE
          INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
          INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
          REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

          REAL(wp) :: GM(NRR*NPHI)
          REAL(wp) :: XXM(NX),YYM(NY)
          REAL(wp) :: RRM(NRR),PHIM(NPHI)

          REAL(wp) :: DR(NX*NY)
          INTEGER :: NX,NY,NZ,NRR,NPHI

! local
          REAL(wp) :: XX,YY  !,RXY                  !!!!!!!!!
          REAL(wp) :: RRR,PHI,HX,HY  !,HPHI     !HZ
          INTEGER :: KK,I,I1,I2,I3,I21,IERR,EI
          INTEGER :: J1,J2,J3,J4

          HX=PAR(18)
          HY=PAR(21)

             DO 6 I2=1,NPHI
                PHI=PHIM(I2)
             DO 6 I1=1,NRR
                RRR=RRM(I1)

                I21=(I2-1)*NRR+I1

                XX=RRR*COS(PHI)
                YY=RRR*SIN(PHI)

                CALL INDNW2(XXM,YYM,NX,NY,HX,HY,XX,YY,J1,J2,J4)

                CALL REGN2(XXM,YYM,NX,NY,XX,YY,EI)

                IF(EI .GT. 0) THEN
                   DR(J4)=GM(I21)
                ELSE
                   DR(J4)=0.0
                ENDIF

6            CONTINUE
          RETURN
          END SUBROUTINE POLAR_DECART2D_2
!------------------------------------------

      SUBROUTINE POLAR_DECART_3D(GM,XXM,YYM,ZZM,NX,NY,NZ,RRM,TETAM,PHIM, &
                              NRR,NTETA,NPHI,DR)
!
!        USE PARAMTR
          IMPLICIT NONE
          INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
          INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
          REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

          REAL(wp) :: GM(NRR*NTETA*NPHI)
          REAL(wp) :: XXM(NX),YYM(NY),ZZM(NZ)
          REAL(wp) :: RRM(NRR),TETAM(NTETA),PHIM(NPHI)

          REAL(wp) :: DR(NX*NY*NZ) 
          REAL(wp) :: FUNL3  !,ALINEAR_3D
          INTEGER :: NX,NY,NZ,NRR,NTETA,NPHI

! local
!          REAL(wp) :: GMM(NRR*NTETA*NPHI)
          REAL(wp) :: XX,YY,ZZ
          REAL(wp) :: RRR,TETA,PHI  !,HPP,HTET,HPHI     !!!!!!!!!!!!!!!
          INTEGER :: I1,I2,I3,I21,I321,EI,IERR

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
                  DR(I321)=FUNL3(PHIM,TETAM,RRM,GM,NPHI,NTETA,NRR,PHI,TETA,RRR,IERR)
                ENDIF

                IF(RRR .GE. RRM(NRR)) THEN
                  DR(I321)=0.0
                ENDIF  

6            CONTINUE
          RETURN
          END SUBROUTINE POLAR_DECART_3D
!------------------------------------------

      FUNCTION Burke_3D(AM,XM,YM,ZM,NX,NY,NZ,XNEW,YNEW,ZNEW)

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 

        INTEGER :: NX,NY,NZ
        REAL(wp) :: AM(NX*NY*NZ)
        REAL(wp) :: XM(NX),YM(NY),ZM(NZ)
        REAL(wp) :: XNEW,YNEW,ZNEW
! local
        REAL(wp) :: HX,HY,HZ,SS1,SS2
        REAL(wp) :: WW(8)
        REAL(wp) :: XX,YY,ZZ,BRK
        INTEGER :: I11,I22,I33,I
        INTEGER :: JM(8)
        INTEGER :: P
! functions
        REAL(wp) :: Burke_3D,DIST_3D

        P=4
        HX=XM(2)-XM(1)
        HY=YM(2)-YM(1)
        HZ=ZM(2)-ZM(1)

        I11=INT((XNEW-XM(1))/HX)+1
        I22=INT((YNEW-YM(1))/HY)+1
        I33=INT((ZNEW-ZM(1))/HZ)+1
        
        JM(1)=(I33-1)*NX*NY+(I22-1)*NX+I11
        JM(2)=(I33-1)*NX*NY+(I22-1)*NX+I11+1
        JM(3)=(I33-1)*NX*NY+(I22+1-1)*NX+I11+1
        JM(4)=(I33-1)*NX*NY+(I22+1-1)*NX+I11
        JM(5)=(I33+1-1)*NX*NY+(I22-1)*NX+I11
        JM(6)=(I33+1-1)*NX*NY+(I22-1)*NX+I11+1
        JM(7)=(I33+1-1)*NX*NY+(I22+1-1)*NX+I11+1
        JM(8)=(I33+1-1)*NX*NY+(I22+1-1)*NX+I11

        IF(XNEW .NE. XM(I11) .AND. YNEW .NE. YM(I22) .AND. ZNEW .NE. ZM(I33)) THEN
           WW(1)= DIST_3D(XM,YM,ZM,NX,NY,NZ,XNEW,YNEW,ZNEW,I11,I22,I33,P)
        ELSE
           BRK=AM(JM(1))
           GOTO 8
        ENDIF   

        IF(XNEW .NE. XM(I11+1) .AND. YNEW .NE. YM(I22) .AND. ZNEW .NE. ZM(I33)) THEN
           WW(2)= DIST_3D(XM,YM,ZM,NX,NY,NZ,XNEW,YNEW,ZNEW,I11+1,I22,I33,P)
        ELSE
           BRK=AM(JM(2))
           GOTO 8
        ENDIF   

        IF(XNEW .NE. XM(I11+1) .AND. YNEW .NE. YM(I22+1) .AND. ZNEW .NE. ZM(I33)) THEN
           WW(3)= DIST_3D(XM,YM,ZM,NX,NY,NZ,XNEW,YNEW,ZNEW,I11+1,I22+1,I33,P) 
         ELSE
           BRK=AM(JM(3))
           GOTO 8
        ENDIF 

        IF(XNEW .NE. XM(I11) .AND. YNEW .NE. YM(I22+1) .AND. ZNEW .NE. ZM(I33)) THEN
           WW(4)= DIST_3D(XM,YM,ZM,NX,NY,NZ,XNEW,YNEW,ZNEW,I11,I22+1,I33,P)
         ELSE
           BRK=AM(JM(4))
           GOTO 8
        ENDIF 

        IF(XNEW .NE. XM(I11) .AND. YNEW .NE. YM(I22) .AND. ZNEW .NE. ZM(I33+1)) THEN
           WW(5)= DIST_3D(XM,YM,ZM,NX,NY,NZ,XNEW,YNEW,ZNEW,I11,I22,I33+1,P)
         ELSE
           BRK=AM(JM(5))
           GOTO 8
        ENDIF 

        IF(XNEW .NE. XM(I11+1) .AND. YNEW .NE. YM(I22) .AND. ZNEW .NE. ZM(I33+1)) THEN
           WW(6)= DIST_3D(XM,YM,ZM,NX,NY,NZ,XNEW,YNEW,ZNEW,I11+1,I22,I33+1,P)
         ELSE
           BRK=AM(JM(6))
           GOTO 8
        ENDIF 

        IF(XNEW .NE. XM(I11+1) .AND. YNEW .NE. YM(I22+1) .AND. ZNEW .NE. ZM(I33+1)) THEN
           WW(7)= DIST_3D(XM,YM,ZM,NX,NY,NZ,XNEW,YNEW,ZNEW,I11+1,I22+1,I33+1,P)
         ELSE
           BRK=AM(JM(7))
           GOTO 8
        ENDIF 

        IF(XNEW .NE. XM(I11) .AND. YNEW .NE. YM(I22+1) .AND. ZNEW .NE. ZM(I33+1)) THEN
           WW(8)= DIST_3D(XM,YM,ZM,NX,NY,NZ,XNEW,YNEW,ZNEW,I11,I22+1,I33+1,P)
         ELSE
           BRK=AM(JM(8))
           GOTO 8
        ENDIF 

        SS1=0.0
        SS2=0.0
        DO 4 I=1,8 
           SS1= SS1+WW(I)*AM(JM(I))
           SS2=SS2+WW(I)
4       CONTINUE
           BRK=SS1/SS2
8       CONTINUE           
        Burke_3D=BRK
      RETURN
      END FUNCTION Burke_3D
!---------------------------------------------
      FUNCTION DIST_3D(XM,YM,ZM,NX,NY,NZ,XX,YY,ZZ,I1,I2,I3,P)
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 

        INTEGER :: NX,NY,NZ
        REAL(wp) :: XM(NX),YM(NY),ZM(NZ)
        REAL(wp) :: XX,YY,ZZ
        INTEGER :: I1,I2,I3,P
! functions
        REAL(wp) :: DIST_3D

        DIST_3D=1./((XX-XM(I1))**2+(YY-YM(I2))**2+(ZZ-ZM(I3))**2)**(P/2)
      RETURN
      END FUNCTION DIST_3D
!---------------------------------------------------------------

      FUNCTION Burke_2D(AM,XM,YM,NX,NY,XNEW,YNEW)

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 

        INTEGER :: NX,NY
        REAL(wp) :: AM(NX*NY)
        REAL(wp) :: XM(NX),YM(NY)
        REAL(wp) :: XNEW,YNEW
! local
        REAL(wp) :: HX,HY,SS1,SS2
        REAL(wp) :: WW(8)
        REAL(wp) :: XX,YY,BRK
        INTEGER :: I11,I22,I33,I
        INTEGER :: JM(8)
        INTEGER :: P
!        REAL(wp) :: P

! functions
        REAL(wp) :: Burke_2D,DIST_2D

        P=3
        HX=XM(2)-XM(1)
        HY=YM(2)-YM(1)

        I11=INT((XNEW-XM(1))/HX)+1
        I22=INT((YNEW-YM(1))/HY)+1

        JM(1)=(I22-1)*NX+I11
        JM(2)=(I22-1)*NX+I11+1
        JM(3)=(I22+1-1)*NX+I11+1
        JM(4)=(I22+1-1)*NX+I11

        IF(XNEW .NE. XM(I11) .AND. YNEW .NE. YM(I22) ) THEN
           WW(1)= DIST_2D(XM,YM,NX,NY,XNEW,YNEW,I11,I22,P)
        ELSE
           BRK=AM(JM(1))
           GOTO 8
        ENDIF   

        IF(XNEW .NE. XM(I11+1) .AND. YNEW .NE. YM(I22) ) THEN
           WW(2)= DIST_2D(XM,YM,NX,NY,XNEW,YNEW,I11+1,I22,P)
        ELSE
           BRK=AM(JM(2))
           GOTO 8
        ENDIF   

        IF(XNEW .NE. XM(I11+1) .AND. YNEW .NE. YM(I22+1) ) THEN
           WW(3)= DIST_2D(XM,YM,NX,NY,XNEW,YNEW,I11+1,I22+1,P) 
         ELSE
           BRK=AM(JM(3))
           GOTO 8
        ENDIF 

        IF(XNEW .NE. XM(I11) .AND. YNEW .NE. YM(I22+1) ) THEN
           WW(4)= DIST_2D(XM,YM,NX,NY,XNEW,YNEW,I11,I22+1,P)
         ELSE
           BRK=AM(JM(4))
           GOTO 8
        ENDIF 
        SS1=0.0
        SS2=0.0
        DO 4 I=1,4 
           SS1= SS1+WW(I)*AM(JM(I))
           SS2=SS2+WW(I)
4       CONTINUE
           BRK=SS1/SS2
8       CONTINUE           
        Burke_2D=BRK
      RETURN
      END FUNCTION Burke_2D
!---------------------------------------------
      FUNCTION DIST_2D(XM,YM,NX,NY,XX,YY,I1,I2,P)
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 

        INTEGER :: NX,NY
        REAL(wp) :: XM(NX),YM(NY)
        REAL(wp) :: XX,YY
        INTEGER :: I1,I2,P
! functions
        REAL(wp) :: DIST_2D

        DIST_2D=1./((XX-XM(I1))**2+(YY-YM(I2))**2)**(P/2)
      RETURN
      END FUNCTION DIST_2D
!-----------------------------------------------------
      SUBROUTINE POLAR_DECART3D_22(GM,XM,YM,ZM,NX,NY,NZ,RRM,TETAM,PHIM, &
                                  NRR,NTETA,NPHI,NTEN,DR)

        USE PARAMTR
          IMPLICIT NONE
          INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
          INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
          REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

          INTEGER :: NX,NY,NZ,NRR,NTETA,NPHI,NTEN
          REAL(wp) :: GM(NRR*NTETA*NPHI,NTEN*NTEN)
          REAL(wp) :: XM(NX),YM(NY),ZM(NZ)
          REAL(wp) :: RRM(NRR),TETAM(NTETA),PHIM(NPHI)

          REAL(wp) :: DR(NX*NY*NZ,NTEN*NTEN)
!          REAL(wp) :: FUNL3,ALINEAR_3D


! local
!          REAL(wp) :: GMM(NRR*NTETA*NPHI)
          REAL(wp) :: XX,YY,ZZ  !,RXY                  !!!!!!!!!
          REAL(wp) :: RRR,TETA,PHI,HX,HY,HZ !,HPP,HTET,HPHI     !!!!!!!!!!!!!!!
          INTEGER :: KK,I,I1,I2,I3,I21,I321,IERR,EI
          INTEGER :: J1,J2,J3,J4
!          HPP =PAR(35)
!          HTET=PAR(15)
!          HPHI=PAR(12)

          HX=PAR(18)
          HY=PAR(21)
          HZ=PAR(24)

          DO 10 KK=1,NTEN*NTEN

             DO 6 I3=1,NRR
                RRR=RRM(I3)
             DO 6 I2=1,NTETA
                TETA=TETAM(I2)
             DO 6 I1=1,NPHI
                PHI=PHIM(I1)
                I321=(I3-1)*NTETA*NPHI+(I2-1)*NPHI +I1

                XX=RRR*SIN(TETA)*COS(PHI)
                YY=RRR*SIN(TETA)*SIN(PHI)
                ZZ=RRR*COS(TETA)

                CALL INDNW3(XM,YM,ZM,NX,NY,NZ,      &
                        HX,HY,HZ,XX,YY,ZZ,&
                        J1,J2,J3,J4)

                CALL REGN(XM,YM,ZM,NX,NY,NZ,XX,YY,ZZ,EI)

                IF(EI .GT. 0) THEN
                   DR(J4,KK)=GM(I321,KK)
                ENDIF
6            CONTINUE
!          write(*,*) 'AAAAAAAAAAAAAAAAAAAAAAAA'
!          write(*,*) DR
10        CONTINUE
          RETURN
          END SUBROUTINE POLAR_DECART3D_22
!------------------------------------------
      SUBROUTINE SPRINT_3D(DR,XM,YM,ZM,NX,NY,NZ,IND,KSEC,NTEN,KK,KAN)
! print of some section  

        IMPLICIT NONE 
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single 

        INTEGER ::IND,KSEC,NX,NY,NZ,NTEN,KK,KAN
        REAL(wp) :: XM(NX),YM(NY),ZM(NZ)
        REAL(wp) :: DR(NX*NY*NZ,NTEN*NTEN)

        REAL(wp), ALLOCATABLE :: DR1(:) !,DR2(:)
        INTEGER :: I1,I2,I3,I21,I321

        ALLOCATE(DR1(NX*NY))  !,DR2(NX*NY*NZ))

!  IND=-1  -- section XZ
!  IND= 0  -- section XY
!  IND= 1  -- section YZ

!        write(*,*) 'AAAAAAAAAAAAA',KK
!         DO 3 I3=1,NZ
!         DO 3 I2=1,NY
!         DO 3 I1=1,NX
!            I21=(I2-1)*NX+I1
!            I321=(I3-1)*NX*NY+I21
!            DR2(I321)=DR(I321,KK)
!3        CONTINUE   

      IF(IND) 2,4,6

 2    CONTINUE
      DO 10 I3=1,NZ
      DO 10 I2=KSEC,KSEC
      DO 10 I1=1,NX
         I21=(I3-1)*NX+I1
         I321=(I3-1)*NY*NX+(I2-1)*NX+I1
         DR1(I21)=DR(I321,KK)  !DR2(I321)
 10   CONTINUE
!      CALL GNUFORM_2D (NX,NZ,DR1,KAN)
      CALL GNUFORM_M (XM,ZM,NX,NZ,DR1,KAN)

      GOTO 30

 4    CONTINUE
      DO 12 I3=KSEC,KSEC
      DO 12 I2=1,NY
      DO 12 I1=1,NX
         I21=(I2-1)*NX+I1
         I321=(I3-1)*NY*NX+(I2-1)*NX+I1
         DR1(I21)=DR(I321,KK)  !DR2(I321)
 12   CONTINUE
!      CALL GNUFORM_2D (NX,NY,DR1,KAN)
      CALL GNUFORM_M (XM,YM,NX,NY,DR1,KAN)
      GOTO 30

 6    CONTINUE
      DO 14 I3=1,NZ
      DO 14 I2=1,NY
      DO 14 I1=KSEC,KSEC
         I21=(I3-1)*NY+I2
         I321=(I3-1)*NY*NX+(I2-1)*NX+I1
         DR1(I21)=DR(I321,KK)  !DR2(I321)
 14   CONTINUE
!      CALL GNUFORM_2D (NY,NZ,DR1,KAN)
      CALL GNUFORM_M (YM,ZM,NY,NZ,DR1,KAN)
 30   CONTINUE
!      CLOSE(KAN)
      DEALLOCATE(DR1)  !,DR2)
      RETURN
      END  SUBROUTINE SPRINT_3D
!------------------------------------------
      SUBROUTINE GNUFORMVEC_X (YM,ZM,NX,NY,NZ,KNX,GC,NTEN,KAN1,KAN2)
          IMPLICIT NONE
          INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
          INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
!          REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

          INTEGER :: NX,NY,NZ,KNX,NTEN,KAN1,KAN2
          REAL(wp) :: XM(NX),YM(NY),ZM(NZ)
          REAL(wp) :: GC(NX*NY*NZ,NTEN*NTEN)

!      DIMENSION XM(*),ZM(*)
!      REAL GC(NX*NY*NZ,3)
!local 
          INTEGER :: K1,K2,I,I3,I2,I1,I123
          REAL(wp) :: ZZ,YY,XX,AA,BB,PHZ
          REAL(wp) :: PHAZE(NX)
      
          K1=1
          K2=2
          DO 10 I3=1,NZ
             ZZ=ZM(I3)
          DO 10 I2=1,NY
             YY=YM(I2)
          DO 10 I1=KNX,KNX    !1,NX
!           XX=XM(I1)
             I123= (I3-1)*NY*NX+(I2-1)*NX+I1

             WRITE(KAN1,100) YY,ZZ,GC(I123,K1),GC(I123,K2)
             WRITE(KAN1,100)
10        CONTINUE 

          DO 14 I3=1,NZ
!             ZZ=ZM(I3)

          DO 12 I2=20,20  !1,NY
!             YY=YM(I2)
          DO 12 I1=KNX,KNX        !1,NX
             I123= (I3-1)*NY*NX+(I2-1)*NX+I1

             AA=GC(I123,K1)
             BB=GC(I123,K2)

             IF(AA .NE. 0.0 .AND. BB .NE. 0.0) THEN  
               PHZ=ATAN2(BB,AA)
            ELSE
               PHZ=0.
            ENDIF
12       CONTINUE    
            PHAZE(I3)=PHZ
14       CONTINUE 
            WRITE(KAN2,110) (PHAZE(I),I=1,NZ)
            WRITE(KAN2,*)
100   FORMAT(4E15.3)
110   FORMAT(E10.3)
      RETURN
      END SUBROUTINE GNUFORMVEC_X
!---------------------------------------------------------
      SUBROUTINE GNUFORMVEC_Y (XM,ZM,NX,NY,NZ,KNY,GC,NTEN,KAN1,KAN2)
          IMPLICIT NONE
          INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
          INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
!          REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

          INTEGER :: NX,NY,NZ,KNY,NTEN,KAN1,KAN2
          REAL(wp) :: XM(NX),YM(NY),ZM(NZ)
          REAL(wp) :: GC(NX*NY*NZ,NTEN*NTEN)

!      DIMENSION XM(*),ZM(*)
!      REAL GC(NX*NY*NZ,3)
!local 
          INTEGER :: K1,K2,I,I3,I2,I1,I123
          REAL(wp) :: ZZ,YY,XX,AA,BB,PHZ
          REAL(wp) :: PHAZE(NX)

          K1=1
          K2=2
          DO 10 I3=1,NZ
             ZZ=ZM(I3)
          DO 10 I2=KNY,KNY   !1,NY
!         YY=YM(I2)
          DO 10 I1=1,NX
             XX=XM(I1)
             I123= (I3-1)*NY*NX+(I2-1)*NX+I1

             WRITE(KAN1,100) XX,ZZ,GC(I123,K1),GC(I123,K2)
             WRITE(KAN1,*)
10        CONTINUE 

          DO 14 I3=1,NZ
 !            ZZ=ZM(I3)

          DO 12 I1=20,20  !1,NX
!             XX=XM(I1)
          DO 12 I2=KNY,KNY        !1,NY
             I123= (I3-1)*NY*NX+(I2-1)*NX+I1
             AA=GC(I123,K1)
             BB=GC(I123,K2)
!                   PHAZE(I1)=ATAN2(BB,AA)
             IF(AA .NE. 0.0 .AND. BB .NE. 0.0) THEN  
               PHZ=ATAN2(BB,AA)
            ELSE
               PHZ=0.
            ENDIF
12       CONTINUE    
            PHAZE(I3)=PHZ
14       CONTINUE 
            WRITE(KAN2,110) (PHAZE(I),I=1,NZ)
            WRITE(KAN2,*)

100   FORMAT(4E15.3)
110   FORMAT(E10.3)
      RETURN
      END SUBROUTINE GNUFORMVEC_Y
!----------------------------------------------------------

      SUBROUTINE GNUFORMVEC_Z(XM,YM,NX,NY,NZ,KNZ,GC,NTEN,KAN1,KAN2)
          IMPLICIT NONE
          INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
          INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
!          REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

          INTEGER :: NX,NY,NZ,KNZ,NTEN,KAN1,KAN2
          REAL(wp) :: XM(NX),YM(NY),ZM(NZ)
          REAL(wp) :: GC(NX*NY*NZ,NTEN*NTEN)
!      DIMENSION XM(*),YM(*)
!      REAL GC(NX*NY*NZ,3)
!local 
          INTEGER :: K1,K2,I,I3,I2,I1,I123
          REAL(wp) :: ZZ,YY,XX,AA,BB,PHZ
          REAL(wp) :: PHAZE(NY)

          K1=1
          K2=2
          DO 10 I3=KNZ,KNZ    !1,NZ
          DO 10 I2=1,NY
             YY=YM(I2)
          DO 10 I1=1,NX
             XX=XM(I1)
             I123= (I3-1)*NY*NX+(I2-1)*NX+I1

             WRITE(KAN1,100) XX,YY,GC(I123,K1),GC(I123,K2)
             WRITE(KAN1,*)
10       CONTINUE 

         DO 14 I2=1,NY
!            YY=YM(I2)

         DO 12 I3=KNZ,KNZ
         DO 12 I1=20,20 !1,NX
!            XX=XM(I1)

            I123= (I3-1)*NY*NX+(I2-1)*NX+I1
            AA=GC(I123,K1)
            BB=GC(I123,K2)
!            write(*,*) 'AA_BB', AA,BB

            IF(AA .NE. 0.0 .AND. BB .NE. 0.0) THEN  
               PHZ=ATAN2(BB,AA)
            ELSE
               PHZ=0.
            ENDIF
12       CONTINUE    
            PHAZE(I2)=PHZ
14       CONTINUE 
            WRITE(KAN2,110) (PHAZE(I),I=1,NY)
            WRITE(KAN2,*)

100      FORMAT(4E15.3)
110      FORMAT(E10.3)
      RETURN
      END SUBROUTINE GNUFORMVEC_Z
!--------------------------------------------------

      SUBROUTINE POLAR_DECART(GM,XXM,YYM,ZZM,NX,NY,NZ,RRM,TETAM,PHIM, &
                              NRR,NTETA,NPHI,NTEN,DR)
!
!        USE PARAMTR
          IMPLICIT NONE
          INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
          INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
          REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

          INTEGER :: NX,NY,NZ,NRR,NTETA,NPHI,NTEN
          REAL(wp) :: GM(NRR*NTETA*NPHI,NTEN*NTEN)
          REAL(wp) :: XXM(NX),YYM(NY),ZZM(NZ)
          REAL(wp) :: RRM(NRR),TETAM(NTETA),PHIM(NPHI)
          REAL(wp) :: DR(NX*NY*NZ,NTEN*NTEN) 
! local
          REAL(wp) :: GMM(NRR*NTETA*NPHI)
          REAL(wp) :: XX,YY,ZZ
          REAL(wp) :: RRR,TETA,PHI  
          INTEGER :: I,I1,I2,I3,I21,I321,EI,IERR
          INTEGER :: K1,K2,KK
!functions
          REAL(wp) :: FUNL3,Burke_3D  

!          write(*,*) 'aaaaaaaaa'
!          write(*,*) NTEN


      DO 10 K2=1,NTEN
      DO 10 K1=1,NTEN
         KK=(K2-1)*NTEN+K1

         DO 2 I=1,NRR*NTETA*NPHI
            GMM(I)=0.
2        CONTINUE   

         DO 4 I3=1,NRR
         DO 4 I2=1,NTETA
         DO 4 I1=1,NPHI

            I21=(I2-1)*NPHI+I1
            I321=(I3-1)*NTETA*NPHI+I21
            GMM(I321)=GM(I321,KK)
4        CONTINUE   

         DO 6 I3=1,NZ  
            ZZ=ZZM(I3)
         DO 6 I2=1,NY
            YY=YYM(I2)
         DO 6 I1=1,NX
            XX=XXM(I1)

            I321=(I3-1)*NX*NY+(I2-1)*NX +I1

            CALL POLAR(XX,YY,ZZ,TETA,PHI)
!            CALL POLAR2(XX,YY,ZZ,TETA,PHI)

            RRR=SQRT(XX**2+YY**2+ZZ**2)
            CALL REGN(PHIM,TETAM,RRM,NPHI,NTETA,NRR,PHI,TETA,RRR,EI)
            IF(EI .GT. 0) THEN
            DR(I321,KK)=FUNL3(PHIM,TETAM,RRM,GMM,NPHI,NTETA,NRR,PHI,TETA,RRR,IERR)
!            DR(I321,KK)= Burke_3D(GMM,PHIM,TETAM,RRM,NPHI,NTETA,NRR,PHI,TETA,RRR)
            ENDIF
            IF(RRR .GE. RRM(NRR)) THEN
               DR(I321,KK)=0.0
            ENDIF

6        CONTINUE
10    CONTINUE    
      RETURN
      END SUBROUTINE POLAR_DECART
!---------------------------------------------------

      FUNCTION FUNL1(XM,GM,KK,N1,X1,IERR)
!-       1D line interpolation

          IMPLICIT NONE
          INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
          INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 

          INTEGER :: KK,N1,IERR
          REAL(wp) :: X1,XM(N1), GM(N1)
! local
          INTEGER :: K1,I1,I11,J1,J2
          REAL(wp) :: A17,R1,HX,HX1,FY1
          CHARACTER*8 ITEX
! function 
          REAL(wp) :: FUNL1

          ITEX='FUNL1'
          A17=9.E17
!
      IERR=0
      R1=(XM(N1)-XM(1))*1.E-5
      IF(X1-XM(1)+R1)  2,10,10
   2  IERR=1
      CALL ERRPRN(ITEX,IERR)
      RETURN
  10  IF(XM(N1)+R1-X1) 4,17,17
   4  IERR=2
      CALL ERRPRN(ITEX,IERR)
      RETURN
  17  CONTINUE
!
      K1=N1-1
      DO 21 I1=1,K1
      I11=I1
      IF(X1-XM(I1+1)) 27,21,21
  21  CONTINUE
!
  27  CONTINUE
!
      IF(GM(I11)+A17) 28,29,29
  28  IERR=5
      CALL ERRPRN(ITEX,IERR)
!
  29  IF(GM(I11+1)+A17) 30,42,42
  30  IERR=6
      CALL ERRPRN(ITEX,IERR)
  42  CONTINUE
!
      J1=(KK-1)*N1+I11
      J2=J1+1
      HX1=X1-XM(I11)
!     HX2=XM(I11+1)-X1
      HX=XM(I11+1)-XM(I11)
!
      FY1=GM(J1)+(GM(J2)-GM(J1))/HX*HX1
!
      FUNL1=FY1
      RETURN
      END FUNCTION FUNL1
!------------------------------------------------
      SUBROUTINE Tensor_Product(AM,IA,JA,BM,IB,JB,CM)
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
            CM(IJ)=AM(IJ1)*BM(IJ2)   
 4          CONTINUE   
 8       CONTINUE   

!      WRITE(*,*) 'Tensor product'
!      DO 10 I=1,IA*IB
!         WRITE(*,*) (CM((I-1)*JA*JB+J),J=1,JA*JB)
!         WRITE(*,*) "-----------------------------------"
! 10   CONTINUE
      RETURN
      END SUBROUTINE Tensor_Product  
!------------------------------------------------------

      FUNCTION ACOBI_SP(NN,AA,BB,ZZ)
! procedure calculate Yacobi polinome at point  Z
! by summation Hiper-Geometrical series
! look Beitmen p.172
! Z=COS(BET),BET- an angle of rotation over new X axis
! look Nikiforov p.122.

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        
        INTEGER :: NN
        REAL(wp) :: AA,BB,ZZ
!local
!        INTEGER :: I
!        REAL(wp) :: RN,RI,P,S,S1,X
!functions
        REAL(sp) :: ACOBI_SP
!----------------------------
        REAL(sp) :: A,B,Z
        INTEGER :: N,I
        REAL(sp) :: RN,RI,P,S,S1,X

        N=REAL(NN)
        A=REAL(AA)
        B=REAL(BB)
        Z=REAL(ZZ)


      RN=N
      P=1.
      DO 1 I=1,N
         RI=I
         P=P*(RN+A+1.-RI)/RI
1     CONTINUE
      S1=1.
      S=1.
      X=(1.-Z)/2.
      DO 2 I=1,N
         RI=I
         S1=S1*(-RN+RI-1.)*(RN+A+B+RI)/(A+RI)*X/RI
         S=S+S1
2     CONTINUE
      ACOBI_SP=P*S
      RETURN
      END FUNCTION ACOBI_SP
!-----------------------------------------------
      FUNCTION GNLMP_SP(N,L,M,PP1)
!  Вычисление функции U(lmn)  формула (20)
! статья Мат. Мод. 2005
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        INTEGER :: N,L,M
        REAL(wp) :: PP1
!local 
        INTEGER :: KK,LM2 
        REAL(sp) :: GGC,PP,PPP,PP2,PP32,ZZ,A,B,ACOB
!functions
        REAL(wp) :: GNLMP,GCOF,ACOBI
        REAL(sp) :: ACOBI_SP,GNLMP_SP
        
        PP=REAL(PP1)
        KK=MOD(L+M,2)

        IF(KK .EQ. 0) THEN 
         GGC=REAL(GCOF(N,L,M))   
         PPP=PP**M
         PP2=1-PP**2
         IF(PP2 .GT. 0.) THEN 
            PP32=PP2*SQRT(PP2)
         ELSE
            PP32=0.
         ENDIF   
         ZZ=1.-2.*PP**2
         LM2=(L-M)/2.
         A=M
         B=3./2.
!         ACOB=ACOBI(N+LM2,A,B,ZZ)
         ACOB=ACOBI_SP(N+LM2,A,B,ZZ)

         write(*,*) 'GNLMP', ACOB

         GNLMP_SP=GGC*PPP*PP32*ACOB

      ELSE
         GNLMP_SP=0.
      ENDIF   
      RETURN
      END FUNCTION GNLMP_SP
!-----------------------------------------------------
      SUBROUTINE INT_GNLMP(PPM,UPPM,NPP,N1,L1,K1,N2,L2,K2, &
                           HUPP,GG12)

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 

        INTEGER :: NPP,N1,L1,K1,N2,L2,K2
        REAL(wp) :: PPM(NPP),UPPM(NPP)
        REAL(wp) :: HUPP,GG12
!local
        REAL(wp) :: UG1M(NPP),UG2M(NPP)
        REAL(wp) :: GG1M(NPP),GG2M(NPP)
        INTEGER :: I1,IERR
        REAL(wp) :: PP,UPP,S,WW1
!functions
        REAL(wp) :: FUNL1,GNLMP


!      DIMENSION PPM(NPP),UPPM(NPP)
!      DIMENSION UG1M(NPP),UG2M(NPP)
!      DIMENSION GG1M(NPP),GG2M(NPP)

      DO 4 I1=1,NPP
         UPP=UPPM(I1)
!         PP=PPM(I1)
!         GG1M(I1)=GNLMP(N1,L1,K1,PP)
!         GG2M(I1)=GNLMP(N2,L2,K2,PP)
         UG1M(I1)=GNLMP(N1,L1,K1,UPP)
         UG2M(I1)=GNLMP(N2,L2,K2,UPP)
 4    CONTINUE

!      CALL GNUFORM (NPP,1,1,GG1M,20) 

! цикл 6 нужен если PPM - неравномерная сетка
!      DO 6 I1=1,NPP
!         UPP=UPPM(I1)
!!         UG1M(I1)= FUNL1(UPPM,GG1M,1,NPP,UPP,IERR)
!!         UG2M(I1)= FUNL1(UPPM,GG2M,1,NPP,UPP,IERR)

!         UG1M(I1)= FUNL1(PPM,GG1M,1,NPP,UPP,IERR)
!         UG2M(I1)= FUNL1(PPM,GG2M,1,NPP,UPP,IERR)
! 6    CONTINUE   

!      CALL GNUFORM (NPP,1,1,UG1M,22)

      S=0.
      DO 8 I1=1,NPP
         WW1=1.0
         IF(I1 .EQ.1 .OR. I1 .EQ. NPP) WW1=0.5
         S=S+UG1M(I1)*UG2M(I1)*WW1
 8    CONTINUE   
      GG12= S*HUPP
      RETURN
      END SUBROUTINE INT_GNLMP
!---------------------------------------------------------
      FUNCTION ACOBI(N,A,B,Z)
! procedure calculate Yacobi polinome at point  Z
! by summation Hiper-Geometrical series
! look Beitmen p.172
! Z=COS(BET),BET- an angle of rotation over new X axis
! look Nikiforov p.122.

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        
        INTEGER :: N
        REAL(wp) :: A,B,Z
!local
        INTEGER :: I
        REAL(wp) :: RN,RI,P,S,S1,X
!functions
        REAL(wp) :: ACOBI

!        REAL(wp) :: A,B,Z
!      INTEGER RN,RI

      RN=N
      P=1.
      DO 1 I=1,N
         RI=I
         P=P*(RN+A+1.-RI)/RI
1     CONTINUE
      S1=1.
      S=1.
      X=(1.-Z)/2.
      DO 2 I=1,N
         RI=I
         S1=S1*(-RN+RI-1.)*(RN+A+B+RI)/(A+RI)*X/RI
         S=S+S1
2     CONTINUE
      ACOBI=P*S
      RETURN
      END FUNCTION ACOBI
!-----------------------------------------------
       SUBROUTINE INT_Proj_GNLM_P(PROJ,PPM,NPP,N1,L1,M1,     &
                 I1,I2,I3,NJ,ALPHA,BETA,GAMMA,PAR,YYRE,YYIM)
!   integration over p variable of 
!   g(p,al,bet,gam)*g(n,l,m)

        USE PARAM_INPUT

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

       
        INTEGER :: NPP,N1,L1,M1
        INTEGER :: I1,I2,I3,NJ
        REAL(wp) :: ALPHA,BETA,GAMMA
        REAL(wp) :: YYRE,YYIM
        REAL(wp) :: PPM(NPP)
        REAL(wp) :: PROJ(NPP*NJ)
        REAL(wp) :: PAR(50)
!local
        REAL(wp) :: PROJ3(NPP) !,PROJ2(NPP),
        REAL(wp) :: PP
        REAL(wp) :: WW1,WWLMRE,WWLMIM
        REAL(wp) :: S1,S2,HUPP
        INTEGER :: NALPHA,NBETA,NGAMMA
        INTEGER :: I4,I321,I4321
!functions
!        REAL(wp) :: FUNL1
        REAL(wp) :: WWC,VVS

      NALPHA=PARAM(14)
      NBETA =PARAM(15)
      NGAMMA=PARAM(16)
      HUPP  =PAR(35)

      DO 6 I4=1,NPP
         I321=(I3-1)*NBETA*NALPHA+(I2-1)*NALPHA+I1
         I4321=(I4-1)*NGAMMA*NBETA*NALPHA + I321
         PROJ3(I4)=PROJ(I4321)
 6    CONTINUE   

!       в векторной версии массив PPM -неравномнрная сетка
!       массив UPPM - равномерная
!       проекции задавались на неравномерной сетке 
!       цикл 8 делает пересчёт не равномерную сетку

!      DO 8 I4=1,NPP
!         UPP=UPPM(I4)
!         PROJ3(I4)= FUNL1(PPM,PROJ2,1,NPP,UPP,IERR)
! 8    CONTINUE   

      S1=0.
      S2=0.
      DO 10 I4=1,NPP
         PP=PPM(I4)

         WWLMRE= WWC(N1,L1,M1,ALPHA,BETA,GAMMA,PP)
         WWLMIM= VVS(N1,L1,M1,ALPHA,BETA,GAMMA,PP)
!         CALL WWCS(N1,L1,M1,ALPHA,BETA,GAMMA,UPP,WWLMRE,WWLMIM)

         WW1=1.0
         IF(I4 .EQ.1 .OR. I4 .EQ. NPP) WW1=0.5
         S1=S1 + WWLMRE * PROJ3(I4)*WW1
         S2=S2 + WWLMIM * PROJ3(I4)*WW1

 10   CONTINUE   
      YYRE= S1*HUPP
      YYIM= S2*HUPP

!      write(*,*) '2222222222',YYRE,YYIM
      RETURN
      END SUBROUTINE INT_Proj_GNLM_P
!------------------------------------------------------------

      FUNCTION INT_GNLMP_GNLMP(N,L,M,N1,L1,M1,PPM,NPP)
! интегрирование произведения двух функций U(lmn,p)*U(l'm'n',p) по переменной P
!  для процедуры Matrix_WW_ 
        USE PARAMTR
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        INTEGER :: N,L,M, N1,L1,M1,NPP
        REAL(wp) :: PPM(NPP)
!        REAL(wp) :: PAR(50)
!local
        REAL(wp) :: PP,S,WW1,HPP
        INTEGER :: I
!functions
        REAL(wp) :: GNLMP, INT_GNLMP_GNLMP

        HPP  =PAR(35)
        S=0.
        DO 4 I=1,NPP
           PP=PPM(I)
           WW1=1.0
           IF(I .EQ.1 .OR. I .EQ. NPP) WW1=0.5
           S=S+GNLMP(N,L,M,PP)*GNLMP(N1,L1,M1,PP)*WW1
4       CONTINUE

           INT_GNLMP_GNLMP=S*HPP
        RETURN
        END FUNCTION INT_GNLMP_GNLMP
!----------------------------------------------------
      SUBROUTINE WW_WW_TEN(L1,M1,N1,L2,M2,N2, &
              PPM,UPPM,NPP,HUPP, &
              NBETA,NGAMMA,NALPHA,                  &
              ALPHAM,BETAM,GAMMAM,PAR,VLM,NI,NJ,NTEN)

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 

        INTEGER :: L1,M1,N1,L2,M2,N2,NPP
        INTEGER :: NALPHA,NBETA,NGAMMA
        INTEGER :: NI,NJ,NTEN
        REAL(wp) :: PPM(NPP),UPPM(NPP)
        REAL(wp) :: ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
!        INTEGER :: NPAR(100)  !!!!!
        REAL(wp) :: PAR(50)   !!!!!!!
        REAL(wp) :: VLM(4)
        REAL(wp) :: HUPP,PP
!local
        INTEGER :: I1,K1,I2,K2
        REAL(wp) :: S1,S2,S3,S4,GG12
        REAL(wp) :: VLMINT(4)
!functions
        REAL(wp) :: INT_GNLMP_GNLMP


!      DIMENSION ALPHAM(*),BETAM(*),GAMMAM(*)
!      DIMENSION PPM(*),UPPM(*)
!      DIMENSION NPAR(*),PAR(*)
!      DIMENSION VLMINT(4)      
!      DIMENSION VLM(4)      

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

      GG12= INT_GNLMP_GNLMP(N1,L1,K1,N2,L2,K2,PPM,NPP)

!      CALL INT_VLMM_TEN(L1,M1,K1,L2,M2,K2,ALPHAM,BETAM,GAMMAM, &
!        NBETA,NGAMMA,NALPHA,PAR,VLMINT,NI,NJ,NTEN,PP)

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
      END SUBROUTINE WW_WW_TEN
!--------------------------------------------------------

      FUNCTION INT_WWC_WWC_ANGLES(N,L,M,N1,L1,M1,ALPHAM,BETAM,GAMMAM,PPM, &
                                  NALPHA,NBETA,NGAMMA,NPP,NI,NJ)
! интегрирование произведения двух функций WWC(...)*WWC(...) по переменной Al,Bet,Gam
! 
        USE PARAMTR
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        INTEGER :: N,L,M, N1,L1,M1
        INTEGER :: NPP,NALPHA,NBETA,NGAMMA,NI,NJ
        REAL(wp) :: PPM(NPP)
        REAL(wp) :: ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
!        REAL(wp) :: PAR(50)
!local
        REAL(wp) :: PP,S,WW1,WW2,WW3
        REAL(wp) :: HPP,HALPHA,HBETA,HGAMMA
        REAL(wp) :: ALPHA,BETA,GAMMA
        REAL(wp) :: S1,S2,S3
        REAL(wp) :: WUnit(3)
        INTEGER :: I1,I2,I3,I23
!functions
        REAL(wp) :: INT_WWC_WWC, INT_WWC_WWC_ANGLES,INT_WWC_WWC_P

        HALPHA=PAR(38)
        HBETA =PAR(39)
        HGAMMA=PAR(40)
        HPP  =PAR(35)

        S3=0.
        DO 14 I3=1,NGAMMA
           GAMMA=GAMMAM(I3)

        S2=0.   
        DO 12 I2=1,NBETA
           BETA=BETAM(I2)

           WUnit(1)=-SIN(BETA)*COS(GAMMA)
           WUnit(2)= SIN(BETA)*SIN(GAMMA)
           WUnit(3)= COS(BETA)

! интегрирование по alpha отлично от нуля только при M=-M1
! справедливо для векторной томографии
           S1=0.
           DO 10 I1=1,NALPHA
              ALPHA=ALPHAM(I1)
!              WUnit(1)=COS(ALPHA)*SIN(BETA)
!              WUnit(2)=SIN(ALPHA)*SIN(BETA)
!              WUnit(3)=COS(BETA)
              WW1=1.0
              IF(I1 .EQ.1 .OR. I1 .EQ. NALPHA) WW1=0.5
              S1=S1+INT_WWC_WWC_P(N,L,M,N1,L1,M1,ALPHA,BETA,GAMMA,PPM,NPP)*WW1
10         CONTINUE
              
              S2=S1*HALPHA
              WW2=1.0
              IF(I2 .EQ.1 .OR. I2 .EQ. NBETA) WW2=0.5
              S2=S2*SIN(BETA)*WW2
12      CONTINUE
              S3=S2*HBETA
              WW3=1.0
              IF(I3 .EQ.1 .OR. I3 .EQ. NGAMMA) WW3=0.5
              S3=S3*HGAMMA*WW3
14      CONTINUE
        INT_WWC_WWC_ANGLES=S3
        RETURN
        END FUNCTION INT_WWC_WWC_ANGLES
!----------------------------------------------------

      FUNCTION INT_WWC_WWC_P(N,L,M,N1,L1,M1,ALPHA,BETA,GAMMA,PPM,NPP)
! интегрирование произведения двух функций WWC(...)*WWC(...) по переменной P
!   FUNCTION WWC(N,L,M,ALPHA,BETA,GAMMA,PP)
        USE PARAMTR
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        INTEGER :: N,L,M, N1,L1,M1,NPP
        REAL(wp) :: PPM(NPP)
        REAL(wp) :: ALPHA,BETA,GAMMA
!        REAL(wp) :: PAR(50)
!local
        REAL(wp) :: PP,S,WW1,HPP
        INTEGER :: I
!functions
        REAL(wp) :: WWC, INT_WWC_WWC_P

        HPP  =PAR(35)
        S=0.
        DO 4 I=1,NPP
           PP=PPM(I)
           WW1=1.0
           IF(I .EQ.1 .OR. I .EQ. NPP) WW1=0.5
!           S=S + WWC(N,L,M,ALPHA,BETA,GAMMA,PP)  &
!                 *WWC(N1,L1,M1,ALPHA,BETA,GAMMA,PP)*WW1
4       CONTINUE

           INT_WWC_WWC_P=S*HPP
        RETURN
        END FUNCTION INT_WWC_WWC_P
!----------------------------------------------------

      FUNCTION INT_VVS_VVS_P(N,L,M,N1,L1,M1,ALPHA,BETA,GAMMA,PPM,NPP)
! интегрирование произведения двух функций VVS(...)*VVS(...) по переменной P
!   FUNCTION VVS(N,L,M,ALPHA,BETA,GAMMA,PP)
        USE PARAMTR
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        INTEGER :: N,L,M, N1,L1,M1,NPP
        REAL(wp) :: PPM(NPP)
        REAL(wp) :: ALPHA,BETA,GAMMA
!        REAL(wp) :: PAR(50)
!local
        REAL(wp) :: PP,S,WW1,HPP
        INTEGER :: I
!functions
        REAL(wp) :: VVS, INT_VVS_VVS_P

        HPP  =PAR(35)
        S=0.
        DO 4 I=1,NPP
           PP=PPM(I)
           WW1=1.0
           IF(I .EQ.1 .OR. I .EQ. NPP) WW1=0.5
!           S=S + VVS(N,L,M,ALPHA,BETA,GAMMA,PP)  &
!                 *VVS(N1,L1,M1,ALPHA,BETA,GAMMA,PP)*WW1
4       CONTINUE

           INT_VVS_VVS_P=S*HPP
        RETURN
        END FUNCTION INT_VVS_VVS_P
!----------------------------------------------------
      FUNCTION INT_WWC_VVS_P(N,L,M,N1,L1,M1,ALPHA,BETA,GAMMA,PPM,NPP)
! интегрирование произведения двух функций WWC(...)*VVS(...) по переменной P
!   FUNCTION VVS(N,L,M,ALPHA,BETA,GAMMA,PP)
        USE PARAMTR
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        INTEGER :: N,L,M, N1,L1,M1,NPP
        REAL(wp) :: PPM(NPP)
        REAL(wp) :: ALPHA,BETA,GAMMA
!        REAL(wp) :: PAR(50)
!local
        REAL(wp) :: PP,S,WW1,HPP
        INTEGER :: I
!functions
        REAL(wp) :: WWC,VVS, INT_WWC_VVS_P

        HPP  =PAR(35)
        S=0.
        DO 4 I=1,NPP
           PP=PPM(I)
           WW1=1.0
           IF(I .EQ.1 .OR. I .EQ. NPP) WW1=0.5
!           S=S + WWC(N,L,M,ALPHA,BETA,GAMMA,PP)  &
!                 *VVS(N1,L1,M1,ALPHA,BETA,GAMMA,PP)*WW1
4       CONTINUE

           INT_WWC_VVS_P=S*HPP
        RETURN
        END FUNCTION INT_WWC_VVS_P
!----------------------------------------------------
      SUBROUTINE WWCS(N,L,M,ALPHA,BETA,GAMMA,PP,WWLMRE,WWLMIM)
! culculation  W(nlm;alha,beta,gamma) -functions
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 

        INTEGER :: N,L,M
        REAL(wp) :: ALPHA,BETA,GAMMA,PP,WWLMRE,WWLMIM
! local
        INTEGER ::I,M1
        REAL(wp) :: S1,S2,GGLM
        REAL(wp) :: VLMRE,VLMIM
! functions
        REAL(wp) :: GNLMP
        REAL(wp) :: WWC,VVS
        

!         write(*,*) 'WWCS', N,L,M,PP

        S1=0.
        S2=0.
        DO 10 I=1,L+1
           M1=I-1
!!           CALL VLMM_RE_IM(L,M,M1,ALPHA,BETA,GAMMA,VLMRE,VLMIM)
!           VLMRE=WWC(N,L,M1,ALPHA,BETA,GAMMA,PP)
!           VLMIM=VVS(N,L,M1,ALPHA,BETA,GAMMA,PP)
           GGLM= GNLMP(N,L,M1,PP)

!         write(*,*) 'GGLN',GGLM

           S1=S1+GGLM*VLMRE
           S2=S2+GGLM*VLMIM

!           S1=S1+VLMRE  !*GGLM
!           S2=S2+VLMIM  !*GGLM

10      CONTINUE
      WWLMRE=S1
      WWLMIM=S2

 !     write(*,*) 'WWCS',WWLMRE,WWLMIM
      RETURN
      END SUBROUTINE WWCS
!------------------------------------------------
