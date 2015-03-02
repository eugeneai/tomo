
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
        INTEGER :: NNMAX,NTEN,NTEN2
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
        INTEGER :: JJ,KJ

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

                  GMRE(NRR*NTETA*NPHI,NTEN*NTEN),            &
                  DGMRE(NX*NY*NZ,NTEN*NTEN),      &

                  DR(NX*NY*NZ),DI(NX*NY*NZ),   &
                  ROOTB(SMAXIMUM,LMAXIMUM),    &
                  DGMIM(NX*NY*NZ,NTEN),           &
                  GMIM(NX*NY*NZ,NTEN),             &
                  DR1(NX*NY*NZ),DI1(NX*NY*NZ),   &
                   ACOEF((LMAX+1)*SMAX),          &
                  RPART2(2*(LMAX+1)**2*SMAX),          &
!                  RPARTRE((LMAX+1)**2*SMAX),          &
!                  RPARTIM((LMAX+1)**2*SMAX),          &
!                  RPART((LMAX+1)**2*SMAX),          &
                  ACOFRE((LMAX+1)**2*SMAX),          &
                  ACOFIM((LMAX+1)**2*SMAX),          &
                  ACOF2(2*(LMAX+1)**2*SMAX),          &
                  ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA), &
!                  DR2(NX*NY*NZ),DI2(NX*NY*NZ),   &
                  DREX((LMAX+1)**2*SMAX**2),   &
                  DRIX((LMAX+1)*SMAX),   &
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
                  BES(NB),                     &
!                  SVZKARE(NITER),              & 
!                  SVZKAIM(NITER),              & 
!                  PRFURRE(NUU*NVV*(LMAX+1)**2),   &
!                  PRFURIM(NUU*NVV*(LMAX+1)**2),   &
!                  PRFURC(NUU*NVV*(LMAX+1)**2),   &

!                  FURPRJ3DRE(NUU*NVV*NJ),   &
!                  FURPRJ3DIM(NUU*NVV*NJ),   &
                  SMATRE2(4*(LMAX+1)**2*SMAX*(LMAX+1)**2*SMAX),       &
!                  SMATRE((LMAX+1)**2*(LMAX+1)**2),       &
!                  SMATIM((LMAX+1)**2*(LMAX+1)**2),       &
!                  SMATC((LMAX+1)**2*SMAX*(LMAX+1)**2*SMAX),       &  !!!!!!!!!!

                  SMAT((LMAX+1)**2),       &
                  SMA1(LMAX+1),       &
                  SMA2(2*(LMAX+1)**2),       &
                  PMATC((LMAX+1)**2*NPP),       &
                  PMATRE((LMAX+1)**2*NPP),       &
                  PMATIM((LMAX+1)**2*NPP),       &
                  AM(4*L3*(L3+2),4*L3*(L3+2)),                    &
                  BMODRE(NPP*NTETA*NPHI),            &
                  BMODIM(NPP*NTETA*NPHI),            &     
                  SMODRE(NTETA*NPHI),                &
                  SMODIM(NTETA*NPHI),                 &
                  BMODC1(NPP*NPHI),            &
                  BMODRE1(NPP*NPHI),                 &
                  BMODIM1(NPP*NPHI),                 &
                  BMODC(NX*NY*NZ,3),                 &
                  BMOD(NPP),BMOD1(NPP), &
                  MLMM(3),MLMM1(3),  &
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

 !   computation grids in a fourier space

!      CALL NUSET_3D(NX,NY,NZ,HX,HY,HZ,XNU,YNU,ZNU,RNU, &
!                    HNUX,HNUY,HNUZ,IERR)


!---------------------------------------------------------
!      CALL SIZE_MODELS(NDAT1,NDAT2,NDAT3)
!       write(*,*) 'NDAT1=',NDAT1,'NDAT2',NDAT2,'NDAT3',NDAT3


!      CALL INPUT_MODELS(NDAT1,NDAT2,NDAT3)  
!----------------------------------

!      L1=2
!      M1=2
!      M2=-1
!      RL1=2.
!      RM1=0.
!      RM2=-2.
!      TET=PI/3.
      
!      XX=5.0
!      DAX= PLGNDR_NORM(L1,M1,0.0)
!      DAX= PLGNDR(L1,M1,0.0)
!      DRR=-SQRT(1.-XX**2)
!      DRR=-1.5*(5.0*XX**2-1.)*(SQRT(1.-XX**2))
!      DRR=(1.-XX**2)/8.
!      DAX=  BessJ(XX,M1,NB,ERR)
!      write(*,*) 'DAX,DRR 11111',DAX  !,DRR

!      TET=PI/3.
!      PHI=PI/3.
!       CYY= D_WIGNER(L1,M1,M2,PHI,TET,0.0)
!      write(*,*) 'CYY=',CYY
!      DAX= PLM1M2(L1,M1,M2,PI-TET)
!      DAX=DAX*(-1)**(L1+M1)
!      DAX= PLM1M2(L1,M1,M2,TET)
!      DRR=SQRT(3./2.)*SIN(TET)*COS(TET)
!      DRR=COS(TET)
!      DRR=(1.0-COS(TET))**2/4.
!      DRR=(1.0+COS(TET))/2.
!      DRR=0.5*SQRT(1.5)*SIN(TET)**2
!     DRR=-SQRT(1.5)*SIN(TET)*COS(TET)
!      DRR=-SIN(TET)*(1.0-COS(TET))/2.
!      write(*,*) 'DAX,DRR',DAX,DRR

!      L1=2
!      CALL INDUXIS2(L1)

!      write(*,*) '11111111111',YM
!         DO 6 I3=1,NZ  
!            ZZ=ZM(I3)
!         DO 6 I2=1,NY
!            YY=YM(I2)
!         DO 6 I1=1,NX
!            XX=XM(I1)
!            I21=(I2-1)*NX +I1
!            I321=(I3-1)*NX*NY+(I2-1)*NX +I1
!            DI1(I321)= EXP(-(XX**2+YY**2+ZZ**2)/0.1)
!6        CONTINUE   
!            KNZ=NZ/2
!         CALL SPRINT_3D(DI1,XM,YM,ZM,NX,NY,NZ,-1,KNZ,KAN(26))
!         CALL SPRINT_3D(DI1,XM,YM,ZM,NX,NY,NZ, 0,KNZ,KAN(23))
!         CALL SPRINT_3D(DI1,XM,YM,ZM,NX,NY,NZ, 1,KNZ,KAN(20))
!         CALL GNUFORM(NX,NY,DI1,KAN(20))
!         CALL GNUFORM_M(XM,YM,NX,NY,DI1,KAN(20))
!      SUBROUTINE GNUFORM (NX,NY,GC,KAN)
!         CALL SPRINT_1(DR,NX,NY,NZ,18,1,1,20)

! Proverka indeksov matrisi      
!      LL3=3

!      DO 4 I7=1,3
!         WRITE(*,*) 'I=',I7
!      DO 2 L=1,LL3
!      DO 2 M1=1,L+1
!         M=M1-1
!         IA=(I7-1)*(LL3**2+2*LL3)+ L**2+M
!         WRITE(*,*) 'IA=',IA
! 2    CONTINUE   
!      DO 3 L=1,LL3
!      DO 3 M1=1,L
!         IB=(I7-1)*(LL3**2+2*LL3)+ L**2+L+M1
!         WRITE(*,*) 'IB=',IB
! 3    CONTINUE    
! 4    CONTINUE 




!=================================================================

      CALL XXMM(LL,NNMAX,XXM,NTEN)  

      write(*,*) 'solution exact XXM:'
      write(*,117) (XXM(I,1), I=1, LN)
      write(*,117) (XXM(I,2), I=1, LN)
      write(*,117) (XXM(I,3), I=1, LN)
      write(*,117) (XXM(I,4), I=1, LN)
      write(*,117) (XXM(I,5), I=1, LN)
      write(*,117) (XXM(I,6), I=1, LN)
      write(*,117) (XXM(I,7), I=1, LN)
      write(*,117) (XXM(I,8), I=1, LN)
      write(*,117) (XXM(I,9), I=1, LN)

      CALL Vec_Line(XXM,XMVEC,LN,NTEN)

      CALL GNUFORM_0 (NTEN*NTEN*LN,XMVEC,KAN(30))   !! plot cof_exc.gnu
      write(*,*) 'solution exact XMVEC'
      write(*,117) (XMVEC(I), I=1, NTEN*NTEN*LN)

      CALL Summa_Harm_Vec(LL,NNMAX,RRM,TETM,PHIM, &
                         NRR,NTETA,NPHI,NTEN,XXM,GMRE)

      CALL POLAR_DECART(GMRE,XM,YM,ZM,NX,NY,NZ,RRM,TETM,PHIM, &
                              NRR,NTETA,NPHI,NTEN,DGMRE)

!     write(*,*) 'aaaaaaaaaa', DGMRE


      KNZ=NZ/2
      KNY=NY/3
      KNX=NX/2
      KKK=1
!  IND=-1  -- section XZ
!  IND= 0  -- section XY
!  IND= 1  -- section YZ
      IND=-1
      CALL SPRINT_3D(DGMRE,XM,YM,ZM,NX,NY,NZ,IND,KNY,NTEN,KKK,KAN(26))
      IND=0
      CALL SPRINT_3D(DGMRE,XM,YM,ZM,NX,NY,NZ,IND,KNZ,NTEN,KKK,KAN(23))
      IND=1
      CALL SPRINT_3D(DGMRE,XM,YM,ZM,NX,NY,NZ,IND,KNX,NTEN,KKK,KAN(20))

      CALL GNUFORMVEC_Z(XM,YM,NX,NY,NZ,KNZ,DGMRE,NTEN,KAN(46),KAN(58))

      CALL GNUFORMVEC_Y(XM,ZM,NX,NY,NZ,KNY,DGMRE,NTEN,KAN(47),KAN(56))

      CALL GNUFORMVEC_X(YM,ZM,NX,NY,NZ,KNX,DGMRE,NTEN,KAN(51),KAN(54))


!       GOTO 777

!       CALL ProjNum_Vec(DGMRE,PRJ3DIM,UM,VM,WM,  &
!              ALPHAM,BETAM,GAMMAM,NALPHA,NBETA,NGAMMA, &
!              XM,YM,ZM,NX,NY,NZ,NU,NV,NW,NTEN)

       CALL ProjNum_TEN(DGMRE,PRJ3DIM,UM,VM,WM,  &
              ALPHAM,BETAM,GAMMAM,NALPHA,NBETA,NGAMMA, &
              XM,YM,ZM,NX,NY,NZ,NU,NV,NW,NTEN)

       CALL Proj_Summa_W_Harm(LL,NJ,NNMAX,PPM,ALPHAM,BETAM, &
              GAMMAM,NPP,NALPHA,NBETA,NGAMMA,XXM,NTEN,PRJ3DRE)

!       write(*,*) '11111111',PRJ3DIM

       KJ=15

       DO 11 I4=1,NPP
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

          write(*,*) 'alpha_n=',ALPHA*180/PI
          write(*,*) 'beta_n=',BETA*180/PI
          write(*,*) 'gamma_n=',GAMMA*180/PI
!          CALL GNUFORM1(UM,NU,PROJIM1D,KAN(44)) 
!          CALL GNUFORM1(UM,NU,PROJRE1D,KAN(45)) 
          CALL GNUFORM (NU,1,PROJIM1D,KAN(44)) 
          CALL GNUFORM (NU,1,PROJRE1D,KAN(45)) 
!                  SUBROUTINE GNUFORM (NX,NY,GC,KAN)

!          HALPHA=PAR(38)
!          HBETA =PAR(39)
!          HGAMMA=PAR(40)

!       write(*,*) '11111111',PRJ3DIM
!          CALL RightPart_WW_VEC(PRJ3DIM,LL,NJ,NNMAX,  &
!                PPM,PPM,ALPHAM,BETAM,GAMMAM,        &
!                NPP,NALPHA,NBETA,NGAMMA,             &
!                NTEN,PAR,YMVEC)

          CALL RightPart_WW_TEN(PRJ3DIM,LL,NJ,NNMAX,  &
                PPM,PPM,ALPHAM,BETAM,GAMMAM,        &
                NPP,NALPHA,NBETA,NGAMMA,             &
                NTEN,PAR,YMVEC)

          CALL Vec_Line(YMVEC,YMVEC3,LN,NTEN)

       write(*,*) 'right part:YMVEC(1)'
       write(*,117) (YMVEC(I,1), I=1, LN)
       write(*,*) 'right part:YMVEC(2)'
       write(*,117) (YMVEC(I,2), I=1, LN)
       write(*,*) 'right part:YMVEC(3)'
       write(*,117) (YMVEC(I,3), I=1, LN)



       write(*,*) 'right part:YMVEC3'
       write(*,117) (YMVEC3(I), I=1, NTEN*NTEN*LN)

      NTEN2=NTEN*NTEN
      DO 81 NII=1,NTEN2
      DO 81 NJJ=1,NTEN2     !NII

!      CALL IALBET(NII,3,I1,I2)
!      CALL IALBET(NJJ,3,J1,J2)
!      write(*,*) '',I1,I2,J1,J2

      CALL Matrix_WW_TEN(LL,NNMAX,PPM,RRM,ALPHAM,BETAM,GAMMAM, &
                     NPP,NALPHA,NBETA,NGAMMA,AMVEC,PAR,NII,NJJ)

      CALL Sort_Matr(AMVEC,LN,AMVEC3,NTEN,NII,NJJ)

 81   CONTINUE

!      write(*,*) 'Matrix AMVEC3'
!      DO 97 I=1,NTEN*NTEN*LN
!      write(*,117) (AMVEC3(I,J), J=1, NTEN*NTEN*LN)
! 97   CONTINUE
      write(*,*) 'Matrix AMVEC3'
      DO 97 I=1,LN
      write(*,117) (AMVEC3(I,J), J=1,LN)
 97   CONTINUE
!-----------------------------------------------

! initial values for the solution
      DO 13 I=1,NTEN2*LN
         XMVEC0(I)=0.
 13   CONTINUE 

      CALL XYANG_1(AMVEC3,YMVEC3,XMVEC0,NTEN2*LN, NTEN2*LN,NITER)
!... am(n1,n2) - matrix 
!....gm(n1) - right part 
! ...fm1(n2) - initial estimation (input) and solution (output)
!-   niter -the number of iterations
      CALL Vec_Line_Inv(XMVEC0,XXM2,LN,NTEN)

      write(*,*) 'solution reconstr. XMVEC0'
      write(*,117) (XMVEC0(I), I=1, NTEN2*LN)

      CALL GNUFORM_0 (NTEN2*LN,XMVEC0,KAN(31)) 


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

      SUBROUTINE Sort_Matr(AAM,LN,DDM,NTEN,NCP1,NCP2)
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 

        INTEGER :: LN,NTEN,NCP1,NCP2
        REAL(wp) :: AAM(LN,LN),DDM(NTEN*NTEN*LN,NTEN*NTEN*LN)
!local 
        INTEGER :: I,J,II,JJ

!      DIMENSION AAM(LN,LN)
!      DIMENSION DDM(3*LN,3*LN)

!      DO 2 I=1,3*LN
!      DO 2 J=1,3*LN
!         DDM(I,J)=0.
! 2    CONTINUE   

      DO 4 J=1,LN
           JJ=(NCP2-1)*LN+J
      DO 4 I=1,LN
           II=(NCP1-1)*LN+I
         DDM(II,JJ)=AAM(I,J)
 4    CONTINUE

      RETURN
      END SUBROUTINE Sort_Matr
!--------------------------------------------------------
      SUBROUTINE Matrix_WW_TEN(LL,NNMAX,PPM,RRM,ALPHAM,BETAM,GAMMAM, &
                     NPP,NALPHA,NBETA,NGAMMA,AM,PAR,NI,NJ)


        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

        INTEGER :: LL,NNMAX,NPP,NALPHA,NBETA,NGAMMA,NI,NJ
!        INTEGER :: NPAR(100)  !!!!!
        REAL(wp) :: PAR(50)   !!!!!!!
        REAL(wp) :: PPM(NPP),ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
        REAL(wp) :: RRM(NPP)   !!! !NRR=NPP
        REAL(wp) :: AM((LL+1)*(LL+1)*NNMAX,(LL+1)*(LL+1)*NNMAX)

!local
        INTEGER :: LMAX,NI1,NI2,NJ1,NJ2
        INTEGER :: N,L,M,LL2
        INTEGER :: N1,L1,L11,L22,M1
        INTEGER :: IIC1,IIC,IIS1,IIS
        REAL(wp) :: VLM(4)
        REAL(wp) :: HALPHA,HBETA,HGAMMA,HPP

!      DIMENSION AM((LL+1)*(LL+1)*NNMAX,(LL+1)*(LL+1)*NNMAX)
!      DIMENSION PPM(NPP),ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
!      DIMENSION UPPM(NPP)
!      DIMENSION NPAR(*),PAR(*)
!      DIMENSION VLM(4)
!      PAR(35)=HPP
!       PAR(38)=HALPHA
!       PAR(39)=HBETA
!       PAR(40)=HGAMMA

      HALPHA=PAR(38)
      HBETA =PAR(39)
      HGAMMA=PAR(40)
      HPP  =PAR(35)
    
      CALL IALBET(NI,3,NI1,NI2)
      CALL IALBET(NJ,3,NJ1,NJ2)

!      write(*,*) 'Matrix',NI1,NI2,NJ1,NJ2


      LMAX=(LL+1)*(LL+1)

!---------------------------------------
!   cos * cos  (1,1)

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

!        write(*,*) 'Matrix_W: 30-20', IIC,IIC1


!      CALL WW_WW_VEC(L11,M1-1,N1,L22,M-1,N,      &
!              PPM,RMM,NPP,HPP,                &
!              NBETA,NGAMMA,NALPHA,              & 
!              ALPHAM,BETAM,GAMMAM,NPAR,PAR,VLM,NCP1,NCP2)

      CALL WW_WW_TEN(L11,M1-1,N1,L22,M-1,N,       &
              PPM,RRM,NPP,HPP,                    &
              NBETA,NGAMMA,NALPHA,                  &
              ALPHAM,BETAM,GAMMAM,PAR,VLM,NI1,NI2,NJ1,NJ2)

      AM(IIC,IIC1)=VLM(1)
 20   CONTINUE
 30   CONTINUE

!---------------------------------------------------
!   sin * sin  (2,2)

      DO 50  N1=1,NNMAX
      DO 50  L1=1,LL+1
         L11=L1-1
!        write(*,*) 'Matrix_W: 50 L1=  ',L1
      DO 50  M1=1,L11
         IIS1=(N1-1)*LMAX + (L11*L11+(M1+1)+L11)
      DO 40  N=1,NNMAX
      DO 40  L=1,LL+1
         L22=L-1
      DO 40  M=1,L22
         IIS=(N-1)*LMAX + (L22*L22+(M+1)+L22)

!      CALL WW_WW_VEC(L11,M1,N1,L22,M,N,      &
!              PPM,RMM,NPP,HPP,             &
!              NBETA,NGAMMA,NALPHA,           &
!              ALPHAM,BETAM,GAMMAM,NPAR,PAR,VLM,NCP1,NCP2)

      CALL WW_WW_TEN(L11,M1,N1,L22,M,N,       &
              PPM,RRM,NPP,HPP,                    &
              NBETA,NGAMMA,NALPHA,                  &
              ALPHAM,BETAM,GAMMAM,PAR,VLM,NI1,NI2,NJ1,NJ2)

        AM(IIS,IIS1)=VLM(4)
 40   CONTINUE
 50   CONTINUE
!-----------------------------------------------

!  cos * sin  (1,2)

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

!      CALL WW_WW_VEC(L11,M1-1,N1,L22,M,N,      &
!              PPM,RMM,NPP,HPP,               & 
!              NBETA,NGAMMA,NALPHA,             &
!              ALPHAM,BETAM,GAMMAM,NPAR,PAR,VLM,NCP1,NCP2)

      CALL WW_WW_TEN(L11,M1-1,N1,L22,M,N,       &
              PPM,RRM,NPP,HPP,                    &
              NBETA,NGAMMA,NALPHA,                  &
              ALPHAM,BETAM,GAMMAM,PAR,VLM,NI1,NI2,NJ1,NJ2)

         AM(IIS,IIC1)=VLM(3)
 60   CONTINUE
 70   CONTINUE
!------------------------------------------

!  sin * cos   (2,1)

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

!      CALL WW_WW_VEC(L11,M1,N1,L22,M-1,N,       &
!              PPM,RMM,NPP,HPP,                &
!              NBETA,NGAMMA,NALPHA,              &
!              ALPHAM,BETAM,GAMMAM,NPAR,PAR,VLM,NCP1,NCP2)

      CALL WW_WW_TEN(L11,M1,N1,L22,M-1,N,       &
              PPM,RRM,NPP,HPP,                    &
              NBETA,NGAMMA,NALPHA,                  &
              ALPHAM,BETAM,GAMMAM,PAR,VLM,NI1,NI2,NJ1,NJ2)

         AM(IIC,IIS1)=VLM(2)
 80   CONTINUE
 90   CONTINUE
 112  CONTINUE
      RETURN
      END SUBROUTINE Matrix_WW_TEN
!---------------------------------------------------

      SUBROUTINE WW_WW_TEN(L1,M1,N1,L2,M2,N2, &
              PPM,UPPM,NPP,HUPP, &
              NBETA,NGAMMA,NALPHA,                  &
              ALPHAM,BETAM,GAMMAM,PAR,VLM,NI1,NI2,NJ1,NJ2)


        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 

        INTEGER :: L1,M1,N1,L2,M2,N2,NPP
        INTEGER :: NALPHA,NBETA,NGAMMA
        INTEGER :: NI1,NI2,NJ1,NJ2
        REAL(wp) :: PPM(NPP),UPPM(NPP)
        REAL(wp) :: ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
!        INTEGER :: NPAR(100)  !!!!!
        REAL(wp) :: PAR(50)   !!!!!!!
        REAL(wp) :: VLM(4)
        REAL(wp) :: HUPP
!local
        INTEGER :: I1,K1,I2,K2
        REAL(wp) :: S1,S2,S3,S4,GG12
        REAL(wp) :: VLMINT(4)

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

      CALL INT_GNLMP(PPM,UPPM,NPP,N1,L1,K1,N2,L2,K2,  &
                     HUPP,GG12)

!      CALL INT_VLMM_VEC(L1,M1,K1,L2,M2,K2,ALPHAM,BETAM,GAMMAM, &
!             NBETA,NGAMMA,NALPHA,NPAR,PAR,VLMINT,NCP1,NCP2)


      CALL INT_VLMM_TEN(L1,M1,K1,L2,M2,K2,ALPHAM,BETAM,GAMMAM, &
        NBETA,NGAMMA,NALPHA,PAR,VLMINT,NI1,NI2,NJ1,NJ2)

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

      SUBROUTINE INT_VLMM_TEN(L1,M1,K1,L2,M2,K2,ALPHAM,BETAM,GAMMAM, &
        NBETA,NGAMMA,NALPHA,PAR,VLMINT,NI1,NI2,NJ1,NJ2)

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 

        INTEGER :: L1,M1,K1,L2,M2,K2
        INTEGER :: LL1,MM1,KK1,LL2,MM2,KK2
        INTEGER :: NALPHA,NBETA,NGAMMA
        INTEGER :: NPP,NI1,NI2,NJ1,NJ2
        REAL(wp) :: ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
        REAL(wp) :: VLMINT(4)
!        INTEGER :: NPAR(100)  !!!!!
        REAL(wp) :: PAR(50)   !!!!!!!
!local
        INTEGER :: I1,I2,I3,I23
        REAL(wp) :: S1,S2,S3,S4
        REAL(wp) :: ALPHA,BETA,GAMMA,HALPHA,HBETA,HGAMMA
        REAL(wp) :: WW1,WW2,WW3
        REAL(wp) :: VLCC,VLSC,VLCS,VLSS
        REAL(wp) :: VLMRE1,VLMIM1,VLMRE2,VLMIM2
!        REAL(wp) :: WLMRE1,WLMIM1,WLMRE2,WLMIM2
        REAL(wp) :: WUnit(3)
        REAL(wp) :: BM1(NBETA*NGAMMA,4),BM2(NBETA,4)   !!!!!!!!
!functions
        REAL(wp) :: WLMM_RE_TEN,WLMM_IM_TEN


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
         GAMMA=GAMMAM(I3)
         I23=(I2-1)*NGAMMA+I3


!         GAMMA=GAMMAM(I3)
!         WW3=1.0
!         IF(I3 .EQ.1 .OR. I3 .EQ. NGAMMA) WW3=0.5
         S1=0.
         S2=0.
         S3=0.
         S4=0.
      DO 10 I1=1,NALPHA
         ALPHA=ALPHAM(I1)

!         WUnit(1)=COS(ALPHA)*SIN(BETA)
!         WUnit(2)=SIN(ALPHA)*SIN(BETA)
!         WUnit(3)=COS(BETA)

         WUnit(1)=-SIN(BETA)*COS(GAMMA)
         WUnit(2)= SIN(BETA)*SIN(GAMMA)
         WUnit(3)= COS(BETA)

         WW1=1.0
         IF(I1 .EQ.1 .OR. I1 .EQ. NALPHA) WW1=0.5
! correspond to i,i'
         VLMRE1=WLMM_RE_TEN(L1,M1,K1,ALPHA,BETA,GAMMA,NI1,NI2)
         VLMIM1=WLMM_IM_TEN(L1,M1,K1,ALPHA,BETA,GAMMA,NI1,NI2)

! correspond to j,j'
         VLMRE2=WLMM_RE_TEN(L2,M2,K2,ALPHA,BETA,GAMMA,NJ1,NJ2)
         VLMIM2=WLMM_IM_TEN(L2,M2,K2,ALPHA,BETA,GAMMA,NJ1,NJ2)

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
         I23=(I3-1)*NGAMMA+I3
!         S1=0.
!         S2=0.
!         S3=0.
!         S4=0.
!      DO 14 I3=1,NGAMMA
!         I23=(I2-1)*NGAMMA+I3

         WW3=1.0
         IF(I3 .EQ.1 .OR. I3 .EQ. NGAMMA) WW3=0.5
         S1=S1+BM1(I23,1)*SIN(BETA)*WW3
         S2=S2+BM1(I23,2)*SIN(BETA)*WW3
         S3=S3+BM1(I23,3)*SIN(BETA)*WW3
         S4=S4+BM1(I23,4)*SIN(BETA)*WW3
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
      END SUBROUTINE INT_VLMM_TEN
!--------------------------------------------------------


      SUBROUTINE INT_VLMM_TEN_1(L1,M1,K1,L2,M2,K2,ALPHAM,BETAM,GAMMAM, &
        NBETA,NGAMMA,NALPHA,PAR,VLMINT,NI1,NI2,NJ1,NJ2)

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 

        INTEGER :: L1,M1,K1,L2,M2,K2
        INTEGER :: LL1,MM1,KK1,LL2,MM2,KK2
        INTEGER :: NALPHA,NBETA,NGAMMA
        INTEGER :: NPP,NI1,NI2,NJ1,NJ2
        REAL(wp) :: ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
        REAL(wp) :: VLMINT(4)
!        INTEGER :: NPAR(100)  !!!!!
        REAL(wp) :: PAR(50)   !!!!!!!
!local
        INTEGER :: I1,I2,I3,I23
        REAL(wp) :: S1,S2,S3,S4
        REAL(wp) :: ALPHA,BETA,GAMMA,HALPHA,HBETA,HGAMMA
        REAL(wp) :: WW1,WW2,WW3
        REAL(wp) :: VLCC,VLSC,VLCS,VLSS
        REAL(wp) :: VLMRE1,VLMIM1,VLMRE2,VLMIM2
!        REAL(wp) :: WLMRE1,WLMIM1,WLMRE2,WLMIM2
        REAL(wp) :: WUnit(3)
        REAL(wp) :: BM1(NBETA*NGAMMA,4),BM2(NBETA,4)   !!!!!!!!
!functions
        REAL(wp) :: WLMM_RE_TEN,WLMM_IM_TEN


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
      
      DO 12 I3=1,NGAMMA
         GAMMA=GAMMAM(I3)
      DO 12 I2=1,NBETA
         BETA=BETAM(I2)

!      DO 12 I3=1,NGAMMA
!         I23=(I2-1)*NGAMMA+I3
         I23=(I3-1)*NBETA+I2

!         GAMMA=GAMMAM(I3)
!         WW3=1.0
!         IF(I3 .EQ.1 .OR. I3 .EQ. NGAMMA) WW3=0.5
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
! correspond to i,i'
         VLMRE1=WLMM_RE_TEN(L1,M1,K1,ALPHA,BETA,GAMMA,NI1,NI2)
         VLMIM1=WLMM_IM_TEN(L1,M1,K1,ALPHA,BETA,GAMMA,NI1,NI2)

! correspond to j,j'
         VLMRE2=WLMM_RE_TEN(L2,M2,K2,ALPHA,BETA,GAMMA,NJ1,NJ2)
         VLMIM2=WLMM_IM_TEN(L2,M2,K2,ALPHA,BETA,GAMMA,NJ1,NJ2)

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

      DO 16 I3=1,NGAMMA
         S1=0.
         S2=0.
         S3=0.
         S4=0.
      DO 14 I2=1,NBETA
         BETA=BETAM(I2)
         I23=(I3-1)*NBETA+I2
!         S1=0.
!         S2=0.
!         S3=0.
!         S4=0.
!      DO 14 I3=1,NGAMMA
!         I23=(I2-1)*NGAMMA+I3

         WW2=1.0
         IF(I2 .EQ.1 .OR. I2 .EQ. NBETA) WW2=0.5
         S1=S1+BM1(I23,1)*SIN(BETA)*WW2
         S2=S2+BM1(I23,2)*SIN(BETA)*WW2
         S3=S3+BM1(I23,3)*SIN(BETA)*WW2
         S4=S4+BM1(I23,4)*SIN(BETA)*WW2
 14   CONTINUE

         BM2(I3,1)=S1*HBETA
         BM2(I3,2)=S2*HBETA
         BM2(I3,3)=S3*HBETA
         BM2(I3,4)=S4*HBETA
 16   CONTINUE

         VLMINT(1)=0.
         VLMINT(2)=0.
         VLMINT(3)=0.
         VLMINT(4)=0.
      DO 18 I3=1,NGAMMA
!         GAMMA=GAMMAM(I3)
         WW3=1.0
         IF(I3 .EQ.1 .OR. I3 .EQ. NGAMMA) WW3=0.5

         VLMINT(1)=VLMINT(1)+BM2(I3,1)*HGAMMA*WW3
         VLMINT(2)=VLMINT(2)+BM2(I3,2)*HGAMMA*WW3
         VLMINT(3)=VLMINT(3)+BM2(I3,3)*HGAMMA*WW3
         VLMINT(4)=VLMINT(4)+BM2(I3,4)*HGAMMA*WW3

 18   CONTINUE
      RETURN
      END SUBROUTINE INT_VLMM_TEN_1
!--------------------------------------------------------

      FUNCTION WLMM_RE_TEN(L,M,M1,ALPHA,BETA,GAMMA,NI1,NI2)

!          Computer the W harmonics
!          VLMRE - real part

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

        INTEGER :: L,M,M1,NI1,NI2 
        REAL(wp) :: ALPHA,BETA,GAMMA
! local
        INTEGER ::KK
        REAL(wp) :: DM,SIGNUM,CC1,CC2
        REAL(wp) :: PLM1,PLM2
        REAL(wp) :: VLMRE
        REAL(wp) :: WUnit(3) 
! functions
        REAL(wp) :: PLM1M2,WLMM_RE_TEN 

!         WUnit(1)=COS(ALPHA)*SIN(BETA)
!         WUnit(2)=SIN(ALPHA)*SIN(BETA)
!         WUnit(3)=COS(BETA)
         WUnit(1)=-SIN(BETA)*COS(GAMMA)
         WUnit(2)= SIN(BETA)*SIN(GAMMA)
         WUnit(3)= COS(BETA)

        DM=1.
        IF(M1 .EQ. 0) DM=0.5

!      KK=MOD((M-M1),2)
        KK=MOD(M1,2)
        SIGNUM=1.
        IF(KK .NE. 0) SIGNUM=-1.

        CC1=COS(M*ALPHA + M1*GAMMA)
        CC2=COS(M*ALPHA - M1*GAMMA)

!        SS1= SIN(M*ALPHA + M1*GAMMA)
!        SS2= SIN(M*ALPHA - M1*GAMMA)

        PLM1 = PLM1M2(L,M, M1,BETA)
        PLM2 = PLM1M2(L,M,-M1,BETA)

        VLMRE=DM*(CC1*PLM1 + SIGNUM*CC2*PLM2)
        WLMM_RE_TEN = VLMRE*WUnit(NI1)*WUnit(NI2)

      RETURN
      END FUNCTION WLMM_RE_TEN
!----------------------------------------------
      FUNCTION WLMM_IM_TEN(L,M,M1,ALPHA,BETA,GAMMA,NJ1,NJ2)

!          Computer the V harmonics
!          VLMIM - image part
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

        INTEGER :: L,M,M1,NJ1,NJ2
        REAL(wp) :: ALPHA,BETA,GAMMA
! local
        INTEGER ::KK
        REAL(wp) :: DM,SIGNUM,SS1,SS2
        REAL(wp) :: PLM1,PLM2
        REAL(wp) :: VLMIM
        REAL(wp) :: WUnit(3) 
! functions
        REAL(wp) :: PLM1M2,WLMM_IM_TEN 

!         WUnit(1)=COS(ALPHA)*SIN(BETA)
!         WUnit(2)=SIN(ALPHA)*SIN(BETA)
!         WUnit(3)=COS(BETA)

         WUnit(1)=-SIN(BETA)*COS(GAMMA)
         WUnit(2)= SIN(BETA)*SIN(GAMMA)
         WUnit(3)= COS(BETA)

        DM=1.
        IF(M1 .EQ. 0) DM=0.5

!      KK=MOD((M-M1),2)
        KK=MOD(M1,2)
        SIGNUM=1.
        IF(KK .NE. 0) SIGNUM=-1.

        SS1= SIN(M*ALPHA + M1*GAMMA)
        SS2= SIN(M*ALPHA - M1*GAMMA)

        PLM1 = PLM1M2(L,M, M1,BETA)
        PLM2 = PLM1M2(L,M,-M1,BETA)

        VLMIM=DM*(SS1*PLM1 + SIGNUM*SS2*PLM2)
        WLMM_IM_TEN=VLMIM*WUnit(NJ1)*WUnit(NJ2)
      RETURN
      END FUNCTION WLMM_IM_TEN
!----------------------------------------------

      SUBROUTINE Matrix_WW_VEC(LL,NNMAX,PPM,RMM,ALPHAM,BETAM,GAMMAM, &
                     NPP,NALPHA,NBETA,NGAMMA,NPAR,PAR,AM,NCP1,NCP2)


        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

        INTEGER :: LL,NNMAX,NPP,NALPHA,NBETA,NGAMMA,NCP1,NCP2
        INTEGER :: NPAR(100)  !!!!!
        REAL(wp) :: PAR(50)   !!!!!!!
        REAL(wp) :: PPM(NPP),ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
        REAL(wp) :: RMM(NPP)   !!! !NRR=NPP
        REAL(wp) :: AM((LL+1)*(LL+1)*NNMAX,(LL+1)*(LL+1)*NNMAX)

!local
        INTEGER :: LMAX
        INTEGER :: N,L,M,LL2
        INTEGER :: N1,L1,L11,L22,M1
        INTEGER :: IIC1,IIC,IIS1,IIS
        REAL(wp) :: VLM(4)
        REAL(wp) :: HALPHA,HBETA,HGAMMA,HPP

!      DIMENSION AM((LL+1)*(LL+1)*NNMAX,(LL+1)*(LL+1)*NNMAX)
!      DIMENSION PPM(NPP),ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
!      DIMENSION UPPM(NPP)
!      DIMENSION NPAR(*),PAR(*)
!      DIMENSION VLM(4)
!      PAR(35)=HPP
!       PAR(38)=HALPHA
!       PAR(39)=HBETA
!       PAR(40)=HGAMMA

      HALPHA=PAR(38)
      HBETA =PAR(39)
      HGAMMA=PAR(40)
      HPP  =PAR(35)
    

!      write(*,*) 'Matrix_WW',HALPHA,HBETA,HGAMMA

      LMAX=(LL+1)*(LL+1)


!   cos * cos  (1,1)

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

!        write(*,*) 'Matrix_W: 30-20', IIC,IIC1


      CALL WW_WW_VEC(L11,M1-1,N1,L22,M-1,N,      &
              PPM,RMM,NPP,HPP,                &
              NBETA,NGAMMA,NALPHA,              & 
              ALPHAM,BETAM,GAMMAM,NPAR,PAR,VLM,NCP1,NCP2)

      AM(IIC,IIC1)=VLM(1)
 20   CONTINUE
 30   CONTINUE
!---------------------------------------------------

!   sin * sin  (2,2)

      DO 50  N1=1,NNMAX
      DO 50  L1=1,LL+1
         L11=L1-1
!        write(*,*) 'Matrix_W: 50 L1=  ',L1
      DO 50  M1=1,L11
         IIS1=(N1-1)*LMAX + (L11*L11+(M1+1)+L11)
      DO 40  N=1,NNMAX
      DO 40  L=1,LL+1
         L22=L-1
      DO 40  M=1,L22
         IIS=(N-1)*LMAX + (L22*L22+(M+1)+L22)

      CALL WW_WW_VEC(L11,M1,N1,L22,M,N,      &
              PPM,RMM,NPP,HPP,             &
              NBETA,NGAMMA,NALPHA,           &
              ALPHAM,BETAM,GAMMAM,NPAR,PAR,VLM,NCP1,NCP2)

        AM(IIS,IIS1)=VLM(4)
 40   CONTINUE
 50   CONTINUE
!-----------------------------------------------

!  cos * sin  (1,2)

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

      CALL WW_WW_VEC(L11,M1-1,N1,L22,M,N,      &
              PPM,RMM,NPP,HPP,               & 
              NBETA,NGAMMA,NALPHA,             &
              ALPHAM,BETAM,GAMMAM,NPAR,PAR,VLM,NCP1,NCP2)

         AM(IIS,IIC1)=VLM(3)
 60   CONTINUE
 70   CONTINUE
!------------------------------------------

!  sin * cos   (2,1)

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

      CALL WW_WW_VEC(L11,M1,N1,L22,M-1,N,       &
              PPM,RMM,NPP,HPP,                &
              NBETA,NGAMMA,NALPHA,              &
              ALPHAM,BETAM,GAMMAM,NPAR,PAR,VLM,NCP1,NCP2)

         AM(IIC,IIS1)=VLM(2)
 80   CONTINUE
 90   CONTINUE
 112  CONTINUE
      RETURN
      END SUBROUTINE Matrix_WW_VEC
!---------------------------------------------------

      SUBROUTINE WW_WW_VEC(L1,M1,N1,L2,M2,N2,       &
              PPM,UPPM,NPP,HUPP,                    &
              NBETA,NGAMMA,NALPHA,                  &
              ALPHAM,BETAM,GAMMAM,NPAR,PAR,VLM,NCP1,NCP2)


        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 

        INTEGER :: L1,M1,N1,L2,M2,N2,NPP
        INTEGER :: NALPHA,NBETA,NGAMMA
        INTEGER :: NCP1,NCP2
        REAL(wp) :: PPM(NPP),UPPM(NPP)
        REAL(wp) :: ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
        INTEGER :: NPAR(100)  !!!!!
        REAL(wp) :: PAR(50)   !!!!!!!
        REAL(wp) :: VLM(4)
        REAL(wp) :: HUPP
!local
        INTEGER :: I1,K1,I2,K2
        REAL(wp) :: S1,S2,S3,S4,GG12
        REAL(wp) :: VLMINT(4)

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

      CALL INT_GNLMP(PPM,UPPM,NPP,N1,L1,K1,N2,L2,K2,  &
                     HUPP,GG12)


      CALL INT_VLMM_VEC(L1,M1,K1,L2,M2,K2,ALPHAM,BETAM,GAMMAM, &
             NBETA,NGAMMA,NALPHA,NPAR,PAR,VLMINT,NCP1,NCP2)

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
      END SUBROUTINE WW_WW_VEC
!--------------------------------------------------------
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
         PP=PPM(I1)
         GG1M(I1)=GNLMP(N1,L1,K1,PP)
         GG2M(I1)=GNLMP(N2,L2,K2,PP)
 4    CONTINUE

!      CALL GNUFORM (NPP,1,1,GG1M,20) 

      DO 6 I1=1,NPP
         UPP=UPPM(I1)
!         UG1M(I1)= FUNL1(UPPM,GG1M,1,NPP,UPP,IERR)
!         UG2M(I1)= FUNL1(UPPM,GG2M,1,NPP,UPP,IERR)
         UG1M(I1)= FUNL1(PPM,GG1M,1,NPP,UPP,IERR)
         UG2M(I1)= FUNL1(PPM,GG2M,1,NPP,UPP,IERR)
 6    CONTINUE   

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
         VLMRE1=WLMM_RE(L1,M1,K1,ALPHA,BETA,GAMMA,NCP1)
! correspond to j
         VLMRE2=WLMM_RE(L2,M2,K2,ALPHA,BETA,GAMMA,NCP2)

         VLMIM1=WLMM_IM(L1,M1,K1,ALPHA,BETA,GAMMA,NCP1)
         VLMIM2=WLMM_IM(L2,M2,K2,ALPHA,BETA,GAMMA,NCP2)

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

      FUNCTION WLMM_RE(L,M,M1,ALPHA,BETA,GAMMA,NIJ)

!          Computer the W harmonics
!          VLMRE - real part

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

        INTEGER :: L,M,M1,NIJ
        REAL(wp) :: ALPHA,BETA,GAMMA
! local
        INTEGER ::KK
        REAL(wp) :: DM,SIGNUM,CC1,CC2
        REAL(wp) :: PLM1,PLM2
        REAL(wp) :: VLMRE
        REAL(wp) :: WUnit(3) 
! functions
        REAL(wp) :: PLM1M2,WLMM_RE 

         WUnit(1)=COS(ALPHA)*SIN(BETA)
         WUnit(2)=SIN(ALPHA)*SIN(BETA)
         WUnit(3)=COS(BETA)

        DM=1.
        IF(M1 .EQ. 0) DM=0.5

!      KK=MOD((M-M1),2)
        KK=MOD(M1,2)
        SIGNUM=1.
        IF(KK .NE. 0) SIGNUM=-1.

        CC1=COS(M*ALPHA + M1*GAMMA)
        CC2=COS(M*ALPHA - M1*GAMMA)

!        SS1= SIN(M*ALPHA + M1*GAMMA)
!        SS2= SIN(M*ALPHA - M1*GAMMA)

        PLM1 = PLM1M2(L,M, M1,BETA)
        PLM2 = PLM1M2(L,M,-M1,BETA)

        VLMRE=DM*(CC1*PLM1 + SIGNUM*CC2*PLM2)
        WLMM_RE=VLMRE*WUnit(NIJ)
      RETURN
      END FUNCTION WLMM_RE
!----------------------------------------------

      FUNCTION WLMM_IM(L,M,M1,ALPHA,BETA,GAMMA,NIJ)

!          Computer the V harmonics
!          VLMIM - image part
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

        INTEGER :: L,M,M1,NIJ
        REAL(wp) :: ALPHA,BETA,GAMMA
! local
        INTEGER ::KK
        REAL(wp) :: DM,SIGNUM,SS1,SS2
        REAL(wp) :: PLM1,PLM2
        REAL(wp) :: VLMIM
        REAL(wp) :: WUnit(3) 
! functions
        REAL(wp) :: PLM1M2,WLMM_IM 

         WUnit(1)=COS(ALPHA)*SIN(BETA)
         WUnit(2)=SIN(ALPHA)*SIN(BETA)
         WUnit(3)=COS(BETA)

        DM=1.
        IF(M1 .EQ. 0) DM=0.5

!      KK=MOD((M-M1),2)
        KK=MOD(M1,2)
        SIGNUM=1.
        IF(KK .NE. 0) SIGNUM=-1.

        SS1= SIN(M*ALPHA + M1*GAMMA)
        SS2= SIN(M*ALPHA - M1*GAMMA)

        PLM1 = PLM1M2(L,M, M1,BETA)
        PLM2 = PLM1M2(L,M,-M1,BETA)

        VLMIM=DM*(SS1*PLM1 + SIGNUM*SS2*PLM2)
        WLMM_IM=VLMIM*WUnit(NIJ)
      RETURN
      END FUNCTION WLMM_IM
!----------------------------------------------

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
        REAL(wp) :: YM((LL+1)*(LL+1)*NNMAX,NTEN*NTEN)
        REAL(wp) :: PPM(NPP),UPPM(NPP)
        REAL(wp) :: ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
        REAL(wp) :: PAR(50)

!      DIMENSION PROJ(NPP*NALPHA*NBETA*NGAMMA)

!      DIMENSION YM((LL+1)*(LL+1)*NNMAX,3)
!      DIMENSION PPM(NPP),UPPM(NPP)
!      DIMENSION ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
!      DIMENSION NPAR(*),PAR(*)

!local 
        REAL(wp) :: BM1(NBETA*NGAMMA,NTEN*NTEN),BM2(NBETA,NTEN*NTEN)
        REAL(wp) :: WUnit(3)
        INTEGER :: LMAX,N,L,L22,M,IIC,IIS
        INTEGER :: I,J,IJ,I1,I2,I3,I23
        REAL(wp) :: S1(NTEN),S2(NTEN,NTEN),S3,WW1,WW2,WW3
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

          DO 12 I2=1,NBETA
            BETA=BETAM(I2)

          DO 12 I3=1,NGAMMA
             GAMMA=GAMMAM(I3)
             I23=(I2-1)*NGAMMA+I3

!            DO 3 J=1,NTEN
!               S1(J)=0.
!            DO 3 I=1,NTEN
!               IJ=(J-1)*NTEN+I
!               S2(IJ)=0.
!3           CONTINUE   

         DO 3 J=1,NTEN
         DO 3 I=1,NTEN
            IJ=(J-1)*NTEN+I
            S2(I,J)=0.
3        CONTINUE  

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

            CALL INT_Proj_GNLM_P(PROJ,PPM,UPPM,NPP,         &
                             N,L22,M-1,                          & 
                             ALPHA,BETA,GAMMA,                   &
                             NJ,I1,I2,I3,PAR,YYRE,YYIM)

            DO 2 J=1,NTEN
            DO 2 I=1,NTEN   
               IJ=(J-1)*NTEN+I
               S2(I,J)=S2(I,J)+YYRE*WW1*WUnit(I)*WUnit(J)
!               S1(J)=1.
2           CONTINUE
!               S2(IJ)=S2(IJ)+S1(J)*WUnit(J)
!!               S2(IJ)=S2(IJ)+S1(J)
!4           CONTINUE

10       CONTINUE
         
         DO 5 J=1,NTEN
         DO 5 I=1,NTEN   
            IJ=(J-1)*NTEN+I
            BM1(I23,IJ)=S2(I,J)*HALPHA
5        CONTINUE

12    CONTINUE   

!            write(*,*) '1111111111',BM1
      DO 16 I2=1,NBETA

         DO 7 J=1,NTEN
         DO 7 I=1,NTEN   
            IJ=(J-1)*NTEN+I
            S2(I,J)=0.
7        CONTINUE

            DO 14 I3=1,NGAMMA

               I23=(I2-1)*NGAMMA+I3
               
               WW3=1.0
               IF(I3 .EQ.1 .OR. I3 .EQ. NGAMMA) WW3=0.5
               DO 9 J=1,NTEN
               DO 9 I=1,NTEN   
                  IJ=(J-1)*NTEN+I
                  S2(I,J)=S2(I,J)+BM1(I23,IJ)*WW3
9              CONTINUE
14          CONTINUE
            DO 11 J=1,NTEN
            DO 11 I=1,NTEN   
               IJ=(J-1)*NTEN+I
               BM2(I2,IJ)=S2(I,J)*HGAMMA
11          CONTINUE
16    CONTINUE


      DO 13 J=1,NTEN
      DO 13 I=1,NTEN   
         IJ=(J-1)*NTEN+I
         S2(I,J)=0.
13    CONTINUE

      DO 18 I2=1,NBETA
         BETA=BETAM(I2)
         WW2=1.0
         IF(I2 .EQ.1 .OR. I2 .EQ. NBETA) WW2=0.5
         DO 17 J=1,NTEN
         DO 17 I=1,NTEN   
            IJ=(J-1)*NTEN+I
            S2(I,J)=S2(I,J)+BM2(I2,IJ)*SIN(BETA)*WW2
17       CONTINUE
18    CONTINUE
         DO 19 J=1,NTEN
         DO 19 I=1,NTEN   
            IJ=(J-1)*NTEN+I
            YM(IIC,IJ)=S2(I,J)*HBETA
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

         DO 31 J=1,NTEN
!            S1(J)=0.
         DO 31 I=1,NTEN
            IJ=(J-1)*NTEN+I
            S2(I,J)=0.
31       CONTINUE   

         WUnit(1)=-SIN(BETA)*COS(GAMMA)
         WUnit(2)= SIN(BETA)*SIN(GAMMA)
         WUnit(3)= COS(BETA)

      DO 20 I1=1,NALPHA
         ALPHA=ALPHAM(I1)
!         WUnit(1)=COS(ALPHA)*SIN(BETA)
!         WUnit(2)=SIN(ALPHA)*SIN(BETA)
!         WUnit(3)=COS(BETA)
         WW1=1.0
         IF(I1 .EQ.1 .OR. I1 .EQ. NALPHA) WW1=0.5

         CALL INT_Proj_GNLM_P(PROJ,PPM,UPPM,NPP,         &
                     N,L22,M,                          & 
                     ALPHA,BETA,GAMMA,                   &
                     NJ,I1,I2,I3,PAR,YYRE,YYIM)
         
         DO 41 J=1,NTEN
         DO 41 I=1,NTEN   
            IJ=(J-1)*NTEN+I
!            S1(J)=S1(J)+YYIM*WW1*WUnit(I)
!!            S1(J)=1.
!21       CONTINUE
            S2(I,J)=S2(I,J)+YYIM*WW1*WUnit(I)*WUnit(J)
!!            S2(IJ)=S2(IJ)+S1(J)
41       CONTINUE

20    CONTINUE
         
         DO 51 J=1,NTEN
         DO 51 I=1,NTEN   
            IJ=(J-1)*NTEN+I
            BM1(I23,IJ)=S2(I,J)*HALPHA
51       CONTINUE

22    CONTINUE   

      DO 26 I2=1,NBETA
         BETA=BETAM(I2)

         DO 71 J=1,NTEN
         DO 71 I=1,NTEN   
            IJ=(J-1)*NTEN+I
            S2(I,J)=0.
71       CONTINUE

      DO 24 I3=1,NGAMMA
         I23=(I2-1)*NGAMMA+I3

         WW3=1.0
         IF(I3 .EQ.1 .OR. I3 .EQ. NGAMMA) WW3=0.5
         DO 91 J=1,NTEN
         DO 91 I=1,NTEN   
            IJ=(J-1)*NTEN+I
            S2(I,J)=S2(I,J)+BM1(I23,IJ)*WW3
91       CONTINUE
24    CONTINUE

         DO 111 J=1,NTEN
         DO 111 I=1,NTEN   
            IJ=(J-1)*NTEN+I
            BM2(I2,IJ)=S2(I,J)*HGAMMA
111      CONTINUE
26    CONTINUE


         DO 131 J=1,NTEN
         DO 131 I=1,NTEN   
            IJ=(J-1)*NTEN+I
            S2(I,J)=0.
131       CONTINUE

      DO 28 I2=1,NBETA
         BETA=BETAM(I2)
         WW2=1.0
         IF(I2 .EQ.1 .OR. I2 .EQ. NBETA) WW2=0.5
         DO 171 J=1,NTEN
         DO 171 I=1,NTEN   
            IJ=(J-1)*NTEN+I
            S2(I,J)=S2(I,J)+BM2(I2,IJ)*SIN(BETA)*WW2
171       CONTINUE

28    CONTINUE
         DO 191 J=1,NTEN
         DO 191 I=1,NTEN   
            IJ=(J-1)*NTEN+I
            YM(IIS,IJ)=S2(I,J)*HBETA
191       CONTINUE
50    CONTINUE
      RETURN
    END SUBROUTINE RightPart_WW_TEN
!--------------------------------------

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
        REAL(wp) :: YM((LL+1)*(LL+1)*NNMAX,NTEN*NTEN)
        REAL(wp) :: PPM(NPP),UPPM(NPP)
        REAL(wp) :: ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
        REAL(wp) :: PAR(50)

!      DIMENSION PROJ(NPP*NALPHA*NBETA*NGAMMA)

!      DIMENSION YM((LL+1)*(LL+1)*NNMAX,3)
!      DIMENSION PPM(NPP),UPPM(NPP)
!      DIMENSION ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
!      DIMENSION NPAR(*),PAR(*)

!local 
        REAL(wp) :: BM1(NBETA*NGAMMA,NTEN*NTEN),BM2(NBETA,NTEN*NTEN)
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

            DO 3 J=1,NTEN
               S1(J)=0.
            DO 3 I=1,NTEN
               IJ=(J-1)*NTEN+I
               S2(IJ)=0.
3           CONTINUE   

         DO 10 I1=1,NALPHA
            ALPHA=ALPHAM(I1)

            WUnit(1)=COS(ALPHA)*SIN(BETA)
            WUnit(2)=SIN(ALPHA)*SIN(BETA)
            WUnit(3)=COS(BETA)

!         WUnit(1)=-SIN(BETA)*COS(GAMMA)
!         WUnit(2)= SIN(BETA)*SIN(GAMMA)
!         WUnit(3)= COS(BETA)

!            write(*,*) 'RightPart_WW_VEC_11',L22,M-1

            WW1=1.0
            IF(I1 .EQ.1 .OR. I1 .EQ. NALPHA) WW1=0.5

            CALL INT_Proj_GNLM_P(PROJ,PPM,UPPM,NPP,         &
                             N,L22,M-1,                          & 
                             ALPHA,BETA,GAMMA,                   &
                             NJ,I1,I2,I3,PAR,YYRE,YYIM)

            DO 4 J=1,NTEN
            DO 2 I=1,NTEN   
               IJ=(J-1)*NTEN+I
               S1(J)=S1(J)+YYRE*WW1*WUnit(I)
!               S1(J)=1.
2           CONTINUE
               S2(IJ)=S2(IJ)+S1(J)*WUnit(J)
!               S2(IJ)=S2(IJ)+S1(J)
4           CONTINUE

10       CONTINUE
         
         DO 5 J=1,NTEN
         DO 5 I=1,NTEN   
            IJ=(J-1)*NTEN+I
            BM1(I23,IJ)=S2(IJ)*HALPHA
5        CONTINUE

12    CONTINUE   

!            write(*,*) '1111111111',BM1
      DO 16 I3=1,NGAMMA

         DO 7 J=1,NTEN
         DO 7 I=1,NTEN   
            IJ=(J-1)*NTEN+I
            S2(IJ)=0.
7        CONTINUE

            DO 14 I2=1,NBETA
               BETA=BETAM(I2)
               I23=(I3-1)*NBETA+I2
               
               WW2=1.0
               IF(I2 .EQ.1 .OR. I2 .EQ. NBETA) WW2=0.5
               DO 9 J=1,NTEN
               DO 9 I=1,NTEN   
                  IJ=(J-1)*NTEN+I
                  S2(IJ)=S2(IJ)+BM1(I23,IJ)*SIN(BETA)*WW2
9              CONTINUE
14          CONTINUE
            DO 11 J=1,NTEN
            DO 11 I=1,NTEN   
               IJ=(J-1)*NTEN+I
               BM2(I3,IJ)=S2(IJ)*HBETA
11          CONTINUE
16    CONTINUE


      DO 13 J=1,NTEN
      DO 13 I=1,NTEN   
         IJ=(J-1)*NTEN+I
         S2(IJ)=0.
13    CONTINUE

      DO 18 I3=1,NGAMMA
         WW3=1.0
         IF(I3 .EQ.1 .OR. I3 .EQ. NGAMMA) WW3=0.5
         DO 17 J=1,NTEN
         DO 17 I=1,NTEN   
            IJ=(J-1)*NTEN+I
            S2(IJ)=S2(IJ)+BM2(I3,IJ)*WW3
17       CONTINUE
18    CONTINUE
         DO 19 J=1,NTEN
         DO 19 I=1,NTEN   
            IJ=(J-1)*NTEN+I
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

         DO 31 J=1,NTEN
            S1(J)=0.
         DO 31 I=1,NTEN
            IJ=(J-1)*NTEN+I
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
         
         DO 41 J=1,NTEN
         DO 21 I=1,NTEN   
            IJ=(J-1)*NTEN+I
            S1(J)=S1(J)+YYIM*WW1*WUnit(I)
!            S1(J)=1.
21       CONTINUE
            S2(IJ)=S2(IJ)+S1(J)*WUnit(J)
!            S2(IJ)=S2(IJ)+S1(J)
41       CONTINUE

20    CONTINUE
         
         DO 51 J=1,NTEN
         DO 51 I=1,NTEN   
            IJ=(J-1)*NTEN+I
            BM1(I23,IJ)=S2(IJ)*HALPHA
51       CONTINUE

22    CONTINUE   

      DO 26 I3=1,NGAMMA

         DO 71 J=1,NTEN
         DO 71 I=1,NTEN   
            IJ=(J-1)*NTEN+I
            S2(IJ)=0.
71       CONTINUE

      DO 24 I2=1,NBETA
         BETA=BETAM(I2)
         I23=(I3-1)*NBETA+I2

         WW2=1.0
         IF(I2 .EQ.1 .OR. I2 .EQ. NBETA) WW2=0.5
         DO 91 J=1,NTEN
         DO 91 I=1,NTEN   
            IJ=(J-1)*NTEN+I
            S2(IJ)=S2(IJ)+BM1(I23,IJ)*SIN(BETA)*WW2
91       CONTINUE
24    CONTINUE

         DO 111 J=1,NTEN
         DO 111 I=1,NTEN   
            IJ=(J-1)*NTEN+I
            BM2(I3,IJ)=S2(IJ)*HBETA
111      CONTINUE
26    CONTINUE


         DO 131 J=1,NTEN
         DO 131 I=1,NTEN   
            IJ=(J-1)*NTEN+I
            S2(IJ)=0.
131       CONTINUE

      DO 28 I3=1,NGAMMA
         WW3=1.0
         IF(I3 .EQ.1 .OR. I3 .EQ. NGAMMA) WW3=0.5
         DO 171 J=1,NTEN
         DO 171 I=1,NTEN   
            IJ=(J-1)*NTEN+I
            S2(IJ)=S2(IJ)+BM2(I3,IJ)*WW3
171       CONTINUE

28    CONTINUE
         DO 191 J=1,NTEN
         DO 191 I=1,NTEN   
            IJ=(J-1)*NTEN+I
            YM(IIS,IJ)=S2(IJ)*HGAMMA
191       CONTINUE
50    CONTINUE
      RETURN
    END SUBROUTINE RightPart_WW_VEC
!--------------------------------------

       SUBROUTINE INT_Proj_GNLM_P(PROJ,PPM,UPPM,NPP,  &
                        N1,L1,M1,                              &
                        ALPHA,BETA,GAMMA,                      &
                        NJ,I1,I2,I3,PAR,YYRE,YYIM)
!   integration over p variable of 
!   g(p,al,bet,gam)*g(n,l,m)
!      
        USE PARAM_INPUT

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

       
        INTEGER :: NPP,N1,L1,M1,NJ,I1,I2,I3
        REAL(wp) :: ALPHA,BETA,GAMMA
        REAL(wp) :: YYRE,YYIM
        REAL(wp) :: PPM(NPP),UPPM(NPP)
        REAL(wp) :: PROJ(NPP*NJ)
        REAL(wp) :: PAR(50)
!local
        REAL(wp) :: PROJ2(NPP),PROJ3(NPP)
        REAL(wp) :: PP,UPP,HUPP
        REAL(wp) :: WW1,WWLMRE,WWLMIM
        REAL(wp) :: S1,S2
        INTEGER :: NALPHA,NBETA,NGAMMA
        INTEGER :: I4,I321,I4321,IERR
!functions
        REAL(wp) :: FUNL1

!      DIMENSION PPM(NPP),UPPM(NPP)
!      DIMENSION GG1M(NPP),GG2m(NPP)
!      DIMENSION PROJ(NPP*NJ)

!      DIMENSION PROJ2(NPP),PROJ3(NPP)
!      DIMENSION NPAR(*),PAR(*)

      NALPHA=PARAM(14)
      NBETA =PARAM(15)
      NGAMMA=PARAM(16)
      HUPP  =PAR(35)


!      write(*,*) 'INT_Proj_GNLM_P',L1,M1

      DO 6 I4=1,NPP
!         PP=PPM(I4)
         I321=(I3-1)*NBETA*NALPHA+(I2-1)*NALPHA+I1
         I4321=(I4-1)*NGAMMA*NBETA*NALPHA + I321

!         PROJ2(I4)=PROJ(I4321)
         PROJ3(I4)=PROJ(I4321)
 6    CONTINUE   
!         write(*,*) '222222222',PROJ

!      DO 8 I4=1,NPP
!         UPP=UPPM(I4)
!!         PROJ3(I4)= FUNL1(UPPM,PROJ2,1,NPP,UPP,IERR)
!         PROJ3(I4)= FUNL1(PPM,PROJ2,1,NPP,UPP,IERR)
! 8    CONTINUE   

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

!      write(*,*) '2222222222',YYRE,YYIM
      RETURN
      END SUBROUTINE INT_Proj_GNLM_P
!------------------------------------------------------------
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
        
!        write(*,*) 'WWCS',L,M

        S1=0.
        S2=0.
        DO 10 I=1,L+1
           M1=I-1
           CALL VLMM_RE(L,M,M1,ALPHA,BETA,GAMMA,VLMRE,VLMIM)
           GGLM= GNLMP(N,L,M1,PP)

           S1=S1+VLMRE*GGLM
           S2=S2+VLMIM*GGLM
10      CONTINUE
      WWLMRE=S1
      WWLMIM=S2

!      write(*,*) '3333333333',WWLMRE,WWLMIM
      RETURN
      END SUBROUTINE WWCS
!------------------------------------------------
      SUBROUTINE VLMM_RE(L,M,M1,ALPHA,BETA,GAMMA,VLMRE,VLMIM)

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

        VLMRE=DM*(CC1*PLM1 + SIGNUM*CC2*PLM2)
        VLMIM=DM*(SS1*PLM1 + SIGNUM*SS2*PLM2)

!      VLMRE=DM*(CC1*PLM1*(1+SIG2) + CC2*PLM2*(1+SIG1))
!      VLMIM=DM*(SS1*PLM1*(1+SIG2) + SS2*PLM2*(1+SIG1))

      RETURN
      END SUBROUTINE VLMM_RE
!----------------------------------------------

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
        REAL(wp) :: WUnit(3)
        INTEGER :: II,I1,I2,I3,I321
        INTEGER :: I4,I4321,I,J,IJ
        INTEGER :: N,L,L22,M,IIC,IIS,LMAX
        REAL(wp) :: ALPHA,BETA,GAMMA
        REAL(wp) :: WWLMRE,WWLMIM
        REAL(wp) :: S(NTEN),PP,S1

!      REAL PPM(NPP),ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
!      REAL GM(NPP*NALPHA*NBETA*NGAMMA)
!      REAL XXM((LL+1)*(LL+1)*NNMAX,3)
!      REAL WWLMRE(3),WWLMIM(3)

!      COMMON/PRIN/IPRNT
!      PARAMETER (PI = 3.1415926536)

      LMAX=(LL+1)*(LL+1)

!      DO 10 I4=1,NPP
!         PP=PPM(I4)
!      DO 10 I3=1,NGAMMA
!         GAMMA=GAMMAM(I3)
!      DO 10 I2=1,NBETA
!         BETA=BETAM(I2)
!      DO 10 I1=1,NALPHA
!         ALPHA=ALPHAM(I1)


      DO 10 II=1,NJ

         CALL IALBETGAM(II,NALPHA,NBETA,I1,I2,I3)

         ALPHA=ALPHAM(I1)
         BETA =BETAM(I2)
         GAMMA=GAMMAM(I3)
         I321=(I3-1)*NBETA*NALPHA+(I2-1)*NALPHA+I1

!         WUnit(1)=COS(ALPHA)*SIN(BETA)
!         WUnit(2)=SIN(ALPHA)*SIN(BETA)
!         WUnit(3)=COS(BETA)

         WUnit(1)=-SIN(BETA)*COS(GAMMA)
         WUnit(2)= SIN(BETA)*SIN(GAMMA)
         WUnit(3)= COS(BETA)


      DO 10 I4=1,NPP
         PP=PPM(I4)

         I4321=(I4-1)*NGAMMA*NBETA*NALPHA + I321

!         I4321=(II-1)*NPP + I4

        S1=0.
        DO 4  N=1,NNMAX
        DO 4  L=1,LL+1
           L22=L-1
        DO 4  M=1,L22+1
           IIC=(N-1)*LMAX + (L22*L22+M)

           CALL WWCS(N,L22,M-1,ALPHA,BETA,GAMMA,PP,WWLMRE,WWLMIM)

           DO 5 J=1,NTEN
              S(J)=0.
           DO 3 I=1,NTEN   
              IJ=(J-1)*NTEN+I
              S(J)=S(J)+XXM(IIC,IJ)*WWLMRE*WUnit(I)
3          CONTINUE
              S1=S1+S(J)*WUnit(J)
5          CONTINUE
!           S1=S1+(XXM(IIC,1)*WUnit(1)+    &
!                  XXM(IIC,2)*WUnit(2)+    &
!                  XXM(IIC,3)*WUnit(3))*WWLMRE
4       CONTINUE

        DO 6 N=1,NNMAX
        DO 6 L=1,LL+1
           L22=L-1
        DO 6 M=1,L22
           IIS=(N-1)*LMAX + (L22*L22+(M+1)+L22)

           CALL WWCS(N,L22,M,ALPHA,BETA,GAMMA,PP,WWLMRE,WWLMIM)

!c           CALL YLM_RE (L,M,TETA,PHI,YLMRE,YLMIM)

           DO 9 J=1,NTEN
              S(J)=0.
           DO 7 I=1,NTEN   
              IJ=(J-1)*NTEN+I
              S(J)=S(J)+XXM(IIS,IJ)*WWLMIM*WUnit(I)
7          CONTINUE
              S1=S1+S(J)*WUnit(J)
9          CONTINUE

!           S1=S1+(XXM(IIS,1)*WUnit(1)+   &
!                  XXM(IIS,2)*WUnit(2)+   &
!                  XXM(IIS,3)*WUnit(3))*WWLMIM
 6      CONTINUE
        GM(I4321)=S1
 10   CONTINUE
      RETURN
      END SUBROUTINE Proj_Summa_W_Harm
!------------------------------------------------------
      SUBROUTINE ProjNum_TEN(DR,PROJ,UM,VM,WM,  &
              ALPHAM,BETAM,GAMMAM,NALPHA,NBETA,NGAMMA, &
              XM,YM,ZM,NX,NY,NZ,NU,NV,NW,NTEN)
!
        USE PARAM_INPUT

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 
      
        INTEGER :: NX,NY,NZ,NU,NV,NW,NTEN
        INTEGER :: NALPHA,NBETA,NGAMMA
        REAL(wp) :: DR(NX*NY*NZ,NTEN*NTEN)
        REAL(wp) :: ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
        REAL(wp) ::  UM(NU),VM(NV),WM(NW)
        REAL(wp) :: XM(NX),YM(NY),ZM(NZ)
        REAL(wp) :: PROJ(NU*NALPHA*NBETA*NGAMMA)

!local 
        REAL(wp) :: DDR1(NX*NY*NZ,NTEN), DDR(NX*NY*NZ)
        REAL(wp) :: PROJ2(NU)
        REAL(wp) :: WUnit(3)
        REAL(wp) :: ALPHA,BETA,GAMMA 
        REAL(wp) :: CALPHA,SALPHA,CBETA,SBETA,SGAMMA,CGAMMA
        REAL(wp) :: UU,VV,WW,HW,QQQ,WW0,S,SS
        REAL(wp) :: XX,YY,ZZ
        REAL(wp) :: RRR,RRE,WW1,EI
        REAL(wp) :: AA11,AA12,AA13,AA21,AA22,AA23,AA31,AA32,AA33
        INTEGER :: NJ,ITEN
        INTEGER :: I,I4,I1,I2,I3,I321,I4321
        INTEGER :: J,IJ,J1,J3,J321
        INTEGER :: IERR  ! ?????????????????????
        CHARACTER*8 ITEX
!functions
        REAL(wp) :: FUNL3


!      DIMENSION  ALPHAM(*),BETAM(*),GAMMAM(*)
!      DIMENSION  BKAPPAM(*)
!      DIMENSION  UM(*),VM(*),WM(*),XM(*),YM(*),ZM(*)
!      DIMENSION  PROJ(*),PROJ2(*)
!      DIMENSION DR(NX*NY*NZ,3)

!      DIMENSION DDR(NX*NY*NZ)
!      DIMENSION WUnit(3)
!      DIMENSION  NPAR(*),PAR(*)
!      CHARACTER*8 ITEX


!      NALPHA = NPAR(9)
!      NBETA  = NPAR(10)
!      NGAMMA = NPAR(18)
!      NX    = NPAR(12)
!      NY    = NPAR(13)
!      NZ    = NPAR(14)
!      NU    = NPAR(4)
!      NV    = NPAR(5)
!      NW    = NPAR(6)

!      DD  = PAR(36)
!      RCC = PAR(42)
!      UM1 = PAR(43)

      NJ = NALPHA*NBETA*NGAMMA     !NPAR(11)
      RRE= PARAM(5) 
      
      IERR=0
      ITEX='ProjNum_Vec'
      IF(NW.LT.100) GOTO 2
      IERR=1
      CALL ERRPRN(ITEX,IERR)
      RETURN
2     CONTINUE

!      IERR=0
!      ITEX='ProjNum_Vec'
!      IF(NTEN .EQ. 3) GOTO 3
!      IERR=2
!      CALL ERRPRN(ITEX,IERR)
!      RETURN
!3     CONTINUE

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

         WUnit(1)=-SBETA*CGAMMA
         WUnit(2)= SBETA*SGAMMA
         WUnit(3)= CBETA

!      DO 6 I3=1,NZ
!      DO 6 I2=1,NY
!      DO 6 I1=1,NX
!         J321=(I3-1)*NX*NY+(I2-1)*NX +I1

!         DDR(J321)= DR(J321,1)*WUnit(1) &
!                 + DR(J321,2)*WUnit(2) &
!                 + DR(J321,3)*WUnit(3)
! 6    CONTINUE   
      DO 4 I3=1,NZ
      DO 4 I2=1,NY
      DO 4 I1=1,NX
         J321=(I3-1)*NX*NY+(I2-1)*NX +I1
         DDR(J321)=0.
4     CONTINUE

      DO 30 I3=1,NZ
      DO 30 I2=1,NY
      DO 30 I1=1,NX
         J321=(I3-1)*NX*NY+(I2-1)*NX +I1
         DO 24 J=1,NTEN   
            DDR1(J321,J)=0
         DO 22 I=1,NTEN
            IJ=(J-1)*NTEN+I
            DDR1(J321,J)=DDR1(J321,J)+DR(J321,IJ)*WUnit(I)
22       CONTINUE
            DDR(J321)= DDR(J321)+DDR1(J321,J)*WUnit(J) 
24       CONTINUE
30    CONTINUE


!      DO 65 J1=1,NV
!         VV=VM(J1)
!         QQQ=RRE**2-VV**2
!         IF(QQQ .LT. 0) QQQ=0.
            VV=0.                       !!!!!!!!!!!!!!!!!!!
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

!      -1
!--   T  (AL,BET,GAM)-------------------------------------------

!            XX=  AA11*UU + AA21*VV + AA31*WW
!            YY=  AA12*UU + AA22*VV + AA32*WW
!            ZZ=  AA13*UU + AA23*VV + AA33*WW


            XX=  AA11*UU + AA12*VV + AA13*WW
            YY=  AA21*UU + AA22*VV + AA23*WW
            ZZ=  AA31*UU + AA32*VV + AA33*WW

            CALL REGN(XM,YM,ZM,NX,NY,NZ,XX,YY,ZZ,EI)

            IF(EI .LT. 0.) GO TO 55

            S=FUNL3(XM,YM,ZM,DDR,NX,NY,NZ,XX,YY,ZZ,IERR)

            WW1=1.
            IF(J3.EQ.1.OR.J3.EQ.NW) WW1=0.5
            SS=SS+S*WW1

 55   CONTINUE
      
            PROJ2(J1)=SS*HW 
   
 65   CONTINUE

!      DO 75 J1=1,NV
      DO 75 J1=1,NU
         I4321=(J1-1)*NGAMMA*NBETA*NALPHA + I321
         PROJ(I4321)=PROJ2(J1)       

  75  CONTINUE
  85  CONTINUE
      RETURN
      END SUBROUTINE ProjNum_TEN
!------------------------------------------


      SUBROUTINE ProjNum_Vec(DR,PROJ,UM,VM,WM,  &
              ALPHAM,BETAM,GAMMAM,NALPHA,NBETA,NGAMMA, &
              XM,YM,ZM,NX,NY,NZ,NU,NV,NW,NTEN)
!
        USE PARAM_INPUT

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 


       
        INTEGER :: NX,NY,NZ,NU,NV,NW,NTEN
        INTEGER :: NALPHA,NBETA,NGAMMA
        REAL(wp) :: DR(NX*NY*NZ,NTEN*NTEN)
        REAL(wp) :: ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
        REAL(wp) ::  UM(NU),VM(NV),WM(NW)
        REAL(wp) :: XM(NX),YM(NY),ZM(NZ)
        REAL(wp) :: PROJ(NU*NALPHA*NBETA*NGAMMA)

!local 
        REAL(wp) :: DDR1(NX*NY*NZ,NTEN), DDR(NX*NY*NZ)
        REAL(wp) :: PROJ2(NU)
        REAL(wp) :: WUnit(3)
        REAL(wp) :: ALPHA,BETA,GAMMA 
        REAL(wp) :: CALPHA,SALPHA,CBETA,SBETA,SGAMMA,CGAMMA
        REAL(wp) :: UU,VV,WW,HW,QQQ,WW0,S,SS
        REAL(wp) :: XX,YY,ZZ
        REAL(wp) :: RRR,RRE,WW1,EI
        REAL(wp) :: AA11,AA12,AA13,AA21,AA22,AA23,AA31,AA32,AA33
        INTEGER :: NJ,ITEN
        INTEGER :: I,I4,I1,I2,I3,I321,I4321
        INTEGER :: J,IJ,J1,J3,J321
        INTEGER :: IERR  ! ?????????????????????
        CHARACTER*8 ITEX
!functions
        REAL(wp) :: FUNL3


!      DIMENSION  ALPHAM(*),BETAM(*),GAMMAM(*)
!      DIMENSION  BKAPPAM(*)
!      DIMENSION  UM(*),VM(*),WM(*),XM(*),YM(*),ZM(*)
!      DIMENSION  PROJ(*),PROJ2(*)
!      DIMENSION DR(NX*NY*NZ,3)

!      DIMENSION DDR(NX*NY*NZ)
!      DIMENSION WUnit(3)
!      DIMENSION  NPAR(*),PAR(*)
!      CHARACTER*8 ITEX


!      NALPHA = NPAR(9)
!      NBETA  = NPAR(10)
!      NGAMMA = NPAR(18)
!      NX    = NPAR(12)
!      NY    = NPAR(13)
!      NZ    = NPAR(14)
!      NU    = NPAR(4)
!      NV    = NPAR(5)
!      NW    = NPAR(6)

!      DD  = PAR(36)
!      RCC = PAR(42)
!      UM1 = PAR(43)

      NJ = NALPHA*NBETA*NGAMMA     !NPAR(11)
      RRE= PARAM(5) 
      
      IERR=0
      ITEX='ProjNum_Vec'
      IF(NW.LT.100) GOTO 2
      IERR=1
      CALL ERRPRN(ITEX,IERR)
      RETURN
2     CONTINUE

!      IERR=0
!      ITEX='ProjNum_Vec'
!      IF(NTEN .EQ. 3) GOTO 3
!      IERR=2
!      CALL ERRPRN(ITEX,IERR)
!      RETURN
!3     CONTINUE

      DO 85 I4=1,NJ 
!      CALL IALBET(I4,NALPHA,I1,I2)
         CALL IALBETGAM(I4,NALPHA,NBETA,I1,I2,I3)
         ALPHA=ALPHAM(I1)
         BETA =BETAM(I2)
         GAMMA=GAMMAM(I3)
         I321=(I3-1)*NBETA*NALPHA+(I2-1)*NALPHA+I1
!         I321=(I2-1)*NALPHA+I1

!         GAMMA=0.    !!!!!!!!!!!!!!!!!!!!!!

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

         WUnit(1)= CALPHA*SBETA  
         WUnit(2)= SALPHA*SBETA 
         WUnit(3)= CBETA        

!         WUnit(1)=-SBETA*CGAMMA
!         WUnit(2)= SBETA*SGAMMA
!         WUnit(3)= CBETA


      DO 4 I3=1,NZ
      DO 4 I2=1,NY
      DO 4 I1=1,NX
         J321=(I3-1)*NX*NY+(I2-1)*NX +I1
         DDR(J321)=0.
4     CONTINUE

      DO 30 I3=1,NZ
      DO 30 I2=1,NY
      DO 30 I1=1,NX
         J321=(I3-1)*NX*NY+(I2-1)*NX +I1
         DO 24 J=1,NTEN   
            DDR1(J321,J)=0
         DO 22 I=1,NTEN
            IJ=(J-1)*NTEN+I
            DDR1(J321,J)=DDR1(J321,J)+DR(J321,IJ)*WUnit(I)
22       CONTINUE
            DDR(J321)= DDR(J321)+DDR1(J321,J)*WUnit(J) 
24       CONTINUE
30    CONTINUE


!      DO 65 J1=1,NV
!         VV=VM(J1)
!         QQQ=RRE**2-VV**2
!         IF(QQQ .LT. 0) QQQ=0.
            VV=0.                       !!!!!!!!!!!!!!!!!!!
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

!      -1
!--   T  (AL,BET,GAM)-------------------------------------------

            XX=  AA11*UU + AA21*VV + AA31*WW
            YY=  AA12*UU + AA22*VV + AA32*WW
            ZZ=  AA13*UU + AA23*VV + AA33*WW


!            XX=  AA11*UU + AA12*VV + AA13*WW
!            YY=  AA21*UU + AA22*VV + AA23*WW
!            ZZ=  AA31*UU + AA32*VV + AA33*WW

            CALL REGN(XM,YM,ZM,NX,NY,NZ,XX,YY,ZZ,EI)

            IF(EI .LT. 0.) GO TO 55

            S=FUNL3(XM,YM,ZM,DDR,NX,NY,NZ,XX,YY,ZZ,IERR)

            WW1=1.
            IF(J3.EQ.1.OR.J3.EQ.NW) WW1=0.5
            SS=SS+S*WW1

 55   CONTINUE
      
            PROJ2(J1)=SS*HW 
   
 65   CONTINUE

!      DO 75 J1=1,NV
      DO 75 J1=1,NU
         I4321=(J1-1)*NGAMMA*NBETA*NALPHA + I321
         PROJ(I4321)=PROJ2(J1)       

  75  CONTINUE
  85  CONTINUE
      RETURN
      END SUBROUTINE ProjNum_Vec
!------------------------------------------

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

!      REAL A,B
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
        REAL(wp) :: GNLMP,ACOBI,GCOF

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
      END FUNCTION GNLMP
!-----------------------------------------------------
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
        REAL(wp) :: XM((LL+1)*(LL+1)*NNMAX,NTEN*NTEN)

!local 
        INTEGER ::LMAX,N,L,L22,M,IIC,IIS
        INTEGER :: I,J,IJ,MU,NU
        REAL(wp) :: AAC(9),AAS(9)
! functions
!        REAL(wp) :: AACOS1,AACOS2,AACOS3,AACOS4,AACOS5
!        REAL(wp) :: AASIN1,AASIN2,AASIN3,AASIN4,AASIN5

        IF(NTEN .GT. 3) THEN  
           WRITE(*,*) 'Procedure XXMM: NTEN > 3'
           STOP
        ENDIF

        LMAX=(LL+1)*(LL+1)

        DO 10  N=1,NNMAX
        DO 10  L=1,LL+1
           L22=L-1
        DO 10  M=1,L22+1
           IIC=(N-1)*LMAX + (L22*L22+M)
!         write(*,*) 'XXMM-cos',L,M,IIC
           CALL ACOS(L22,M-1,N,AAC)
!           AAC1=AACOS1(L22,M-1,N)
!           AAC2=AACOS2(L22,M-1,N)
!           AAC3=AACOS3(L22,M-1,N)
!           AAC4=AACOS4(L22,M-1,N)
!           AAC5=AACOS5(L22,M-1,N)
!           AAC6=AACOS6(L22,M-1,N)
!           AAC7=AACOS7(L22,M-1,N)
!           AAC8=AACOS8(L22,M-1,N)
!           AAC9=AACOS9(L22,M-1,N)
           DO 2 J=1,NTEN
           DO 2 I=1,NTEN
              IJ=(J-1)*NTEN+I
              MU=LMAX*(IJ-1)*NNMAX+IIC
!              XM(IIC,IJ)=AAC(IJ)
              XM(MU,IJ)=AAC(IJ)
2          CONTINUE

!           XM(IIC,1)=AAC(1)
!           XM(IIC,2)=AAC(2)
!           XM(IIC,3)=AAC(3)
!           XM(IIC,4)=AAC(4)
!           XM(IIC,5)=AAC(5)
!           XM(IIC,6)=AAC(6)
!           XM(IIC,7)=AAC(7)
!           XM(IIC,8)=AAC(8)
!           XM(IIC,9)=AAC(9)
           
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
!           AAS1=AASIN1(L22,M,N)
!           AAS2=AASIN2(L22,M,N)
!           AAS3=AASIN3(L22,M,N)
!           AAS4=AASIN4(L22,M,N)
!           AAS5=AASIN5(L22,M,N)
!           AAS6=AASIN6(L22,M,N)
!           AAS7=AASIN7(L22,M,N)
!           AAS8=AASIN8(L22,M,N)
!           AAS9=AASIN9(L22,M,N)
           DO 4 J=1,NTEN
           DO 4 I=1,NTEN
              IJ=(J-1)*NTEN+I
              NU=LMAX*(IJ-1)*NNMAX+IIS
              XM(NU,IJ)=AAS(IJ)
4          CONTINUE
!           XM(IIS,1)=AAS(1)
!           XM(IIS,2)=AAS(2)
!           XM(IIS,3)=AAS(3)
!           XM(IIS,4)=AAS(4)
!           XM(IIS,5)=AAS(5)
!           XM(IIS,6)=AAS(6)
!           XM(IIS,7)=AAS(7)
!           XM(IIS,8)=AAS(8)
!           XM(IIS,9)=AAS(9)

!           write(*,*) 'XXMM:IIS'
!           write(*,*) IIS,L22,M
 20     CONTINUE
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

        AACOS(4)=0.5   
        IF(L .EQ.0 .AND. M .EQ.0) AACOS(4)=0.5 

        AACOS(5)=0.   
        IF(L .EQ.0 .AND. M .EQ.0) AACOS(5)=0. 

        AACOS(6)=0.   
        IF(L .EQ.0 .AND. M .EQ.0) AACOS(6)=0. 

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

      SUBROUTINE Vec_Line(XM,YM,LD,NTEN)
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        INTEGER:: LD,NTEN
        REAL(wp):: XM(LD,NTEN*NTEN), YM(NTEN*NTEN*LD)
! local 
        INTEGER:: K,I,II

        DO 4 K=1,NTEN*NTEN
        DO 4 I=1,LD
           II=(K-1)*LD+I
           YM(II)=XM(I,K)      
!           IF(II .GT. NTEN2*LD*0.6) THEN   !можно будет убрать
!              YM(II)=0.0                 !можно будет убрать  
!           ENDIF                          !можно будет убрать  
4       CONTINUE
      RETURN
      END SUBROUTINE Vec_Line
!---------------------------------------------
      SUBROUTINE Vec_Line_Inv(YM,XM,LD,NTEN)
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 

        INTEGER:: LD,NTEN
        REAL(wp):: XM(LD,NTEN*NTEN), YM(NTEN*NTEN*LD)
! local 
        INTEGER:: K,I,II

        DO 4 K=1,NTEN*NTEN
        DO 4 I=1,LD
           II=(K-1)*LD+I
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



      SUBROUTINE Matrix_WW_VEC_1(LL,PPM,ALPHAM,NPP,NALPHA,AM)

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 

        INTEGER :: LL,NPP,NALPHA
!        REAL(wp) :: AM(4*LL*(LL+2),4*LL*(LL+2))
        INTEGER :: AM(4*LL*(LL+2),4*LL*(LL+2))
        REAL(wp) :: PPM(NPP),ALPHAM(NALPHA)
! local 
        INTEGER :: N1,L1,M1,N,L,M
       INTEGER :: IIC,IIC1,IIS,IIS1
        INTEGER :: LMAX

!      DIMENSION UPPM(NPP)
!      DIMENSION NPAR(*),PAR(*)
!      DIMENSION VLM(4)


!      HALPHA=PAR(38)
!      HBETA =PAR(39)
!      HGAMMA=PAR(40)
!      HUPP  =PAR(41)
    

!      write(*,*) 'Matrix_WW',HALPHA,HBETA,HGAMMA

      LMAX=LL*(LL+2)


!   cos * cos  (1,1)

      DO 30  N1=1,2
      DO 30  L1=1,LL
      DO 30  M1=-L1,L1  
         IIC1=(N1-1)*LMAX + (L1**2+L1+M1)

      DO 20  N=1,2
      DO 20  L=1,LL
      DO 20  M=-L,L
         IIC=(N-1)*LMAX + (L**2+L+M)

!        write(*,*) 'Matrix_W: 30-20', IIC,IIC1


!      CALL WW_WW_VEC(L11,M1-1,N1,L22,M-1,N,   &
!             PPM,UPPM,NPP,HUPP,              &
!             NBETA,NGAMMA,NALPHA,            &
!             ALPHAM,BETAM,GAMMAM,NPAR,PAR,VLM,NCP1,NCP2)
!      AM(IIC,IIC1)=VLM(1)
         AM(IIC,IIC1)=1.

 20   CONTINUE
 30   CONTINUE
!c---------------------------------------------------

!   sin * sin  (2,2)

      DO 50  N1=1,2
      DO 50  L1=1,LL
      DO 50  M1=-L1,L1
         IIS1=2*LMAX+(N1-1)*LMAX + L1**2+L1+M1
      DO 40  N=1,2
      DO 40  L=1,LL
      DO 40  M=-L,L
         IIS=2*LMAX+(N-1)*LMAX + L**2+L+M

!      CALL WW_WW_VEC(L11,M1,N1,L22,M,N,   &
!             PPM,UPPM,NPP,HUPP,           &  
!             NBETA,NGAMMA,NALPHA,         &
!             ALPHAM,BETAM,GAMMAM,NPAR,PAR,VLM,NCP1,NCP2)
!        AM(IIS,IIS1)=VLM(4)
        AM(IIS,IIS1)=4. 
 40   CONTINUE
 50   CONTINUE
!-----------------------------------------------

!  cos * sin  (2,1)

      DO 70  N1=1,2
      DO 70  L1=1,LL
      DO 70  M1=-L1,L1
         IIC1=(N1-1)*LMAX + (L1**2+L1+M1)
      DO 60  N=1,2
      DO 60  L=1,LL
      DO 60  M=-L,L
         IIS=2*LMAX+(N-1)*LMAX + (L**2+L+M)

!      CALL WW_WW_VEC(L11,M1-1,N1,L22,M,N,  &
!             PPM,UPPM,NPP,HUPP,            &
!             NBETA,NGAMMA,NALPHA,          &
!             ALPHAM,BETAM,GAMMAM,NPAR,PAR,VLM,NCP1,NCP2)
!         AM(IIS,IIC1)=VLM(3)
         AM(IIS,IIC1)=3.
 60   CONTINUE
 70   CONTINUE
!------------------------------------------

!  sin * cos   (1,2)

      DO 90  N1=1,2
      DO 90  L1=1,LL
      DO 90  M1=-L1,L1
         IIS1=2*LMAX+(N1-1)*LMAX + (L1**2+L1+M1)

      DO 80  N=1,2
      DO 80  L=1,LL
      DO 80  M=-L,L
         IIC=(N-1)*LMAX + (L**2+L+M)

!      CALL WW_WW_VEC(L11,M1,N1,L22,M-1,N, &
!             PPM,UPPM,NPP,HUPP,           &  
!             NBETA,NGAMMA,NALPHA,         &  
!             ALPHAM,BETAM,GAMMAM,NPAR,PAR,VLM,NCP1,NCP2)
!         AM(IIC,IIS1)=VLM(2)
         AM(IIC,IIS1)=2.
 80   CONTINUE
 90   CONTINUE
 112  CONTINUE
      RETURN
      END SUBROUTINE Matrix_WW_VEC_1
!-----------------------------------------------------------


      SUBROUTINE INDUXIS(LMAX)
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        INTEGER :: LMAX
        INTEGER :: AM1(4*LMAX*(LMAX+2)),AM2(4*LMAX*(LMAX+2))
        INTEGER :: L,M,SIGMA,J
        INTEGER :: NU1,NU2
        J=0
        DO 4 SIGMA=1,2
        DO 4 L=1,LMAX   
        DO 4 M=-L,L   
           J=J+1
           NU1=(SIGMA-1)*LMAX*(LMAX+2)+L**2+L+M
           AM1(J)=NU1
           NU2=2*LMAX*(LMAX+2)+(SIGMA-1)*LMAX*(LMAX+2)+L**2+L+M
           AM2(J)=NU2
           WRITE(*,*) 'INDUXIS: AM1,AM2=', AM1(J),AM2(J)
4       CONTINUE

!      WRITE(*,*) 'INDUXIS: AM1,AM2=', AM1,AM2
      RETURN
      END SUBROUTINE INDUXIS  
!---------------------------------------------

      SUBROUTINE INDUXIS2(LMAX)
! Проверка индексов для матрицы в задаче 
! для Бессилового поля  
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        INTEGER :: LMAX,LLL
        INTEGER :: AM1(2*LMAX*(LMAX+2)),AM2(2*LMAX*(LMAX+2))
        INTEGER :: L,M,L1,M1,I,J
        INTEGER :: MU,MU1,NU,NU1,NU2
        INTEGER :: AM(2*LMAX*(LMAX+2),2*LMAX*(LMAX+2))

        LLL=2*LMAX*(LMAX+2)

        DO 6 L1=1,LMAX  
        DO 6 M1=-L1,L1
           NU=L1**2+L1+M1
           NU1=LMAX*(LMAX+2)+L1**2+L1+M1
        DO 6 L=1,LMAX  
        DO 6 M=-L,L
           MU=L**2+L+M
           MU1=LMAX*(LMAX+2)+L**2+L+M
           AM(MU,NU)=1
           AM(MU,NU1)=3
           AM(MU1,NU)=2
           AM(MU1,NU1)=4

!           WRITE(*,*) 'MU,NU,MU1,NU1', MU1 !,NU,MU1,NU1
6       CONTINUE 
           WRITE(*,*) 'INDUXIS2'
           WRITE(*,100), ((AM(I,J),I=1,LLL),J=1,LLL)

100   FORMAT(1X,16I3)
      RETURN
      END SUBROUTINE INDUXIS2  
!---------------------------------------------
      FUNCTION D_M1_M2(M1,M2)
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        INTEGER :: M1,M2
!local  
        INTEGER :: M11,M22        
! functions
        INTEGER :: D_M1_M2

        M11=ABS(M1)
        M22=ABS(M2)

        IF(M1 .GE. 0 .AND. M2 .GE. 0) THEN
           D_M1_M2=1
           RETURN
        ENDIF

        IF(M1 .GE. 0 .AND. M2 .LT. 0) THEN
           D_M1_M2=(-1)**M22
           RETURN
        ENDIF

        IF(M1 .LT. 0 .AND. M2 .GE. 0) THEN
           D_M1_M2=(-1)**M11
           RETURN
        ENDIF

        IF(M1 .LT. 0 .AND. M2 .LT. 0) THEN
           D_M1_M2=(-1)**(M11+M22)
           RETURN
        ENDIF
      RETURN
      END FUNCTION D_M1_M2
!-------------------------------------------

      SUBROUTINE M_LM_DECART_POINT (L,M,KK,XX,YY,ZZ,ALPHA,BETA,CLMM,MMRE,MMIM)
!  Vectors M_LM in Decart coordinate system at point (x,y,z)

!        USE PARAMTR
!        USE PARAM_INPUT

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

        INTEGER :: L,M
        REAL(wp) :: KK  
        REAL(wp) :: XX,YY,ZZ
        REAL(wp) :: ALPHA,BETA
        REAL(wp) :: MMRE(3),MMIM(3)
        COMPLEX(wp) :: CLMM(3)

!local 
        COMPLEX(wp) :: CLMM1(3),CLMM2(3)    
        INTEGER :: I
!        INTEGER :: II3,II2,II1
        REAL(wp) :: UU,VV,WW,SIGNW,TETA1,PHI1,RRR1 
        REAL(wp) :: RRR,TETA,PHI,SIGNZ
        REAL(wp) ::BB(3,3),BI(3,3),AA(3,3),BB1(3,3)
        REAL(wp) :: CALPHA,SALPHA,CBETA,SBETA
!---------------------------
!  M_LM  .NE. 0 if  L=1,2....
!----------------------------


!         SIGNZ=SIGN(1.0,ZZ)
!         RRR=SQRT(XX**2+YY**2+ZZ**2)

         CALL POLAR(XX,YY,ZZ,TETA,PHI)

!         IF(TETA .EQ. PI .OR. TETA .EQ. 0.) THEN
!            TETA=PI*(1.0+SIGNZ/360.)
!         ENDIF   

        write(*,*) '_POINT0000', TETA*180./PI,PHI*180./PI
! the matrix of transformation polar to decart 
!         CALL MATRIX_B(TETA,PHI,BB,BI)

!         CALL M_LM (L,M,KK,RRR,TETA,PHI,CLMM1)

! decart components in Lab. system
!         CLMM(1)= CLMM1(1)*BB(1,1)+CLMM1(2)*BB(1,2)+CLMM1(3)*BB(1,3) 
!         CLMM(2)= CLMM1(1)*BB(2,1)+CLMM1(2)*BB(2,2)+CLMM1(3)*BB(2,3)  
!         CLMM(3)= CLMM1(1)*BB(3,1)+CLMM1(2)*BB(3,2)+CLMM1(3)*BB(3,3)  
!--------------------------------------------------

         CALPHA = COS(ALPHA)
         SALPHA = SIN(ALPHA)
         CBETA  = COS(BETA)
         SBETA  = SIN(BETA)

!--  inverse  matrix AA(AL,BET,GAMM=0)----------------------------
         AA(1,1)= CALPHA*CBETA
         AA(1,2)= SALPHA*CBETA
         AA(1,3)=-SBETA
         AA(2,1)=-SALPHA
         AA(2,2)= CALPHA
         AA(2,3)= 0.0
         AA(3,1)= CALPHA*SBETA
         AA(3,2)= SALPHA*SBETA
         AA(3,3)= CBETA
!-----------------------------------------------------------------
!         UU= AA(1,1)*XX + AA(1,2)*YY + AA(1,3)*ZZ  
!         VV= AA(2,1)*XX + AA(2,2)*YY + AA(2,3)*ZZ
!         WW= AA(3,1)*XX + AA(3,2)*YY + AA(3,3)*ZZ

!  direct
         UU= AA(1,1)*XX + AA(2,1)*YY + AA(3,1)*ZZ  
         VV= AA(1,2)*XX + AA(2,2)*YY + AA(3,2)*ZZ
         WW= AA(1,3)*XX + AA(2,3)*YY + AA(3,3)*ZZ


!         write(*,*) '_POINT', UU,VV,WW
!         SIGNZ=SIGN(1.0,WW)
         RRR1=SQRT(UU**2+VV**2+WW**2)
         CALL POLAR(UU,VV,WW,TETA1,PHI1)

!         IF(TETA1 .EQ. PI .OR. TETA1 .EQ. 0.) THEN
!            TETA1=PI*(1.0+SIGNZ/360.)
!         ENDIF   

         write(*,*) '_POINT', TETA1*180./PI,PHI1*180./PI

! the matrix of transformation polar to decart 
         CALL MATRIX_B(TETA1,PHI1,BB1,BI)

         CALL M_LM (L,M,KK,RRR1,TETA1,PHI1,CLMM2)

! decart components rotated in S^prime system
         CLMM(1)= CLMM2(1)*BB1(1,1)+CLMM2(2)*BB1(1,2)+CLMM2(3)*BB1(1,3) 
         CLMM(2)= CLMM2(1)*BB1(2,1)+CLMM2(2)*BB1(2,2)+CLMM2(3)*BB1(2,3)  
         CLMM(3)= CLMM2(1)*BB1(3,1)+CLMM2(2)*BB1(3,2)+CLMM2(3)*BB1(3,3)  
!--------------------------------------------------

         DO 4 I =1,3
            MMRE(I)=REAL(CLMM(I))
            MMIM(I)=AIMAG(CLMM(I))
4        CONTINUE
6     CONTINUE

      RETURN
      END SUBROUTINE M_LM_DECART_POINT
!-------------------------------------------------

      SUBROUTINE ProjExact_M_LM_1(L,M,KK,NB,ERR,ALPHA,BETA,    &
                                PPM,NPP,PSIM,NPSI,BMODRE1,BMODIM1)
!  Exact Ray transform  of the vectors M_LM 

!      MODULE PARAMTR
      IMPLICIT NONE
      INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
      INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single 
      REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

      INTEGER :: L,M,NPP,NPSI
      INTEGER :: NB,ERR 
      REAL(wp) :: KK  !,BES(NB)  
      REAL(wp) :: ALPHA,BETA
      REAL(wp) :: PPM(NPP),PSIM(NPSI)
      REAL(wp) :: BMODRE1(NPP*NPSI),BMODIM1(NPP*NPSI)
!      REAL(wp), PARAMETER :: AL=0.0_wp  ! for 
! local
      INTEGER :: I,M1,M11,I1,I2,I21,MM
      REAL(wp) :: PNRM,PLG,PLM1,PLM2
      REAL(wp) :: X,RR,PSI,BES,ALMTET
      COMPLEX(wp) :: DWIG,SS
      COMPLEX(wp) :: BMODC1(NPP*NPSI)
      COMPLEX(wp),PARAMETER :: ci = (0.,1.)
! functions
!      REAL(wp) ::PNORM, PLGNDR, BessJ,PLM1M2
      REAL(wp) ::BessJ,PLM1M2,PNORM,PLGNDR
      COMPLEX(wp) :: CCEXP

! complex exp definition
      CCEXP(X)=CMPLX(COS(X),SIN(X))
      MM=ABS(M)

      DO 8 I2=1,NPSI
         PSI=PSIM(I2)
      DO 8 I1=1,NPP
         RR=PPM(I1)
         I21=(I2-1)*NPP+I1

         SS=(0.0,0.0)
         DO 6 M1=1,L
!            M1=ABS(I)
!            M11=I
!            PNRM=PNORM(L,M1)
            PLG= PLGNDR(L,M1,0.0)*PNORM(L,M1)
            BES= BessJ(KK*RR,M1,NB,ERR)/KK
            ALMTET= M1 * PLG * BES

            PLM1=(-1)**(L+M+1)*PLM1M2(L,M,M1,PI-BETA)*CCEXP(-M1*PSI)
            PLM2=PLM1M2(L,M,M1,BETA)*CCEXP(M1*PSI)

            SS=SS + ci**M1 * (PLM1+PLM2)*ALMTET

6        CONTINUE

!         SS=  CCEXP(M*ALPHA) * SS
!         write(*,*) 'ProExact_M_LM: SS=',SS

!         BMODC1(I21)=ci**(L+1) * SQRT(PI/2.) * CCEXP(M*ALPHA) * SS
         BMODC1(I21)=(-1)**M*ci**L * SQRT(PI/2.) * CCEXP(M*ALPHA) * SS
         BMODRE1(I21)=REAL(BMODC1(I21))/ SQRT(L*(L+1.0))
         BMODIM1(I21)=AIMAG(BMODC1(I21))/ SQRT(L*(L+1.0))

8     CONTINUE         

      RETURN
      END SUBROUTINE ProjExact_M_LM_1
!------------------------------------------------------------

      SUBROUTINE ProjExact_M_LM(L,M,KK,NB,ERR,ALPHA,BETA,    &
                                PPM,NPP,PSIM,NPSI,BMODRE1,BMODIM1)
!  Exact Ray transform  of the vectors M_LM 

!      MODULE PARAMTR
      IMPLICIT NONE
      INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
      INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single 
      REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

      INTEGER :: L,M,NPP,NPSI
      INTEGER :: NB,ERR 
      REAL(wp) :: KK  !,BES(NB)  
      REAL(wp) :: ALPHA,BETA,GAMMA
      REAL(wp) :: PPM(NPP),PSIM(NPSI)
      REAL(wp) :: BMODRE1(NPP*NPSI),BMODIM1(NPP*NPSI)
!      REAL(wp), PARAMETER :: AL=0.0_wp  ! for 
! local
      INTEGER :: I,M1,M11,I1,I2,I21,LL
      REAL(wp) :: PNRM,PLG,PLAG
      REAL(wp) :: X,RR,PSI,BES,ALMTET
      REAL(wp) :: RL,RM,RM11
      COMPLEX(wp) :: DWIG,SS
      COMPLEX(wp) :: BMODC1(NPP*NPSI)
      COMPLEX(wp),PARAMETER :: ci = (0.,1.)
! functions
      REAL(wp) ::djmn
      REAL(wp) ::BessJ,PLM1M2,PNORM,PLGNDR,PLGNDR_NORM
      COMPLEX(wp) :: CCEXP

! complex exp definition
      CCEXP(X)=CMPLX(COS(X),SIN(X))

!      write(*,*) ' ProjExact_M_LM', ALPHA,BETA
      RL=L
      RM=M
      GAMMA=0.  !PI/3.

      DO 8 I2=1,NPSI
         PSI=PSIM(I2)
      DO 8 I1=1,NPP
         RR=PPM(I1)
         I21=(I2-1)*NPP+I1

         SS=(0.0,0.0)
         DO 6 I=-L,L
            M1=ABS(I)
            M11=I
            RM11=M11
            LL=MOD(L+M1,2)
! при нечётных L+M1  P_L^M в нуле равны 0 
            IF(LL .EQ. 0) THEN
!               PNRM=PNORM(L,M1)
!               PLAG=PLGNDR(L,M11,0.0)*PNRM
               PLG= PLGNDR_NORM(L,M11,0.0)

               ALMTET= M11 * PLG     
!               DWIG= djmn(RL,RM,RM11,BETA) * CCEXP(M11*PSI+M*ALPHA+(M11+L)*PI/2.-PI/2.)  
               DWIG= PLM1M2(L,M,M11,BETA) * CCEXP(M11*PSI+M*ALPHA+(M11+GAMMA-L)*PI/2.+PI/2.)  
               BES= BessJ(KK*RR,M11,NB,ERR)/KK
               SS=SS + DWIG * ALMTET * BES 
            ENDIF
6        CONTINUE
!         BMODC1(I21)=-SS /SQRT(2*PI)  !*(L*(L+1))
         BMODC1(I21)=-SS /SQRT(8*PI)
         BMODRE1(I21)=REAL(BMODC1(I21))    
         BMODIM1(I21)=AIMAG(BMODC1(I21))   

8     CONTINUE         

      RETURN
      END SUBROUTINE ProjExact_M_LM
!------------------------------------------------------------

      SUBROUTINE ProjExact_M_LM_2(L,M,KK,NB,ERR,ALPHA,BETA,    &
                                PPM,NPP,PSIM,NPSI,BMODRE1,BMODIM1)
!  Exact Ray transform  of the vectors M_LM 

!      MODULE PARAMTR
      IMPLICIT NONE
      INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
      INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single 
      REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

      INTEGER :: L,M,NPP,NPSI
      INTEGER :: NB,ERR 
      REAL(wp) :: KK  !,BES(NB)  
      REAL(wp) :: ALPHA,BETA
      REAL(wp) :: PPM(NPP),PSIM(NPSI)
      REAL(wp) :: BMODRE1(NPP*NPSI),BMODIM1(NPP*NPSI)
!      REAL(wp), PARAMETER :: AL=0.0_wp  ! for 
! local
      INTEGER :: I,MM,M1,M11,I1,I2,I21,K,J,LL
      REAL(wp) :: PNRM,PLG
      REAL(wp) :: X,RR,PSI,BES,ALMTET
      REAL(wp) :: MALPHA,MPSI
      REAL(wp) :: DWIG,SS,SS1
!      COMPLEX(wp) :: BMODC1(NPP*NPSI)
!      COMPLEX(wp),PARAMETER :: ci = (0.,1.)
! functions
!      REAL(wp) ::PNORM, PLGNDR, BessJ,PLM1M2
      REAL(wp) ::BessJ,PLM1M2,PLGNDR_NORM,IMprimeRO !, PNORM

      COMPLEX(wp) :: CCEXP

! complex exp definition
!      CCEXP(X)=CMPLX(COS(X),SIN(X))


      MALPHA=M*ALPHA+L*PI/2.

!---------------------------

      DO 8 I2=1,NPSI
         PSI=PSIM(I2)
      DO 8 I1=1,NPP
         RR=PPM(I1)
         I21=(I2-1)*NPP+I1

         SS=0.0
         DO 6 I=-L,L
            M1=ABS(I)
            M11=I
            MPSI=M11*(PSI+PI/2.)
            LL=MOD(L+M1,2)
            IF(LL .EQ. 0) THEN
!               PNRM=PNORM(L,M1)                 !NOR
               PLG= PLGNDR_NORM(L,M11,0.0)
               ALMTET= M11 * PLG 

               DWIG= PLM1M2(L,M,M11,BETA) * COS(MALPHA+MPSI) * IMprimeRO(M11,RR)   
               SS=SS + DWIG * ALMTET
            ENDIF 
6        CONTINUE

         BMODRE1(I21)=-SS/SQRT(8*PI)
!         BMODIM1(I21)=AIMAG(BMODC1(I21))/ SQRT(L*(L+1.0)) 
8     CONTINUE         

      RETURN
      END SUBROUTINE ProjExact_M_LM_2
!------------------------------------------------------------
      FUNCTION IMprimeRO(M,RR)
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
!        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

        INTEGER :: M
        REAL(wp) :: RR
        REAL(wp) :: IMprimeRO
! local 
        INTEGER :: LL
        REAL(wp) :: SS

        LL=MOD(M,2)
        IF(LL .EQ. 0) THEN
           SS=1.
        ELSE   
           IF(M .LT. 0) THEN 
              SS=-1.
           ENDIF
        ENDIF   

!        IMprimeRO=SS/RR
        IMprimeRO=SS
      RETURN
      END FUNCTION IMprimeRO
!------------------------------------------------------------
      SUBROUTINE Y0LM (L,M,TETA,PHI,CLMM )
!  Vector spherical harmonics   Y^(0)_(LM) 

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

        INTEGER :: L,M
        REAL(wp) :: TETA,PHI
! vector spherical harmonics
        COMPLEX(wp) :: CLMM(3)

!local 
        COMPLEX(wp),PARAMETER :: ci = (0.,1.)
        COMPLEX(wp) :: CKLM,YLM    !  ,CYY,SS
        REAL(wp) :: COST, SINT, COSF, SINF
        REAL(wp) :: R3,SIGNUM
        INTEGER :: MM,JJ,I
! spherical unit vectors
        REAL(wp) :: AM(3,3)
!        REAL(wp) :: rUnit(3), tetaUnit(3), fiUnit(3)
! functions
      COMPLEX(wp) :: CYLM

!=======================
! initialization 
!=======================
      MM=ABS(M)

      COST=COS(TETA)
      SINT=SIN(TETA)
      COSF=COS(PHI)
      SINF=SIN(PHI)

!=======================
! matrix  M(x,y,z <-- r,\theta,\phi)
!=======================

      AM(1,1) = COSF *SINT
      AM(2,1) = SINF *SINT
      AM(3,1) = COST
      
      AM(1,2) = COSF * COST
      AM(2,2) = SINF * COST
      AM(3,2) =-SINT

      AM(1,3) =-SINF
      AM(2,3) = COSF
      AM(3,3) = 0. 

!=======================
! unit vectors
!=======================

!      rUnit(1) =     COSF *SINT
!      rUnit(2) =     SINF *SINT
!      rUnit(3) =           COST
!      tetaUnit(1) = COSF * COST
!      tetaUnit(2) = SINF * COST
!      tetaUnit(3) =       -SINT
!      fiUnit(1)  =  -SINF
!      fiUnit(2)  =   COSF
!      fiUnit(3)  =   0. 
!-------------------------------------
!      write(*,*) 'Y0LM: J,M1',L,M

!  derivatives over teta
     CALL DYLM_TETA (L,MM,TETA,PHI,CKLM)
!------------------------------------
!  spherical harmonics Y(LM)

      YLM= CYLM (L,MM,TETA,PHI)

!-------------------------------------
         IF(L .EQ. 0) THEN
            PRINT*, 'Error,Y0LM: L must not be = 0! STOP'
            STOP
         ENDIF   

         R3=1./SQRT(FLOAT(L*(L+1)))

!
         IF(TETA .EQ. PI .OR. TETA .EQ. 0.) THEN
            TETA= PI/360.   !  CLMM(2)=(0.0,0.0) 
         ENDIF

         CLMM(1)=(0.0,0.0)                    ! radial component
         CLMM(2)=-R3*YLM *FLOAT(MM)/SIN(TETA)    !teta Component
         CLMM(3)=-ci*R3*CKLM                 !fi Componnent

!         write(*,*) 'Y0LM: CLMM(3)=',CLMM(3)

!   Y^(0)(LM) coincides with X(LM)
!   article B.C. Brock (002217r.pdf)
!            CLMM(I) = R3*(-tetaUnit(I)*YLM *FLOAT(MM)/SINT   &
!                    -fiUnit(I)*CKLM*ci)
! Varshalovich (english) p.215 (39)
! Varshalovich (russian) p.188 (39)
         IF( M .LT. 0) THEN 
               JJ=MOD(MM+1,2)
               IF(JJ .EQ. 0) THEN 
                  SIGNUM=1.
               ELSE
                  SIGNUM=-1.
               ENDIF
            DO 10 I=1,3
! this is harmonics Y^0(LM)
               CLMM(I)= SIGNUM*CONJG(CLMM(I))   
! this is harmonic C(LM) (Hansen, Morse - Feshbach) 
!               CLMM(I)=-ci* SIGNUM*CONJG(CLMM(I))    
 10         CONTINUE   
         ENDIF
      RETURN
      END SUBROUTINE Y0LM
!-----------------------------------------------------

      SUBROUTINE Rotation_YLM
! Rotation of the  Scalar Spherical harmonics 

!        USE PARAMTR
!        USE PARAM_INPUT

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

        INTEGER :: L,M
        REAL(wp) :: ALPHA,BETA 

!local 
        COMPLEX(wp) :: CLMM(3)    
        INTEGER :: J321,I,J
!        INTEGER :: II3,II2,II1
        INTEGER :: M1,I1
        REAL(wp) :: XX,YY,ZZ,UU,VV,WW
        REAL(wp) :: CALPHA,SALPHA,CBETA,SBETA 
        REAL(wp) :: RRR,TET,PHI,SIGNW,SIGNZ
        REAL(wp) :: RRR1,TET1,PHI1,X
        COMPLEX(wp) :: CYY,CYY1,CYY2,SS,SS1
        COMPLEX(wp) :: MLMM1(3),MLMM(3)
        REAL(wp) :: AA(3,3)
        REAL(wp) :: RL,RM,RM1
        COMPLEX(wp),PARAMETER :: ci = (0.,1.)
! functions
        COMPLEX(wp) :: D_WIGNER,CCEXP,CYLM
        REAL(wp) :: PLM1M2,djmn
! complex exp definition
        CCEXP(X)=CMPLX(COS(X),SIN(X))

        L=3
        M=-2
        RL=L
        RM=M
        TET1=PI/4.
        PHI1=PI/4.
        XX=SIN(TET1)*COS(PHI1)
        YY=SIN(TET1)*SIN(PHI1)
        ZZ=COS(TET1)
        
        ALPHA=0.0
        BETA=PI/4.0
!        GAMMA=0.0
        
!         SIGNZ=SIGN(1.0,ZZ)
!         RRR1=SQRT(XX**2+YY**2+ZZ**2)

!         CALL POLAR(XX,YY,ZZ,TETA1,PHI1)
!         write(*,*) 'POLAR111111', PHI1*180/PI, TETA1*180/PI
         CYY= CYLM (L,M,TET1,PHI1)
!         CALL Y0LM (L,M,TETA1,PHI1,MLMM1)
!         CALL M_LM (L,M,KK,RRR1,TETA1,PHI1,MLMM1)
!---------------------------------------------------
         CALPHA = COS(ALPHA)
         SALPHA = SIN(ALPHA)
         CBETA  = COS(BETA)
         SBETA  = SIN(BETA)

!--  direct  matrix AA(AL,BET,GAMM=0)----------------------------
         AA(1,1)= CALPHA*CBETA
         AA(1,2)= SALPHA*CBETA
         AA(1,3)=-SBETA
         AA(2,1)=-SALPHA
         AA(2,2)= CALPHA
         AA(2,3)= 0.0
         AA(3,1)= CALPHA*SBETA
         AA(3,2)= SALPHA*SBETA
         AA(3,3)= CBETA

!         UU= AA(1,1)*XX + AA(1,2)*YY + AA(1,3)*ZZ  
!         VV= AA(2,1)*XX + AA(2,2)*YY + AA(2,3)*ZZ
!         WW= AA(3,1)*XX + AA(3,2)*YY + AA(3,3)*ZZ

         UU= AA(1,1)*XX + AA(2,1)*YY + AA(3,1)*ZZ  
         VV= AA(1,2)*XX + AA(2,2)*YY + AA(3,2)*ZZ
         WW= AA(1,3)*XX + AA(2,3)*YY + AA(3,3)*ZZ

         CALL POLAR(UU,VV,WW,TET,PHI)
!         IF(TET .EQ. PI .OR. TET .EQ. 0.) THEN
!            TET=PI*(1.0+SIGNW/180.)
!         ENDIF   
      SS=(0.0,0.0)   
      DO 12 M1=-L,L
         RM1=M1
         CYY1=CYLM (L,M1,TET,PHI)
         CYY2=PLM1M2(L,M1,M,BETA) * CCEXP(-M1*ALPHA)
!         CYY2= djmn (RL,RM,RM1,BETA)*CCEXP(-M*ALPHA)
!         CYY2= D_WIGNER(L,M,M1,ALPHA,BETA,0.0)
!         CALL Y0LM (L,M1,TETA,PHI,CLMM)
!         CALL M_LM (L,M1,KK,RRR,TETA,PHI,CLMM)
         SS=SS+CYY1*CYY2
12    CONTINUE
!         SS1=1.+SIN(BETA)/SQRT(2.)
         write(*,*) 'Rotation:',CYY,SS
!         write(*,*) 'Rotation:',SS,SS1
            
      RETURN
      END SUBROUTINE Rotation_YLM
!-------------------------------------------------

      SUBROUTINE Rotation_M_LM_Direct
! Rotation of the  M_LM vectors

!        USE PARAMTR
!        USE PARAM_INPUT

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

        INTEGER :: L,M
        INTEGER :: NX,NY,NZ
        REAL(wp) :: ALPHA,BETA  !,GAMMA
!        REAL(wp) :: XM(NX),YM(NY),ZM(NZ)
!        REAL(wp) :: MLMRE(NX*NY*NZ,3),MLMIM(NX*NY*NZ,3)

!local 
        COMPLEX(wp) :: CLMM(3),CLMM1(3)     
        INTEGER :: I,J
        INTEGER :: M1,I1,I2
        REAL(wp) :: XX,YY,ZZ,UU,VV,WW
        REAL(wp) :: CALPHA,SALPHA,CBETA,SBETA 
        REAL(wp) :: CTET,STET,CPHI,SPHI 
        REAL(wp) :: CTET1,STET1,CPHI1,SPHI1 
        REAL(wp) :: RRR,TET,PHI,SIGNW,SIGNZ
        REAL(wp) :: TET1,PHI1,X
        COMPLEX(wp) :: CYY
        COMPLEX(wp) :: MLMM1(3),MLMM(3),MLMM2(3),MLMM3(3),MLMM4(3)
        REAL(wp) :: AA(3,3),BB1(3,3),BB(3,3)
        REAL(wp) :: KK,RR
! functions
        COMPLEX(wp) :: CCEXP  !,CYLM,D_WIGNER
        REAL(wp) :: PLM1M2
! complex exp definition
        CCEXP(X)=CMPLX(COS(X),SIN(X))

!------------------------------------------------------
! ALPHA must not be  = PHI1  !!!!!!!!!!!
!-----------------------------------------------------

        L=1
        M=0

        KK=10.0
        RR=1.0

        ALPHA=PI/3.1    !ALPHA must not be  = PHI1  !!!!!!!!!!!
        BETA=PI/3.0

        TET1=PI/3.
        PHI1=PI/3.
        
         CTET1 = COS(TET1)
         STET1 = SIN(TET1)
         CPHI1 = COS(PHI1)
         SPHI1 = SIN(PHI1)

         BB1(1,1)= STET1*CPHI1
         BB1(1,2)= CTET1*CPHI1
         BB1(1,3)=-SPHI1
         BB1(2,1)= STET1*SPHI1
         BB1(2,2)= CTET1*SPHI1
         BB1(2,3)= CPHI1
         BB1(3,1)= CTET1
         BB1(3,2)=-STET1
         BB1(3,3)= 0.

!--------------------------------------------
!computer the right part
!------------------------------------------
      DO 2 J=1,3   
         MLMM(J)=(0.0,0.0)
2     CONTINUE

      DO 12 M1=-L,L

!         CALL Y0LM(L,M1,TET1,PHI1,CLMM)
         CALL M_LM (L,M1,KK,RR,TET1,PHI1,CLMM )
         CYY=PLM1M2(L,M1,M,BETA) * CCEXP(-M1*ALPHA)

         DO 4 J =1,3
            MLMM(J)= MLMM(J) + CLMM(J)*CYY
4        CONTINUE

12    CONTINUE
! decart components
         MLMM1(1)= MLMM(1)*BB1(1,1)+MLMM(2)*BB1(1,2)+MLMM(3)*BB1(1,3) 
         MLMM1(2)= MLMM(1)*BB1(2,1)+MLMM(2)*BB1(2,2)+MLMM(3)*BB1(2,3)  
         MLMM1(3)= MLMM(1)*BB1(3,1)+MLMM(2)*BB1(3,2)+MLMM(3)*BB1(3,3)  

!-----------------------------------------------
! compute the left  part
!-----------------------------------------------
        XX=SIN(TET1)*COS(PHI1)
        YY=SIN(TET1)*SIN(PHI1)
        ZZ=COS(TET1)
        
        CALPHA = COS(ALPHA)
        SALPHA = SIN(ALPHA)
        CBETA  = COS(BETA)
        SBETA  = SIN(BETA)

!--  inverse  matrix AA(AL,BET,GAMM=0)----------------------------
         AA(1,1)= CALPHA*CBETA
         AA(1,2)= SALPHA*CBETA
         AA(1,3)=-SBETA
         AA(2,1)=-SALPHA
         AA(2,2)= CALPHA
         AA(2,3)= 0.0
         AA(3,1)= CALPHA*SBETA
         AA(3,2)= SALPHA*SBETA
         AA(3,3)= CBETA

         UU= AA(1,1)*XX + AA(1,2)*YY + AA(1,3)*ZZ  
         VV= AA(2,1)*XX + AA(2,2)*YY + AA(2,3)*ZZ
         WW= AA(3,1)*XX + AA(3,2)*YY + AA(3,3)*ZZ
!  direct
!         UU= AA(1,1)*XX + AA(2,1)*YY + AA(3,1)*ZZ  
!         VV= AA(1,2)*XX + AA(2,2)*YY + AA(3,2)*ZZ
!         WW= AA(1,3)*XX + AA(2,3)*YY + AA(3,3)*ZZ

         CALL POLAR(UU,VV,WW,TET,PHI)
!         write(*,*) 'TET_PHI', TET*180./PI,PHI*180./PI

!         IF(TET .EQ. PI .OR. TET .EQ. 0.) THEN
!            TET=TET+PI/180.
!         ENDIF   

         CTET = COS(TET)
         STET = SIN(TET)
         CPHI = COS(PHI)
         SPHI = SIN(PHI)

         BB(1,1)= STET*CPHI
         BB(1,2)= CTET*CPHI
         BB(1,3)=-SPHI
         BB(2,1)= STET*SPHI
         BB(2,2)= CTET*SPHI
         BB(2,3)= CPHI
         BB(3,1)= CTET
         BB(3,2)=-STET
         BB(3,3)= 0.

!  M_LM  .NE. 0 if  L=1,2....
!----------------------------

!        CALL Y0LM(L,M,TET,PHI,MLMM2)
        CALL M_LM (L,M,KK,RR,TET,PHI,MLMM2 )
! decart components
        MLMM3(1)= MLMM2(1)*BB(1,1)+MLMM2(2)*BB(1,2)+MLMM2(3)*BB(1,3) 
        MLMM3(2)= MLMM2(1)*BB(2,1)+MLMM2(2)*BB(2,2)+MLMM2(3)*BB(2,3)  
        MLMM3(3)= MLMM2(1)*BB(3,1)+MLMM2(2)*BB(3,2)+MLMM2(3)*BB(3,3)  

! component transformation
      MLMM4(1)= AA(1,1)*MLMM3(1) + AA(2,1)*MLMM3(2) + AA(3,1)*MLMM3(3)  
      MLMM4(2)= AA(1,2)*MLMM3(1) + AA(2,2)*MLMM3(2) + AA(3,2)*MLMM3(3)
      MLMM4(3)= AA(1,3)*MLMM3(1) + AA(2,3)*MLMM3(2) + AA(3,3)*MLMM3(3) 

!---------------------------------------------------

         DO 20 J=1,3
            write(*,*) 'Rotation:',MLMM1(J), MLMM4(J)
!            write(*,*) 'Rotation:',MLMM(J),MLMM2(J)
20       CONTINUE
            
      RETURN
      END SUBROUTINE Rotation_M_LM_Direct
!-------------------------------------------------

      SUBROUTINE Fourier_Exact_M_LM (L,M,KK,XNU,YNU,ZNU,NX,NY,NZ,DGRE,DGIM)
!  
        USE FREQUENCY
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

        INTEGER :: L,M
        INTEGER :: NX,NY,NZ
        REAL(wp) :: KK
        REAL(wp) :: XNU(NX),YNU(NY),ZNU(NZ)
        REAL(wp) :: DGRE(NX*NY*NZ,3),DGIM(NX*NY*NZ,3)
!local 
        INTEGER :: I,II1,II2,II3,J321
        REAL(wp) :: RRR,TETA,PHI,SIGNZ
        REAL(wp) :: XX,YY,ZZ,NOR
        REAL(wp) :: SNAIX,SNAIY,SNAIZ
        COMPLEX(wp) :: MLMM(3),MLMM1(3)
!        COMPLEX(wp) :: CLMM(3)
        COMPLEX(wp),PARAMETER :: ci = (0.,1.)
        
        REAL(wp) :: DeltaF

!        KK= FREQ(1)  ! SNAIX  
        SNAIX=FREQ(1) 
!        SNAIY=FREQ(2) 
!        SNAIZ=FREQ(3) 
        
        NOR=SQRT(FLOAT(L*(L+1)))/(4.0*PI)/KK**2

!        write(*,*) 'Fourier_Exact_M_LM', NX,NY,NZ

!        KK=5.0
        DO 6 II3=1,NZ
           ZZ=ZNU(II3)  
        DO 6 II2=1,NY
           YY=YNU(II2)  
        DO 6 II1=1,NX
           XX=XNU(II1)  

           J321=(II3-1)*NX*NY+(II2-1)*NX +II1

           SIGNZ=SIGN(1.0,ZZ)
           RRR=SQRT(XX**2+YY**2+ZZ**2)

           CALL POLAR(XX,YY,ZZ,TETA,PHI)

           IF(TETA .EQ. 0. .OR. TETA .EQ. PI) THEN
              TETA=PI*(1.0+SIGNZ/180.)
           ENDIF

           CALL  Y0LM (L,M,TETA,PHI,MLMM) 
  
        DO 4 I=1,3   

           MLMM(I)=MLMM(I)*ci**(L-1)
           IF(RRR .LT. SNAIX) THEN
              DGRE(J321,I)=REAL(MLMM(I)) * NOR * DeltaF(RRR-KK)
              DGIM(J321,I)=AIMAG(MLMM(I)) * NOR * DeltaF(RRR-KK)
           ELSE
              DGRE(J321,I)=0.0
              DGIM(J321,I)=0.0
           ENDIF
   
4       CONTINUE
6       CONTINUE
        RETURN
        END SUBROUTINE Fourier_Exact_M_LM
!---------------------------------------------------------------
        FUNCTION DeltaF(X)
          IMPLICIT NONE
          INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
          INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single 
          REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

          REAL(wp) :: X,EPS
! function
          REAL(wp) :: DeltaF
          
          EPS=0.54

          DeltaF=1.0/(EPS*SQRT(PI))*EXP(-(X/EPS)**2)

        RETURN
        END FUNCTION DeltaF
!--------------------------------------------------------
      SUBROUTINE FFT3D_VEC(DGMRE,DGMIM,M1,M2,M3,N1,N2,N3,H1,H2,H3,ISIGNUM,IERR)
!..... three dimensional fourier transform
!..... FRM ,FIM   - input massives
!......QR, QI   - worker massives
!......N1=2**M1+1 - the number points on first coordinate
!......N2=2**M2+1 - the number points on second coordinate
!......N3=2**M3+1 - the number points on third coordinate
!......H1 , H2,H3 - steps of the 
!......SIGNUM=-1  - direct transform
!......SIGNUM=+1  - invers transform                

      IMPLICIT NONE
      INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
      INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single 
      REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

      INTEGER :: M1,M2,M3,N1,N2,N3,IERR,ISIGNUM
      REAL(wp) ::  DGMRE(N1*N2*N3,3), DGMIM(N1*N2*N3,3)
      REAL(wp) :: H1,H2,H3
! local
      REAL(wp) :: FRM(N1*N2*N3),FIM(N1*N2*N3)
      INTEGER :: MM(3)
      INTEGER :: I1,I2,I3,I21,I321,I4
      INTEGER :: N11,N22,N33 
      REAL(wp) :: HN1,HN2,HN3,DELTA
      CHARACTER*8 ITEX
      ITEX='FFT3D_VEC'

          DO 14 I4=1,3
             DO 8 I3=1,N3
             DO 8 I2=1,N2
             DO 8 I1=1,N1
                I21=(I2-1)*N1+I1
                I321=(I3-1)*N1*N2+I21
                FRM(I321)=DGMRE(I321,I4)
                FIM(I321)=DGMIM(I321,I4)
 8           CONTINUE   
!             ISIGNUM=-1
             CALL FFT3D(FRM,FIM,M1,M2,M3,N1,N2,N3,H1,H2,H3,ISIGNUM,IERR)
             DO 10 I3=1,N3
             DO 10 I2=1,N2
             DO 10 I1=1,N1
                I21=(I2-1)*N1+I1
                I321=(I3-1)*N1*N2+I21
                DGMRE(I321,I4)= FRM(I321)
                DGMIM(I321,I4)= FIM(I321)
10           CONTINUE   
14        CONTINUE
      RETURN
      END SUBROUTINE FFT3D_VEC
!--------------------------------------------------------------

        SUBROUTINE MATRIX_B(TET,PHI,BD,BI)
! BD(3,3) direct:   (r,tet,phi) => (x,y,z)
! BI(3,3) inverse: (x,y,z) => (r,tet,phi) 
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single
      
        REAL(wp) :: BD(3,3),BI(3,3)
        REAL(wp) :: TET,PHI
! local
        REAL(wp) :: CTET,STET,CPHI,SPHI 
           
         CTET = COS(TET)
         STET = SIN(TET)
         CPHI = COS(PHI)
         SPHI = SIN(PHI)
! direct matrix
         BD(1,1)= STET*CPHI
         BD(1,2)= CTET*CPHI
         BD(1,3)=-SPHI
         BD(2,1)= STET*SPHI
         BD(2,2)= CTET*SPHI
         BD(2,3)= CPHI
         BD(3,1)= CTET
         BD(3,2)=-STET
         BD(3,3)= 0.

! inverse  matrix
         BI(1,1)= STET*CPHI
         BI(1,2)= STET*SPHI
         BI(1,3)= CTET
         BI(2,1)= CTET*CPHI
         BI(2,2)= CTET*SPHI
         BI(2,3)=-STET
         BI(3,1)=-SPHI
         BI(3,2)= CPHI
         BI(3,3)= 0.

        RETURN
        END SUBROUTINE MATRIX_B
!----------------------------------------------------------

      SUBROUTINE M_LM_DECART (L,M,KK,XM,YM,ZM,NX,NY,NZ,MLMM,MLMRE,MLMIM)
!  Vectors M_LM in Decart coordinate system

!        USE PARAMTR
!        USE PARAM_INPUT

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

        INTEGER :: L,M
        INTEGER :: NX,NY,NZ
        REAL(wp) :: KK  
        REAL(wp) :: XM(NX),YM(NY),ZM(NZ)
        REAL(wp) :: MLMRE(NX*NY*NZ,3),MLMIM(NX*NY*NZ,3)
! vector spherical harmonics
        COMPLEX(wp) :: MLMM(NX*NY*NZ,3)

!local 
        COMPLEX(wp) :: CLMM(3),CLMM1(3)    
        INTEGER :: J321,I
        INTEGER :: II3,II2,II1
        REAL(wp) :: XX,YY,ZZ
        REAL(wp) :: RRR,TETA,PHI,SIGNZ
        REAL(wp) ::BB(3,3),BI(3,3)
!---------------------------
!  M_LM  .NE. 0 if  L=1,2....
!----------------------------

      DO 6 II3=1,NZ
         ZZ=ZM(II3)
      DO 6 II2=1,NY
         YY=YM(II2)
      DO 6 II1=1,NX
         XX=XM(II1)

         J321=(II3-1)*NX*NY+(II2-1)*NX +II1

         SIGNZ=SIGN(1.0,ZZ)
         RRR=SQRT(XX**2+YY**2+ZZ**2)

         CALL POLAR(XX,YY,ZZ,TETA,PHI)

!         IF(TETA .EQ. PI .OR. TETA .EQ. 0.) THEN
!            TETA=PI*(1.0+SIGNZ/360.)
!         ENDIF   
! the matrix of transformation polar to decart 
         CALL MATRIX_B(TETA,PHI,BB,BI)

         CALL M_LM (L,M,KK,RRR,TETA,PHI,CLMM)

!         write(*,*) 'M_LM_DECART:M_LMM1=', CLMM


! decart components
         CLMM1(1)= CLMM(1)*BB(1,1)+CLMM(2)*BB(1,2)+CLMM(3)*BB(1,3) 
         CLMM1(2)= CLMM(1)*BB(2,1)+CLMM(2)*BB(2,2)+CLMM(3)*BB(2,3)  
         CLMM1(3)= CLMM(1)*BB(3,1)+CLMM(2)*BB(3,2)+CLMM(3)*BB(3,3)  

!         write(*,*) 'M_LM_DECART: CLMM1:', CLMM1


         DO 4 I =1,3
            MLMM(J321,I)=CLMM1(I)
            MLMRE(J321,I)=REAL(MLMM(J321,I))
            MLMIM(J321,I)=AIMAG(MLMM(J321,I))
4        CONTINUE
6     CONTINUE

      RETURN
      END SUBROUTINE M_LM_DECART
!-------------------------------------------------

      SUBROUTINE POLAR_DECART2D(GM,XXM,YYM,NX,NY,RRM,PHIM,NRR,NPHI,DR)
!
        USE PARAMTR
          IMPLICIT NONE
          INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
          INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
          REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

          INTEGER :: NX,NY,NRR,NPHI
          REAL(wp) :: GM(NRR*NPHI)
          REAL(wp) :: XXM(NX),YYM(NY) !,ZZM(NZ)
          REAL(wp) :: RRM(NRR),PHIM(NPHI)
          REAL(wp) :: DR(NX*NY)
          REAL(wp) :: FUNL3,ALINEAR_3D,Burke_2D
! local
!          REAL(wp) :: GMM(NRR*NPHI)
          REAL(wp) :: XX,YY   !,ZZ,RXY                  !!!!!!!!!
          REAL(wp) :: RR2,PHI   !,HPP,HPHI  !,HTET     !!!!!!!!!!!!!!!
          INTEGER :: I,I1,I2,I3,I21,IERR,EI
! functions 
          REAL(wp) :: FUNL2

!          HPP =PAR(35)
!          HTET=PAR(15)
!          HPHI=PAR(12)

             DO 6 I2=1,NY
                YY=YYM(I2)
             DO 6 I1=1,NX
                XX=XXM(I1)

                I21=(I2-1)*NX +I1

!                RR2=SQRT(XX**2+YY**2)
!                PHI=ATAN2(YY,XX)
!                CALL POLAR_2D(XX,YY,RR2,PHI)
                RR2=SQRT(XX**2+YY**2)
                CALL POLAR2D(XX,YY,PHI)

!                write(*,*) 'RR2-PHI=',NRR,NPHI
                CALL REGN2(PHIM,RRM,NPHI,NRR,PHI,RR2,EI)

!                write(*,*) 'RR2-PHI=',EI

                IF(EI .GT. 0) THEN
                   DR(I21)=FUNL2(GM,RRM,PHIM,NRR,NPHI,RR2,PHI,IERR)
!                   DR(I21)= Burke_2D(GM,RRM,PHIM,NRR,NPHI,RR2,PHI)
                ELSE
                   DR(I21)=0.0
                ENDIF
6            CONTINUE
          RETURN
          END SUBROUTINE POLAR_DECART2D
!------------------------------------------

      SUBROUTINE POLAR_DECART3D_2(GM,XXM,YYM,ZZM,NX,NY,NZ,RRM,TETAM,PHIM, &
                                  NRR,NTETA,NPHI,DR)

        USE PARAMTR
          IMPLICIT NONE
          INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
          INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
          REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

          REAL(wp) :: GM(NRR*NTETA*NPHI)
          REAL(sp) :: XXM(NX),YYM(NY),ZZM(NZ)
          REAL(sp) :: RRM(NRR),TETAM(NTETA),PHIM(NPHI)

          REAL(wp) :: DR(NX*NY*NZ)
          REAL(wp) :: FUNL3,ALINEAR_3D
          INTEGER :: NX,NY,NZ,NRR,NTETA,NPHI

! local
          REAL(wp) :: GMM(NRR*NTETA*NPHI)
          REAL(sp) :: XX,YY,ZZ,RXY                  !!!!!!!!!
          REAL(sp) :: RRR,TETA,PHI,HX,HY,HZ,HPP,HTET,HPHI     !!!!!!!!!!!!!!!
          INTEGER :: KK,I,I1,I2,I3,I21,I321,IERR,EI
          INTEGER :: J1,J2,J3,J4
!          HPP =PAR(35)
!          HTET=PAR(15)
!          HPHI=PAR(12)

          HX=PAR(18)
          HY=PAR(21)
          HZ=PAR(24)

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

                CALL INDNW3(XXM,YYM,ZZM,NX,NY,NZ,      &
                        HX,HY,HZ,XX,YY,ZZ,&
                        J1,J2,J3,J4)

                CALL REGN(XXM,YYM,ZZM,NX,NY,NZ,XX,YY,ZZ,EI)

                IF(EI .GT. 0) THEN
                   DR(J4)=GM(I321)
                ENDIF

6            CONTINUE
          RETURN
          END SUBROUTINE POLAR_DECART3D_2
!------------------------------------------

      SUBROUTINE From3D_To_2D(PRJ3D,PRJ2D,JJ,NJ,NGAMMA,NPP)
! select 2D projection from 3D massive projections

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 

        INTEGER :: JJ,NJ,NGAMMA,NPP
        REAL(wp) :: PRJ3D(NJ*NGAMMA*NPP)
        REAL(wp) :: PRJ2D(NGAMMA*NPP)
!local
        INTEGER :: I3,J2,J1,J21,I321

         DO 85 I3=JJ,JJ      
            DO 75 J2=1,NGAMMA
            DO 75 J1=1,NPP
               J21=(J2-1)*NPP+J1
               I321=(I3-1)*NPP*NGAMMA + J21
               PRJ2D(J21)=PRJ3D(I321)       
75          CONTINUE
85       CONTINUE

        RETURN
        END SUBROUTINE From3D_To_2D
!---------------------------------------------------------
      SUBROUTINE Model(XXM,YYM,ZZM,NX,NY,NZ,GM)

! simple model for vector field
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp

        INTEGER :: NX,NY,NZ
        REAL(wp) ::  XXM(NX),YYM(NY),ZZM(NZ)
        REAL(wp) ::  GM(NX*NY*NZ,3)
! local
        INTEGER :: I1,I2,I3,I21,I321
        REAL(wp) :: XX,YY,ZZ,RR2

      DO 10 I3=1,NZ
         ZZ=ZZM(I3)

      DO 10 I2=1,NY
         YY=YYM(I2)

      DO 10 I1=1,NX
         XX=XXM(I1)

         RR2=XX**2+YY**2+ZZ**2
         I21=(I2-1)*NX+I1
         I321=(I3-1)*NY*NX+I21
 
         GM(I321,1)=0.0
         GM(I321,2)=1.0   !EXP(-(RR2/0.3))
         GM(I321,3)=1.0   !EXP(-(RR2/0.2))
 10   CONTINUE
      RETURN
      END SUBROUTINE Model
!-------------------------------------------------------------

      SUBROUTINE ProjNum_Vec_Decart_2(DR,PROJ,UM,VM,WM,ALPHA,BETA,GAMMAM, &
                                 XM,YM,ZM,NX,NY,NZ,NU,NV,NW,NGAMMA,IERR)

        USE PARAM_INPUT
!        USE CHANALS

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp

        INTEGER :: NX,NY,NZ,NU,NV,NW,NGAMMA,IERR
        REAL(wp) :: DR(NX*NY*NZ,3)
        REAL(wp) :: XM(NX),YM(NY),ZM(NZ)
        REAL(wp) :: PROJ(NU*NV)
        REAL(wp) :: UM(NU),VM(NV),WM(NW)  
        REAL(wp) :: GAMMAM(NGAMMA)
        REAL(wp) :: ALPHA,BETA
! local
        REAL(wp) :: DDR(NX*NY*NZ)
        REAL(wp) :: WUnit(3)
        REAL(wp) :: DD,RRE,RCC  !,BKAPMAX
        REAL(wp) :: GAMMA
        REAL(wp) :: CALPHA,SALPHA,CBETA,SBETA,CGAMMA,SGAMMA
        REAL(wp) :: UU,VV,WW,WW0,HW,WW1,QQQ
        REAL(wp) :: XX,YY,ZZ
        REAL(wp) :: S1,S2,S3,SS,RR2,RRR
        REAL(wp) :: AA(3,3)
        REAL(wp) :: DR1(NX*NY*NZ),DR2(NX*NY*NZ),DR3(NX*NY*NZ)
        REAL(wp) :: DRR3(NX*NY*NZ),DRR1(NX*NY*NZ),DRR2(NX*NY*NZ)
        INTEGER :: NJ,EI
        INTEGER :: I1,I2,I3,I4,I21,I421 !, I321
        INTEGER :: II1,II2,II3,I321
        INTEGER :: J1,J2,J21,J3,J321
        CHARACTER*8 ITEX

        REAL(wp) :: PROJ2(NU*NV)
!functions
        REAL(wp) :: FUNL3

!        MUU    = PARAM(1)
        RRE    = PARAM(5)
!        NALPHA = PARAM(14)
!        NBETA  = PARAM(15)
!        NGAMMA = PARAM(16)
!        NJ     = NALPHA*NBETA  !*NGAMMA
!        NU = 2**MUU+1
!        NV = NU
!        NW = NU

!        ALLOCATE(ALPHAM(NALPHA), BETAM(NBETA), GAMMAM(NGAMMA))

!        DD  = 1.5
!        RCC = 0.0
!        BKAPMAX=ASIN(RRE/DD)

        IERR=0
        ITEX='Procedure: "ProjNum_Vec_Decart" '
        IF(NW.LT.100) GO TO 20
        IERR=1
        CALL ERRPRN(ITEX,IERR)
        RETURN
!
20      CONTINUE
!        ALPHA=ALPHA+PI/2. 

         CALPHA = COS(ALPHA)
         SALPHA = SIN(ALPHA)
         CBETA  = COS(BETA)
         SBETA  = SIN(BETA)

!--  inverse  matrix AA(AL,BET,GAMM=0)----------------------------
         AA(1,1)= CALPHA*CBETA
         AA(1,2)= SALPHA*CBETA
         AA(1,3)=-SBETA
         AA(2,1)=-SALPHA
         AA(2,2)= CALPHA
         AA(2,3)= 0.0
         AA(3,1)= CALPHA*SBETA
         AA(3,2)= SALPHA*SBETA
         AA(3,3)= CBETA
!-----------------------------------------------------------------


         WUnit(1)=CALPHA*SBETA
         WUnit(2)=SALPHA*SBETA
         WUnit(3)=CBETA
!         WUnit(1)=-SBETA*CALPHA
!         WUnit(2)= SBETA*SALPHA
!         WUnit(3)= CBETA

      DO 6 II3=1,NZ
      DO 6 II2=1,NY
      DO 6 II1=1,NX
         I321=(II3-1)*NX*NY+(II2-1)*NX +II1

         DR1(I321)= DR(I321,1)
         DR2(I321)= DR(I321,2)
         DR3(I321)= DR(I321,3)
 6    CONTINUE   

      DO 65 J3=1,NZ
         ZZ=ZM(J3)

      DO 65 J2=1,NY
         YY=YM(J2)

      DO 65 J1=1,NX
         XX=XM(J1)

         J21=(J2-1)*NX+J1
         J321=(J3-1)*NX*NY+J21


!         RR2=SQRT(VV**2+UU**2)
!         QQQ=RRE**2-RR2**2
!         IF(QQQ .LT. 0.0) QQQ=0.

!         WW0=SQRT(QQQ)
!         HW =2.*WW0/(NW-1)

!         CALL POLAR2D(UU,VV,GAMMA)
!         CGAMMA = COS(GAMMA)
!         SGAMMA = SIN(GAMMA)

!         S=0.
!         SS=0.
!      DO 55 J3=1,NW
!         WW =-WW0+(J3-1)*HW
!         WWDD=1.+WW/DD


         UU= AA(1,1)*XX + AA(1,2)*YY + AA(1,3)*ZZ  
         VV= AA(2,1)*XX + AA(2,2)*YY + AA(2,3)*ZZ
         WW= AA(3,1)*XX + AA(3,2)*YY + AA(3,3)*ZZ

!         XX= AA(1,1)*UU + AA(2,1)*VV + AA(3,1)*WW  
!         YY= AA(1,2)*UU + AA(2,2)*VV + AA(3,2)*WW
!         ZZ= AA(1,3)*UU + AA(2,3)*VV + AA(3,3)*WW

!      XX=  (CALPHA*CBETA*CGAMMA - SALPHA*SGAMMA)*UU +  &
!!    %    (-CALPHA*CBETA*SGAMMA - SALPHA*CGAMMA)*VV+
!           CALPHA*SBETA*WW

!      YY=  (SALPHA*CBETA*CGAMMA + CALPHA*SGAMMA)*UU + &
!!     %    (-SALPHA*CBETA*SGAMMA + CALPHA*CGAMMA)*VV+
!           SALPHA*SBETA*WW

!      ZZ=  -SBETA*CGAMMA*UU +  &
!!     %     SBETA*SGAMMA*VV +
!           CBETA*WW

!         RRR=SQRT(XX**2+YY**2+ZZ**2)

         CALL REGN(UM,VM,WM,NU,NV,NW,UU,VV,WW,EI)

         IF(EI .LT. 0.) GO TO 65

         DRR1(J321)=FUNL3(UM,VM,WM,DR1,NU,NV,NW,UU,VV,WW,IERR)
         DRR2(J321)=FUNL3(UM,VM,WM,DR2,NU,NV,NW,UU,VV,WW,IERR)
         DRR3(J321)=FUNL3(UM,VM,WM,DR3,NU,NV,NW,UU,VV,WW,IERR)

65    CONTINUE

      DO 75 J3=1,NW
      DO 75 J2=1,NV
      DO 75 J1=1,NU
         J21=(J2-1)*NU+J1
         J321=(J3-1)*NU*NV+J21

         DR(J321,1)= AA(1,1)* DRR1(J321) + AA(2,1)* DRR2(J321) + AA(3,1)* DRR3(J321)  
         DR(J321,2)= AA(1,2)* DRR1(J321) + AA(2,2)* DRR2(J321) + AA(3,2)* DRR3(J321)
         DR(J321,3)= AA(1,3)* DRR1(J321) + AA(2,3)* DRR2(J321) + AA(3,3)* DRR3(J321)

75    CONTINUE

      DO 85 J2=1,NV
         VV=VM(J2)
      DO 85 J1=1,NU
         UU=UM(J1)
         J21=(J2-1)*NU+J1

         RR2=SQRT(VV**2+UU**2)
         QQQ=RRE**2-RR2**2
         IF(QQQ .LT. 0.0) QQQ=0.

         WW0=SQRT(QQQ)
         HW =2.*WW0/(NW-1)

         SS=0.
      DO 80 J3=1,NW
         WW1=1.
         IF(J3.EQ.1.OR.J3.EQ.NW) WW1=0.5
         SS=SS+WW1*DR(J321,3)
80    CONTINUE

      IF(RR2 .LT. RRE) THEN
         PROJ(J21)=SS*HW 
      ELSE
         PROJ(J21)=0. 
      ENDIF

85   CONTINUE
      RETURN
      END SUBROUTINE ProjNum_Vec_Decart_2
!--------------------------------------------

      SUBROUTINE ProjNum_Vec_D(DR,PROJ,UM,VM,WM,ALPHAM,BETAM,GAMMAM, &
                                    XM,YM,ZM,NX,NY,NZ,NU,NV,NW,NJ,  &
                                    NALPHA,NBETA,NGAMMA,IERR)

        USE PARAM_INPUT

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp

        INTEGER :: NX,NY,NZ,NU,NV,NW,NJ
        INTEGER :: NALPHA,NBETA,NGAMMA,IERR
        REAL(wp) :: DR(NX*NY*NZ,3)
        REAL(wp) :: XM(NX),YM(NY),ZM(NZ)
        REAL(wp) :: PROJ(NU*NV*NJ)
        REAL(wp) :: UM(NU),VM(NV),WM(NW)  
        REAL(wp) :: ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
! local
        REAL(wp) :: ALPHA,BETA  !,GAMMA
        INTEGER :: I1,I2,I4,I421 !, I321
        INTEGER :: J1,J2,J21
        CHARACTER*8 ITEX
        REAL(wp) :: PROJ2(NU*NV)
!functions
        REAL(wp) :: FUNL3

        IERR=0
        ITEX='Procedure: "ProjNum_Vec" '
        IF(NW.LT.100) GO TO 20
        IERR=1
        CALL ERRPRN(ITEX,IERR)
        RETURN
!
20      CONTINUE
!
      DO 85 I4=1,NJ
         CALL IALBET(I4,NALPHA,I1,I2)

         ALPHA=ALPHAM(I1)
         BETA =BETAM(I2)

      CALL ProjNum_Vec_Decart(DR,PROJ2,UM,VM,WM,ALPHA,BETA,GAMMAM, &
                              XM,YM,ZM,NX,NY,NZ,NU,NV,NW,NGAMMA,IERR)

      DO 75 J2=1,NV
      DO 75 J1=1,NU
         J21=(J2-1)*NU+J1
         I421=(I4-1)*NU*NV + J21
         PROJ(I421)=PROJ2(J21)       
  75  CONTINUE
  85  CONTINUE
      RETURN
      END SUBROUTINE ProjNum_Vec_D
!--------------------------------------------

      SUBROUTINE ProjNum_Vec_Decart(DR,PROJ,UM,VM,WM,ALPHA,BETA,GAMMAM, &
                                 XM,YM,ZM,NX,NY,NZ,NU,NV,NW,NGAMMA,IERR)

        USE PARAM_INPUT
!        USE CHANALS

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp

        INTEGER :: NX,NY,NZ,NU,NV,NW,NGAMMA,IERR
        REAL(wp) :: DR(NX*NY*NZ,3)
        REAL(wp) :: XM(NX),YM(NY),ZM(NZ)
        REAL(wp) :: PROJ(NU*NV)
        REAL(wp) :: UM(NU),VM(NV),WM(NW)  
        REAL(wp) :: GAMMAM(NGAMMA)
        REAL(wp) :: ALPHA,BETA
! local
        REAL(wp) :: DDR(NX*NY*NZ)
        REAL(wp) :: WUnit(3)
        REAL(wp) :: DD,RRE,RCC  !,BKAPMAX
        REAL(wp) :: GAMMA
        REAL(wp) :: CALPHA,SALPHA,CBETA,SBETA,CGAMMA,SGAMMA
        REAL(wp) :: UU,VV,WW,WW0,HW,WW1,QQQ
        REAL(wp) :: XX,YY,ZZ
        REAL(wp) :: S,SS,RR2,RRR
        REAL(wp) :: AA(3,3)
        INTEGER :: NJ,EI
        INTEGER :: I1,I2,I3,I4,I21,I421 !, I321
        INTEGER :: II1,II2,II3
        INTEGER :: J1,J2,J21,J3,J321
        CHARACTER*8 ITEX

        REAL(wp) :: PROJ2(NU*NV)
!functions
        REAL(wp) :: FUNL3,Burke_3D

!        MUU    = PARAM(1)
        RRE    = PARAM(5)
!        NALPHA = PARAM(14)
!        NBETA  = PARAM(15)
!        NGAMMA = PARAM(16)
!        NJ     = NALPHA*NBETA  !*NGAMMA
!        NU = 2**MUU+1
!        NV = NU
!        NW = NU

!        ALLOCATE(ALPHAM(NALPHA), BETAM(NBETA), GAMMAM(NGAMMA))

!        DD  = 1.5
!        RCC = 0.0
!        BKAPMAX=ASIN(RRE/DD)

        IERR=0
        ITEX='Procedure: "ProjNum_Vec_Decart" '
        IF(NW.LT.100) GO TO 20
        IERR=1
        CALL ERRPRN(ITEX,IERR)
        RETURN
!
20      CONTINUE
!        ALPHA=ALPHA+PI/2. 

         CALPHA = COS(ALPHA)
         SALPHA = SIN(ALPHA)
         CBETA  = COS(BETA)
         SBETA  = SIN(BETA)

!--  direct  matrix AA(AL,BET,GAMM=0)----------------------------
         AA(1,1)= CALPHA*CBETA
         AA(1,2)= SALPHA*CBETA
         AA(1,3)=-SBETA
         AA(2,1)=-SALPHA
         AA(2,2)= CALPHA
         AA(2,3)= 0.0
         AA(3,1)= CALPHA*SBETA
         AA(3,2)= SALPHA*SBETA
         AA(3,3)= CBETA
!-----------------------------------------------------------------


         WUnit(1)=CALPHA*SBETA
         WUnit(2)=SALPHA*SBETA
         WUnit(3)=CBETA
!         WUnit(1)=-SBETA*CALPHA
!         WUnit(2)= SBETA*SALPHA
!         WUnit(3)= CBETA

      DO 6 II3=1,NZ
      DO 6 II2=1,NY
      DO 6 II1=1,NX
         J321=(II3-1)*NX*NY+(II2-1)*NX +II1

         DDR(J321)= DR(J321,1)*WUnit(1)   &
                  + DR(J321,2)*WUnit(2)   &
                  + DR(J321,3)*WUnit(3)
 6    CONTINUE   

      DO 65 J2=1,NV
         VV=VM(J2)

      DO 65 J1=1,NU
         UU=UM(J1)

         J21=(J2-1)*NU+J1

         RR2=SQRT(VV**2+UU**2)
         QQQ=RRE**2-RR2**2
         IF(QQQ .LT. 0.0) QQQ=0.

         WW0=SQRT(QQQ)
         HW =2.*WW0/(NW-1)

!         CALL POLAR2D(UU,VV,GAMMA)
!         CGAMMA = COS(GAMMA)
!         SGAMMA = SIN(GAMMA)

         S=0.
         SS=0.
      DO 55 J3=1,NW
         WW =-WW0+(J3-1)*HW

!         XX= AA(1,1)*UU + AA(1,2)*VV + AA(1,3)*WW  
!         YY= AA(2,1)*UU + AA(2,2)*VV + AA(2,3)*WW
!         ZZ= AA(3,1)*UU + AA(3,2)*VV + AA(3,3)*WW

         XX= AA(1,1)*UU + AA(2,1)*VV + AA(3,1)*WW  
         YY= AA(1,2)*UU + AA(2,2)*VV + AA(3,2)*WW
         ZZ= AA(1,3)*UU + AA(2,3)*VV + AA(3,3)*WW

!      XX=  (CALPHA*CBETA*CGAMMA - SALPHA*SGAMMA)*UU +  &
!!    %    (-CALPHA*CBETA*SGAMMA - SALPHA*CGAMMA)*VV+
!           CALPHA*SBETA*WW

!      YY=  (SALPHA*CBETA*CGAMMA + CALPHA*SGAMMA)*UU + &
!!     %    (-SALPHA*CBETA*SGAMMA + CALPHA*CGAMMA)*VV+
!           SALPHA*SBETA*WW

!      ZZ=  -SBETA*CGAMMA*UU +  &
!!     %     SBETA*SGAMMA*VV +
!           CBETA*WW

         RRR=SQRT(XX**2+YY**2+ZZ**2)

         CALL REGN(XM,YM,ZM,NX,NY,NZ,XX,YY,ZZ,EI)

         IF(EI .LT. 0.) GO TO 55

!         S=FUNL3(XM,YM,ZM,DDR,NX,NY,NZ,XX,YY,ZZ,IERR)
         S= Burke_3D(DDR,XM,YM,ZM,NX,NY,NZ,XX,YY,ZZ)

         WW1=1.
         IF(J3.EQ.1.OR.J3.EQ.NW) WW1=0.5
         SS=SS+S*WW1
55    CONTINUE

      IF(RR2 .LT. RRE) THEN
         PROJ(J21)=SS*HW 
      ELSE
         PROJ(J21)=0. 
      ENDIF

 65   CONTINUE
      RETURN
      END SUBROUTINE ProjNum_Vec_Decart
!--------------------------------------------

      SUBROUTINE ProjNum_Vec_P(DR,PROJ,PPM,NPP,WM,NW,ALPHAM,BETAM,GAMMAM,    &
                             NALPHA,NBETA,NGAMMA,XM,YM,ZM,NX,NY,NZ,IERR)

        USE PARAM_INPUT

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp

!
        INTEGER :: NPP,NALPHA,NBETA,NGAMMA
        INTEGER :: NX,NY,NZ,NW,IERR
        REAL(wp) :: DR(NX*NY*NZ,3)
        REAL(wp) :: PPM(NPP)
        REAL(wp) :: XM(NX),YM(NY),ZM(NZ),WM(NW)
        REAL(wp) :: ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)
        REAL(wp) :: PROJ(NPP*NGAMMA*NALPHA*NBETA)
! local
        REAL(wp), ALLOCATABLE :: PROJ2(:)
        REAL(wp) :: ALPHA,BETA,GAMMA
        INTEGER :: NJ,EI
        INTEGER :: I1,I2,I4,I421 !, I3,I21,I321
        INTEGER :: J1,J2,J21  !,J3,J321
        CHARACTER*8 ITEX
!functions

        NJ = NALPHA*NBETA  !*NGAMMA

        ALLOCATE(PROJ2(NPP*NGAMMA))

        IERR=0
        ITEX='Procedure: "ProjNum_Vec_P" '
        IF(NW.LT.100) GO TO 20
        IERR=1
        CALL ERRPRN(ITEX,IERR)
        RETURN
!
20      CONTINUE
!
      DO 85 I4=1,NJ
         CALL IALBET(I4,NALPHA,I1,I2)
!         CALL IALBETGAM(I4,NALPHA,NBETA,I1,I2,I3)
         ALPHA=ALPHAM(I1)
         BETA =BETAM(I2)

!         write(*,*) 'ProjNum_Vec_P', ALPHA, BETA

      CALL ProjNum_Vec_Polar(DR,PROJ2,PPM,NPP,WM,NW,ALPHA,BETA,    &
                             GAMMAM,NGAMMA,XM,YM,ZM,NX,NY,NZ,IERR)

      DO 75 J2=1,NGAMMA
      DO 75 J1=1,NPP
         J21=(J2-1)*NPP+J1
         I421=(I4-1)*NPP*NGAMMA + J21
         PROJ(I421)=PROJ2(J21)       
  75  CONTINUE
  85  CONTINUE
      RETURN
      END SUBROUTINE ProjNum_Vec_P
!--------------------------------------------

      SUBROUTINE ProjNum_Vec_Polar(DR,PROJ2,PPM,NPP,WM,NW,ALPHA,BETA,      &
                                   GAMMAM,NGAMMA,XM,YM,ZM,NX,NY,NZ,IERR)

        USE PARAM_INPUT
        USE CHANALS

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp
!
        INTEGER :: NX,NY,NZ,IERR
        INTEGER :: NW,NPP,NGAMMA
        REAL(wp) :: DR(NX*NY*NZ,3)
        REAL(wp) :: PROJ2(NPP*NGAMMA)
        REAL(wp) :: PPM(NPP)
        REAL(wp) :: WM(NW)
        REAL(wp) :: GAMMAM(NGAMMA)
        REAL(wp) :: ALPHA,BETA
        REAL(wp) :: XM(NX),YM(NY),ZM(NZ)

! local
        REAL(wp) :: DDR1(NX*NY*NZ),DDR2(NX*NY*NZ),DDR3(NX*NY*NZ)
        REAL(wp) :: WUnit(3)
        REAL(wp) :: RRE,RCC
        REAL(wp) :: GAMMA
        REAL(wp) :: CALPHA,SALPHA,CBETA,SBETA,CGAMMA,SGAMMA
        REAL(wp) :: PP,UU,VV,WW,WW0,HW,WW1,QQQ
        REAL(wp) :: XX,YY,ZZ
        REAL(wp) :: S,S1,S2,S3,SS,RRR
        REAL(wp) :: AA(3,3)
!        INTEGER :: MUU,NU,NV
        INTEGER :: NALPHA,NBETA
        INTEGER :: EI   !,NJ
        INTEGER :: I1,I2,I3,I4,I21,I421 !, I321
        INTEGER :: II1,II2,II3
        INTEGER :: J1,J2,J21,J3,J321,IND,KSEC
        CHARACTER*8 ITEX
!functions
        REAL(wp) :: FUNL3

        RRE    = PARAM(5)
!        NALPHA = PARAM(14)
!        NBETA  = PARAM(15)
!        NJ     = NALPHA*NBETA  !*NGAMMA

        RCC = 0.0

        IERR=0
        ITEX='Procedure: "ProjNum_Vec_Polar" '
        IF(NW.LT.100) GO TO 20
        IERR=1
        CALL ERRPRN(ITEX,IERR)
        RETURN
!
20      CONTINUE
!
!      DO 85 I4=1,NJ
!         CALL IALBET(I4,NALPHA,I1,I2)
!!         CALL IALBETGAM(I4,NALPHA,NBETA,I1,I2,I3)
!         ALPHA=ALPHAM(I1)
!         BETA =BETAM(I2)

         CALPHA = COS(ALPHA)
         SALPHA = SIN(ALPHA)
         CBETA  = COS(BETA)
         SBETA  = SIN(BETA)
!--  inverse  matrix AA(AL,BET,GAMM=0)----------------------------
         AA(1,1)= CALPHA*CBETA
         AA(1,2)= SALPHA*CBETA
         AA(1,3)=-SBETA
         AA(2,1)=-SALPHA
         AA(2,2)= CALPHA
         AA(2,3)= 0.0
         AA(3,1)= CALPHA*SBETA
         AA(3,2)= SALPHA*SBETA
         AA(3,3)= CBETA
!-----------------------------------------------------------------

         WUnit(1)=CALPHA*SBETA
         WUnit(2)=SALPHA*SBETA
         WUnit(3)=CBETA
!         WUnit(1)=-SBETA*CGAMMA
!         WUnit(2)= SBETA*SGAMMA
!         WUnit(3)= CBETA

      DO 6 II3=1,NZ
      DO 6 II2=1,NY
      DO 6 II1=1,NX
         J321=(II3-1)*NX*NY+(II2-1)*NX +II1

         DDR1(J321)= DR(J321,1)  !*WUnit(1)   &
         DDR2(J321)= DR(J321,2)  !*WUnit(2)   &
         DDR3(J321)= DR(J321,3) !         + DR(J321,3)*WUnit(3)
 6    CONTINUE   

!      IND=0
!      KSEC=NZ/2
!      CALL SPRINT_3D(IND,KSEC,DDR,NX,NY,NZ,KAN(80))
!  direct matrix
!         AA(1,1)= CALPHA*CBETA*CGAMMA - SALPHA*SGAMMA
!         AA(1,2)= SALPHA*CBETA*CGAMMA + CALPHA*SGAMMA
!         AA(1,3)=-SBETA*CGAMMA
!         AA(2,1)=-CALPHA*CBETA*SGAMMA - SALPHA*CGAMMA
!         AA(2,2)=-SALPHA*CBETA*SGAMMA + CALPHA*CGAMMA
!         AA(2,3)= SBETA*SGAMMA
!         AA(3,1)= CALPHA*SBETA
!         AA(3,2)= SALPHA*SBETA
!         AA(3,3)= CBETA

      DO 65 J2=1,NGAMMA
         GAMMA=GAMMAM(J2)

         CGAMMA = COS(GAMMA)
         SGAMMA = SIN(GAMMA)

      DO 65 J1=1,NPP
         J21=(J2-1)*NPP+J1
         PP=PPM(J1)

         UU=PP*CGAMMA
         VV=PP*SGAMMA

         QQQ=RRE**2-PP**2
         IF(QQQ .LT. 0.0) QQQ=0.

         WW0=SQRT(QQQ)
         HW =2.*WW0/(NW-1)

         S1=0.0; S2=0.0; S3=0.0
         SS=0.
      DO 55 J3=1,NW
         WW =-WW0+(J3-1)*HW
!         WWDD=1.+WW/DD

!--   AA (AL,BET,GAM)-------------------------------------------
         XX= AA(1,1)*UU + AA(1,2)*VV + AA(1,3)*WW  
         YY= AA(2,1)*UU + AA(2,2)*VV + AA(2,3)*WW
         ZZ= AA(3,1)*UU + AA(3,2)*VV + AA(3,3)*WW

!         XX= AA(1,1)*UU + AA(3,1)*WW  
!         YY= AA(1,2)*UU + AA(3,2)*WW
!         ZZ= AA(1,3)*UU + AA(3,3)*WW

!      XX=  (CALPHA*CBETA*CGAMMA - SALPHA*SGAMMA)*UU +  &
!!    %    (-CALPHA*CBETA*SGAMMA - SALPHA*CGAMMA)*VV+
!           CALPHA*SBETA*WW

!      YY=  (SALPHA*CBETA*CGAMMA + CALPHA*SGAMMA)*UU + &
!!     %    (-SALPHA*CBETA*SGAMMA + CALPHA*CGAMMA)*VV+
!           SALPHA*SBETA*WW

!      ZZ=  -SBETA*CGAMMA*UU +  &
!!     %     SBETA*SGAMMA*VV +
!           CBETA*WW

         RRR=SQRT(UU**2+VV**2)

         CALL REGN(XM,YM,ZM,NX,NY,NZ,XX,YY,ZZ,EI)

         IF(EI .LT. 0) GO TO 55

         S1=FUNL3(XM,YM,ZM,DDR1,NX,NY,NZ,XX,YY,ZZ,IERR)
         S2=FUNL3(XM,YM,ZM,DDR2,NX,NY,NZ,XX,YY,ZZ,IERR)
         S3=FUNL3(XM,YM,ZM,DDR3,NX,NY,NZ,XX,YY,ZZ,IERR)

         S=S1*WUnit(1)+S2*WUnit(2)+S3*WUnit(3)

         WW1=1.
         IF(J3 .EQ. 1 .OR. J3 .EQ. NW) WW1=0.5
         SS=SS+S*WW1
55    CONTINUE

      IF(RRR .GE. RCC) THEN
         PROJ2(J21)=SS*HW 
      ELSE
         PROJ2(J21)=0. 
      ENDIF
   
 65   CONTINUE

!      DO 75 J2=1,NGAMMA
!      DO 75 J1=1,NPP
!         J21=(J2-1)*NPP+J1
!         I421=(I4-1)*NPP*NGAMMA + J21
!         PROJ(I421)=PROJ2(J21)       
!  75  CONTINUE
!  85  CONTINUE
      RETURN
      END SUBROUTINE ProjNum_Vec_Polar
!--------------------------------------------

      SUBROUTINE RayTransform_M_LM(L,M,KK,RR,UM,VM,WM,NU,NV,NW,ALPHA,BETA,CLMM,PROJRE2D)

!  Ray transform of vectors M_LM

        USE PARAM_INPUT
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp


        INTEGER :: L,M,NU,NV,NW
        REAL(wp) :: UM(NU),VM(NV),WM(NW) 
        REAL(wp) :: ALPHA,BETA 
        REAL(wp) :: KK,RR
        REAL(wp) :: PROJRE2D(NU*NV)
        COMPLEX(wp) :: CLMM(NU*NV)
! local
        INTEGER :: I1,I2,I3,I21,I321,J
        COMPLEX(wp) :: MLMM1(3),MLMM2(3),MLMM(3),SS1,SS2
        REAL(wp) :: UNIT(3)
        REAL(wp) :: UU,VV,WW,WW0,HW,WW1,QQQ
        REAL(wp) :: XX,YY,ZZ,RRR,RR2,RRE 
        REAL(wp) :: AA(3,3),BB(3,3)
        REAL(wp) :: CALPHA,CBETA,SALPHA,SBETA
        REAL(wp) :: CTET,STET,CPHI,SPHI 
        REAL(wp) :: TETA,PHI,SIGNZ 

        RRE    = PARAM(5)

         CALPHA = COS(ALPHA)
         SALPHA = SIN(ALPHA)
         CBETA  = COS(BETA)
         SBETA  = SIN(BETA)

!  direct matrix
         AA(1,1)= CALPHA*CBETA
         AA(1,2)= SALPHA*CBETA
         AA(1,3)=-SBETA
         AA(2,1)=-SALPHA
         AA(2,2)= CALPHA
         AA(2,3)= 0.
         AA(3,1)= CALPHA*SBETA
         AA(3,2)= SALPHA*SBETA
         AA(3,3)= CBETA

         Unit(1)  = AA(3,1)  !SIN(BETA)*COS(ALPHA)
         Unit(2)  = AA(3,2)  !SIN(BETA)*SIN(ALPHA)
         Unit(3)  = AA(3,3)  !COS(BETA) 

!         Unit(1)  = AA(1,3)  
!         Unit(2)  = AA(2,3)  
!         Unit(3)  = AA(3,3)  


         DO 16 I2=1,NV
            VV=VM(I2)
         DO 16 I1=1,NU
            UU=UM(I1)
            I21=(I2-1)*NU+I1

            RR2=SQRT(UU**2+VV**2)

            QQQ=RRE**2-RR2**2
            IF(QQQ .LT. 0.0) QQQ=0.

            WW0=SQRT(QQQ)
            HW =2.*WW0/(NW-1)

            SS1=(0.0,0.0)
            DO 14 I3=1,NW
               WW =-WW0+(I3-1)*HW
               WW1=1.
               IF(I3.EQ.1.OR.I3.EQ.NW) WW1=0.5
   
               XX=AA(1,1)*UU + AA(2,1)*VV + AA(3,1)*WW
               YY=AA(1,2)*UU + AA(2,2)*VV + AA(3,2)*WW
               ZZ=AA(1,3)*UU + AA(2,3)*VV + AA(3,3)*WW

!               XX=AA(1,1)*UU + AA(1,2)*VV + AA(1,3)*WW
!               YY=AA(2,1)*UU + AA(2,2)*VV + AA(2,3)*WW
!               ZZ=AA(3,1)*UU + AA(3,2)*VV + AA(3,3)*WW

               RRR=SQRT(XX**2+YY**2+ZZ**2)
               SIGNZ=SIGN(1.0,ZZ)
               CALL POLAR(XX,YY,ZZ,TETA,PHI)
! in order to avoide TETA=0 in DYLM_TETA procedure
               IF(TETA .EQ. PI .OR. TETA .EQ. 0.) THEN 
                  TETA=PI*(1.+SIGNZ/180.)
               ENDIF   

            CALL M_LM (L,M,KK,RRR,TETA,PHI,MLMM1)

!            write(*,*) 'RayTransform:M_LMM1=', MLMM1

            CTET = COS(TETA)
            STET = SIN(TETA)
            CPHI = COS(PHI)
            SPHI = SIN(PHI)

            BB(1,1)= STET*CPHI
            BB(1,2)= CTET*CPHI
            BB(1,3)=-SPHI
            BB(2,1)= STET*SPHI
            BB(2,2)= CTET*SPHI
            BB(2,3)= CPHI
            BB(3,1)= CTET
            BB(3,2)=-STET
            BB(3,3)= 0.

! decart components
            MLMM(1)= MLMM1(1)*BB(1,1)+MLMM1(2)*BB(1,2)+MLMM1(3)*BB(1,3)
            MLMM(2)= MLMM1(1)*BB(2,1)+MLMM1(2)*BB(2,2)+MLMM1(3)*BB(2,3)
            MLMM(3)= MLMM1(1)*BB(3,1)+MLMM1(2)*BB(3,2)+MLMM1(3)*BB(3,3)

            SS2=(0.0,0.0)
            DO 12 J=1,3
               SS2=SS2+MLMM(J)*Unit(J)
12          CONTINUE
               SS1=SS1+SS2
14          CONTINUE

               CLMM(I21)=SS1*HW*WW1    
               PROJRE2D(I21)=REAL(CLMM(I21))  
!                  write(*,*) 'SS1=',CLMM
16       CONTINUE  
      RETURN
      END SUBROUTINE RayTransform_M_LM
!-------------------------------------------------

      SUBROUTINE RayTransform_M_LM_REAL(L,M,KK,RR,UM,VM,WM,NU,NV,NW,ALPHA,BETA,CLMM,PROJRE2D)

!  Ray transform of vectors M_LM

        USE PARAM_INPUT
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp


        INTEGER :: L,M,NU,NV,NW
        REAL(wp) :: UM(NU),VM(NV),WM(NW) 
        REAL(wp) :: ALPHA,BETA 
        REAL(wp) :: KK,RR
        REAL(wp) :: PROJRE2D(NU*NV)
        COMPLEX(wp) :: CLMM(NU*NV)
! local
        INTEGER :: I1,I2,I3,I21,I321,J
        COMPLEX(wp) :: MLMM1(3)
        REAL(wp) :: MLMM(3),MLMMRE(3),SS1,SS2
        REAL(wp) :: UNIT(3)
        REAL(wp) :: UU,VV,WW,WW0,HW,WW1,QQQ
        REAL(wp) :: XX,YY,ZZ,RRR,RR2,RRE 
        REAL(wp) :: AA(3,3),BB(3,3)
        REAL(wp) :: CALPHA,CBETA,SALPHA,SBETA
        REAL(wp) :: CTET,STET,CPHI,SPHI 
        REAL(wp) :: TETA,PHI,SIGNZ 

        RRE    = PARAM(5)

         CALPHA = COS(ALPHA)
         SALPHA = SIN(ALPHA)
         CBETA  = COS(BETA)
         SBETA  = SIN(BETA)

!  direct matrix
         AA(1,1)= CALPHA*CBETA
         AA(1,2)= SALPHA*CBETA
         AA(1,3)=-SBETA
         AA(2,1)=-SALPHA
         AA(2,2)= CALPHA
         AA(2,3)= 0.
         AA(3,1)= CALPHA*SBETA
         AA(3,2)= SALPHA*SBETA
         AA(3,3)= CBETA

         Unit(1)  = AA(3,1)  !SIN(BETA)*COS(ALPHA)
         Unit(2)  = AA(3,2)  !SIN(BETA)*SIN(ALPHA)
         Unit(3)  = AA(3,3)  !COS(BETA) 

         DO 16 I2=1,NV
            VV=VM(I2)
         DO 16 I1=1,NU
            UU=UM(I1)
            I21=(I2-1)*NU+I1

            RR2=SQRT(UU**2+VV**2)

            QQQ=RRE**2-RR2**2
            IF(QQQ .LT. 0.0) QQQ=0.

            WW0=SQRT(QQQ)
            HW =2.*WW0/(NW-1)

            SS1=0.0
            DO 14 I3=1,NW
               WW =-WW0+(I3-1)*HW
               WW1=1.
               IF(I3.EQ.1.OR.I3.EQ.NW) WW1=0.5
   
               XX=AA(1,1)*UU + AA(2,1)*VV + AA(3,1)*WW
               YY=AA(1,2)*UU + AA(2,2)*VV + AA(3,2)*WW
               ZZ=AA(1,3)*UU + AA(2,3)*VV + AA(3,3)*WW

!               XX=AA(1,1)*UU + AA(1,2)*VV + AA(1,3)*WW
!               YY=AA(2,1)*UU + AA(2,2)*VV + AA(2,3)*WW
!               ZZ=AA(3,1)*UU + AA(3,2)*VV + AA(3,3)*WW

               RRR=SQRT(XX**2+YY**2+ZZ**2)
               SIGNZ=SIGN(1.0,ZZ)
               CALL POLAR(XX,YY,ZZ,TETA,PHI)
! in order to avoide TETA=0 in DYLM_TETA procedure
               IF(TETA .EQ. PI .OR. TETA .EQ. 0.) THEN 
                  TETA=PI*(1.+SIGNZ/180.)
               ENDIF   

            CALL M_LM (L,M,KK,RRR,TETA,PHI,MLMM1)

            MLMMRE(1)=REAL(MLMM1(1))
            MLMMRE(2)=REAL(MLMM1(2))
            MLMMRE(3)=REAL(MLMM1(3))

            CTET = COS(TETA)
            STET = SIN(TETA)
            CPHI = COS(PHI)
            SPHI = SIN(PHI)

            BB(1,1)= STET*CPHI
            BB(1,2)= CTET*CPHI
            BB(1,3)=-SPHI
            BB(2,1)= STET*SPHI
            BB(2,2)= CTET*SPHI
            BB(2,3)= CPHI
            BB(3,1)= CTET
            BB(3,2)=-STET
            BB(3,3)= 0.

! decart components 
            MLMM(1)= MLMMRE(1)*BB(1,1)+MLMMRE(2)*BB(1,2)+MLMMRE(3)*BB(1,3)
            MLMM(2)= MLMMRE(1)*BB(2,1)+MLMMRE(2)*BB(2,2)+MLMMRE(3)*BB(2,3)
            MLMM(3)= MLMMRE(1)*BB(3,1)+MLMMRE(2)*BB(3,2)+MLMMRE(3)*BB(3,3)


            SS2=0.0
            DO 12 J=1,3
               SS2=SS2+MLMM(J)*Unit(J)
12          CONTINUE
               SS1=SS1+SS2
14          CONTINUE

!               CLMM(I21)=SS1*HW*WW1    
               PROJRE2D(I21)=SS1*HW*WW1    
!                  write(*,*) 'SS1=',CLMM
16       CONTINUE  
      RETURN
      END SUBROUTINE RayTransform_M_LM_REAL
!-------------------------------------------------



      SUBROUTINE RayTransform_M_LM_Polar(L,M,KK,PPM,GAMMAM,NPP,NGAMMA,WM,NW,   &
                                   ALPHA,BETA,CLMM,BMODRE)

!  Ray transform of vectors M_LM

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp


        INTEGER :: L,M,NPP,NGAMMA,NW
        REAL(wp) :: PPM(NPP),GAMMAM(NGAMMA),WM(NW) 
        REAL(wp) :: ALPHA,BETA 
        REAL(wp) :: KK
        REAL(wp) :: BMODRE(NPP*NGAMMA)
        COMPLEX(wp) :: CLMM(NPP*NGAMMA)
! local
        INTEGER :: I1,I2,I3,I21,I321,J
        COMPLEX(wp) :: MLMM1(3),MLMM2(3),MLMM(3),SS1,SS2
        REAL(wp) :: UNIT(3)
        REAL(wp) :: UU,VV,WW,RP,GAMMA,TETA,PHI
        REAL(wp) :: XX,YY,ZZ,RRR
        REAL(wp) :: AA(3,3)
        REAL(wp) :: CALPHA,CBETA,SALPHA,SBETA

         CALPHA = COS(ALPHA)
         SALPHA = SIN(ALPHA)
         CBETA  = COS(BETA)
         SBETA  = SIN(BETA)

!  direct matrix
         AA(1,1)= CALPHA*CBETA
         AA(1,2)= SALPHA*CBETA
         AA(1,3)=-SBETA
         AA(2,1)=-SALPHA
         AA(2,2)= CALPHA
         AA(2,3)= 0.
         AA(3,1)= CALPHA*SBETA
         AA(3,2)= SALPHA*SBETA
         AA(3,3)= CBETA

         Unit(1)  = SIN(BETA)*COS(ALPHA)
         Unit(2)  = SIN(BETA)*SIN(ALPHA)
         Unit(3)  = COS(BETA) 

         DO 16 I2=1,NGAMMA
            GAMMA=GAMMAM(I2)
         DO 16 I1=1,NPP
            RP=PPM(I1)

            I21=(I2-1)*NPP+I1

            UU=RP*COS(GAMMA)
            VV=RP*SIN(GAMMA)

            SS1=(0.0,0.0)
            DO 14 I3=1,NW
               WW=WM(I3)
            
               I321=(I3-1)*NPP*NGAMMA+I21

               XX=AA(1,1)*UU + AA(2,1)*VV + AA(3,1)*WW
               YY=AA(1,2)*UU + AA(2,2)*VV + AA(3,2)*WW
               ZZ=AA(1,3)*UU + AA(2,3)*VV + AA(3,3)*WW
            
               RRR=SQRT(XX**2+YY**2+ZZ**2)
            
               CALL POLAR(XX,YY,ZZ,TETA,PHI)

            SS2=(0.0,0.0)
            DO 12 J=1,3
               CALL M_LM (L,M,KK,RRR,TETA,PHI,MLMM1)
               CALL M_LM (L,M,KK,RRR,PI-TETA,PI+PHI,MLMM2)
               SS2=SS2+(MLMM1(J)+MLMM2(J))*Unit(J)
!               SS2=SS2+MLMM1(J)*Unit(J)
12          CONTINUE
               SS1=SS1+SS2
14          CONTINUE
               CLMM(I21)=SS1  
               BMODRE(I21)=REAL(SS1)
16       CONTINUE  
      RETURN
      END SUBROUTINE RayTransform_M_LM_Polar
!-------------------------------------------------

      SUBROUTINE Y_1LM (L,M,TETA,PHI,CLMM )

!  Vector spherical harmonics  P_LM  = Y^(-1)(LM)

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp

        INTEGER :: L,M
        REAL(wp) :: TETA,PHI
! vector spherical harmonics
        COMPLEX(wp) :: CLMM(3)

!local 
        COMPLEX(wp),PARAMETER :: ci = (0.,1.)
        COMPLEX(wp) :: CKLM,YLM
        REAL(wp) :: COST, SINT, COSF, SINF
        REAL(wp) :: R3,SIGNUM
        INTEGER :: I
! spherical unit vectors
        REAL(wp) :: AM(3,3)
! functions
      COMPLEX(wp) :: CYLM

!=======================
! initialization 
!=======================
!      MM=ABS(M)

      COST=COS(TETA)
      SINT=SIN(TETA)
      COSF=COS(PHI)
      SINF=SIN(PHI)

!=======================
! matrix  M(x,y,z   <-- r,\theta,\phi)
!=======================

      AM(1,1) = COSF *SINT
      AM(2,1) = SINF *SINT
      AM(3,1) = COST
      
      AM(1,2) = COSF * COST
      AM(2,2) = SINF * COST
      AM(3,2) =-SINT

      AM(1,3) =  -SINF
      AM(2,3) =   COSF
      AM(3,3) =   0. 
 
!------------------------------------
!  spherical harmonics Y(LM)

      YLM= CYLM (L,M,TETA,PHI)
!-------------------------------------
!         DO 6 I =1,3
!            CLMM(I) =rUnit(I)*YLM 
! 6          CONTINUE

      CLMM(1) =YLM 
      CLMM(2) =(0.0,0.0) 
      CLMM(3) =(0.0,0.0) 

      RETURN
      END SUBROUTINE Y_1LM
!--------------------------------------------------
      FUNCTION DerSpherBessel(L,X)
! computer the drivative of function X*Z_L(X),  
! Z_L - any spher. Bessell function.
        USE PARAMTR
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

        INTEGER :: L
        REAL(wp) :: X
        REAL(wp) :: DerSpherBessel
! local
        REAL(wp) :: BES(QMAX)
        REAL(wp) :: SHR1,SHR2,DER
        INTEGER :: ERR

        CALL SpherBessJ(X,QMAX,BES,ERR) 
        SHR1= BES(L+1) 
        SHR2= BES(L+2) 

        DER=SHR1*(L+1)-X*SHR2

        DerSpherBessel=DER

        RETURN
      END FUNCTION DerSpherBessel
!-------------------------------------------------

      SUBROUTINE N_LM (L,M,KK,RR,TETA,PHI,NLMM )
!  Vectors N_LM

        USE PARAMTR
!        USE PARAM_INPUT

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

        INTEGER :: L,M
        REAL(wp) :: KK,RR,TETA,PHI
! vector spherical harmonics
        COMPLEX(wp) :: NLMM(3)

!local 
        COMPLEX(wp),PARAMETER :: ci = (0.,1.)
        COMPLEX(wp) :: CKLM,YLM    !  ,CYY,SS
        REAL(wp) :: COST, SINT, COSF, SINF
        REAL(wp) :: R3,SIGNUM
        REAL(wp) :: DerSpherBessel
        INTEGER :: MM,JJ,I
! spherical unit vectors
!        REAL(wp) :: rUnit(3), tetaUnit(3), fiUnit(3)
! functions
        COMPLEX(wp) :: CLMM1(3),CLMM2(3)

        REAL(wp) :: BES(QMAX)
        REAL(wp) :: SHR1,DER
        INTEGER :: ERR

         R3=SQRT(FLOAT(L*(L+1)))

!        R0=PARAM(5)
!---------------------------
!  N_LM  .NE. 0 if  L=1,2....
!----------------------------
        CALL SpherBessJ(KK*RR,QMAX,BES,ERR) 
        SHR1= BES(L+1)/(KK*RR) 

        CALL Y_1LM (L,M,TETA,PHI,CLMM1 )

        DER = DerSpherBessel(L,KK*RR)
        CALL Y1LM  (L,M,TETA,PHI,CLMM2 )

!        DO 6 I =1,3
!           NLMM(I)=ci*SHR1*CLMM1(I)*SQRT(FLOAT(L*(L+1))) + &
!                   ci*DER/(KK*RR)*CLMM2(I)
!6       CONTINUE   

        DO 6 I =1,3
           NLMM(I)=R3**2*SHR1*CLMM1(I) + R3*DER/(KK*RR)*CLMM2(I)
6       CONTINUE  

      RETURN
      END SUBROUTINE N_LM
!-------------------------------------------------

      SUBROUTINE M_LM (L,M,KK,RR,TETA,PHI,MLMM )
!  Vectors M_LM

        USE PARAMTR
!        USE PARAM_INPUT

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
!        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

        INTEGER :: L,M
        REAL(wp) :: KK,RR,TETA,PHI
! vector spherical harmonics
        COMPLEX(wp) :: MLMM(3)

!local 
        COMPLEX(wp),PARAMETER :: ci = (0.,1.)
        COMPLEX(wp) :: CLMM(3)   
        REAL(wp) :: BES(QMAX)
        REAL(wp) :: SHR1
        INTEGER :: MM,JJ,I,ERR

!---------------------------
!  M_LM  .NE. 0 if  L=1,2....
!----------------------------
        CALL SpherBessJ(KK*RR,QMAX,BES,ERR) 
        SHR1= BES(L+1) 

        CALL  Y0LM (L,M,TETA,PHI,CLMM )

        DO 6 I =1,3
           MLMM(I)=SHR1*CLMM(I)*(-ci)   ! for C_{lm}=-ci*Y_{lm}^{0}
6       CONTINUE   

      RETURN
      END SUBROUTINE M_LM
!-------------------------------------------------


      SUBROUTINE Y1LM (L,M,TETA,PHI,CLMM )
!  Vector spherical harmonics  B_LM  = Y^(1)(LM)

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

        INTEGER :: L,M
        REAL(wp) :: TETA,PHI
! vector spherical harmonics
        COMPLEX(wp) :: CLMM(3)

!local 
        COMPLEX(wp),PARAMETER :: ci = (0.,1.)
        COMPLEX(wp) :: CKLM,YLM   ! ,CYY,SS
        REAL(wp) :: COST, SINT, COSF, SINF
        REAL(wp) :: R3,SIGNUM
       INTEGER :: MM,JJ,I
! spherical unit vectors
        REAL(wp) :: AM(3,3)
! functions
      COMPLEX(wp) :: CYLM
      
!=======================
! initialization 
!=======================
      MM=ABS(M)

      COST=COS(TETA)
      SINT=SIN(TETA)
      COSF=COS(PHI)
      SINF=SIN(PHI)

!=======================
! matrix  M(x,y,z   <-- r,\theta,\phi)
!=======================

      AM(1,1) = COSF *SINT
      AM(2,1) = SINF *SINT
      AM(3,1) = COST
      
      AM(1,2) = COSF * COST
      AM(2,2) = SINF * COST
      AM(3,2) =-SINT

      AM(1,3) =-SINF
      AM(2,3) = COSF
      AM(3,3) = 0. 

!=======================
! unit vectors
!=======================
!      rUnit(1) =     COSF *SINT
!      rUnit(2) =     SINF *SINT
!      rUnit(3) =           COST
!      tetaUnit(1) = COSF * COST
!      tetaUnit(2) = SINF * COST
!      tetaUnit(3) =       -SINT
!      fiUnit(1)  =  -SINF
!      fiUnit(2)  =   COSF
!      fiUnit(3)  =   0. 
!-----------------------------------
 
!  derivatives over teta
      CALL DYLM_TETA (L,MM,TETA,PHI,CKLM)
!------------------------------------
!  spherical harmonics Y(LM)

      YLM= CYLM (L,MM,TETA,PHI)
!-------------------------------------
         IF(L .EQ. 0) THEN
            PRINT*, 'Error,Y1LM: L must not be = 0! STOP'
            STOP
         ENDIF   

         R3=1./SQRT(FLOAT(L*(L+1)))

!         IF(SINT .EQ. 0.) THEN
!            PRINT*, 'Error,Y1LM:TETA must not be = PI or 0! STOP'
!            STOP
!         ENDIF   

         IF(TETA .EQ. PI .OR. TETA .EQ. 0.) THEN
            TETA= PI/360.   !  CLMM(2)=(0.0,0.0) 
         ENDIF

         CLMM(1)=(0.0,0.0)                      ! radial component
         CLMM(2)=R3*CKLM                        !teta Component
         CLMM(3)=ci*R3*FLOAT(MM)/SIN(TETA)*YLM  !fi Componnent

!         DO 6 I =1,3
!            CLMM(I) =R3*(tetaUnit(I)*CKLM +  &
!                     fiUnit(I)*YLM *FLOAT(MM)*ci/SINT)
!!            CLMM(I) = R3*(-tetaUnit(I)*FLOAT(MM)*YLM/SINT  &
!!                     -fiUnit(I)*CKLM *FLOAT(MM)*ci)
!6          CONTINUE

         IF( M .LT. 0) THEN 
               JJ=MOD(MM,2)
               IF(JJ .EQ. 0) THEN 
                  SIGNUM=1.
               ELSE
                  SIGNUM=-1.
               ENDIF
            DO 10 I=1,3
! this is harmonics Y^1(LM)
               CLMM(I)= SIGNUM*CONJG(CLMM(I))   
! this is harmonic C(LM) (Hansen, Morse-Feshbach) 
!               CLMM(I)=-ci* SIGNUM*CONJG(CLMM(I))    
10          CONTINUE   
         ENDIF
! orthogonality property         
!         SS=(0.0,0.0)
!         DO 12 I=1,3
!            SS=SS+rUnit(I)*CLMM(I)
!12       CONTINUE  
!            CYY=SS
!      write(*,*) 'Y1LM:CYY=',CYY     

      RETURN
      END SUBROUTINE Y1LM
 !----------------------------------------------------          

      SUBROUTINE DYLM_TETA (L,M,TETA,PHI,CKLM)
!      Calculate the derivative of the spherical harmonics 
!      over TETA 
!      CKLM - (real,image) part
!      the relation between (l,m) and lm is given by 
!      lm = l*(l+1) + m +1 
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
       REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

          INTEGER :: L,M
          REAL(wp) :: TETA,PHI
          COMPLEX(wp) :: CKLM
!functions
          COMPLEX(wp) :: CYLM,CCEXP
          REAL(wp) :: PLGNDR,PNORM
! local
          COMPLEX(wp) ::  PH1
          REAL(wp) :: X,Y,SIGNUM,PNN,PLAG,PLAG1 
          INTEGER :: MM,KK

      CCEXP(X)=CMPLX(COS(X),SIN(X))

      MM=ABS(M)
      KK=MOD(MM,2)

! L'hopitals rule is used
      IF(TETA .EQ. PI .OR. TETA .EQ. 0.) THEN
        TETA=PI/360.    ! CKLM=(0.0,0.0)
      ENDIF   

      X=COS(TETA) 
      Y=SIN(TETA)

!            PNN=PNORM(L,MM)
      PNN=PNORM(L,MM)/SQRT(2.*PI)
      PLAG=PLGNDR(L,MM,X)
      PLAG1=PLGNDR(L+1,MM,X)

!      CKLM=PNN/Y*((L-MM+1)*PLAG1-(L+1)*X*PLAG)* &
!          CCEXP(MM*PHI)

      IF(M .LT. 0) THEN 
           IF(KK .EQ. 1) THEN
              SIGNUM=-1.
           ELSE
              SIGNUM=1.
           ENDIF  
         CKLM=SIGNUM*PNN/Y*((L-MM+1)*PLAG1-(L+1)*X*PLAG)* &
          CONJG(CCEXP(MM*PHI))
      ELSE
         CKLM=PNN/Y*((L-MM+1)*PLAG1-(L+1)*X*PLAG)* &
             CCEXP(MM*PHI)
      ENDIF
      RETURN
      END SUBROUTINE DYLM_TETA
!-----------------------------------------------------------
      SUBROUTINE DYLM_PHI (L,M,TETA,PHI,CKLM)
!      Calculate the derivative of the spherical harmonics 
!      over PHI 
!      CKLM - (real,image) part
!      the relation between (l,m) and lm is given by 
!      lm = l*(l+1) + m +1 
      
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
!        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

          INTEGER :: L,M
          REAL(wp) :: TETA,PHI
          COMPLEX(wp) :: CKLM
!functions
          COMPLEX(wp) :: CYLM
!local 
          COMPLEX(wp),PARAMETER :: ci = (0.,1.)

      CKLM=ci*M*CYLM (L,M,TETA,PHI)
      RETURN
      END SUBROUTINE DYLM_PHI

!----------------------------------------------------------------

      SUBROUTINE FourProjec(XNU,YNU,NX,NY,L,M1,M2,LMAX,N,GN,R0,  &
                            NJ,NALPHA,NBETA,NGAMMA,ALPHAM,BETAM,GAMMAM, &
                            FPROJ)
! Fourier projections
! N- is N_th root of spher. bessel function
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

          INTEGER :: NX,NY,NZ
          INTEGER :: L,M1,M2,LMAX,N
          INTEGER :: NJ,NALPHA,NBETA,NGAMMA
          REAL(wp) :: XNU(NX),YNU(NY)
          REAL(wp) :: HNUX,HNUY
          REAL(wp) :: GN,R0
          COMPLEX(wp) :: FPROJ(NX*NY*NJ)
          REAL(wp) :: ALPHAM(NALPHA),BETAM(NBETA),GAMMAM(NGAMMA)          
! local 
          INTEGER :: I1,I2,I21,I4,I421
          INTEGER :: J1,J2,J3,J321
          REAL(wp) :: XX,YY,RRR
          REAL(wp) :: TETA,PHI
          REAL(wp) :: FCSWF,PLG
          REAL(wp) :: ALPHA,BETA,GAMMA
          COMPLEX(wp) :: WIG,SS,CYY
          REAL(wp) :: Y
! functions 
          COMPLEX(wp) :: FCoefSWF,D_WIGNER,ylm
          REAL(wp) :: PLGNDR_NORM
          COMPLEX(wp) :: CCEXP,CYLM

!          CCEXP(Y)=CMPLX(COS(Y),SIN(Y))

        TETA=PI/2.0
  
        DO 6 I4=1,NJ
           CALL IALBETGAM(I4,NALPHA,NBETA,J1,J2,J3)

!           J321=(J3-1)*NBETA*NALPHA+(J2-1)*NALPHA+J1
           ALPHA=ALPHAM(J1)
           BETA =BETAM(J2)
           GAMMA=GAMMAM(J3)

        DO 6 I2=1,NY
           YY=YNU(I2)
        DO 6 I1=1,NX
           XX=XNU(I1)
           I21=(I2-1)*NX +I1

           I421=(I4-1)*NX*NY+I21

           RRR=SQRT(XX**2+YY**2)

           CALL POLAR2D(XX,YY,PHI)

        SS=0.0   
        DO 4 L=0,LMAX   

           FCSWF= FCoefSWF(L,N,RRR,GN,R0)

        DO 4 M2=-L,L   
        DO 4 M1=-L,L
   
!           WIG= D_WIGNER(L,M1,M2,ALPHA,BETA,GAMMA)
           WIG= D_WIGNER(L,M2,M1,ALPHA,BETA,GAMMA)

!          CYY=ylm(L, M2, TETA, PHI)  
          CYY= CYLM (L,M2,TETA,PHI)

          SS=SS+FCSWF*WIG*CYY

4       CONTINUE
        
        FPROJ(I421)=SS   

6       CONTINUE

      RETURN
      END SUBROUTINE FourProjec
!----------------------------------------------------------------
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

      FUNCTION FCoefSWF(L,N,KK,GN,R0)
!  Exact 3D Fourier transform of the Spherical Wave Functions 
!  (with zero boundary condition on a sphere radius of R0) 
!  GN -is the N-th root of the L harmonics  
!  R0 - radius of the sphere
!  KK - current variable in a Fourier space 

        USE PARAMTR    ! QMAX, is described in PARAM module 
!                        QMAX must be >= L1

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 
        complex(wp), parameter :: eye = (0.0_wp,1.0_wp)  !SQRT(-1)

          REAL(wp) :: GN,KK,R0
          INTEGER :: L,N
! function
          COMPLEX(wp) :: ylm,FCoefSWF
! local
          COMPLEX(wp) :: CRAB
          REAL(wp) :: RAB,SHR1,SHR2,GNN
          REAL(wp) :: BES(QMAX)
!          COMPLEX(wp) :: YY
          INTEGER :: ERR
          
          GNN=GN/R0

!          CALL ROOTS_BESSEL_2(ROOTB,SMAXIMUM,LMAXIMUM)

!          GN=ROOTB(N,L+1)  ! N -is number of the root of the L-order harmonic

          CALL SpherBessJ(KK*R0,QMAX,BES,ERR) 
          SHR1= BES(L+1) 

          CALL SpherBessJ(GNN*R0,QMAX,BES,ERR) 
          SHR2= BES(L+2)    

          IF(KK .NE. GNN) THEN
             RAB= -GNN*R0**2*SHR1*SHR2/(KK**2-GNN**2)
          ELSE
             RAB=R0**3/2.0*SHR2**2
          ENDIF

!          YY=ylm(L, M, TETA, PHI)  
          FCoefSWF=4.0*PI*eye**L*RAB  
        RETURN
        END FUNCTION FCoefSWF
!--------------------------------------------

      FUNCTION Fourier3DSWF(L,M,N,RR,TETA,PHI,GN,R0)
!  Exact 3D Fourier transform of the Spherical Wave Functions 
!  (with zero boundary condition on a sphere radius of R0) 
!  GN -is the N-th root of the L harmonics  
!  R0 - radius of the sphere
!  RR - current variable in a Fourier space 
!!  Q - может потребоваться для задачи в неограниченной области

        USE PARAMTR    ! QMAX, is described in PARAM module 
!                        QMAX must be >= L1

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 
        complex(wp), parameter :: eye = (0.0_wp,1.0_wp)  !SQRT(-1)

          REAL(wp) :: GN,RR,TETA,PHI,R0
          INTEGER :: L,M,N
! function
          COMPLEX(wp) :: ylm,Fourier3DSWF
! local
          COMPLEX(wp) :: CRAB
          REAL(wp) :: RAB,SHR1,SHR2
          REAL(wp) :: BES(QMAX)
          COMPLEX(wp) :: YY
          INTEGER :: ERR

!          CALL ROOTS_BESSEL_2(ROOTB,SMAXIMUM,LMAXIMUM)

!          GN=ROOTB(N,L+1)  ! N -is number of the root of the L-order harmonic

          CALL SpherBessJ(RR*R0,QMAX,BES,ERR) 
          SHR1= BES(L+1) 

          CALL SpherBessJ(GN*R0,QMAX,BES,ERR) 
          SHR2= BES(L+2)    

          IF(RR .NE. GN) THEN
             RAB= -GN*R0**2*SHR1*SHR2/(RR**2-GN**2)
          ELSE
             RAB=R0**3/2.0*SHR2**2
          ENDIF

          YY=ylm(L, M, TETA, PHI)  
          Fourier3DSWF=4.0*PI*eye**L*YY*RAB  
        RETURN
        END FUNCTION Fourier3DSWF
!--------------------------------------------

        FUNCTION SWF(L1,M1,Q,RR,TETA,PHI)
! Spherical wave-functions \SWF=j_L1(Q*RR)*Y_L1_M1(TETA,PHI)

          USE PARAMTR    ! QMAX, are described in PARAM module 
!                        QMAX must be >= L1
          IMPLICIT NONE
          INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
          INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=double 
!          REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

          REAL(wp) :: Q,RR,TETA,PHI
          INTEGER :: L1,M1
! functions
          COMPLEX(wp) :: ylm,SWF
! local 
          REAL(wp) :: BES(QMAX)
          COMPLEX(wp) :: YY
          INTEGER :: ERR

             CALL SpherBessJ(Q*RR,QMAX,BES,ERR) 

             YY=ylm(L1, M1, TETA, PHI)

             SWF= BES(L1+1) * YY   

        RETURN
        END FUNCTION SWF
!--------------------------------------------

      FUNCTION D_WIGNER(L,M1,M2,ALPHA,BETA,GAMMA)
! Wigner D-functions

      IMPLICIT NONE
      INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
      INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single 
      REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 
     
      INTEGER :: L,M1,M2
      REAL(wp) :: ALPHA,BETA,GAMMA
! functions
      COMPLEX(wp) :: D_WIGNER
      REAL(wp) :: PLM1M2
! local
      COMPLEX(wp) :: CCEXP
      REAL(wp) :: X

      CCEXP(X)=CMPLX(COS(X),-SIN(X))

      D_WIGNER=CCEXP(M1*ALPHA+M2*GAMMA) * PLM1M2(L,M1,M2,BETA) ! &
!              * SQRT((2.0*L+1.)/8.0*PI)

      RETURN
      END FUNCTION D_WIGNER

!-------------------------------------------------------------

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



      SUBROUTINE ROOTS_BESSEL_2(ROOTB,SMAXIMUM,LMAXIMUM)

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single 
!        INTEGER :: LMAX,SMAX
        INTEGER :: SMAXIMUM,LMAXIMUM
!        REAL(wp) :: ROOTBB(SMAXIMUM,LMAX+1)

! local
        REAL(wp) ::ROOTB(SMAXIMUM,LMAXIMUM) 
        INTEGER ::I,J,JI

!        SMAXIMUM=20
!        LMAXIMUM=10
        
!        DO 2 J=1,LMAX
!        DO 2 I=1,SMAX
!           JI=(J-1)*SMAX+I
 
        IF(LMAXIMUM .GT. 11 .OR. SMAXIMUM .GT. 20) THEN
         WRITE(*,*) 'Procedure ROOTS_BESSEL_2:Error','LMAX=',LMAXIMUM,'SMAX=',SMAXIMUM 
         STOP
        ENDIF
   
!  j_0
        ROOTB(1,1)= 3.14159_wp
        ROOTB(2,1)= 6.28319_wp
        ROOTB(3,1)= 9.42478_wp
        ROOTB(4,1)= 12.56637_wp
        ROOTB(5,1)= 15.70796_wp
        ROOTB(6,1)= 18.84955_wp
        ROOTB(7,1)= 21.99115_wp
        ROOTB(8,1)= 25.13274_wp
        ROOTB(9,1)= 28.27433_wp
        ROOTB(10,1)= 31.41593_wp
        ROOTB(11,1)= 34.55752_wp
        ROOTB(12,1)= 37.69911_wp
        ROOTB(13,1)= 40.84070_wp
        ROOTB(14,1)= 43.98230_wp
        ROOTB(15,1)= 47.12389_wp
        ROOTB(16,1)= 50.26548_wp
        ROOTB(17,1)= 53.40708_wp
        ROOTB(18,1)= 56.54867_wp
        ROOTB(19,1)= 59.69026_wp
        ROOTB(20,1)= 62.83185_wp

! j_1
        ROOTB(1,2)= 4.49341_wp
        ROOTB(2,2)= 7.72525_wp
        ROOTB(3,2)= 10.90412_wp
        ROOTB(4,2)= 14.06619_wp
        ROOTB(5,2)= 17.22076_wp
        ROOTB(6,2)= 20.37130_wp
        ROOTB(7,2)= 23.51945_wp
        ROOTB(8,2)= 26.66605_wp
        ROOTB(9,2)= 29.81160_wp
        ROOTB(10,2)= 32.95639_wp
        ROOTB(11,2)= 36.10062_wp
        ROOTB(12,2)= 39.24443_wp
        ROOTB(13,2)= 42.38791_wp
        ROOTB(14,2)= 45.53113_wp
        ROOTB(15,2)= 48.67414_wp
        ROOTB(16,2)= 51.81698_wp
        ROOTB(17,2)= 54.95968_wp
        ROOTB(18,2)= 58.10225_wp
        ROOTB(19,2)= 61.24473_wp
        ROOTB(20,2)= 64.38712_wp
!j_2
        ROOTB(1,3)= 5.76346_wp
        ROOTB(2,3)= 9.09501_wp
        ROOTB(3,3)= 12.32294_wp
        ROOTB(4,3)= 15.51460_wp
        ROOTB(5,3)= 18.68904_wp
        ROOTB(6,3)= 21.85387_wp
        ROOTB(7,3)= 25.01280_wp
        ROOTB(8,3)= 28.16783_wp
        ROOTB(9,3)= 31.32014_wp
        ROOTB(10,3)= 34.47049_wp
        ROOTB(11,3)= 37.61937_wp
        ROOTB(12,3)= 40.76712_wp
        ROOTB(13,3)= 43.91398_wp
        ROOTB(14,3)= 47.06014_wp
        ROOTB(15,3)= 50.20573_wp
        ROOTB(16,3)= 53.35084_wp
        ROOTB(17,3)= 56.49557_wp
        ROOTB(18,3)= 59.63996_wp
        ROOTB(19,3)= 62.78407_wp
        ROOTB(20,3)= 65.92794_wp

!j_3
        ROOTB(1,4)= 6.98793_wp
        ROOTB(2,4)= 10.41712_wp
        ROOTB(3,4)= 13.69802_wp
        ROOTB(4,4)= 16.92362_wp
        ROOTB(5,4)= 20.12181_wp
        ROOTB(6,4)= 23.30425_wp
        ROOTB(7,4)= 26.47676_wp
        ROOTB(8,4)= 29.64260_wp
        ROOTB(9,4)= 32.80373_wp
        ROOTB(10,4)= 35.96141_wp
        ROOTB(11,4)= 39.11647_wp
        ROOTB(12,4)= 42.26951_wp
        ROOTB(13,4)= 45.42096_wp
        ROOTB(14,4)= 48.57113_wp
        ROOTB(15,4)= 51.72025_wp
        ROOTB(16,4)= 54.86850_wp
        ROOTB(17,4)= 58.01603_wp
        ROOTB(18,4)= 61.16295_wp
        ROOTB(19,4)= 64.30934_wp
        ROOTB(20,4)= 67.45528_wp
!j_4
        ROOTB(1,5)= 8.18256_wp
        ROOTB(2,5)= 11.70491_wp
        ROOTB(3,5)= 15.03966_wp
        ROOTB(4,5)= 18.30126_wp
        ROOTB(5,5)= 21.52542_wp
        ROOTB(6,5)= 24.72757_wp
        ROOTB(7,5)= 27.91558_wp
        ROOTB(8,5)= 31.09393_wp
        ROOTB(9,5)= 34.26539_wp
        ROOTB(10,5)= 37.43174_wp
        ROOTB(11,5)= 40.59419_wp
        ROOTB(12,5)= 43.75361_wp
        ROOTB(13,5)= 46.91061_wp
        ROOTB(14,5)= 50.06565_wp
        ROOTB(15,5)= 53.21910_wp
        ROOTB(16,5)= 56.37121_wp
        ROOTB(17,5)= 59.52220_wp
        ROOTB(18,5)= 62.67225_wp
        ROOTB(19,5)= 65.82148_wp
        ROOTB(20,5)= 68.97001_wp
!j_5
        ROOTB(1,6)= 9.35581_wp
        ROOTB(2,6)= 12.96653_wp
        ROOTB(3,6)= 16.35471_wp
        ROOTB(4,6)= 19.65315_wp
        ROOTB(5,6)= 22.90455_wp
        ROOTB(6,6)= 26.12775_wp
        ROOTB(7,6)= 29.33256_wp
        ROOTB(8,6)= 32.52466_wp
        ROOTB(9,6)= 35.70758_wp
        ROOTB(10,6)= 38.88363_wp
        ROOTB(11,6)= 42.05442_wp
        ROOTB(12,6)= 45.22107_wp
        ROOTB(13,6)= 48.38440_wp
        ROOTB(14,6)= 51.54505_wp
        ROOTB(15,6)= 54.70348_wp
        ROOTB(16,6)= 57.86006_wp
        ROOTB(17,6)= 61.01508_wp
        ROOTB(18,6)= 64.16878_wp
        ROOTB(19,6)= 67.32133_wp
        ROOTB(20,6)= 70.47290_wp
!j_6
        ROOTB(1,7)= 10.51284_wp
        ROOTB(2,7)= 14.20739_wp
        ROOTB(3,7)= 17.64797_wp
        ROOTB(4,7)= 20.98346_wp
        ROOTB(5,7)= 24.26277_wp
        ROOTB(6,7)= 27.50787_wp
        ROOTB(7,7)= 30.73038_wp
        ROOTB(8,7)= 33.93711_wp
        ROOTB(9,7)= 37.13233_wp
        ROOTB(10,7)= 40.31889_wp
        ROOTB(11,7)= 43.49876_wp
        ROOTB(12,7)= 46.67333_wp
        ROOTB(13,7)= 49.84366_wp
        ROOTB(14,7)= 53.01050_wp
        ROOTB(15,7)= 56.17448_wp
        ROOTB(16,7)= 59.33604_wp
        ROOTB(17,7)= 62.49557_wp
        ROOTB(18,7)= 65.65336_wp
        ROOTB(19,7)= 68.80966_wp
        ROOTB(20,7)= 71.96465_wp
!j_7
        ROOTB(1,8)= 11.65703_wp
        ROOTB(2,8)= 15.43129_wp
        ROOTB(3,8)= 18.92300_wp
        ROOTB(4,8)= 22.29535_wp
        ROOTB(5,8)= 25.60286_wp
        ROOTB(6,8)= 28.87037_wp
        ROOTB(7,8)= 32.11120_wp
        ROOTB(8,8)= 36.33319_wp
        ROOTB(9,8)= 38.54136_wp
        ROOTB(10,8)= 41.73905_wp
        ROOTB(11,8)= 44.92859_wp
        ROOTB(12,8)= 48.11165_wp
        ROOTB(13,8)= 51.28949_wp
        ROOTB(14,8)= 54.46304_wp
        ROOTB(15,8)= 57.63302_wp
        ROOTB(16,8)= 60.80001_wp
        ROOTB(17,8)= 63.96446_wp
        ROOTB(18,8)= 67.12673_wp
        ROOTB(19,8)= 70.28713_wp
        ROOTB(20,8)= 73.44590_wp

!j_8
        ROOTB(1,9)= 12.79078_wp
        ROOTB(2,9)= 16.64100_wp
        ROOTB(3,9)= 20.18247_wp
        ROOTB(4,9)= 23.59127_wp
        ROOTB(5,9)= 26.92704_wp
        ROOTB(6,9)= 30.21726_wp
        ROOTB(7,9)= 33.47680_wp
        ROOTB(8,9)= 36.71453_wp
        ROOTB(9,9)= 39.93613_wp
        ROOTB(10,9)= 43.14543_wp
        ROOTB(11,9)= 46.34511_wp
        ROOTB(12,9)= 49.53712_wp
        ROOTB(13,9)= 52.72290_wp
        ROOTB(14,9)= 55.90356_wp
        ROOTB(15,9)= 59.07995_wp
        ROOTB(16,9)= 62.25274_wp
        ROOTB(17,9)= 65.42247_wp
        ROOTB(18,9)= 68.58956_wp
        ROOTB(19,9)= 71.75438_wp
        ROOTB(20,9)= 74.91722_wp

!j_9
        ROOTB(1,10)= 13.91582_wp
        ROOTB(2,10)= 17.83864_wp
        ROOTB(3,10)= 21.42849_wp
        ROOTB(4,10)= 24.87321_wp
        ROOTB(5,10)= 28.23713_wp
        ROOTB(6,10)= 31.55019_wp
        ROOTB(7,10)= 34.82870_wp
        ROOTB(8,10)= 38.08248_wp
        ROOTB(9,10)= 41.31786_wp
        ROOTB(10,10)= 44.53914_wp
        ROOTB(11,10)= 47.74935_wp
        ROOTB(12,10)= 50.95067_wp
        ROOTB(13,10)= 54.14477_wp
        ROOTB(14,10)= 57.33289_wp
        ROOTB(15,10)= 60.51602_wp
        ROOTB(16,10)= 63.69493_wp
        ROOTB(17,10)= 66.87024_wp
        ROOTB(18,10)= 70.04245_wp
        ROOTB(19,10)= 73.21197_wp
        ROOTB(20,10)= 76.37914_wp

!j_10
        ROOTB(1,11)= 15.03347_wp
        ROOTB(2,11)= 19.02585_wp
        ROOTB(3,11)= 22.66272_wp
        ROOTB(4,11)= 26.14277_wp
        ROOTB(5,11)= 29.53463_wp
        ROOTB(6,11)= 32.87053_wp
        ROOTB(7,11)= 36.16816_wp
        ROOTB(8,11)= 39.43821_wp
        ROOTB(9,11)= 42.68765_wp
        ROOTB(10,11)= 45.92120_wp
        ROOTB(11,11)= 49.14222_wp
        ROOTB(12,11)= 52.35316_wp
        ROOTB(13,11)= 55.55587_wp
        ROOTB(14,11)= 58.75175_wp
        ROOTB(15,11)= 61.94190_wp
        ROOTB(16,11)= 65.12721_wp
        ROOTB(17,11)= 68.30836_wp
        ROOTB(18,11)= 71.48594_wp
        ROOTB(19,11)= 74.66040_wp
        ROOTB(20,11)= 77.83215_wp

!        DO 2 J=1,LMAX+1
!        DO 2 I=1,SMAX
!!           JI=(J-1)*SMAX+I
!           ROOTBB(I,J)=ROOTB(I,J)
!2       CONTINUE         
           
      RETURN  
      END SUBROUTINE  ROOTS_BESSEL_2
!-------------------------------------------------------------

      FUNCTION S_Bessel_J(XXX,L1,ERR)
! calculate spherical Bessel functions
! BES(1) -> J_(1/2),  BES(2) -> J_(3/2) ans so on
      IMPLICIT NONE
      INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
      INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single 
      REAL(wp) :: XXX
      INTEGER :: L1,ERR 
      REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 
      REAL(wp), PARAMETER :: AL=0.5_wp
! local
!      REAL(wp) :: BES(L)
      INTEGER :: I,L
      REAL(wp) :: S_Bessel_J, XX

      REAL(wp), allocatable:: BES(:)
      
      L=L1+1

      ALLOCATE (BES(L))

      XX=ABS(XXX)

      IF(XX .EQ. 0.0_wp) XX=0.0000001_wp

      CALL rjbesl(XX,AL,L,BES,ERR)

!      DO 2 I=L,L 
         BES(L) = BES(L)*SQRT(PI/(2.0_wp*XX))
!         S_Bessel_J=BES(I)
!2     CONTINUE   
         S_Bessel_J=BES(L)

      IF(XXX .LT. 0.0_wp) THEN
!      DO 4 I=L,L 
         BES(L) = BES(L)*(-1)**(L-1)
         S_Bessel_J=BES(L)
!4     CONTINUE  
      ENDIF   
      RETURN
      END FUNCTION S_Bessel_J
!------------------------------------------------------------

      SUBROUTINE NM_SPHER(JJ,NN,MM) 
! computer NN and MM by use of single index JJ    
!        USE PARAMTR    
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0) !working precision=double 
        INTEGER :: JJ,NN,MM
        REAL(wp) RJJ,RNN

        RJJ=JJ

        NN=CEILING(SQRT(RJJ))-1

        RNN=NN

        MM=INT((RJJ-RNN**2)/2.0)

      RETURN
      END SUBROUTINE NM_SPHER
!--------------------------------------------------------------

      SUBROUTINE XYANG_1(am,gm,fm1,n1,n2,niter)
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

      CC=1.0_wp
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

      IF(q2 .NE. 0.0_wp) GO TO 12
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
      END SUBROUTINE XYANG_1  
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
      FUNCTION  BessJ(X,N,NB,ERR)
!   calculate  Bessel functions J_N(X)
!   max(N) = NB-1
      IMPLICIT NONE
      INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
      INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single 
      INTEGER :: N,NB,ERR 
      REAL(wp) :: X
      REAL(wp) :: BessJ
      REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 
      REAL(wp), PARAMETER :: AL=0.0_wp
! local
      INTEGER :: I,NN
      REAL(wp) :: XX
      REAL(wp) :: BES(NB)

      XX=ABS(X)
      NN=ABS(N)

      CALL rjbesl(XX,AL,NB,BES,ERR)

      IF(X .LT. 0.0 .OR. N .LT. 0) THEN
         DO 4 I=1,NB 
            BES(I) = BES(I)*(-1)**(I-1)
4        CONTINUE  
      ENDIF   
      BessJ=BES(NN+1)
      RETURN
      END FUNCTION BessJ
!------------------------------------------------------------  
      SUBROUTINE SpherBessJ(XXX,NB,BES,ERR)
! calculate spherical Bessel functions
      IMPLICIT NONE
      INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
      INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single 
      INTEGER :: NB,ERR
      REAL(wp) :: XXX,BES(NB)
      REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 
      REAL(wp), PARAMETER :: AL=0.5_wp
      INTEGER :: I
      REAL(wp) :: XX,NN

      XX=ABS(XXX)
      

      IF(XX .EQ. 0.0_wp) XX=0.0000001_wp

      CALL rjbesl(XX,AL,NB,BES,ERR)

      DO 2 I=1,NB 
         BES(I) = BES(I)*SQRT(PI/(2.0_wp*XX))
2     CONTINUE   
      IF(XXX .LT. 0.0_wp) THEN
         DO 4 I=1,NB 
            BES(I) = BES(I)*(-1)**(I-1)
4        CONTINUE  
      ENDIF   
      RETURN
      END SUBROUTINE SpherBessJ
!------------------------------------------------------------

      FUNCTION ylm (l, m, thrad, phirad)
        implicit none
!
! Computes the spherical harmonic Y_lm (theta,phi) using the
! reduced rotation matrix d^l_{m 0} (theta) and using the
! external function fac10(n) = factorial(n)/10**n
!
! input: angular momentum quantum numbers l, m (integers)
!        angles theta and phi (radian)
!
! Reference: D.M. Brink and G.R. Satchler, Angular Momentum,
!            second edition, Oxford University Press, p.22 and p. 145
!
        integer, parameter  :: wp = kind(1.0d0)  ! working precision = double 
!
        real(wp), parameter :: pi = 3.14159265358979323846_wp    ! Schaum's Math handbook
        complex(wp), parameter :: eye = (0.0_wp,1.0_wp)
!
!   formal arguments
!
        integer, intent(in)  :: l, m
        real(wp), intent(in) :: thrad, phirad
!
!   local variables
!
        integer :: itmin1, itmin2, itmin, itmax1, itmax2, itmax, it, iphase, &
             ia, ib, ic
        real(wp) :: sqrt_fac, sumt, denom, term, dlm0, const, cosb2, sinb2
        complex(wp) :: ylm, exphi
!
!   external function
!
        real(wp), external :: fac10
!
!  program starts here
!  first calculate d^l_{m 0} (theta)
!
        cosb2 = cos(thrad/2.0_wp)
        sinb2 = sin(thrad/2.0_wp)
!
! determine lower and upper limits for summation index it; these
! are derived from the requirement that all factorials n! in the
! denominator are restricted to values with n >=0.
!
        itmin1 = 0
        itmin2 = m
        itmin = max(itmin1,itmin2)
        itmax1 = l+m
        itmax2 = l
        itmax = min(itmax1,itmax2)
!  write (6,'(10X,A,2I6)') ' itmin, itmax = ', itmin, itmax
        sqrt_fac = sqrt( fac10(l+m) * fac10(l-m) * fac10(l) * fac10(l) )
!
        sumt = 0.0_wp
        do it = itmin, itmax
           iphase = (-1)**it
           ia = l + m - it
           ib = l - it
           ic = it - m
!     write (6,'(10X,A,5I6)') ' it, iphase, ia, ib, ic  = ', it, iphase, ia, ib, ic
           denom = fac10(ia) * fac10(ib) * fac10(it) * fac10(ic)
           term = iphase * cosb2**(ia+ib) * sinb2**(it+ic) / denom
           sumt = sumt + term
        end do
        dlm0 = sqrt_fac * sumt
!
!  now compute Y_{l m} (theta,phi) from d^l_{m 0} (theta)
!
        const = sqrt( (2.0_wp *l + 1.0_wp) / (4.0_wp * pi) )
        exphi = exp( eye * m * phirad )
        ylm = const * exphi * dlm0
!
      return
      END FUNCTION ylm
!----------------------------------------------------------
      function fac10 (n)
        implicit none
! function fac10(n) calculates factorial(n)/10**n
! input: integer n >= 0 (you may want to check this
!        in the program calling this function)

        integer, parameter :: wp = kind(1.0d0)  ! working precision = double 

!      formal arguments
        integer, intent(in) :: n

!      local variables
        integer :: i
        real(wp) :: fac10, q
! 
        if (n == 0) then
           fac10 = 1.0_wp
        else
           fac10 = 1.0_wp
           q = 1.0_wp
           do i = 1, n
              fac10 = fac10 * q / 10.0_wp
              q = q + 1.0_wp
           end do
        endif
!
        return
      end function fac10
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

      subroutine rjbesl ( x, alpha, nb, b, ncalc )

! RJBESL calculates J Bessel function with non-integer orders.
!
!  Discussion:
!
!    This routine calculates Bessel functions J sub(N+ALPHA) (X)
!    for non-negative argument X, and non-negative order N+ALPHA.
!

!  Parameters:
!
!    Input, double precision X, the argument for which the
!    J_s are to be calculated.
!
!    Input, double precision ALPHA, the fractional part of order for which
!    the J's or exponentially scaled J's (J*exp(X)) are to be calculated.  
!    0 <= ALPHA < 1.0.
!
!    Input, integer NB, the number of functions to be calculated.
!    0 < NB.  The first function calculated is of order ALPHA, and the
!    last is of order (NB - 1 + ALPHA).
!
!    Output, double precision B(NB).  If RJBESL terminates normally, with
!    NCALC = NB, then B contains the functions J/ALPHA/(X) through 
!    J/NB-1+ALPHA/(X), or the corresponding exponentially scaled functions.
!
!    Output, integer NCALC, error indicator.  If NCALC = NB, then all the 
!    requested values were calculated to the desired accuracy.  
!
!  Local Parameters:
!
!    IT, the number of bits in the mantissa of a working precision
!    variable.
!
!    NSIG, the decimal significance desired.  Should be set to
!    INT(LOG10(2)*IT+1).  Setting NSIG lower will result
!    in decreased accuracy while setting NSIG higher will
!    increase CPU time without increasing accuracy.  The
!    truncation error is limited to a relative error of
!    T=.5*10**(-NSIG).
!
!    ENTEN = 10.0**K, where K is the largest integer such that
!    ENTEN is machine-representable in working precision
!
!    ENSIG = 10.0**NSIG
!
!    RTNSIG = 10.0**(-K) for the smallest integer K such that
!    K .GE. NSIG/4
!
!    ENMTEN, the smallest ABS(X) such that X/4 does not underflow
!
!    XLARGE, the upper limit on the magnitude of X.  If ABS(X)=N,
!    then at least N iterations of the backward recursion
!    will be executed.  The value of 10000.0 is used on
!    every machine.
!
      implicit none

      integer nb

      double precision alpha
      double precision alpem
      double precision alp2em
      double precision b(nb)
      double precision capp
      double precision capq
      double precision eighth
      double precision em
      double precision en
      double precision enmten
      double precision ensig
      double precision enten
      double precision fact(25)
      double precision four
      double precision gnu
      double precision half
      double precision halfx
      integer i
      integer j
      integer k
      integer l
      integer m
      integer magx
      integer n
      integer nbmx
      integer ncalc
      integer nend
      integer nstart
      double precision one
      double precision one30
      double precision p
      double precision pi2
      double precision plast
      double precision pold
      double precision psave
      double precision psavel
      double precision r8_gamma
      double precision rtnsig
      double precision s
      double precision sum
      double precision t
      double precision t1
      double precision tempa
      double precision tempb
      double precision tempc
      double precision test
      double precision three
      double precision three5
      double precision tover
      double precision two
      double precision twofiv
      double precision twopi1
      double precision twopi2
      double precision x
      double precision xc
      double precision xin
      double precision xk
      double precision xlarge
      double precision xm
      double precision vcos
      double precision vsin
      double precision z
      double precision zero
!
!  Mathematical constants
!
!   PI2    - 2 / PI
!   TWOPI1 - first few significant digits of 2 * PI
!   TWOPI2 - (2*PI - TWOPI) to working precision, i.e.,
!            TWOPI1 + TWOPI2 = 2 * PI to extra precision.
!
      data pi2 / 0.636619772367581343075535d0 /
      data twopi1 / 6.28125d0 /
      data twopi2 / 1.935307179586476925286767d-3 /
      data zero /0.0d0 /
      data eighth / 0.125d0 /
      data half / 0.5d0 /
      data one / 1.0d0/
      data two /2.0d0 / 
      data three / 3.0d0 /
      data four / 4.0d0 /
      data twofiv /25.0d0/
      data one30 /130.0d0 /
      data three5 / 35.0d0/
!
!  Machine-dependent parameters
!
      data enten /1.0d38 /
      data ensig / 1.0d17 /
      data rtnsig / 1.0d-4/
      data enmten /1.2d-37 /
      data xlarge / 1.0d4/
!
!  Factorial(N)
!
      data fact /  1.0d0,    &
      1.0d0,    &
      2.0d0,    &
      6.0d0,    &
      24.0d0,   &
      1.2d2,    &
      7.2d2,    & 
      5.04d3,   & 
      4.032d4,  &
      3.6288d5,3.6288d6,3.99168d7,4.790016d8,6.2270208d9,  &
      8.71782912d10,1.307674368d12,2.0922789888d13,3.55687428096d14, &
      6.402373705728d15,1.21645100408832d17,2.43290200817664d18,     &
      5.109094217170944d19,1.12400072777760768d21,   &
      2.585201673888497664d22,   &
      6.2044840173323943936d23/
!
!  Check for out of range arguments.
!
      magx = int ( x )

      if (  0 .lt. nb .and.   & 
       zero .le. x .and.      &
       x .le. xlarge .and.    &
       zero .le. alpha .and.  &
       alpha .lt. one ) then  
!
!  Initialize result array to zero.
!
        ncalc = nb
        do i = 1, nb
          b(i) = zero
        end do
!
!  Branch to use 2-term ascending series for small X and asymptotic
!  form for large X when NB is not too large.
!
        if ( x .lt. rtnsig ) then
!
!  Two-term ascending series for small X.
!
          tempa = one
          alpem = one + alpha
          halfx = zero

          if ( enmten .lt. x ) then
            halfx = half * x
          end if

          if ( alpha .ne. zero ) then
            tempa = halfx**alpha / ( alpha * r8_gamma ( alpha ) )
          end if

          tempb = zero

          if ( one .lt. x + one ) then
            tempb = -halfx * halfx
          end if

          b(1) = tempa + tempa * tempb / alpem

          if ( x .ne. zero .and. b(1) .eq. zero ) then
            ncalc = 0
          end if

          if ( nb .ne. 1 ) then

            if ( x .le. zero ) then

              do n = 2, nb
                b(n) = zero
              end do
!
!  Calculate higher order functions.
!
            else

              tempc = halfx
              tover = ( enmten + enmten ) / x

              if ( tempb .ne. zero ) then
                tover = enmten / tempb
              end if

              do n = 2, nb

                tempa = tempa / alpem
                alpem = alpem + one
                tempa = tempa * tempc

                if ( tempa .le. tover * alpem ) then
                  tempa = zero
                end if

                b(n) = tempa + tempa * tempb / alpem

                if ( b(n) .eq. zero .and. n .lt. ncalc ) then
                  ncalc = n - 1
                end if

              end do

            end if
          end if
!
!  Asymptotic series for 21 < X.
!
        else if ( twofiv .lt. x .and. nb .le. magx + 1 ) then

          xc = sqrt ( pi2 / x )
          xin = ( eighth / x )**2
          m = 11

          if ( x .ge. three5 ) then
            m = 8
          end if

          if ( x .ge. one30 ) then
            m = 4
          end if

          xm = four * dble ( m )
!
!  Argument reduction for SIN and COS routines.
!
          t = aint ( x / ( twopi1 + twopi2 ) + half )
          z = ( ( x - t * twopi1 ) - t * twopi2 )   & 
           - ( alpha + half ) / pi2
          vsin = sin ( z )
          vcos = cos ( z )
          gnu = alpha + alpha

          do i = 1, 2

            s = ( ( xm - one ) - gnu ) * ( ( xm - one ) + gnu )   & 
             * xin * half
            t = ( gnu - ( xm - three ) ) * ( gnu + ( xm - three ) )
            capp = s * t / fact(2*m+1)
            t1 = ( gnu - ( xm + one ) ) * ( gnu + ( xm + one ) )
            capq = s * t1 / fact(2*m+2)
            xk = xm
            k = m + m
            t1 = t

            do j = 2, m
              xk = xk - four
              s = ( ( xk - one ) - gnu ) * ( ( xk - one ) + gnu )
              t = ( gnu - ( xk - three ) ) * ( gnu + ( xk - three ) )
              capp = ( capp + one / fact(k-1) ) * s * t * xin
              capq = ( capq + one / fact(k) ) * s * t1 * xin
              k = k - 2
              t1 = t
            end do

            capp = capp + one
            capq = ( capq + one ) * ( gnu * gnu - one ) * ( eighth / x )
            b(i) = xc * ( capp * vcos - capq * vsin )

            if ( nb .eq. 1 ) then
              return
            end if

            t = vsin
            vsin = -vcos
            vcos = t
            gnu = gnu + two

          end do
!
!  If 2 < NB, compute J(X,ORDER+I)  I = 2, NB-1.
!
          if ( 2 .lt. nb ) then
            gnu = alpha + alpha + two
            do j = 3, nb
              b(j) = gnu * b(j-1) / x - b(j-2)
              gnu = gnu + two
            end do
          end if
!
!  Use recurrence to generate results.  First initialize the
!  calculation of P's.
!
        else

          nbmx = nb - magx
          n = magx + 1
          en = dble ( n + n ) + ( alpha + alpha )
          plast = one
          p = en / x
!
!  Calculate general significance test.
!
          test = ensig + ensig
!
!  Calculate P's until N = NB-1.  Check for possible overflow.
!
          if ( 3 .le. nbmx ) then

            tover = enten / ensig
            nstart = magx + 2
            nend = nb - 1
            en = dble ( nstart + nstart ) - two + ( alpha + alpha )

            do k = nstart, nend

              n = k
              en = en + two
              pold = plast
              plast = p
              p = en * plast / x - pold
!
!  To avoid overflow, divide P's by TOVER.  Calculate P's until
!  1 < ABS(P).
!
              if ( tover .lt. p ) then

                tover = enten
                p = p / tover
                plast = plast / tover
                psave = p
                psavel = plast
                nstart = n + 1

  100           continue

                n = n + 1
                en = en + two
                pold = plast
                plast = p
                p = en * plast / x - pold

                if ( p .le. one ) then
                  go to 100
                end if

                tempb = en / x
!
!  Calculate backward test and find NCALC, the highest N such that
!  the test is passed.
!
                test = pold * plast     & 
                 * ( half - half / ( tempb * tempb ) )
                test = test / ensig
                p = plast * tover
                n = n - 1
                en = en - two
                nend = min ( nb, n )

                do l = nstart, nend
                  pold = psavel
                  psavel = psave
                  psave = en * psavel / x - pold
                  if ( test .lt. psave * psavel ) then
                    ncalc = l - 1
                    go to 190
                  end if
                end do

                ncalc = nend
                go to 190

              end if

            end do

            n = nend
            en = dble ( n + n ) + ( alpha + alpha )
!
!  Calculate special significance test for 2 < NBMX.
!
            test = max ( test, sqrt ( plast * ensig ) * sqrt ( p + p ) )

          end if
!
!  Calculate P's until significance test passes.
!
  140     continue

          n = n + 1
          en = en + two
          pold = plast
          plast = p
          p = en * plast / x - pold
          if ( p .lt. test ) then
            go to 140
          end if
!
!  Initialize the backward recursion and the normalization sum.
!
  190     continue

          n = n + 1
          en = en + two
          tempb = zero
          tempa = one / p
          m = 2 * n - 4 * ( n / 2 )
          sum = zero
          em = dble ( n / 2 )
          alpem = ( em - one ) + alpha
          alp2em = ( em + em ) + alpha

          if ( m .ne. 0 ) then
            sum = tempa * alpem * alp2em / em
          end if

          nend = n - nb
!
!  Recur backward via difference equation, calculating (but not
!  storing) B(N), until N = NB.
!
          if ( 0 .lt. nend ) then

            do l = 1, nend

              n = n - 1
              en = en - two
              tempc = tempb
              tempb = tempa
              tempa = ( en * tempb ) / x - tempc
              m = 2 - m

              if ( m .ne. 0 ) then
                em = em - one
                alp2em = ( em + em ) + alpha
                if ( n .eq. 1 ) then
                  go to 210
                end if
                alpem = ( em - one ) + alpha
                if ( alpem .eq. zero ) then
                  alpem = one
                end if
                sum = ( sum + tempa * alp2em ) * alpem / em
              end if

            end do

          end if
!
!  Store B(NB).
!
  210     continue

          b(n) = tempa

          if ( nend .ge. 0 ) then

            if ( nb .le. 1 ) then

              alp2em = alpha
              if ( alpha + one .eq. one ) then
                alp2em = one
              end if
              sum = sum + b(1) * alp2em
              go to 250

            else
!
!  Calculate and store B(NB-1).
!
              n = n - 1
              en = en - two
              b(n) = ( en * tempa ) / x - tempb

              if ( n .eq. 1 ) then
                go to 240
              end if

              m = 2 - m

              if ( m .ne. 0 ) then
                em = em - one
                alp2em = ( em + em ) + alpha
                alpem = ( em - one ) + alpha
                if ( alpem .eq. zero ) then
                  alpem = one
                end if
                sum = ( sum + b(n) * alp2em ) * alpem / em
              end if

            end if

          end if

          nend = n - 2
!
!  Calculate via difference equation and store B(N), until N = 2.
!
          if ( nend .ne. 0 ) then

            do l = 1, nend
              n = n - 1
              en = en - two
              b(n) = ( en * b(n+1) ) / x - b(n+2)
              m = 2 - m
              if ( m .ne. 0 ) then
                em = em - one
                alp2em = ( em + em ) + alpha
                alpem = ( em - one ) + alpha
                if ( alpem .eq. zero ) then
                  alpem = one
                end if
                sum = ( sum + b(n) * alp2em ) * alpem / em
              end if
            end do

          end if
!
!  Calculate B(1).
!
          b(1) = two * ( alpha + one ) * b(2) / x - b(3)

  240     continue

          em = em - one
          alp2em = ( em + em ) + alpha

          if ( alp2em .eq. zero ) then
            alp2em = one
          end if

          sum = sum + b(1) * alp2em
!
!  Normalize.  Divide all B(N) by sum.
!
  250     continue

          if ( alpha + one .ne. one ) then
            sum = sum * r8_gamma ( alpha ) * ( x * half )**( -alpha )
          end if

          tempa = enmten
          if ( one .lt. sum ) then
            tempa = tempa * sum
          end if

          do n = 1, nb
            if ( abs ( b(n) ) .lt. tempa ) then
              b(n) = zero
            end if
            b(n) = b(n) / sum
          end do

        end if
!
!  Error return: X, NB, or ALPHA is out of range.
!
      else
        b(1) = zero
        ncalc = min ( nb, 0 ) - 1
      end if

      return
      end  subroutine rjbesl
!--------------------------------------------------------- 
       function r8_gamma ( x )

! R8_GAMMA evaluates Gamma(X) for a real argument.
!
!  Discussion:
!    This routine calculates the GAMMA function for a real argument X.
!    Computation is based on an algorithm outlined in reference 1.
!  Parameters:
!
!    Input, double precision X, the argument of the function.
!
!    Output, double precision R8_GAMMA, the value of the function.
!
      implicit none

      double precision r8_gamma
      double precision c(7)
      double precision eps
      double precision fact
      double precision half
      integer i
      integer n
      double precision one
      double precision p(8)
      logical parity
      double precision pi
      double precision q(8)
      double precision res
      double precision sqrtpi
      double precision sum
      double precision twelve
      double precision two
      double precision x
      double precision xbig
      double precision xden
      double precision xinf
      double precision xminin
      double precision xnum
      double precision y
      double precision y1
      double precision ysq
      double precision z
      double precision zero
!
!  Mathematical constants
!
      data one /1.0D+00 /
      data half /0.5D+00/
      data twelve /12.0D+00/
      data two /2.0D+00 /
      data zero /0.0D+00/
      data sqrtpi /0.9189385332046727417803297D+00/
      data pi /3.1415926535897932384626434D+00/
!
!  Machine dependent parameters
!
      data xbig / 171.624D+00 /
      data xminin / 2.23D-308 /
      data eps /2.22D-16/
      data xinf /1.79D+308/
!
!  Numerator and denominator coefficients for rational minimax
!  approximation over (1,2).
!
      data p/-1.71618513886549492533811d+00,  &
       2.47656508055759199108314d+01,  &
      -3.79804256470945635097577d+02,  &
       6.29331155312818442661052d+02,  &
       8.66966202790413211295064d+02,  &
      -3.14512729688483675254357d+04,  &
      -3.61444134186911729807069d+04,  &
       6.64561438202405440627855d+04/

      data q/-3.08402300119738975254353d+01,  &
       3.15350626979604161529144d+02,  &
      -1.01515636749021914166146d+03,  &
      -3.10777167157231109440444d+03,  &
       2.25381184209801510330112d+04,  &
       4.75584627752788110767815d+03,  &
      -1.34659959864969306392456d+05,  &
      -1.15132259675553483497211d+05/
!
!  Coefficients for minimax approximation over (12, INF).
!
      data c/-1.910444077728D-03,      &
       8.4171387781295D-04,            &
      -5.952379913043012D-04,          &
       7.93650793500350248D-04,        &
      -2.777777777777681622553D-03,    &
       8.333333333333333331554247D-02, &
       5.7083835261D-03/

      parity = .false.
      fact = one
      n = 0
      y = x
!
!  Argument is negative.
!
      if ( y .le. zero ) then

        y = - x
        y1 = aint ( y )
        res = y - y1

        if ( res .ne. zero ) then

          if ( y1 .ne. aint ( y1 * half ) * two ) then
            parity = .true.
          end if

          fact = - pi / sin ( pi * res )
          y = y + one

        else

          res = xinf
          r8_gamma = res
          return

        end if

      end if
!
!  Argument is positive.
!
      if ( y .lt. eps ) then
!
!  Argument < EPS.
!
        if ( xminin .le. y ) then
          res = one / y
        else
          res = xinf
          r8_gamma = res
          return
        end if

      else if ( y .lt. twelve ) then

        y1 = y
!
!  0.0 < argument < 1.0.
!
        if ( y .lt. one ) then

          z = y
          y = y + one
!
!  1.0 < argument < 12.0.
!  Reduce argument if necessary.
!
        else

          n = int ( y ) - 1
          y = y - dble ( n )
          z = y - one

        end if
!
!  Evaluate approximation for 1.0 < argument < 2.0.
!
        xnum = zero
        xden = one
        do i = 1, 8
          xnum = ( xnum + p(i) ) * z
          xden = xden * z + q(i)
        end do

        res = xnum / xden + one
!
!  Adjust result for case  0.0 < argument < 1.0.
!
        if ( y1 .lt. y ) then

          res = res / y1
!
!  Adjust result for case 2.0 < argument < 12.0.
!
        else if ( y .lt. y1 ) then

          do i = 1, n
            res = res * y
            y = y + one
          end do

        end if

      else
!
!  Evaluate for 12.0 <= argument.
!
        if ( y .le. xbig ) then

          ysq = y * y
          sum = c(7)
          do i = 1, 6
            sum = sum / ysq + c(i)
          end do
          sum = sum / y - y + sqrtpi
          sum = sum + ( y - half ) * log ( y )
          res = exp ( sum )

        else

          res = xinf
          r8_gamma = res
          return

        end if

      end if
!
!  Final adjustments and return.
!
      if ( parity ) then
        res = - res
      end if

      if ( fact .ne. one ) then
        res = fact / res
      end if

      r8_gamma = res

      return
      end function r8_gamma
!-----------------------------------------------
        SUBROUTINE VPV(XM,YM,NN)
!   Y=X+Y
        IMPLICIT NONE
      INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
      INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single 
        INTEGER :: NN
        REAL(wp) :: XM(NN),YM(NN)
        INTEGER :: I,J,IJ

        DO J = 1, NN
          YM(J) = XM(J) + YM(J)
        END DO  

      RETURN
      END SUBROUTINE VPV
!---------------------------------------
      SUBROUTINE MATRIX(SMAT,XM,NX,NB,AX)
! calculate matrix for spherical Bessel functions
      IMPLICIT NONE
      INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
      INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single 
      REAL(wp) :: SMAT(NX*NB),XM(NX)
      REAL(wp) :: AX
      INTEGER :: NX,NB 

      INTEGER :: I1,I2,J,ERR  
      REAL(wp) :: XX,BES(NB)

        DO 6 I1=1,NX
           XX=XM(I1)
           CALL SpherBessJ(AX*XX,NB,BES,ERR)
        DO 5 I2=1,NB 
           J=(I2-1)*NX + I1
           SMAT(J)=BES(I2)
5       CONTINUE   
6       CONTINUE   

      RETURN
      END SUBROUTINE MATRIX
!--------------------------------------------------------------
        SUBROUTINE VMV(XM,YM,NN)
!   Y=X-Y
        IMPLICIT NONE
      INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
      INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single 
        INTEGER :: NN
        REAL(4) :: XM(NN),YM(NN)
        INTEGER :: I,J,IJ

        DO J = 1, NN
          YM(J) = XM(J) - YM(J)
        END DO  

      RETURN
      END SUBROUTINE VMV
!---------------------------------------
        SUBROUTINE MXV(AM,MM,NN,XM,YM)
!   A*X=Y,  AM(MM,NN)
        IMPLICIT NONE
      INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
      INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single 
        INTEGER :: MM,NN
        REAL(4) :: AM(MM*NN),XM(NN),YM(MM)
        INTEGER :: I,J,IJ

      DO I = 1, MM
        YM(I) = 0.0
        DO J = 1, NN
           IJ=(J-1)*MM+I
          YM(I) = YM(I) + AM(IJ) * XM(J)
        END DO
      END DO

      RETURN
      END SUBROUTINE MXV
!---------------------------------------------------
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

          SIGNY=SIGN(1.0_wp,Y)
          SIGNZ=SIGN(1.0_wp,Z)

          IF(X.NE.0.0_wp)GO TO 30
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

          SIGNY=SIGN(1.0_wp,Y)
          SIGNZ=SIGN(1.0_wp,Z)

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
      SUBROUTINE GNUFORM (NX,NY,GC,KAN)
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

      SUBROUTINE INPUT_MODELS(NDAT1,NDAT2,NDAT3)

        USE MODEL_GAUSS
        USE MODEL_PARALL
        USE MODEL_CYLINDER
        USE CHANALS

        IMPLICIT NONE

        INTEGER(4) I,NDAT1,NDAT2,NDAT3,IERR,ERR_ALLOC 

        ALLOCATE (GA(NDAT1),GAH(NDAT1),GB(NDAT1), &
                  GBH(NDAT1),GC(NDAT1),GCH(NDAT1),&
                  GAMPL(NDAT1),GTET1(NDAT1),GPHI1(NDAT1),&
                  PA(NDAT2),PAH(NDAT2),PB(NDAT2),PBH(NDAT2),&
                  PC(NDAT2),PCH(NDAT2),PAMPL(NDAT2),PTET1(NDAT2),PPHI1(NDAT2), &
                  CA(NDAT3),CAH(NDAT3),CB(NDAT3),CBH(NDAT3),&
                  CC(NDAT3),CCH(NDAT3),CAMPL(NDAT3),CTET1(NDAT3),CPHI1(NDAT3), &
                  STAT=ERR_ALLOC)

        IF(ERR_ALLOC .NE. 0) THEN
           IERR=3
           WRITE(*,*) 'INPUT_MODELS: IERR=',IERR
           GOTO 888
        ENDIF 

        READ(KAN(12),*)
        DO 30 I=1,NDAT1
           READ(KAN(12),110) GA(I),GAH(I),GB(I),GBH(I),GC(I),GCH(I),GAMPL(I),GTET1(I),GPHI1(I)
30      CONTINUE

        READ(KAN(14),*)
        DO 32 I=1,NDAT2
           READ(KAN(14),110) PA(I),PAH(I),PB(I),PBH(I),PC(I),PCH(I),PAMPL(I),PTET1(I),PPHI1(I)
32      CONTINUE

        READ(KAN(15),*)
        DO 34 I=1,NDAT3
           READ(KAN(15),110) CA(I),CAH(I),CB(I),CBH(I),CC(I),CCH(I),CAMPL(I),CTET1(I),CPHI1(I)
34      CONTINUE

110   FORMAT(1X,9F10.4)   
888   CONTINUE
      RETURN
      END SUBROUTINE INPUT_MODELS

!-----------------------------------------------

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
        PHI1=  2.0*PI  
!        RE1 =  RR
!        RE2 =  2.*RR 
!        NU=NX
!        NV=NY
!        NW=NZ

        PP0   = 0.01   !-RR 
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


           ALPHA0=0.0
           ALPHA1=2.0*PI
           HALPHA=(ALPHA1-ALPHA0)/FLOAT(NALPHA)
        DO 32 I=1,NALPHA
32         ALPHAM(I)=ALPHA0+(I-1)*HALPHA

           BETA0=0.0
           BETA1=PI
           HBETA=(BETA1-BETA0)/FLOAT(NBETA)
        DO 34 I=1,NBETA
34         BETAM(I)=BETA0+(I-1)*HBETA

           GAMMA0=0.0
           GAMMA1=2.0*PI
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

      RETURN
      END SUBROUTINE DATA
!----------------------------------------------------------
 
      SUBROUTINE NUSET_3D(NX,NY,NZ,HX,HY,HZ,XNU,YNU,ZNU,RNU, &
                          HNUX,HNUY,HNUZ)
!CL...  grids in a fourier space

        USE FREQUENCY
        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) ! working precision =single 

        INTEGER :: NX,NY,NZ
        REAL(wp) :: XNU(NX),YNU(NY),ZNU(NZ),RNU(NX)
        REAL(wp) :: HX,HY,HZ
        REAL(wp) :: HNUX,HNUY,HNUZ
! local
        INTEGER :: I,IERR
        CHARACTER*8 ITEX
        REAL(wp) :: SNAIX,SNAIY,SNAIZ

        ITEX='NUSET_3D'
        IERR=0

        IF (HX .NE. 0. .AND. HY .NE. 0. .AND. HZ .NE. 0.) GO TO 4
        IERR=1
        CALL ERRPRN(ITEX,IERR)
        GOTO 777
4       IF (NX .NE. 0  .AND. NY .NE. 0  .AND. NZ .NE. 0 ) GO TO 6
        IERR=2
        CALL ERRPRN(ITEX,IERR)
        GOTO 777

   6  CONTINUE

      HNUX=1./HX/FLOAT(NX)
      HNUY=1./HY/FLOAT(NY)
      HNUZ=1./HZ/FLOAT(NZ)
      SNAIX=0.5/HX
      SNAIY=0.5/HY
      SNAIZ=0.5/HZ

      FREQ(1) = SNAIX        
      FREQ(2) = SNAIY        
      FREQ(3) = SNAIZ  
      FREQ(4) = HNUX
      FREQ(5) = HNUY
      FREQ(6) = HNUZ

      write(*,*) 'NUSET_3D: Nyquist frequency=', SNAIX,HX

      DO 8 I=1,NX
         XNU(I)=-SNAIX+(I-1)*HNUX
8    CONTINUE
 
      DO 10 I=1,NY
         YNU(I)=-SNAIY+(I-1)*HNUY
10    CONTINUE

      DO 20 I=1,NZ
         ZNU(I)=-SNAIZ+(I-1)*HNUZ
20    CONTINUE

!      DO 22 I=1,NX/2
!22    RNU(I)=XNU(NX/2+1+I)

      DO 22 I=1,NX
22    RNU(I)=I*HNUX

      write(*,*) 'NUSET_3D: FREQUENCY'
      write(*,*) (XNU(I),I=1,NX) 
!      write(*,*) 'NUSET_3D: FREQUENCY: RNU'
!      write(*,*) (RNU(I),I=1,NX/2) 
!      write(*,*) (RNU(I),I=1,NX) 

777   CONTINUE
      RETURN
      END SUBROUTINE NUSET_3D
!----------------------------------------------------

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

      FUNCTION ALINEAR_3D(AM,XM,YM,ZM,N1,N2,N3,  &
                          DELX,DELY,DELZ,XX1,YY1,ZZ1)

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=singl 

        REAL(wp) :: ALINEAR_3D
        REAL(wp) :: AM(N1*N2*N3)
        REAL(wp) :: FILTR,FILTR3_1_N
        REAL(wp) :: XM(N1),YM(N2),ZM(N3)  
        REAL(wp) :: DELX,DELY,DELZ
        REAL(wp) :: XX1,YY1,ZZ1
        INTEGER*4 N1,N2,N3

        REAL(wp) :: XXX,YYY,ZZZ
        REAL(wp) :: FIL1,FIL2,FIL3
        REAL(wp) :: FILL1,FILL2,FILL3
        REAL(wp) :: SS1,SS2

        INTEGER*4 I1,I2,I3,I321,J1,J2,J3,J4

!        IF(DELX .LT. 1.E-8 ) THEN 
!        write(*,*) 'ALINEAR_3D',DELX,DELY,DELZ
!        ENDIF

        CALL INDNW3(XM,YM,ZM,N1,N2,N3,  &
             DELX,DELY,DELZ,XX1,YY1,ZZ1,&
             J1,J2,J3,J4)

! for first and end points of the gread
!--------------------------------------------------
        SS1=0.
        DO 2 I3=1,N3,N3-1
           ZZZ=ZM(I3)
           FILL3=FILTR3_1_N(N3,J3,ZZZ,ZZ1,DELZ)

        DO 2 I2=1,N2,N2-1
           YYY=YM(I2)
           FILL2=FILTR3_1_N(N2,J2,YYY,YY1,DELY)

        DO 2 I1=1,N1,N1-1
           XXX=XM(I1)
           FILL1=FILTR3_1_N(N1,J1,XXX,XX1,DELX)

           I321=(I3-1)*N2*N1+(I2-1)*N1+I1

           SS1=SS1+AM(I321)*FILL1*FILL2*FILL3
2       CONTINUE
!----------------------------------------------------

        IF(J1 .NE. 1 .AND. J2 .NE. 1 .AND. J3 .NE. 1) THEN
        IF(J1 .NE. N1 .AND. J2 .NE. N2 .AND. J3 .NE. N3) THEN
        SS2=0.
        DO 4 I3=J3-1,J3+1
           ZZZ=ZM(I3)
           FIL3=FILTR(ZZ1-ZZZ,DELZ)

        DO 4 I2=J2-1,J2+1
           YYY=YM(I2)
           FIL2=FILTR(YY1-YYY,DELY)

        DO 4 I1=J1-1,J1+1
           XXX=XM(I1)
           FIL1=FILTR(XX1-XXX,DELX)

           I321=(I3-1)*N2*N1+(I2-1)*N1+I1

           SS2=SS2+AM(I321)*FIL1*FIL2*FIL3

4       CONTINUE
           ALINEAR_3D=SS2 + SS1
        ENDIF   
        ENDIF 
      RETURN
      END FUNCTION ALINEAR_3D
!------------------------------------------


      FUNCTION ALINEAR_1_3D(AM,XM,YM,ZM,N1,N2,N3,  &
                          DELX,DELY,DELZ,XX1,YY1,ZZ1)

        IMPLICIT NONE

        REAL(4) ALINEAR_1_3D,FILTR_0
        REAL(4) AM(N1*N2*N3),XM(N1),YM(N2),ZM(N3)  
        REAL(4) DELX,DELY,DELZ,XX1,YY1,ZZ1
        INTEGER*4 N1,N2,N3

        REAL(4) XXX,YYY,ZZZ,FIL1,FIL2,FIL3
        REAL(4) SS1,SS2
        REAL(4) FILL1,FILL2,FILL3

        INTEGER*4 I1,I2,I3,I321,J1,J2,J3,J4

        CALL INDNW3(XM,YM,ZM,N1,N2,N3,  &
             DELX,DELY,DELZ,XX1,YY1,ZZ1,&
             J1,J2,J3,J4)

        IF(J1 .NE. N1 .AND. J2 .NE. N2 .AND. J3 .NE. N3) THEN

        DO 4 I3=J3,J3+1
           ZZZ=ZM(I3)

           FIL3=FILTR_0(ZZ1-ZZZ,DELZ)

        DO 4 I2=J2,J2+1
           YYY=YM(I2)

           FIL2=FILTR_0(YY1-YYY,DELY)

        DO 4 I1=J1,J1+1
           XXX=XM(I1)

           FIL1=FILTR_0(XX1-XXX,DELX)

           I321=(I3-1)*N2*N1+(I2-1)*N1+I1

           SS2=SS2+AM(I321)*FIL1*FIL2*FIL3

4       CONTINUE
           ALINEAR_1_3D=SS2 
        ENDIF   
      RETURN
      END FUNCTION ALINEAR_1_3D
!--------------------------------------------
      FUNCTION FILTR_0(XX,DELX)
        IMPLICIT NONE 

        REAL(4) FILTR_0
        REAL(4) XX,DELX
        REAL(4) HH,AA,BB

!        HH(AA,BB)=1.0-ABS(AA)/BB

        IF(ABS(XX) .LE. DELX) THEN
           FILTR_0=1.
        ELSE
           FILTR_0=0.
        ENDIF
      RETURN
      END FUNCTION FILTR_0
!---------------------------------------

      FUNCTION FILTR(XX,DELX)

        IMPLICIT NONE 
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=singl 

        REAL(sp) FILTR
        REAL(sp) XX,DELX
        REAL(sp) AA,BB
        REAL(sp) HH

        HH(AA,BB)=1.0-ABS(AA)/BB

        IF(ABS(XX) .LE. DELX) THEN
!        IF(ABS(XX) .LT. DELX) THEN
           FILTR=HH(XX,DELX)
        ELSE
           FILTR=0.
        ENDIF
      RETURN
      END FUNCTION FILTR
!------------------------------------------------

      FUNCTION FILTR3_1_N(N,J,XXX,XX1,DELX)

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=singl 

        REAL(sp) FILTR3_1_N
        REAL(sp) XXX,XX1,DELX
        INTEGER*4 N,J

        IF(J .EQ. 1) FILTR3_1_N= (XXX-XX1)/DELX
        IF(J .EQ. N) FILTR3_1_N= (XX1-XXX)/DELX

      RETURN
      END FUNCTION FILTR3_1_N

!--------------------------------------------------------

      SUBROUTINE IFITET(JJ,N1,N2,I1,I2)
! calculation I2, I1-- tet and fi indices by number of projection JJ
!  JJ=(I2-1)*N1+I1

        INTEGER :: JJ,N1,N2,I1,I2
        INTEGER :: MODUL,IC

        MODUL=MOD(JJ,N1)
        IC=1
        IF(MODUL.EQ.0) IC=0
        I2=INT(JJ/N1)+IC
        I1=JJ-(I2-1)*N1
      RETURN
      END SUBROUTINE IFITET
!-------------------------------------------------

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

!      FUNCTION SIGNUM(SA)
      
!        IMPLICIT NONE 
!        REAL SA
!        INTEGER*4 NN
!        REAL(4) SIGNUM

!!        NN=ABS(AMOD(SA,2.0))
!!        SIGNUM=1.
!!        IF(NN .EQ.1) SIGNUM=-1.

!        IF(SA .NE. 0.) THEN 
!           SIGNUM=ABS(SA)/SA
!        ELSE 
!           SIGNUM=0.
!        ENDIF
!      RETURN
!      END FUNCTION SIGNUM
!--------------------------------------------------------



      FUNCTION SpherBessY(NN,XX)

      IMPLICIT NONE

      INTEGER*4 NN
      REAL(8) XX
      REAL(8) SpherBessY

      REAL(8) X,FACT
      REAL(8) FY0,FY1
      REAL(8) SY,SY0,SY1
      REAL(8) HX
      INTEGER*4 I

!      FJ0(X)=SIN(X)/X
!      FJ1(X)=SIN(X)/X**2 - COS(X)/X
!!      FJ2(X)=(3/X**3-1/X)*SIN(X) - 3/X**2*COS(X)
      FY0(X)=-COS(X)/X
      FY1(X)=-COS(X)/X**2 - SIN(X)/X

      HX=0.001

      SELECT CASE(NN)
         CASE(0)
            IF(XX .NE. 0.) THEN
 !           SJ=FJ0(XX)
            SY=FY0(XX)
         ELSE
 !           SJ=1.
            SY=FACT(NN,HX)
         ENDIF   

         CASE(1)
            IF(XX .NE. 0. ) THEN
 !           SJ=FJ1(XX)
            SY=FY1(XX)
         ELSE
 !           SJ=0.
            SY=FACT(NN,HX)
         ENDIF  
 
         CASE(2:)
            IF(XX .NE. 0. ) THEN
 !              SJ0=FJ0(XX)
 !              SJ1=FJ1(XX)

               SY0=FY0(XX)
               SY1=FY1(XX)
            ELSE
!               SJ0=1.
!               SJ1=0.
               
               SY0=FACT(0,HX)
               SY1=FACT(1,HX)
            ENDIF
            
          DO 4 I=2,NN  
            IF(XX .NE. 0. ) THEN
!               SJ=(2*(I-1)+1)/XX*SJ1-SJ0
               SY=(2*(I-1)+1)/XX*SY1-SY0
            ELSE
!               SJ=0.
               SY=FACT(I,HX)
            ENDIF
!               SJ0=SJ1
!               SJ1=SJ

               SY0=SY1
               SY1=SY
4         CONTINUE
      ENDSELECT      

      SpherBessY=SY

      RETURN
      END FUNCTION SpherBessY
!-------------------------------------------------------

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
      SUBROUTINE RayTransform1_M_LM(L,M,KK,PPM,NPP,TETA,PHI,ALPHA,BETA,CLMM)

!  Ray transform of vectors M_LM

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp


        INTEGER :: L,M,NPP
        REAL(wp) :: TETA,PHI,ALPHA,BETA 
        REAL(wp) :: PPM(NPP) 
        REAL(wp) :: KK
        COMPLEX(wp) :: CLMM
! local
        INTEGER :: I,J
        COMPLEX(wp) :: MLMM1(3),MLMM2(3),MLMM(3),SS1,SS2
        REAL(wp) :: RR,UNIT(3)

        Unit(1)  = SIN(BETA)*COS(ALPHA)
        Unit(2)  = SIN(BETA)*SIN(ALPHA)
        Unit(3)  = COS(BETA) 

        SS1=(0.0,0.0)
        DO 4 I=1,NPP
           RR=PPM(I)
           SS2=(0.0,0.0)
           DO 2 J=1,3
              CALL M_LM (L,M,KK,RR,TETA,PHI,MLMM1)
              CALL M_LM (L,M,KK,RR,PI-TETA,PI+PHI,MLMM2)
              SS2=SS2+(MLMM1(J)+MLMM2(J))*Unit(J)
2          CONTINUE
              SS1=SS1+SS2
4       CONTINUE
        CLMM=SS1      
      RETURN
      END SUBROUTINE RayTransform1_M_LM
!-------------------------------------------------

      function djmn (j, m, n, beta_rad)
        implicit none
! -------------------------------------------------------------------
! input: angular momentum quantum numbers j, m, n (all real)
!        beta_rad = Euler angle beta (radian)
! -------------------------------------------------------------------
        integer, parameter :: wp = kind(1.0d0)  ! working precision = double (portable)
!--------------------------------------------------------------------
!    formal arguments
!--------------------------------------------------------------------
        real(wp), intent(in) :: j, m, n, beta_rad
!--------------------------------------------------------------------
!    local variables
!--------------------------------------------------------------------
        integer :: itmin1, itmin2, itmin, itmax1, itmax2, itmax, ij1, ij2, &
             it, iphase, ia, ib, ic
        real(wp) :: djmn, cosb2, sinb2, sqrt_fac, sumt, denom, term
!--------------------------------------------------------------------
!   external functions
!--------------------------------------------------------------------
        real(wp), external :: fac10           ! computes factorial(n)/10**n
!--------------------------------------------------------------------
!  program starts here
!--------------------------------------------------------------------
        cosb2 = cos(beta_rad/2.0_wp)
        sinb2 = sin(beta_rad/2.0_wp)
!--------------------------------------------------------------------
! determine lower and upper limits for summation index it; these
! are derived from the requirement that all factorials n! in the
! denominator are restricted to values with n >=0.
!--------------------------------------------------------------------
        itmin1 = 0
        itmin2 = nint(m-n)
        itmin = max(itmin1,itmin2)
        itmax1 = nint(j+m)
        itmax2 = nint(j-n)
        itmax = min(itmax1,itmax2)
!  write (6,'(10X,A,2I6)') ' itmin, itmax = ', itmin, itmax
        ij1 = nint(j-m)
        ij2 = nint(j+n)
        sqrt_fac = sqrt( fac10(itmax1) * fac10(ij1) * fac10(ij2) * fac10(itmax2) )
!
        sumt = 0.0_wp
        do it = itmin, itmax
           iphase = (-1)**it
           ia = itmax1 - it
           ib = itmax2 - it
           ic = it + nint(n-m)
!     write (6,'(10X,A,5I6)') ' it, iphase, ia, ib, ic  = ', it, iphase, ia, ib, ic
           denom = fac10(ia) * fac10(ib) * fac10(it) * fac10(ic)
           term = iphase * cosb2**(ia+ib) * sinb2**(it+ic) / denom
           sumt = sumt + term
        end do
        djmn = sqrt_fac * sumt
!
        return
      end function djmn
!-------------------------------------------------------------

      SUBROUTINE FFT3D(FRM,FIM,M1,M2,M3,N1,N2,N3,H1,H2,H3,ISIGNUM,IERR)
!CL..... three dimensional fourier transform
!CL..... FRM ,FIM   - input massives
!CL......QR, QI   - worker massives
!CL......N1=2**M1+1 - the number points on first coordinate
!CL......N2=2**M2+1 - the number points on second coordinate
!CL......N3=2**M3+1 - the number points on third coordinate
!CL......H1 , H2,H3 - steps of the 
!CL......SIGNUM=-1  - direct transform
!CL......SIGNUM=+1  - invers transform

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp

      REAL(wp) :: FRM(N1*N2*N3),FIM(N1*N2*N3)
      REAL(wp) :: H1,H2,H3
      INTEGER :: M1,M2,M3,N1,N2,N3,ISIGNUM,IERR

!      REAL(4) QR(N1),QI(N1)
!!    QR ad QI -1D massive; QR(MAX(N1,N2,N3)),QI(MAX(N1,N2,N3))
      INTEGER :: MM(3)
      INTEGER :: I1,I2,I3,I321,N11,N22,N33 

      REAL(wp) :: HN1,HN2,HN3,DELTA

      CHARACTER*8 ITEX
      ITEX='FFT3D'


!      IERR=0
!      DELTA=1./FLOAT(N1*N2*N3)
!      IF(ISIGNUM .EQ. -1) DELTA=1.

      N11=2**M1
      N22=2**M2
      N33=2**M3

      HN1=H1*N11
      HN2=H2*N22
      HN3=H3*N33

      IF(ISIGNUM .EQ.1)  DELTA=H1*H2*H3
      IF(ISIGNUM .EQ.-1) DELTA=1./(HN1*HN2*HN3)

!      IF(ISIGNUM .EQ.-1)  DELTA=1./(HN1*HN2*HN3)
!      IF(ISIGNUM .EQ.1)   DELTA=H1*H2*H3



      MM(1)=M1
      MM(2)=M2
      MM(3)=M3

!      IF(N1.EQ.2**M1) GO TO 4
!      IERR=1
!      CALL ERRPRN(ITEX,IERR)
!      RETURN
!   4  CONTINUE

!      IF(N2.EQ.2**M2) GO TO 6
!      IERR=2
!      CALL ERRPRN(ITEX,IERR)
!      RETURN
!   6  CONTINUE

!      IF(N3.EQ.2**M3) GO TO 8
!      IERR=3
!      CALL ERRPRN(ITEX,IERR)
!      RETURN
!   8  CONTINUE

      CALL SHFT3 (N1,N2,N3,FRM,FIM,IERR)


      CALL NLOGN33(FRM,FIM,MM,N1,N2,N3,ISIGNUM)


      CALL SHFT3 (N1,N2,N3,FRM,FIM,IERR)

      DO 12 I3 = 1, N3
      DO 12 I2 = 1, N2
      DO 12 I1 = 1, N1
      I321 = (I3-1)*N2*N1 + (I2-1)*N1 + I1
      FRM(I321) = FRM(I321)*DELTA
      FIM(I321) = FIM(I321)*DELTA
12    CONTINUE
      RETURN
      END SUBROUTINE FFT3D
!-----------------------------------------------------

      SUBROUTINE NLOGN33(WR,WI,MM,N1,N2,N3,ISIGNUM)
!CL..... 3-D fourier transform based on NLOGN1_1 procedure

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp


        REAL(wp) :: WR(N1*N2*N3),WI(N1*N2*N3)
        INTEGER :: N1,N2,N3,MM(3),ISIGNUM
        INTEGER :: I1,I2,I3,M1,M2,M3,J 

        REAL(wp), ALLOCATABLE :: QR(:),QI(:)

        ALLOCATE(QR(N1),QI(N1))

        M1=MM(1)
        M2=MM(2)
        M3=MM(3)

      DO 12 I3=1,N3
      DO 12 I2=1,N2
      DO 10 I1=1,N1
      J=(I3-1)*N1*N2+(I2-1)*N1+I1   
      QR(I1)=WR(J)
  10  QI(I1)=WI(J)
!      CALL NLOGN1(M1,QR,QI,ISIGNUM)
      CALL NLOGN1_1(M1,QR,QI,ISIGNUM)
      DO 12 I1=1,N1
      J=(I3-1)*N1*N2+(I2-1)*N1+I1  
      WR(J)=QR(I1)
      WI(J)=QI(I1)
  12  CONTINUE
!-------------------------------------
      DO 16 I3=1,N3
      DO 16 I1=1,N1
      DO 14 I2=1,N2
      J=(I3-1)*N1*N2+(I2-1)*N1+I1   
      QR(I2)=WR(J)
  14  QI(I2)=WI(J)
!      CALL NLOGN1(M2,QR,QI,ISIGNUM)
      CALL NLOGN1_1(M2,QR,QI,ISIGNUM)
      DO 16 I2=1,N2
      J=(I3-1)*N1*N2+(I2-1)*N1+I1 
      WR(J)=QR(I2)
      WI(J)=QI(I2)
  16  CONTINUE
!----------------------------------
      DO 20 I1=1,N1
      DO 20 I2=1,N2
      DO 18 I3=1,N3
      J=(I3-1)*N1*N2+(I2-1)*N1+I1 
      QR(I3)=WR(J)
  18  QI(I3)=WI(J)
!      CALL NLOGN1(M3,QR,QI,ISIGNUM)
      CALL NLOGN1_1(M3,QR,QI,ISIGNUM)
      DO 20 I3=1,N3
      J=(I3-1)*N1*N2+(I2-1)*N1+I1 
      WR(J)=QR(I3)
      WI(J)=QI(I3)
  20  CONTINUE

      DEALLOCATE(QR,QI)
      RETURN
      END SUBROUTINE NLOGN33
!--------------------------------------------

      SUBROUTINE NLOGN1_1(N,XR,XI,ISIGNUM)

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp

        REAL(wp) :: XR(2**N+1),XI(2**N+1)
        INTEGER :: N,M(25)

        INTEGER :: I,II,J,JH,L,K,LX,ISIGNUM
        INTEGER :: NBLOCK,LBLOCK,LBHALF,IBLOCK
        INTEGER :: ISTART
        REAL(wp) :: FLX,ALFA,FK,V,COSV,SINV,QR,QI
        REAL(wp) :: HOLDR,HOLDI

        LX=2**N
        FLX=LX

        ALFA=ISIGNUM*2.*PI/FLX
        DO 1 I=1,N
1          M(I)=2**(N-I)
        DO 4 L=1,N
           NBLOCK=2**(L-1)
           LBLOCK=LX/NBLOCK
           LBHALF=LBLOCK/2
           K=0
        DO 4 IBLOCK=1,NBLOCK
           FK=K
           V=ALFA*FK
           COSV=COS(V)
           SINV=SIN(V)
           ISTART=LBLOCK*(IBLOCK-1)
        DO 2 I=1,LBHALF
           J=ISTART+I
           JH=J+LBHALF
           QR=XR(JH)*COSV-XI(JH)*SINV
           QI=XR(JH)*SINV+XI(JH)*COSV
           XR(JH)=XR(J)-QR
           XI(JH)=XI(J)-QI
           XR(J) =XR(J)+QR
           XI(J) =XI(J)+QI
2       CONTINUE
        DO 3 I=2,N
           II=I
           IF(K.LT.M(I)) GO TO 4
3          K=K-M(I)
4          K=K+M(II)
           K=0
        DO 7 J=1,LX
           IF(K.LT.J) GO TO 5
           HOLDR=XR(J)
           HOLDI=XI(J)
           XR(J)=XR(K+1)
           XI(J)=XI(K+1)
           XR(K+1)=HOLDR
           XI(K+1)=HOLDI
5       DO 6 I=1,N
              II=I
              IF(K.LT.M(I)) GO TO 7
6             K=K-M(I)
7             K=K+M(II)
      RETURN
      END SUBROUTINE NLOGN1_1
!-------------------------------------------------------------

      SUBROUTINE SHFT3(N1,N2,N3,SHIM,SHRE,IERR)

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp

        REAL(wp) :: SHIM(N1*N2*N3),SHRE(N1*N2*N3)
        INTEGER :: N1,N2,N3,IERR

        REAL(wp), ALLOCATABLE :: AR(:),AI(:)
        !      REAL(4) AR(N1),AI(N1)
        
        INTEGER :: I1,I2,I3,I321

        ALLOCATE(AR(N1),AI(N1))


        DO 30 I3=1,N3     
        DO 30 I2=1,N2
        DO 10 I1=1,N1  
           I321=(I3-1)*N1*N2+(I2-1)*N1+I1
           AR(I1)=SHRE(I321)
           AI(I1)=SHIM(I321)
10      CONTINUE

           CALL SHFT1(N1,AR,AI,IERR)

        DO 20 I1=1,N1  
            I321=(I3-1)*N1*N2+(I2-1)*N1+I1
            SHRE(I321)=AR(I1)
            SHIM(I321)=AI(I1)
20      CONTINUE
30      CONTINUE
!---------
        DO 50 I3=1,N3     
        DO 50 I1=1,N1
        DO 40 I2=1,N2  
           I321=(I3-1)*N1*N2+(I2-1)*N1+I1
           AR(I2)=SHRE(I321)
           AI(I2)=SHIM(I321)
40      CONTINUE

           CALL SHFT1(N2,AR,AI,IERR)

        DO 42 I2=1,N2  
           I321=(I3-1)*N1*N2+(I2-1)*N1+I1
           SHRE(I321)=AR(I2)
           SHIM(I321)=AI(I2)
42      CONTINUE
50      CONTINUE
!---------
        DO 80 I1=1,N1
        DO 80 I2=1,N2  
        DO 82 I3=1,N3     
           I321=(I3-1)*N1*N2+(I2-1)*N1+I1
           AR(I3)=SHRE(I321)
           AI(I3)=SHIM(I321)
82      CONTINUE

           CALL SHFT1(N3,AR,AI,IERR)

        DO 84 I3=1,N3  
           I321=(I3-1)*N1*N2+(I2-1)*N1+I1
           SHRE(I321)=AR(I3)
           SHIM(I321)=AI(I3)
84      CONTINUE
80      CONTINUE
!--------
      DEALLOCATE(AR,AI)
      RETURN
      END SUBROUTINE SHFT3
!-------------------------------------------------------------
      SUBROUTINE SHFT1(N,AR,AI,IERR)
!CL...  AR,AI - massive for transformation
!CL.... XR,XI - worker massiv

        IMPLICIT NONE
        INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
        INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
        REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp

        REAL(wp) :: AR(N),AI(N)
        INTEGER :: N,IERR

        REAL(wp), ALLOCATABLE :: XR(:),XI(:)
        INTEGER :: NI,NI2,NI21,I,II

        ALLOCATE(XR(N),XI(N))

        NI=N
        NI2=NI/2
        NI21=NI2+1

        DO 10 I=1,NI2
           II=I+NI2
           XI(I)=AI(II)
10         XR(I)=AR(II)
        DO 16 I=NI21,NI
           II=I-NI2
           XI(I)=AI(II)
16         XR(I)=AR(II)
        DO 18 I=1,NI
           AI(I)=XI(I)
18         AR(I)=XR(I)

      DEALLOCATE(XR,XI)
      RETURN
      END SUBROUTINE SHFT1
!--------------------------------------------------------
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
      FUNCTION PSI_FUN(L1,M1,RR,TETA,PHI)
        USE PARAMTR    ! QMAX  are described in PARAM module 
!                        QMAX must be >= L1-1
          IMPLICIT NONE
          INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
          INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
!          REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

          INTEGER :: L1,M1  
          REAL(wp) :: RR,TETA,PHI

! local 
          REAL(wp) :: BES(QMAX)
          COMPLEX(wp) :: YY
          INTEGER :: ERR
!functions
          COMPLEX(wp) :: ylm, PSI_FUN

             CALL SpherBessJ(RR,QMAX,BES,ERR) 

             YY=ylm(L1, M1, TETA, PHI)

             PSI_FUN= BES(L1+1) * YY   

        RETURN
        END FUNCTION PSI_FUN
!--------------------------------------------

        SUBROUTINE Fourier_PSI_FUN_EXACT_POL(L1,M1,KK,RNU,NUU,TETM,NTET,PHIM,NPHI, &
                                            XNU,YNU,ZNU,NX,NY,NZ,PSIRE,PSIIM) 
! Analitical Fourier transform of PSI- functions 
!          USE PARAMTR
          IMPLICIT NONE
          INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
          INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
          REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

          INTEGER :: L1,M1
          INTEGER :: NUU,NTET,NPHI
          INTEGER :: NX,NY,NZ
          REAL(wp) :: KK 
          REAL(wp) :: RNU(NUU),TETM(NTET),PHIM(NPHI)
          REAL(wp) :: XNU(NX),YNU(NY),ZNU(NZ)
          REAL(wp) :: PSIRE(NX*NY*NZ),PSIIM(NX*NY*NZ)
! local
          INTEGER :: I1,I2,I3,I21,I321,IERR
          REAL(wp) :: XX,YY,ZZ,RRR
          REAL(wp) :: HPHI,HTET
          REAL(wp) :: EI
          REAL(wp) :: TETA,PHI,SIGNZ,RPART
          COMPLEX(wp) :: YYC  
          REAL(wp) :: DR(NUU*NTET*NPHI),DI(NUU*NTET*NPHI)
          COMPLEX(wp),PARAMETER :: ci = (0.,1.)

! functions 
          COMPLEX(wp) :: ylm
          REAL(wp) :: DeltaF

          DO 10 I3=1,NUU
             RRR=RNU(I3)

          DO 10 I2=1,NTET
             TETA=TETM(I2)

          DO 10 I1=1,NPHI
             PHI=PHIM(I1)

             I21=(I2-1)*NPHI+I1
             I321=(I3-1)*NTET*NPHI+I21


             YYC=ylm(L1, M1, TETA, PHI)*ci**L1

             RPART= SQRT(PI/2.)/RRR**2 *DeltaF(RRR-KK)

!             write(*,*) 'DR=',DeltaF(RRR-KK)
             DR(I321)=REAL(YYC)*RPART
             DI(I321)=AIMAG(YYC)*RPART
10        CONTINUE

!             write(*,*) 'DR=',DR


          CALL POLAR_DECART_3D(DR,XNU,YNU,ZNU,NX,NY,NZ,RNU,TETM,PHIM, &
                                        NUU,NTET,NPHI,PSIRE)
          CALL POLAR_DECART_3D(DI,XNU,YNU,ZNU,NX,NY,NZ,RNU,TETM,PHIM, &
                                        NUU,NTET,NPHI,PSIIM)
        RETURN
        END SUBROUTINE Fourier_PSI_FUN_EXACT_POL
!------------------------------------------------------------------

        SUBROUTINE Fourier_PSI_FUN_EXACT_DEC(L1,M1,KK,XNU,YNU,ZNU,NX,NY,NZ,PSIRE,PSIIM) 
! Analitical Fourier transform of PSI- functions 
!          USE PARAMTR
          IMPLICIT NONE
          INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
          INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
          REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

          INTEGER :: L1,M1
          INTEGER :: NX,NY,NZ
!          INTEGER :: MX,MY,MZ
          REAL(wp) :: KK !,RR
          REAL(wp) :: XNU(NX),YNU(NY),ZNU(NZ)
          REAL(wp) :: PSIRE(NX*NY*NZ),PSIIM(NX*NY*NZ)
! local
          INTEGER :: I1,I2,I3,I21,I321,IERR
          REAL(wp) :: XX,YY,ZZ,RRR
          REAL(wp) :: SNAI
          REAL(wp) :: TETA,PHI,SIGNZ,RPART
          COMPLEX(wp) :: YYC  
          COMPLEX(wp),PARAMETER :: ci = (0.,1.)

! functions 
          COMPLEX(wp) :: ylm
          REAL(wp) :: DeltaF

          SNAI=XNU(NX)

          DO 10 I3=1,NZ
             ZZ=ZNU(I3)
             SIGNZ=SIGN(1.0,ZZ)

          DO 10 I2=1,NY
             YY=YNU(I2)

          DO 10 I1=1,NX
             XX=XNU(I1)

             RRR=SQRT(XX**2+YY**2+ZZ**2)
             I21=(I2-1)*NX+I1
             I321=(I3-1)*NY*NX+I21

             CALL POLAR(XX,YY,ZZ,TETA,PHI)

             IF(TETA .EQ. 0. .OR. TETA .EQ. PI) THEN
                TETA=PI*(1.0+SIGNZ/180.)
             ENDIF

             YYC=ylm(L1, M1, TETA, PHI)*ci**L1

             RPART= SQRT(PI/2.)/RRR**2*DeltaF(RRR-KK)

             PSIRE(I321)=REAL(YYC)*RPART
             PSIIM(I321)=AIMAG(YYC)*RPART

             IF(RRR .GE. SNAI) THEN
                PSIRE(I321)=0.0
                PSIIM(I321)=0.0
             ENDIF

10        CONTINUE

        RETURN
        END SUBROUTINE Fourier_PSI_FUN_EXACT_DEC
!-----------------------------------------------
        SUBROUTINE Fourier_PSI_FUN_NUMERIC(L1,M1,KK,RR,XM,YM,ZM,NX,NY,NZ, &
                                           MX,MY,MZ,PSIRE,PSIIM) 
! Numerical Fourier transform of PSI- functions 
          USE PARAMTR
          IMPLICIT NONE
          INTEGER, PARAMETER  :: wp = kind(1.0d0) !working precision=double 
          INTEGER, PARAMETER  :: sp = kind(1.0) !working precision=single 
          REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp 

          INTEGER :: L1,M1
          INTEGER :: NX,NY,NZ
          INTEGER :: MX,MY,MZ
          REAL(wp) :: KK,RR
          REAL(wp) :: XM(NX),YM(NY),ZM(NZ)
          REAL(wp) :: PSIRE(NX*NY*NZ),PSIIM(NX*NY*NZ)
! local
          INTEGER :: I1,I2,I3,I21,I321
          REAL(wp) :: XX,YY,ZZ,RR2
          REAL(wp) :: HX,HY,HZ
          REAL(wp) :: TETA,PHI,SIGNZ
          COMPLEX(wp) :: PSIC
          INTEGER :: ISIGNUM,IERR
! functions 
          COMPLEX(wp) :: PSI_FUN

          HX=PAR(18)
          HY=PAR(21)
          HZ=PAR(24)

          DO 10 I3=1,NZ
             ZZ=ZM(I3)
             SIGNZ=SIGN(1.0,ZZ)

          DO 10 I2=1,NY
             YY=YM(I2)

          DO 10 I1=1,NX
             XX=XM(I1)

             RR2=SQRT(XX**2+YY**2+ZZ**2)
             I21=(I2-1)*NX+I1
             I321=(I3-1)*NY*NX+I21

             CALL POLAR(XX,YY,ZZ,TETA,PHI)

             IF(TETA .EQ. 0. .OR. TETA .EQ. PI) THEN
                TETA=PI*(1.0+SIGNZ/180.)
             ENDIF

             PSIC= PSI_FUN(L1,M1,KK*RR2,TETA,PHI)

             PSIRE(I321)=REAL(PSIC)
             PSIIM(I321)=AIMAG(PSIC)

10        CONTINUE
          ISIGNUM=1   
          CALL FFT3D(PSIRE,PSIIM,MX,MY,MZ,NX,NY,NZ,HX,HY,HZ,ISIGNUM,IERR)

        RETURN
        END SUBROUTINE Fourier_PSI_FUN_NUMERIC
!-----------------------------------------------
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

        REAL(wp), ALLOCATABLE :: DR1(:),DR2(:)
        INTEGER :: I1,I2,I3,I21,I321

        ALLOCATE(DR1(NX*NY),DR2(NX*NY*NZ))

!  IND=-1  -- section XZ
!  IND= 0  -- section XY
!  IND= 1  -- section YZ

!        write(*,*) 'AAAAAAAAAAAAA',KK
         DO 3 I3=1,NZ
         DO 3 I2=1,NY
         DO 3 I1=1,NX
            I21=(I2-1)*NX+I1
            I321=(I3-1)*NX*NY+I21
            DR2(I321)=DR(I321,KK)
3        CONTINUE   

      IF(IND) 2,4,6

 2    CONTINUE
      DO 10 I3=1,NZ
      DO 10 I2=KSEC,KSEC
      DO 10 I1=1,NX
         I21=(I3-1)*NX+I1
         I321=(I3-1)*NY*NX+(I2-1)*NX+I1
         DR1(I21)=DR2(I321)
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
         DR1(I21)=DR2(I321)
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
         DR1(I21)=DR2(I321)
 14   CONTINUE
!      CALL GNUFORM_2D (NY,NZ,DR1,KAN)
      CALL GNUFORM_M (YM,ZM,NY,NZ,DR1,KAN)
 30   CONTINUE
!      CLOSE(KAN)
      DEALLOCATE(DR1,DR2)
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
          INTEGER :: KK,I,I1,I2,I3,I21,I321,EI,IERR
!functions
          REAL(wp) :: FUNL3,Burke_3D  

!          write(*,*) 'aaaaaaaaa'
!          write(*,*) NTEN


      DO 10 KK=1,NTEN*NTEN

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
