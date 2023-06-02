C**********************************************************************
C
C                       ****************
C                       * VERSION 2.00 *
C                       ****************
C
C--This program calculates the production cross section of neutral
C  Higgs boson pairs via gg -> H_1 H_2 and qqbar -> HA at hadron
C  colliders at NLO QCD according to the formulae presented in
C
C  T. Plehn, M. Spira and P.M. Zerwas, Nucl. Phys. B479 (1996) 46,
C                                              (E) B531 (1998) 655;
C  S. Dawson, S. Dittmaier and M. Spira, Phys. Rev. D58 (1998) 115012.
C
C  The QCD corrections are only included in the heavy top quark limit,
C  which should be reliable for the SM Higgs and the MSSM Higgs bosons
C  for small tan(beta).
C
C  The program allows to calculate the total production cross sections
C  for the neutral Higgs bosons of the SM and MSSM. The MSSM Higgs
C  sector is implemented in the approximate two-loop RGE approach of
C
C  M. Carena, J. Espinosa, M. Quiros and C.E.M. Wagner, Phys. Lett. B355
C  (1995) 209.
C
C  The parton densities are defined in the subroutine STRUC at the end
C  of the program and can be changed there. The default is the use of
C  the PDFLIB.
C
C--Description of the input file hpair.in with default values:
C  ===========================================================
C
C--PROCESS:  0 = GG --> hh   1 = GG --> Hh  2 = GG --> HH
C            3 = GG --> Ah   4 = GG --> AH  5 = GG --> AA
C            6 = QQ --> Ah   7 = QQ --> AH
C  PROCESS = 0
C
C--QCD corrections: Born - virt - gg - gq - qqbar
C                   0 = no     1 = yes  (should be left as it is)
C  QCDCORR  =  11111
C
C--Model: MODEL = 0 -> SM
C         MODEL = 1 -> MSSM (2-loop)
C         MODEL = 2 -> MSSM (1-loop)
C  MODEL    = 0
C
C--SILH: SILH = 1 -> SILH
C        SILH = 0 -> NON-LINEAR
C  SILH     = 0
C
C--NON-LINEAR (SILH = 0): PARAMETERS
C  C_3      =  1.0D0
C  C_T      =  1.00D0
C  C_TT     =  0.00D0
C  C_G      =  0.00D0
C  C_GG     =  0.00D0
C  C_B      =  0.D0
C
C--SILH (SILH = 1): PARAMETERS
C  BAR C_H  =  0.00D0
C  BAR C_U  =  0.00D0
C  BAR C_6  =  0.00D0
C  BAR C_G  =  0.00012D0
C
C--tan(beta) for MSSM
C  TANBET   = 3.D0
C
C--Parton densities from LHAPDF. Here: MSTW 2008 NLO (to be chosen via
C                                                     PDFPATH AND PDFNAME)
C  NGROUP   = 0
C  NSET     = 0
C  PDFPATH  =  Grids
C  PDFNAME  =  MSTW2008nlo68cl.LHgrid
C
C--LAMBDA = LAMBDA_N0 for alpha_s
C  N0       = 5
C  LAMBDA   = 0.255097D0
C
C--Order of alpha_s: LOOP = 1: LO     LOOP = 2: NLO
C
C  LOOP     = 2
C
C--Scale: mu_R = mu_F = scale * mu_0
C         scale = [BEG,END] with NUMBER points
C
C  BEG      = 1.D0
C  END      = 1.D0
C  NUMBER   = 1
C
C--Energy of hadron collider in GeV
C  ENERGY   = 14000.D0
C
C--Choice pp (0) or ppbar (1) collider
C  PP/PPBAR = 0
C
C--Heavy quark masses for running masses and alpha_s (top quark
C                                                     decoupled!!!)
C  MSB(1)   = 0.190D0
C  MC       = 1.42D0
C  MB       = 4.75D0
C  MT       = 173.2D0
C
C--Higgs mass loop: MA = [MABEG, MAEND] with MULTMA equidistant points
C                   in total (SM: MH = MA)
C  MABEG    = 125.D0
C  MAEND    = 1000.D0
C  MULTMA   = 2
C
C--VEGAS: VPOINT = points for Born term
C         VITMAX = number of iterations
C         VNPRN: 0  = no output of individual iterations
C                1  = full output of individual iterations
C                10 = table output of individual iterations
C  VPOINT   = 100000
C  VITMAX   = 5
C  VNPRN    = 10
C
C--CUT: EPSILON = cut at the endpoints of integrations:
C                 [0,1] -> [EPSILON,1-EPSILON]
C       Leave it as it is, unless you observe instabilities (check the
C       individual VEGAS iterations!). In case of instabilities increase
C       EPSILON until the stability plateau is reached (the result must
C       still be independent of EPSILON!). Instabilities may occur for
C       very large collider energies compared with the bottom/top mass.
C  EPSILON  = 1.D-8
C
C--Parameters:
C  Fermi constant, QED coupling, Z mass, W mass, Weinberg angle
C  GF       = 1.16639D-5
C  ALPHA    = 137.0359895D0
C  MZ       = 91.187D0
C  MW       = 80.33D0
C  SW2      = 0.2315D0
C
C--MSSM parameters: in GeV
C                   MST_L =  left-handed stop mass
C                   MST_R = right-handed stop mass
C                   A_U   = trilinear coupling A_t
C                   A_D   = trilinear coupling A_b
C                   MU_H  = Higgsino mass paramater mu
C  MST_L    = 1000.D0
C  MST_R    = 1000.D0
C  A_U      = 0.D0
C  A_D      = 0.D0
C  MU_H     = 0.D0
C
C--Z decay width:
C
C  GAMZ     = 2.490D0
C
C
C
C  The following papers have to be cited when using this program:
C
C  \bibitem{Plehn:1996wb}
C  T.~Plehn, M.~Spira and P.~M.~Zerwas,
C  %``Pair production of neutral Higgs particles in gluon-gluon
C  %collisions,''
C  Nucl.\ Phys.\ B {\bf 479} (1996) 46
C  [Nucl.\ Phys.\ B {\bf 531} (1998) 655]
C  [hep-ph/9603205].
C  %%CITATION = HEP-PH/9603205;%%
C
C  \bibitem{Dawson:1998py}
C  S.~Dawson, S.~Dittmaier and M.~Spira,
C  %``Neutral Higgs boson pair production at hadron colliders: QCD
C  %corrections,''
C  Phys.\ Rev.\ D {\bf 58} (1998) 115012
C  [hep-ph/9805244].
C  %%CITATION = HEP-PH/9805244;%%
C
C  \bibitem{Grober:2015cwa}
C  R.~Grober, M.~Muhlleitner, M.~Spira and J.~Streicher,
C  %``NLO QCD Corrections to Higgs Pair Production including Dimension-6
C  %Operators,''
C  arXiv:1504.06577 [hep-ph].
C  %%CITATION = ARXIV:1504.06577;%%
C
C  written by Michael Spira, e-mail: michael.spira@psi.ch
C
C  June 12, 2015
C=======================================================
C
      PROGRAM HPAIR
C--PROGRAM FOR NEUTRAL HIGGS PAIR PRODUCTION AT HADRON COLLIDERS
      PARAMETER(NIN=95, NOUT=96)
      PARAMETER(NV=6)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION M1,M2,MZ
      CHARACTER*100 PDFNAME,PATHNAME
      COMMON/SUSYPAR/SGF,ALPHA,AMZ,AMSQ
      COMMON/BREAK/AMUR,AU,AD,AMU
      COMMON/PARAM/S,T,U,M1,M2,MZ
      COMMON/WEINBERG/SW2
      COMMON/SCALE/ZSS,ZPT,ZQQ,QQ,ISUB
      COMMON/TOTPAR/W,GEVPB,GF,GAMZ,GAMA,GAML,GAMH
      COMMON/HMASS/AMA,AML,AMH,AMCH
      COMMON/SUSY/TGBET
      COMMON/QMASS/AMB,AMT
      COMMON/MASSES/XMS,XMC,XMB,XMT,XMZ,XMW
      COMMON/STRANGE/AMSB
      COMMON/PROC/IPROC,IAPP,ITOP,ISFAC
      COMMON/TREE/ITREE,ITRIA,IWRITE
      COMMON/TRILINEAR/FACTRIA,FACT,FACB,FACTT,FACG,FACGG,ISILH
      COMMON/ALSN/LOOP
      COMMON/PARTONIC/WHAT
      COMMON/CUTS/EPS,REPS
      COMMON/SCAL0/V,TAU,AMQ
      COMMON/INOUT/IIN,IOUT
      COMMON/PDFLIB/NGROUP,NSET,SCALFAC
      COMMON/PDFLIB0/ISETERR,IALS
      COMMON/ALS/XLAMBDA,AMC,AMBB,AMTT,N0
      COMMON/VEGPAR/ABSERR,IPOINT,ITERAT,NPRN
      COMMON/VEGOUT/INV
      COMMON/KNIEHL/IKNIEHL
      COMMON/RESULT/RES,ERR,DCHI2,DUM
      COMMON/WIDTHS/SMWDTH,HLWDTH,HHWDTH,AWDTH,HCWDTH
      COMMON/WDPARAM/WGF,WALPH,WMTAU,WMMUON,WMZ,WMW
      COMMON/VIRTCOF/VC1,VC2,VC3,IHA
      COMMON/COLLIDER/ICOLL
      EXTERNAL SIG,SIGRES,SIGHAT,SIGPART
      EXTERNAL QCDV,QCDGG,QCDGQ,QCDQQ
      EXTERNAL SIGBQQ,QSIGQQ,QSIGGQ

      INTEGER :: NUM_ARGS, IX
      CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE :: ARGS
      NUM_ARGS = COMMAND_ARGUMENT_COUNT()
      ALLOCATE(ARGS(NUM_ARGS))
      IF (NUM_ARGS.NE.2) THEN
         PRINT*, "Please provide exactly two arguments: ./run in out"
         STOP
      ENDIF
      DO IX = 1, NUM_ARGS
         CALL GET_COMMAND_ARGUMENT(IX, ARGS(IX))
      END DO
      PRINT*, "Running with: Input=", ARGS(1), "\n", "Output=", ARGS(2)

      OPEN(NIN,FILE=ARGS(1))
      OPEN(NOUT,FILE=ARGS(2))
C     OPEN(NV,FILE='hpair.veg')
      IIN=NIN
      IOUT=NOUT
      INV=NV

      PI = 4*DATAN(1.D0)

C--VERSION OF INPUT FILE: 1 = FULL WITH ALL FLAGS    0 = PUBLIC VERSION
      IVERSION = 1

C--READING INPUT FILE
      IF(IVERSION.EQ.0)THEN
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,100)IPROC
       IPARTON = 0
       IAPP = 0
       ITOP = 0
       IRES = 0
       READ(NIN,100)IQCD
       IQCDB  = INT(IQCD/10000)
       IQCDV  = INT((IQCD-10000*IQCDB)/1000)
       IQCDGG = INT((IQCD-10000*IQCDB-1000*IQCDV)/100)
       IQCDGQ = INT((IQCD-10000*IQCDB-1000*IQCDV-100*IQCDGG)/10)
       IQCDQQ = IQCD-10000*IQCDB-1000*IQCDV-100*IQCDGG-10*IQCDGQ
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,100)ISUSY
       ITREE = 0
       ITRIA = 0
       IWRITE = 0
       READ(NIN,*)
       READ(NIN,100)ISILH
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,101)FACTRIA
       READ(NIN,101)FACT
       READ(NIN,101)FACTT
       READ(NIN,101)FACG
       READ(NIN,101)FACGG
       READ(NIN,101)FACB
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,101)BARCH
       READ(NIN,101)BARCU
       READ(NIN,101)BARC6
       READ(NIN,101)BARCG
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,101)TGBET
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,100)NGROUP
       READ(NIN,100)NSET
       READ(NIN,102)PATHNAME
       READ(NIN,102)PDFNAME
       ISETERR = 0
       IALS = 0
       SCALFAC = 1
       ISUB = 1
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,101)ALSMZ
       ACC = 1.D-8
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,100)LOOP
       NDEC = 0
       ZSS0 = 1
       ZPT0 = 0
       ZQQ0 = 0
       QQ = 100.D0
       READ(NIN,101)SCBEG
       READ(NIN,101)SCEND
       READ(NIN,100)NSCALE
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,101)W
       READ(NIN,100)ICOLL
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,101)AMSB
       READ(NIN,101)AMC
       READ(NIN,101)AMB
       READ(NIN,101)AMT
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,101)AMABEG
       READ(NIN,101)AMAEND
       READ(NIN,100)NMA
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       ABSERR = 0.D0
       READ(NIN,100)IPOINT
       READ(NIN,100)ITERAT
       READ(NIN,100)NPRN
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,101)EPS
       REPS = 1.D-15
       IKNIEHL = 0
       NBER = 18
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       GEVPB = 389379660.D0
       READ(NIN,101)GF
       READ(NIN,101)ALPHA
       READ(NIN,101)MZ
       READ(NIN,101)AMW
       READ(NIN,101)SW2
       READ(NIN,101)AMSQ 
       READ(NIN,101)AMUR 
       READ(NIN,101)AU 
       READ(NIN,101)AD
       READ(NIN,101)AMU
       READ(NIN,*)
       ISFAC = 0
       READ(NIN,101)GAMZ 
c      write(6,*)'Gamma(Z) = ',gamz
       V = 1.D-1
       TAU = 1.D-2
       AM1 = 5.D0
       AM2 = 10.D0
       AMQ = 175.D0
       WHAT = 1000.D0
C      FACTRIA = 1
      ELSE
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,100)IPROC
C--Flags:
C         PARTON: 0 = HADRON CXN  1 = PARTON CXN
C         APPROX: 0 = MASSIVE     1 = LIMIT MQ --> INFINITY
C         TOP: 0 = TOP + BOTTOM   1 = ONLY TOP   2 = ONLY BOTTOM
C         RESONANCE: 0 = FULL     1 = RESONANCE CONTR. (PROC. 0, 5)
C  PARTON   =  0
C  APPROX   =  0
C  TOP      =  0
C  RESONANCE=  0
C
       READ(NIN,100)IPARTON
       READ(NIN,100)IAPP
       READ(NIN,100)ITOP
       READ(NIN,100)IRES

       READ(NIN,100)IQCD
       IQCDB  = INT(IQCD/10000)
       IQCDV  = INT((IQCD-10000*IQCDB)/1000)
       IQCDGG = INT((IQCD-10000*IQCDB-1000*IQCDV)/100)
       IQCDGQ = INT((IQCD-10000*IQCDB-1000*IQCDV-100*IQCDGG)/10)
       IQCDQQ = IQCD-10000*IQCDB-1000*IQCDV-100*IQCDGG-10*IQCDGQ
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,100)ISUSY
C         TR-LEVEL: 0 = full    1 = tree-level couplings and masses
C  TR-LEVEL =  0
C--Loops: TRIANGLE = 0: full
C         TRIANGLE = 1: no triangle loops
C         TRIANGLE = 2: no Z boson contribution
C         TRIANGLE = 3: only s-channel H exchange
C         TRIANGLE = 4: only s-channel h exchange
C
C--Printout: PR(SUSY) = 0: no printout of MSSM couplings
C            PR(SUSY) = 1: print out MSSM couplings
C  IWRITE   =  0
C
       READ(NIN,100)ITREE
       READ(NIN,100)ITRIA
       READ(NIN,100)IWRITE
       READ(NIN,*)
       READ(NIN,100)ISILH
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,101)FACTRIA
       READ(NIN,101)FACT
       READ(NIN,101)FACTT
       READ(NIN,101)FACG
       READ(NIN,101)FACGG
       READ(NIN,101)FACB
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,101)BARCH
       READ(NIN,101)BARCU
       READ(NIN,101)BARC6
       READ(NIN,101)BARCG
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,101)TGBET
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,100)NGROUP
       READ(NIN,100)NSET
       READ(NIN,102)PATHNAME
       READ(NIN,102)PDFNAME
C--Scale factor for parton densities (not partonic cxn!!!):
C                                             Q = SCALFAC * mu_F
C  SCALFAC  =  1.D0
C
       READ(NIN,101)SCALFAC
C--Substitution in integration: 0 = no    1 = yes (should be left
C                                                  as it is)
C  ISUB     =  1
C
       READ(NIN,100)ISUB
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,101)ALSMZ
C--ACC = accuracy of alpha_s evolution across thresholds
C  ACC      =  1.D-8
C
       READ(NIN,101)ACC
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,100)LOOP
       READ(NIN,100)NDEC
C--Scale for parton densities: mu_0 = zss*sqrt(s) + zpt*m_T + zq*q
C
C  ZSS0     = 1
C  ZPT0     = 0
C  ZQQ0     = 0
C  QQ       = 100.D0
C
       READ(NIN,101)ZSS0
       READ(NIN,101)ZPT0
       READ(NIN,101)QQ
       READ(NIN,101)ZQQ0
       READ(NIN,101)SCBEG
       READ(NIN,101)SCEND
       READ(NIN,100)NSCALE
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,101)W
       READ(NIN,100)ICOLL
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,101)AMSB
       READ(NIN,101)AMC
       READ(NIN,101)AMB
       READ(NIN,101)AMT
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,101)AMABEG
       READ(NIN,101)AMAEND
       READ(NIN,100)NMA
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,101)ABSERR
       READ(NIN,100)IPOINT
       READ(NIN,100)ITERAT
       READ(NIN,100)NPRN
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
C--Integration cut (should be left as it is)
C  EPSILON  = 1.D-8
C
       READ(NIN,101)EPS
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,101)REPS
       READ(NIN,100)IKNIEHL
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,100)NBER
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,101)GEVPB
       READ(NIN,101)GF
       READ(NIN,101)ALPHA
       READ(NIN,101)MZ
       READ(NIN,101)AMW
       READ(NIN,101)SW2
       READ(NIN,101)AMSQ 
       READ(NIN,101)AMUR 
       READ(NIN,101)AU 
       READ(NIN,101)AD
       READ(NIN,101)AMU
       READ(NIN,100)ISFAC
       READ(NIN,101)GAMZ 
c      write(6,*)'Gamma(Z) = ',gamz
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,101)V
       READ(NIN,101)TAU
       READ(NIN,101)AM1
       READ(NIN,101)AM2
       READ(NIN,101)AMQ
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,*)
       READ(NIN,101)WHAT
      ENDIF

      WRITE(6,*)
      IF(ISILH.NE.0)THEN
       ALP2 = DSQRT(2.D0)*GF*AMW**2/PI
       FACTRIA = 1 - 3*BARCH/2 + BARC6
       FACT    = 1 - BARCH/2 - BARCU
       FACTT   = -BARCH/2 - 3*BARCU/2
       FACG    = BARCG*4*PI/ALP2
       FACGG   = BARCG*4*PI/ALP2
       FACB    = 0
       WRITE(6,*)'SILH: OMITTING BOTTOM CONTRIBUTION...'
      WRITE(6,*)
      ENDIF

C--CHECK, IF RESONANCE CONTRIBUTION POSSIBLE
      IF(IRES.NE.0.AND.IPROC.NE.0.AND.IPROC.NE.5)THEN
        WRITE(6,*)'NO RESONANCE CONTRIBUTION!!'
        WRITE(6,*)
      ELSEIF(IRES.NE.0.AND.ISUSY.EQ.0)THEN
        WRITE(6,*)'NO RESONANCE CONTRIBUTION!!'
        WRITE(6,*)
      ELSE

C--SETTING FLAGS FOR MSSM VERSION
      IOLD=0
      IF(ISUSY.EQ.2)THEN
       ISUSY=1
       IOLD=1
      ENDIF

      ALPHA=1/ALPHA
C--INITIALIZING COMMON BLOCKS
      SGF=GF
      AMZ=MZ
      XMZ=MZ
      AMS=AMSB
      XMS=AMS
      XMC=AMC
      XMB=AMB
      XMT=AMT
      XMW=AMW
      XMZ=AMZ
      AMBB=AMB
      IF(NDEC.NE.0)THEN
       AMTT=AMT
      ELSE
       AMTT=1.D10
      ENDIF
      N0 = 5
      XLAMBDA = XITLA(LOOP,ALSMZ,ACC)
      XXLAMBDA = XLAMBDA
      NN0 = N0
C--INTIALIZING RANDOM NUMBER GENERATOR
      CALL RSTART(12,34,56,78)
C--INTIALIZING POLYLOGARITHMS
      CALL BERNINI(NBER)

C--INITIALIZING COMMON BLOCKS
      WGF = GF
      WALPH = ALPHA
      WMTAU = 1.7771D0
      WMMUON = 0.105658389D0
      WMZ = AMZ
      WMW = AMW

C--INITIALIZE PDFSET
      CALL PDFSET(PATHNAME,PDFNAME)

      IF(IPROC.GE.0)THEN

C--WRITING HEADER OF OUTPUT FILE
      WRITE(NOUT,*)
      WRITE(NOUT,*)'PDF: ',PDFNAME
      WRITE(NOUT,*)
      WRITE(NOUT,*)'ALPHA_S(M_Z) = ',ALSMZ
      WRITE(NOUT,*)'LAMBDA_5     = ',XLAMBDA
      WRITE(NOUT,*)
      WRITE(NOUT,*)'VEGAS:'
      WRITE(NOUT,*)'======'
      WRITE(NOUT,*)'ABSERR     = ',ABSERR
      WRITE(NOUT,*)'POINTS     = ',IPOINT
      WRITE(NOUT,*)'ITERATIONS = ',ITERAT
      WRITE(NOUT,*)
      WRITE(NOUT,*)'CUTS:'
      WRITE(NOUT,*)'====='
      WRITE(NOUT,*)'EPS        = ',EPS
      WRITE(NOUT,*)'IMAG       = ',REPS
      WRITE(NOUT,*)
      WRITE(NOUT,*)'M_B        = ',AMB,' GEV'
      WRITE(NOUT,*)'M_T        = ',AMT,' GEV'
      WRITE(NOUT,*)
      WRITE(NOUT,*)'C_T        = ',FACT
      WRITE(NOUT,*)'C_B        = ',FACB
      WRITE(NOUT,*)'C_3        = ',FACTRIA
      WRITE(NOUT,*)'C_TT       = ',FACTT
      WRITE(NOUT,*)'C_G        = ',FACG
      WRITE(NOUT,*)'C_GG       = ',FACGG
      WRITE(NOUT,*)
      IF(FACT.NE.1.D0.OR.FACTRIA.NE.1.D0.OR.FACTT.NE.0.D0.OR.
     .   FACG.NE.0.D0.OR.FACGG.NE.0.D0)THEN
       IF(ISILH.NE.0)THEN
        WRITE(NOUT,*)'SILH'
       ELSE
        WRITE(NOUT,*)'NON-LINEAR'
       ENDIF
       WRITE(NOUT,*)
      ENDIF
      IF(IPARTON.LE.0)THEN
       WRITE(NOUT,*)'ENERGY     = ',W,' GEV'
       WRITE(NOUT,*)
       WRITE(NOUT,*)'HADRON CROSS SECTION'
      ELSE
       WRITE(NOUT,*)'ENERGY     = ',WHAT,' GEV'
       WRITE(NOUT,*)
       WRITE(NOUT,*)'PARTON CROSS SECTION'
      ENDIF
      WRITE(NOUT,*)
      IF(ISUSY.EQ.0)THEN
       WRITE(NOUT,*)'GG --> HH'
      ELSE
       IF(IPROC.EQ.0)THEN
        WRITE(NOUT,*)'GG --> hh'
       ELSEIF(IPROC.EQ.1)THEN
        WRITE(NOUT,*)'GG --> Hh'
       ELSEIF(IPROC.EQ.2)THEN
        WRITE(NOUT,*)'GG --> HH'
       ELSEIF(IPROC.EQ.3)THEN
        WRITE(NOUT,*)'GG --> Ah'
       ELSEIF(IPROC.EQ.4)THEN
        WRITE(NOUT,*)'GG --> AH'
       ELSEIF(IPROC.EQ.5)THEN
        WRITE(NOUT,*)'GG --> AA'
       ELSEIF(IPROC.EQ.6)THEN
        WRITE(NOUT,*)'Q QBAR --> Ah'
       ELSEIF(IPROC.EQ.7)THEN
        WRITE(NOUT,*)'Q QBAR --> AH'
       ENDIF
      ENDIF
      IF(IPROC.LE.5)THEN
       WRITE(NOUT,*)'========='
      ELSE
       WRITE(NOUT,*)'============='
      ENDIF
      WRITE(NOUT,*)
      IF(ITRIA.EQ.1) WRITE(NOUT,*)'NO TRIANGLES'
      IF(ITRIA.EQ.2) WRITE(NOUT,*)'NO Z-TRIANGLE'
      IF(ITRIA.EQ.3) WRITE(NOUT,*)'HEAVY TRIANGLE'
      IF(ITRIA.EQ.4) WRITE(NOUT,*)'LIGHT TRIANGLE'
      IF(ITRIA.NE.0) WRITE(NOUT,*)

      DO 9999 ILOOP=1,NMA
       IF(NMA.GT.1) THEN
        AMA=AMABEG+(AMAEND-AMABEG)/(NMA-1.D0)*(ILOOP-1)
       ELSE
        AMA=AMABEG
       ENDIF

C--INTIALIZING ALPHA_S FOR HIGGS DECAYS [ALPHA_S(MZ) = 0.119]
       XLAMBDA = 0.239D0
       N0 = 5
       CALL ALSINI(ACC)

C--CALCULATE SUSY COUPLINGS AND HIGGS WIDTHS
       IF(ISUSY.EQ.0)THEN
        IPROC=2
        CALL SMCP
        CALL HDEC(TGBET)
        GAMA=SMWDTH
        GAML=SMWDTH
        GAMH=SMWDTH
       ELSE
        FACT    = 1
        FACB    = 1
        FACTRIA = 1
        FACTT   = 0
        FACG    = 0
        FACGG   = 0
        IF(IOLD.NE.1)THEN
         CALL SUSYCP(TGBET)
        ELSE
         CALL SUSYCP0(TGBET)
        ENDIF
        CALL HDEC(TGBET)
        GAMA=AWDTH
        GAML=HLWDTH
        GAMH=HHWDTH
       ENDIF

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     do 9999 iii = 0,40,1
c--68% CL
c      iseterr  = -iii 
c--90% CL
c      iseterr  = iii 

c     n00      = 5
c     xlambda0 = 0.255097d0
c     loop     = 2

c     if(iii.eq.0)then
c      n0 = 0
c      n1 = 4
c     else
c      n0 = 0
c      n1 = 4
c     endif
c     do 9999 jjj = n0,n1
c      ials = jjj

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

       iseterr  = 0 
       ials = 0

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

c--mstw2008 68% CL
c      if(jjj.eq.0)then
c       xlambda0 = 0.255097d0
c      elseif(jjj.eq.1)then
c       xlambda0 = 0.234851d0
c      elseif(jjj.eq.2)then
c       xlambda0 = 0.244764d0
c      elseif(jjj.eq.3)then
c       xlambda0 = 0.263604d0
c      elseif(jjj.eq.4)then
c       xlambda0 = 0.272302d0
c      endif

c--mstw2008 90% CL
c      if(jjj.eq.0)then
c       xlambda0 = 0.255097d0
c      elseif(jjj.eq.1)then
c       xlambda0 = 0.205587d0
c      elseif(jjj.eq.2)then
c       xlambda0 = 0.229424d0
c      elseif(jjj.eq.3)then
c       xlambda0 = 0.277527d0
c      elseif(jjj.eq.4)then
c       xlambda0 = 0.301413d0
c      endif

c      xxlambda = xlambda0

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

C--INTIALIZING ALPHA_S FOR CROSS SECTIONS
       XLAMBDA=XITLA(LOOP,ALSMZ,ACC)
       N0=5
       CALL ALSINI(ACC)
       AMAT=AMT0*1.D8
c      write(6,*)'alpha_s(MZ) = ',alphas(amz,loop)

      DO 9999 ISCALE=1,NSCALE
       IF(NSCALE.GT.1) THEN
        SCA=SCBEG+(SCEND-SCBEG)/(NSCALE-1.D0)*(ISCALE-1)
       ELSE
        SCA=SCBEG
       ENDIF
       ZSS = ZSS0*SCA
       ZPT = ZPT0*SCA
       ZQQ = ZQQ0*SCA

C--DEFINING PROCESS: MASSES AND QCD COEFFICIENTS
       IF(IPROC.EQ.0)THEN
C-- gg -> hh
        M1=AML
        M2=AML
        VC1 = 11.D0/2.D0
        VC2 = 4.D0/9.D0
        VC3 = -4.D0/9.D0
        IHA = 0
       ELSEIF(IPROC.EQ.1)THEN
C-- gg -> Hh
        M1=AMH
        M2=AML
        VC1 = 11.D0/2.D0
        VC2 = 4.D0/9.D0
        VC3 = -4.D0/9.D0
        IHA = 0
       ELSEIF(IPROC.EQ.2)THEN
C-- gg -> HH
        M1=AMH
        M2=AMH
        VC1 = 11.D0/2.D0
        VC2 = 4.D0/9.D0
        VC3 = -4.D0/9.D0
        IHA = 0
       ELSEIF(IPROC.EQ.3)THEN
C-- gg -> hA
        M1=AML
        M2=AMA
        VC1 = 6.D0
        VC2 = 2.D0/3.D0
        VC3 = 2.D0/3.D0
        IHA = 1
       ELSEIF(IPROC.EQ.4)THEN
C-- gg -> HA
        M1=AMH
        M2=AMA
        VC1 = 6.D0
        VC2 = 2.D0/3.D0
        VC3 = 2.D0/3.D0
        IHA = 1
       ELSEIF(IPROC.EQ.5)THEN
C-- gg -> AA
        M1=AMA
        M2=AMA
        VC1 = 11.D0/2.D0
        VC2 = -1.D0
        VC3 = -1.D0
        IHA = 0
       ELSEIF(IPROC.EQ.6)THEN
C-- qqbar -> hA
        M1=AMA
        M2=AML
       ELSEIF(IPROC.EQ.7)THEN
C-- qqbar -> HA
        M1=AMA
        M2=AMH
       ENDIF
       IGRAPH=0
       IF(IPROC.LE.5)THEN
        IF(IPARTON.EQ.0)THEN
        NDIM=3
        IF(IRES.EQ.0)THEN

         SIGLO = 0
         DSIGLO = 0
         SIGV = 0
         DSIGV = 0
         SIGGG = 0
         DSIGGG = 0
         SIGGQ = 0
         DSIGGQ = 0
         SIGQQ = 0
         DSIGQQ = 0

         NDIM=3
         IF(IQCDB.NE.0)THEN
C--INTEGRATION OF BORN TERM
          CALL VEGASN(SIG,ABSERR,NDIM,IPOINT,ITERAT,NPRN,IGRAPH)
          SIGLO = RES
          DSIGLO = ERR
         ENDIF
         IF(IQCDV.NE.0)THEN
C--INTEGRATION OF VIRTUAL CORRECTIONS
          CALL VEGASN(QCDV,ABSERR,NDIM,IPOINT,ITERAT,NPRN,IGRAPH)
          SIGV = RES
          DSIGV = ERR
         ENDIF
         NDIM=4
         IF(IQCDGG.NE.0)THEN
C--INTEGRATION OF REAL CORRECTIONS - GG INITIAL STATE
          CALL VEGASN(QCDGG,ABSERR,NDIM,IPOINT,ITERAT,NPRN,IGRAPH)
          SIGGG = RES
          DSIGGG = ERR
         ENDIF
         IF(IQCDGQ.NE.0)THEN
C--INTEGRATION OF REAL CORRECTIONS - GQ INITIAL STATE
          CALL VEGASN(QCDGQ,ABSERR,NDIM,IPOINT,ITERAT,NPRN,IGRAPH)
          SIGGQ = RES
          DSIGGQ = ERR
         ENDIF
         IF(IQCDQQ.NE.0)THEN
C--INTEGRATION OF REAL CORRECTIONS - QQBAR INITIAL STATE
          CALL VEGASN(QCDQQ,ABSERR,NDIM,IPOINT,ITERAT,NPRN,IGRAPH)
          SIGQQ = RES
          DSIGQQ = ERR
         ENDIF
         NDIM=3
        ELSE
         IF(AMH.LE.(M1+M2))THEN
          RES=0
          ERR=0
         ELSE
C--INTEGRATION OF RESONANCE CONTRIBUTION
          CALL VEGASN(SIGRES,ABSERR,NDIM,IPOINT,ITERAT,NPRN,IGRAPH)
          SIGLO = RES
          DSIGLO = ERR
         ENDIF
        ENDIF
        ELSE IF(IPARTON.EQ.-1)THEN
C--INTEGRATION OF HADRONIC FROM PARTONIC CROSS SECTION
         NDIM=3
         CALL VEGASN(SIGPART,ABSERR,NDIM,IPOINT,ITERAT,NPRN,IGRAPH)
         SIGLO = RES
         DSIGLO = ERR
        ELSE
C--INTEGRATION OF PARTONIC CROSS SECTION
         NDIM=1
         CALL VEGASN(SIGHAT,ABSERR,NDIM,IPOINT,ITERAT,NPRN,IGRAPH)
         SIGLO = RES
         DSIGLO = ERR
        ENDIF
       ELSE
         IF(IQCDB.NE.0)THEN
C--INTEGRATION OF BORN TERM FOR QQBAR -> HA
          NDIM=2
          CALL VEGASN(SIGBQQ,ABSERR,NDIM,IPOINT,ITERAT,NPRN,IGRAPH)
          SIGLO = RES
          DSIGLO = ERR
         ENDIF
         IF(IQCDQQ.NE.0)THEN
C--INTEGRATION OF NLO TERM FOR QQBAR -> HA: QQBAR INITIAL STATE
          NDIM=3
          CALL VEGASN(QSIGQQ,ABSERR,NDIM,IPOINT,ITERAT,NPRN,IGRAPH)
          SIGQQ = RES
          DSIGQQ = ERR
         ENDIF
         IF(IQCDGQ.NE.0)THEN
C--INTEGRATION OF NLO TERM FOR QQBAR -> HA: GQ INITIAL STATE
          NDIM=3
          CALL VEGASN(QSIGGQ,ABSERR,NDIM,IPOINT,ITERAT,NPRN,IGRAPH)
          SIGGQ = RES
          DSIGGQ = ERR
         ENDIF
       ENDIF
       IF(ISUSY.NE.0)THEN
        WRITE(NOUT,*)'TAN(BETA)  = ',TGBET
        WRITE(NOUT,*)'M_A = ',AMA,' GEV    GAMMA = ',GAMA,' GEV'
        WRITE(NOUT,*)'M_h = ',AML,' GEV    GAMMA = ',GAML,' GEV'
       ENDIF
       WRITE(NOUT,*)'M_H = ',AMH,' GEV    GAMMA = ',GAMH,' GEV'
       WRITE(NOUT,*)'SCALE = ',SCA
       WRITE(NOUT,*)
       WRITE(NOUT,*)'SIGMA_BORN = (',SIGLO,' +- ',DSIGLO,') PB'

       IF(IPROC.LE.5.AND.IPARTON.EQ.0.AND.IRES.EQ.0)THEN

       SIGNLO = SIGLO + SIGV + SIGGG + SIGGQ + SIGQQ
       DSIGNLO = DSQRT(DSIGLO**2 + DSIGV**2 + DSIGGG**2
     .               + DSIGGQ**2 + DSIGQQ**2)
       WRITE(NOUT,*)'SIGMA_VIRT = (',SIGV,' +- ',DSIGV,') PB'
       WRITE(NOUT,*)'SIGMA_GG   = (',SIGGG,' +- ',DSIGGG,') PB'
       WRITE(NOUT,*)'SIGMA_GQ   = (',SIGGQ,' +- ',DSIGGQ,') PB'
       WRITE(NOUT,*)'SIGMA_QQ   = (',SIGQQ,' +- ',DSIGQQ,') PB'
       WRITE(NOUT,*)'    K_VIRT =  ',SIGV/SIGLO
       WRITE(NOUT,*)'    K_GG   =  ',SIGGG/SIGLO
       WRITE(NOUT,*)'    K_GQ   =  ',SIGGQ/SIGLO
       WRITE(NOUT,*)'    K_QQ   =  ',SIGQQ/SIGLO
       WRITE(NOUT,*)
       WRITE(NOUT,*)'SIGMA_NLO  = (',SIGNLO,' +- ',DSIGNLO,') PB'
       WRITE(NOUT,*)'    K_NLO  =  ',SIGNLO/SIGLO

       ELSEIF(IPROC.GT.5)THEN

       SIGNLO = SIGLO + SIGQQ + SIGGQ
       DSIGNLO = DSQRT(DSIGLO**2 + DSIGQQ**2 + DSIGGQ**2)
       WRITE(NOUT,*)'SIGMA_QQ   = (',SIGQQ,' +- ',DSIGQQ,') PB'
       WRITE(NOUT,*)'SIGMA_GQ   = (',SIGGQ,' +- ',DSIGGQ,') PB'
       WRITE(NOUT,*)'    K_QQ   =  ',SIGQQ/SIGLO
       WRITE(NOUT,*)'    K_GQ   =  ',SIGGQ/SIGLO
       WRITE(NOUT,*)
       WRITE(NOUT,*)'SIGMA_NLO  = (',SIGNLO,' +- ',DSIGNLO,') PB'
       WRITE(NOUT,*)'    K_NLO  =  ',SIGNLO/SIGLO

       ENDIF
       WRITE(NOUT,*)
c      write(50,200)amh,siglo,signlo,alphas(amz,2),iseterr,ials
c      write(6,200)amh,siglo,signlo,alphas(amz,2),iseterr,ials
c      write(6,200)amh,siglo,signlo,alphas(amz,1),iseterr,ials
       call flush(6)
       call flush(50)
       call flush(nout)
9999  CONTINUE
      ELSEIF(IPROC.EQ.-1)THEN
C--SCALAR INTEGRALS
       M1=AM1
       M2=AM2
       CALL VIRTINT
      ENDIF
      ENDIF

100   FORMAT(10X,I30)
101   FORMAT(10X,G30.20)
102   FORMAT(12X,A100)
200   format(1x,f10.0,3(3x,g12.5),3x,i3,3x,i3)
      CLOSE(NIN)
      CLOSE(NOUT)
      CLOSE(NV)
      STOP
      END

      DOUBLE PRECISION FUNCTION QCDV(XX)
C--FUNCTION FOR GG -> HH: VIRTUAL CORRECTIONS
      PARAMETER(N=3)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DOUBLE PRECISION M1,M2,MZ
      DIMENSION XX(N),Y(N)
      COMPLEX*16 A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/SCALE/ZSS,ZPT,ZQQ,QQ,ISUB
      COMMON/ALSN/LOOP
      COMMON/CUTS/EPS,REPS
      COMMON/PARAM/S,T,U,M1,M2,MZ
      COMMON/PROC/IPROC,IAPP,ITOP,ISFAC
      COMMON/HMASS/AMA,AML,AMH,AMCH
      COMMON/SUSY/TGBET
      COMMON/FORM/A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/TOTPAR/W,GEVPB,GF,GAMZ,GAMA,GAML,GAMH
      COMMON/PARTONIC/WHAT
      COMMON/VIRTCOF/VC1,VC2,VC3,IHA
      PI=4.D0*DATAN(1.D0)
      DO 1 I=1,N-1
       Y(I)=EPS+(1.D0-2.D0*EPS)*XX(I)
1     CONTINUE
      Y(N) = XX(N)
      TH=(M1+M2)**2/W**2
      IF(ISUB.EQ.0)THEN
       TAU=TH+(1.D0-TH)*Y(1)
       X=TAU+(1.D0-TAU)*Y(2)
       DJAC=(1.D0-TH)*(1.D0-TAU)*(1.D0-2.D0*EPS)**(N-1)
      ELSE
       YY=-DLOG(TH)*Y(1)
       TAU=DEXP(-YY)
       XY=-DLOG(TAU)*Y(2)
       X=DEXP(-XY)
       DJAC=DLOG(TH)*DLOG(TAU)*TAU*X*(1.D0-2.D0*EPS)**(N-1)
      ENDIF
      S=TAU*W**2
      WHAT = DSQRT(S)
      Q=ZSS*DSQRT(S) + ZQQ*QQ
      SCLOG = DLOG(Q**2/S)
      DUMMY1 = SIGB(1.D0,XX(N))*(PI**2 + VC1 + (33-2*5)/6.D0*SCLOG)
      DUMMY2 = SIG2(1.D0,XX(N)) * VC2
      DUMMY3 = SIG3(1.D0,XX(N)) * VC3
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     dummy1 = 0
c     dummy2 = 0
c     dummy3 = 0
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DUMMY = DUMMY1 + DUMMY2 + DUMMY3
      QCDV = DJAC*DLUMGG(TAU,X,Q)*DUMMY
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION QCDGG(XX)
C--FUNCTION FOR GG -> HH: GG REAL CORRECTIONS
      PARAMETER(N=4)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DOUBLE PRECISION M1,M2,MZ
      DIMENSION XX(N),Y(N)
      COMPLEX*16 A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/SCALE/ZSS,ZPT,ZQQ,QQ,ISUB
      COMMON/ALSN/LOOP
      COMMON/CUTS/EPS,REPS
      COMMON/PARAM/S,T,U,M1,M2,MZ
      COMMON/PROC/IPROC,IAPP,ITOP,ISFAC
      COMMON/HMASS/AMA,AML,AMH,AMCH
      COMMON/SUSY/TGBET
      COMMON/FORM/A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/TOTPAR/W,GEVPB,GF,GAMZ,GAMA,GAML,GAMH
      COMMON/PARTONIC/WHAT
      PRGG(X) = 6*(1/X - 2 + X*(1-X))
      PDGG(X) = 6/(1-X)
      WGG(X) = -11*(1-X)**3/2
      D1(X) = DLOG(1-X)/(1-X)
      WDGG(X) = 1+X**4+(1-X)**4
      RD1(X) = -DLOG(1-X)**2/2
      RPDGG(X) = -6*DLOG(1-X)
      PI=4.D0*DATAN(1.D0)
      DO 1 I=1,N-1
       Y(I)=EPS+(1.D0-2.D0*EPS)*XX(I)
1     CONTINUE
      Y(N) = XX(N)
      TH=(M1+M2)**2/W**2
      IF(ISUB.EQ.0)THEN
       VAU=TH+(1.D0-TH)*Y(1)
       X=VAU+(1.D0-VAU)*Y(2)
       DJAC=(1.D0-TH)*(1.D0-VAU)*(1.D0-2.D0*EPS)**(N-2)
      ELSE
       YY=-DLOG(TH)*Y(1)
       VAU=DEXP(-YY)
       XY=-DLOG(VAU)*Y(2)
       X=DEXP(-XY)
       DJAC=DLOG(TH)*DLOG(VAU)*VAU*X*(1.D0-2.D0*EPS)**(N-2)
      ENDIF
      Z = VAU/X + (1-VAU/X)*Y(3)
      TAU = VAU/Z
      FACTOR = (1-VAU/X)/Z
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     TAU = VAU
c     Z = TH/TAU + (1-TH/TAU)*Y(3)
c     FACTOR = (1-TH/TAU)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      S=TAU*W**2
      WHAT = DSQRT(S)
      Q=ZSS*DSQRT(Z*S) + ZQQ*QQ
      Q0=Q
      SCLOG = DLOG(Q**2/S)
      SCLOG0 = DLOG(Q0**2/Z/S)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     Q0=ZSS*DSQRT(S) + ZQQ*QQ
c     SCLOG0 = DLOG(Q0**2/S)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DLGG = DLUMGG(TAU,X,Q)
      DLGG0 = DLUMGG(VAU,X,Q0)
      DJACZ = FACTOR*(1-2*EPS)
      DSIGZ = SIGB(Z,XX(4))
      WHAT = DSQRT(Z*TAU*W**2)
      DSIG1 = SIGB(1.D0,XX(4))
      WHAT = DSQRT(TAU*W**2)
      REST = VAU/X
      ZZ = Z
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     DSIGZ = SIGB(Z,XX(4))
c     DSIG1 = SIGB(1.D0,XX(4))
c     REST = TH/TAU
c     ZZ = 1.D0
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DUMMY1 = DSIGZ/Z*WGG(Z)*DJACZ*DLGG
      DUMMY2 = DSIGZ/Z*(-Z*PRGG(Z)*SCLOG)*DJACZ*DLGG
      DUMMY21 = DSIG1*(-(33-2*5)/6.D0*SCLOG0)*DLGG0
      DUMMY22 = (DSIGZ*SCLOG*DLGG - DSIG1*SCLOG0*DLGG0*ZZ)
     .       * (-PDGG(Z)) * DJACZ
      DUMMY23 = - DSIG1 * (-RPDGG(REST)*SCLOG0)*DLGG0
      DUMMY3 =(WDGG(Z)/Z*DSIGZ*DLGG-2*DSIG1*DLGG0*ZZ)
     .       * (6*D1(Z)) * DJACZ
      DUMMY31 = -2*DSIG1 * (6*RD1(REST))*DLGG0
      DUMMY2 = DUMMY2 + DUMMY21 + DUMMY22 + DUMMY23
      DUMMY3 = DUMMY3 + DUMMY31
      DUMMY = DUMMY1 + DUMMY2 + DUMMY3
      QCDGG = DJAC*DUMMY
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION QCDGQ(XX)
C--FUNCTION FOR GG -> HH: GQ REAL CORRECTIONS
      PARAMETER(N=4)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DOUBLE PRECISION M1,M2,MZ
      DIMENSION XX(N),Y(N)
      COMPLEX*16 A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/SCALE/ZSS,ZPT,ZQQ,QQ,ISUB
      COMMON/ALSN/LOOP
      COMMON/CUTS/EPS,REPS
      COMMON/PARAM/S,T,U,M1,M2,MZ
      COMMON/PROC/IPROC,IAPP,ITOP,ISFAC
      COMMON/HMASS/AMA,AML,AMH,AMCH
      COMMON/SUSY/TGBET
      COMMON/FORM/A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/TOTPAR/W,GEVPB,GF,GAMZ,GAMA,GAML,GAMH
      COMMON/PARTONIC/WHAT
      PGQ(X) = 4.D0/3.D0*(1 + (1-X)**2)/X
      WGQ(X) = 2.D0/3.D0*X**2-(1-X)**2
      PI=4.D0*DATAN(1.D0)
      DO 1 I=1,N-1
       Y(I)=EPS+(1.D0-2.D0*EPS)*XX(I)
1     CONTINUE
      Y(N) = XX(N)
      TH=(M1+M2)**2/W**2
      IF(ISUB.EQ.0)THEN
       VAU=TH+(1.D0-TH)*Y(1)
       X=VAU+(1.D0-VAU)*Y(2)
       DJAC=(1.D0-TH)*(1.D0-VAU)*(1.D0-2.D0*EPS)**(N-2)
      ELSE
       YY=-DLOG(TH)*Y(1)
       VAU=DEXP(-YY)
       XY=-DLOG(VAU)*Y(2)
       X=DEXP(-XY)
       DJAC=DLOG(TH)*DLOG(VAU)*VAU*X*(1.D0-2.D0*EPS)**(N-2)
      ENDIF
      Z = VAU/X + (1-VAU/X)*Y(3)
      TAU = VAU/Z
      FACTOR = (1-VAU/X)/Z
      S=TAU*W**2
      WHAT = DSQRT(S)
c     Z = TH/TAU + (1-TH/TAU)*Y(3)
c     FACTOR = (1-TH/TAU)
      Q=ZSS*DSQRT(Z*S) + ZQQ*QQ
      SCLOG = DLOG(Q**2/S/(1-Z)**2)
      DJACZ = FACTOR*(1-2*EPS)
      DUMMY = SIGB(Z,XX(4))/Z*(WGQ(Z)-Z/2*PGQ(Z)*SCLOG)*DJACZ
      QCDGQ = DJAC*DLUMGQ(TAU,X,Q)*DUMMY
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION QCDQQ(XX)
C--FUNCTION FOR GG -> HH: QQBAR REAL CORRECTIONS
      PARAMETER(N=4)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DOUBLE PRECISION M1,M2,MZ
      DIMENSION XX(N),Y(N)
      COMPLEX*16 A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/SCALE/ZSS,ZPT,ZQQ,QQ,ISUB
      COMMON/ALSN/LOOP
      COMMON/CUTS/EPS,REPS
      COMMON/PARAM/S,T,U,M1,M2,MZ
      COMMON/PROC/IPROC,IAPP,ITOP,ISFAC
      COMMON/HMASS/AMA,AML,AMH,AMCH
      COMMON/SUSY/TGBET
      COMMON/FORM/A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/TOTPAR/W,GEVPB,GF,GAMZ,GAMA,GAML,GAMH
      COMMON/PARTONIC/WHAT
      WQQ(X) = 32.D0/27.D0*(1-X)**3
      PI=4.D0*DATAN(1.D0)
      DO 1 I=1,N-1
       Y(I)=EPS+(1.D0-2.D0*EPS)*XX(I)
1     CONTINUE
      Y(N) = XX(N)
      TH=(M1+M2)**2/W**2
      IF(ISUB.EQ.0)THEN
       VAU=TH+(1.D0-TH)*Y(1)
       X=VAU+(1.D0-VAU)*Y(2)
       DJAC=(1.D0-TH)*(1.D0-VAU)*(1.D0-2.D0*EPS)**(N-2)
      ELSE
       YY=-DLOG(TH)*Y(1)
       VAU=DEXP(-YY)
       XY=-DLOG(VAU)*Y(2)
       X=DEXP(-XY)
       DJAC=DLOG(TH)*DLOG(VAU)*VAU*X*(1.D0-2.D0*EPS)**(N-2)
      ENDIF
      Z = VAU/X + (1-VAU/X)*Y(3)
      TAU = VAU/Z
      FACTOR = (1-VAU/X)/Z
      S=TAU*W**2
      WHAT = DSQRT(S)
c     Z = TH/TAU + (1-TH/TAU)*Y(3)
c     FACTOR = (1-TH/TAU)
      Q=ZSS*DSQRT(Z*S) + ZQQ*QQ
      DJACZ = FACTOR*(1-2*EPS)
      DUMMY = SIGB(Z,XX(4))/Z*WQQ(Z)*DJACZ
      QCDQQ = DJAC*DLUMQQ(TAU,X,Q)*DUMMY
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION SIGB(Z,XX)
C--FUNCTION FOR GG -> HH: BORN TERM
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DOUBLE PRECISION M1,M2,MZ
      COMPLEX*16 A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/SCALE/ZSS,ZPT,ZQQ,QQ,ISUB
      COMMON/ALSN/LOOP
      COMMON/CUTS/EPS,REPS
      COMMON/PARAM/S,T,U,M1,M2,MZ
      COMMON/PROC/IPROC,IAPP,ITOP,ISFAC
      COMMON/HMASS/AMA,AML,AMH,AMCH
      COMMON/SUSY/TGBET
      COMMON/FORM/A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/TOTPAR/W,GEVPB,GF,GAMZ,GAMA,GAML,GAMH
      COMMON/PARTONIC/WHAT
      PI=4.D0*DATAN(1.D0)
      Y=EPS+(1.D0-2.D0*EPS)*XX
      TH=(M1+M2)**2/W**2
      TAU = WHAT**2/W**2
      VAU = Z*TAU
      DJAC = 1-2*EPS
      V=Y/2.D0
      S=VAU*W**2
      XLAM=DSQRT((S+M1**2-M2**2)**2-4.D0*S*M1**2)
      BET=XLAM/(S+M1**2-M2**2)
      T1=-(S+M1**2-M2**2)*((1.D0+BET)/2.D0-BET*V)
      U1=-(S+M1**2-M2**2)*((1.D0-BET)/2.D0+BET*V)
      T=T1+M1**2
      U=U1+M1**2
      PT2=(T1*U1-S*M1**2)/S
      AM2=(M1+M2)**2/4.D0
      Q=ZSS*DSQRT(S) + ZPT*DSQRT(AM2+PT2) + ZQQ*QQ
      IF(IPROC.EQ.0.OR.IPROC.EQ.2.OR.IPROC.EQ.5)DJAC=DJAC/2.D0
      CALL MATRIX(DMAT)
      ALPS=ALPHAS(Q,LOOP)
      QCD=2.D0*ALPS**2/(4.D0*PI)**2
      SIGB=GEVPB*QCD*DJAC*XLAM/S**2/4096.D0/PI*DMAT * ALPS/PI
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION SIG2(Z,XX)
C--FUNCTION FOR GG -> HH: C2 TERM OF VIRTUAL CORRECTIONS
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DOUBLE PRECISION M1,M2,MZ
      COMPLEX*16 A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/SCALE/ZSS,ZPT,ZQQ,QQ,ISUB
      COMMON/ALSN/LOOP
      COMMON/CUTS/EPS,REPS
      COMMON/PARAM/S,T,U,M1,M2,MZ
      COMMON/PROC/IPROC,IAPP,ITOP,ISFAC
      COMMON/HMASS/AMA,AML,AMH,AMCH
      COMMON/SUSY/TGBET
      COMMON/FORM/A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/TOTPAR/W,GEVPB,GF,GAMZ,GAMA,GAML,GAMH
      COMMON/PARTONIC/WHAT
      PI=4.D0*DATAN(1.D0)
      Y=EPS+(1.D0-2.D0*EPS)*XX
      TH=(M1+M2)**2/W**2
      TAU = WHAT**2/W**2
      DJAC = 1-2*EPS
      V=Y/2.D0
      S=Z*TAU*W**2
      XLAM=DSQRT((S+M1**2-M2**2)**2-4.D0*S*M1**2)
      BET=XLAM/(S+M1**2-M2**2)
      T1=-(S+M1**2-M2**2)*((1.D0+BET)/2.D0-BET*V)
      U1=-(S+M1**2-M2**2)*((1.D0-BET)/2.D0+BET*V)
      T=T1+M1**2
      U=U1+M1**2
      PT2=(T1*U1-S*M1**2)/S
      AM2=(M1+M2)**2/4.D0
      Q=ZSS*DSQRT(S) + ZPT*DSQRT(AM2+PT2) + ZQQ*QQ
      IF(IPROC.EQ.0.OR.IPROC.EQ.2.OR.IPROC.EQ.5)DJAC=DJAC/2.D0
      CALL MATRIX2(DMAT)
      ALPS=ALPHAS(Q,LOOP)
      QCD=2.D0*ALPS**2/(4.D0*PI)**2
      SIG2=GEVPB*QCD*DJAC*XLAM/S**2/4096.D0/PI*DMAT * ALPS/PI
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION SIG3(Z,XX)
C--FUNCTION FOR GG -> HH: C3 TERM OF VIRTUAL CORRECTIONS
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DOUBLE PRECISION M1,M2,MZ
      COMPLEX*16 A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/SCALE/ZSS,ZPT,ZQQ,QQ,ISUB
      COMMON/ALSN/LOOP
      COMMON/CUTS/EPS,REPS
      COMMON/PARAM/S,T,U,M1,M2,MZ
      COMMON/PROC/IPROC,IAPP,ITOP,ISFAC
      COMMON/HMASS/AMA,AML,AMH,AMCH
      COMMON/SUSY/TGBET
      COMMON/FORM/A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/TOTPAR/W,GEVPB,GF,GAMZ,GAMA,GAML,GAMH
      COMMON/PARTONIC/WHAT
      PI=4.D0*DATAN(1.D0)
      Y=EPS+(1.D0-2.D0*EPS)*XX
      TH=(M1+M2)**2/W**2
      TAU = WHAT**2/W**2
      DJAC = 1-2*EPS
      V=Y/2.D0
      S=Z*TAU*W**2
      XLAM=DSQRT((S+M1**2-M2**2)**2-4.D0*S*M1**2)
      BET=XLAM/(S+M1**2-M2**2)
      T1=-(S+M1**2-M2**2)*((1.D0+BET)/2.D0-BET*V)
      U1=-(S+M1**2-M2**2)*((1.D0-BET)/2.D0+BET*V)
      T=T1+M1**2
      U=U1+M1**2
      PT2=(T1*U1-S*M1**2)/S
      AM2=(M1+M2)**2/4.D0
      Q=ZSS*DSQRT(S) + ZPT*DSQRT(AM2+PT2) + ZQQ*QQ
      IF(IPROC.EQ.0.OR.IPROC.EQ.2.OR.IPROC.EQ.5)DJAC=DJAC/2.D0
      CALL MATRIX3(DMAT)
      ALPS=ALPHAS(Q,LOOP)
      QCD=2.D0*ALPS**2/(4.D0*PI)**2
      SIG3=GEVPB*QCD*DJAC*XLAM/S**2/4096.D0/PI*DMAT * ALPS/PI
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION SIGBQQ(XX)
C--FUNCTION FOR QQBAR -> HA: BORN TERM
      PARAMETER(N=2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION M1,M2,MZ
      DIMENSION XX(N),Y(N)
      COMMON/SCALE/ZSS,ZPT,ZQQ,QQ,ISUB
      COMMON/CUTS/EPS,REPS
      COMMON/PARAM/S,T,U,M1,M2,MZ
      COMMON/PROC/IPROC,IAPP,ITOP,ISFAC
      COMMON/HMASS/AMA,AML,AMH,AMCH
      COMMON/TOTPAR/W,GEVPB,GF,GAMZ,GAMA,GAML,GAMH
      COMMON/COUP/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,
     .            GHHH,GLLL,GHLL,GLHH,GHAA,GLAA,GLVV,GHVV,
     .            GLPM,GHPM,GZT,GZB
      PI=4.D0*DATAN(1.D0)
      VV=1.D0/DSQRT(DSQRT(2.D0)*GF)
      DO 1 I=1,N
       Y(I)=EPS+(1.D0-2.D0*EPS)*XX(I)
1     CONTINUE
      TH=(M1+M2)**2/W**2
      IF(ISUB.EQ.0)THEN
       TAU=TH+(1.D0-TH)*Y(1)
       X=TAU+(1.D0-TAU)*Y(2)
       DJAC=(1.D0-TH)*(1.D0-TAU)*(1.D0-2.D0*EPS)**N
      ELSE
       YY=-DLOG(TH)*Y(1)
       TAU=DEXP(-YY)
       XY=-DLOG(TAU)*Y(2)
       X=DEXP(-XY)
       DJAC=DLOG(TH)*DLOG(TAU)*TAU*X*(1.D0-2.D0*EPS)**N
      ENDIF
      S=TAU*W**2
      QSQ = S
      Q = ZSS*DSQRT(S) + ZQQ*QQ
      SIGBQQ=GEVPB*DJAC*DLUMQQZ(TAU,X,Q)*SIG0QQ(QSQ)
      RETURN
      END
  
      DOUBLE PRECISION FUNCTION SIG0QQ(QSQ)
C--FUNCTION FOR QQBAR -> HA: BORN FACTOR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION M1,M2,MZ
      COMMON/PARAM/SS,T,U,M1,M2,MZ
      COMMON/PROC/IPROC,IAPP,ITOP,ISFAC
      COMMON/HMASS/AMA,AML,AMH,AMCH
      COMMON/TOTPAR/W,GEVPB,GF,GAMZ,GAMA,GAML,GAMH
      COMMON/COUP/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,
     .            GHHH,GLLL,GHLL,GLHH,GHAA,GLAA,GLVV,GHVV,
     .            GLPM,GHPM,GZT,GZB
      PI=4.D0*DATAN(1.D0)
      VV=1.D0/DSQRT(DSQRT(2.D0)*GF)
      S=QSQ
      XLAM=DSQRT((S+M1**2-M2**2)**2-4.D0*S*M1**2)
      BET=XLAM/S
      FACZ=1
      IF(ISFAC.EQ.1)THEN
       FACZ=S/MZ**2
      ENDIF
      PROZ = (S-MZ**2)**2 + MZ**2*GAMZ**2*FACZ**2
      DMAT = MZ**4*BET**3*GF**2*S/288.D0/PI/PROZ
      IF(IPROC.EQ.6)THEN
       XI = (GZAL*VV/MZ)**2
      ELSE
       XI = (GZAH*VV/MZ)**2
      ENDIF
      SIG0QQ=DMAT*XI
      RETURN
      END
  
      DOUBLE PRECISION FUNCTION QSIGQQ(XX)
C--FUNCTION FOR QQBAR -> HA: NLO QQBAR TERM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(3)
      DOUBLE PRECISION M1,M2,MZ
      COMMON/PARAM/S,T,U,M1,M2,MZ
      COMMON/SCALE/ZSS,ZPT,ZQQ,QQ,ISUB
      COMMON/TOTPAR/W,GEVPB,GF,GAMZ,GAMA,GAML,GAMH
      COMMON/CUTS/EPS,REPS
      EXTERNAL F0,F1
      PI = 4*DATAN(1.D0)
      ZETA2 = PI**2/6
      TAU = (M1+M2)**2/W**2
      X0 = EPS + (1-2*EPS)*XX(1)
      Y0 = EPS + (1-2*EPS)*XX(2)
      Z0 = EPS + (1-2*EPS)*XX(3)
      IF(ISUB.EQ.0)THEN
       X = TAU + (1-TAU)*X0
       Y = DEXP(DLOG(X)*Y0)
       Z = TAU/X + (1-TAU/X)*Z0
       FAC = -(1-TAU)*DLOG(X)*Y*(1-TAU/X)*(1-2*EPS)**3
       FAC0 = -(1-TAU)*DLOG(X)*Y*(1-2*EPS)**2
      ELSE
       X = DEXP(DLOG(TAU)*X0)
       Y = DEXP(DLOG(X)*Y0)
       Z = TAU/X + (1-TAU/X)*Z0
       FAC = DLOG(TAU)*X*DLOG(X)*Y*(1-TAU/X)*(1-2*EPS)**3
       FAC0 = DLOG(TAU)*X*DLOG(X)*Y*(1-2*EPS)**2
      ENDIF
      QSQ = X*W**2
      SCA = ZSS*DSQRT(QSQ) + ZQQ*QQ
      ALP = ALPHAS(SCA,2)/PI
      PD0 = DLUMQQZ(X,Y,SCA)
      COUP = SIG0QQ(QSQ)*ALP*PD0*FAC0
      W0 = (-2*DLOG(SCA**2/X/W**2) + 8.D0/3.D0*(ZETA2-2))*COUP
      RQQ = (-1-Z)
      QSQ = X*Z*W**2
      SCA = ZSS*DSQRT(QSQ) + ZQQ*QQ
      ALP = ALPHAS(SCA,2)/PI
      PD = DLUMQQZ(X,Y,SCA)
      COUP = SIG0QQ(QSQ)*ALP*PD*FAC
      W1 = 4.D0/3.D0*(-RQQ*DLOG(SCA**2/X/W**2)-2*(1+Z)*DLOG(1-Z))*COUP
      RES0 = QD0(TAU,X,Z,PD,PD0)*FAC
      RES1 = QD1(TAU,X,Z,PD,PD0)*FAC
      W2 = RES0 + RES1
      WW = W0 + W1 + W2
      QSIGQQ = GEVPB*WW
      RETURN
      END

      DOUBLE PRECISION FUNCTION QSIGGQ(XX)
C--FUNCTION FOR QQBAR -> HA: NLO GQ TERM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(3)
      DOUBLE PRECISION M1,M2,MZ
      COMMON/PARAM/S,T,U,M1,M2,MZ
      COMMON/SCALE/ZSS,ZPT,ZQQ,QQ,ISUB
      COMMON/TOTPAR/W,GEVPB,GF,GAMZ,GAMA,GAML,GAMH
      COMMON/CUTS/EPS,REPS
      PI = 4*DATAN(1.D0)
      TAU = (M1+M2)**2/W**2
      X0 = EPS + (1-2*EPS)*XX(1)
      Y0 = EPS + (1-2*EPS)*XX(2)
      Z0 = EPS + (1-2*EPS)*XX(3)
      IF(ISUB.EQ.0)THEN
       X = TAU + (1-TAU)*X0
       Y = DEXP(DLOG(X)*Y0)
       Z = TAU/X + (1-TAU/X)*Z0
       FAC = -(1-TAU)*DLOG(X)*Y*(1-TAU/X)*(1-2*EPS)**3
      ELSE
       X = DEXP(DLOG(TAU)*X0)
       Y = DEXP(DLOG(X)*Y0)
       Z = TAU/X + (1-TAU/X)*Z0
       FAC = DLOG(TAU)*X*DLOG(X)*Y*(1-TAU/X)*(1-2*EPS)**3
      ENDIF
      QSQ = X*Z*W**2
      SCA = ZSS*DSQRT(QSQ) + ZQQ*QQ
      PD = DLUMGQZ(X,Y,SCA)
      COUP = SIG0QQ(QSQ)
      PGQ = (Z**2 + (1-Z)**2)/2
      WW = (-PGQ/2*DLOG(SCA**2/X/W**2/(1-Z)**2)+(1+6*Z-7*Z**2)/8)*COUP
      ALP = ALPHAS(SCA,2)/PI
      QSIGGQ = GEVPB*WW*ALP*PD*FAC
      RETURN
      END

      DOUBLE PRECISION FUNCTION QD0(TAU0,TAU,X,PD,PD0)
C--PLUS DISTRIBUTION 1/(1-X)_+ FOR QQBAR -> HA: NLO QQBAR TERM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      FFZ = F0(TAU,X)*PD
      FF1 = F0(TAU,1.D0)*PD0
      QD0 = 1/(1-X)*(FFZ - FF1) + DLOG(1-TAU0/TAU)/(1-TAU0/TAU)*FF1
      RETURN
      END

      DOUBLE PRECISION FUNCTION QD1(TAU0,TAU,X,PD,PD0)
C--PLUS DISTRIBUTION [LOG(1-X)/(1-X)]_+ FOR QQBAR -> HA: NLO QQBAR TERM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      FFZ = F1(TAU,X)*PD
      FF1 = F1(TAU,1.D0)*PD0
      QD1 = DLOG(1-X)/(1-X)*(FFZ - FF1)
     .   + DLOG(1-TAU0/TAU)**2/2/(1-TAU0/TAU)*FF1
      RETURN
      END

      DOUBLE PRECISION FUNCTION F0(TAU,Z)
C--FUNCTION FOR QD0(TAU0,TAU,X,PD,PD0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION M1,M2,MZ
      COMMON/PARAM/S,T,U,M1,M2,MZ
      COMMON/SCALE/ZSS,ZPT,ZQQ,QQ,ISUB
      COMMON/TOTPAR/W,GEVPB,GF,GAMZ,GAMA,GAML,GAMH
      QSQ = TAU*Z*W**2
      PI = 4*DATAN(1.D0)
      SCA = ZSS*DSQRT(QSQ) + ZQQ*QQ
      ALP = ALPHAS(SCA,2)/PI
      SLG = DLOG(SCA**2/TAU/W**2)
      F0 = -2*4.D0/3.D0*SIG0QQ(QSQ)*SLG*ALP
      RETURN
      END

      DOUBLE PRECISION FUNCTION F1(TAU,Z)
C--FUNCTION FOR QD1(TAU0,TAU,X,PD,PD0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION M1,M2,MZ
      COMMON/PARAM/S,T,U,M1,M2,MZ
      COMMON/SCALE/ZSS,ZPT,ZQQ,QQ,ISUB
      COMMON/TOTPAR/W,GEVPB,GF,GAMZ,GAMA,GAML,GAMH
      QSQ = TAU*Z*W**2
      PI = 4*DATAN(1.D0)
      SCA = ZSS*DSQRT(QSQ) + ZQQ*QQ
      ALP = ALPHAS(SCA,2)/PI
      F1 = 4*4.D0/3.D0*SIG0QQ(QSQ)*ALP
      RETURN
      END

      DOUBLE PRECISION FUNCTION SIGHAT(XX)
C--FUNCTION FOR PARTONIC CROSS SECTION
      PARAMETER(N=1)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DOUBLE PRECISION M1,M2,MZ
      DIMENSION XX(N),Y(N)
      COMPLEX*16 A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/SCALE/ZSS,ZPT,ZQQ,QQ,ISUB
      COMMON/ALSN/LOOP
      COMMON/CUTS/EPS,REPS
      COMMON/PARAM/S,T,U,M1,M2,MZ
      COMMON/PROC/IPROC,IAPP,ITOP,ISFAC
      COMMON/HMASS/AMA,AML,AMH,AMCH
      COMMON/SUSY/TGBET
      COMMON/FORM/A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/TOTPAR/W,GEVPB,GF,GAMZ,GAMA,GAML,GAMH
      COMMON/PARTONIC/WHAT
      PI=4.D0*DATAN(1.D0)
      DO 1 I=1,N
       Y(I)=EPS+(1.D0-2.D0*EPS)*XX(I)
1     CONTINUE
      TH=(M1+M2)**2/W**2
      TAU = WHAT**2/W**2
      DJAC = 1-2*EPS
      V=Y(1)/2.D0
      S=TAU*W**2
      XLAM=DSQRT((S+M1**2-M2**2)**2-4.D0*S*M1**2)
      BET=XLAM/(S+M1**2-M2**2)
      T1=-(S+M1**2-M2**2)*((1.D0+BET)/2.D0-BET*V)
      U1=-(S+M1**2-M2**2)*((1.D0-BET)/2.D0+BET*V)
      T=T1+M1**2
      U=U1+M1**2
      PT2=(T1*U1-S*M1**2)/S
      AM2=(M1+M2)**2/4.D0
      Q=ZSS*DSQRT(S) + ZPT*DSQRT(AM2+PT2) + ZQQ*QQ
      IF(IPROC.EQ.0.OR.IPROC.EQ.2.OR.IPROC.EQ.5)DJAC=DJAC/2.D0
      CALL MATRIX(DMAT)
      ALPS=ALPHAS(Q,LOOP)
      QCD=2.D0*ALPS**2/(4.D0*PI)**2
      SIGHAT=GEVPB*QCD*DJAC*XLAM/S**2/4096.D0/PI*DMAT
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION SIGPART(XX)
C--FUNCTION FOR HADRONIC FROM PARTONIC CROSS SECTION
      PARAMETER(N=3)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DOUBLE PRECISION M1,M2,MZ
      DIMENSION XX(N),Y(N)
      COMPLEX*16 A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/SCALE/ZSS,ZPT,ZQQ,QQ,ISUB
      COMMON/ALSN/LOOP
      COMMON/CUTS/EPS,REPS
      COMMON/PARAM/S,T,U,M1,M2,MZ
      COMMON/PROC/IPROC,IAPP,ITOP,ISFAC
      COMMON/HMASS/AMA,AML,AMH,AMCH
      COMMON/SUSY/TGBET
      COMMON/FORM/A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/TOTPAR/W,GEVPB,GF,GAMZ,GAMA,GAML,GAMH
      COMMON/PARTONIC/WHAT
      PI=4.D0*DATAN(1.D0)
      DO 1 I=1,N-1
       Y(I)=EPS+(1.D0-2.D0*EPS)*XX(I)
1     CONTINUE
      Y(N) = XX(N)
      TH=(M1+M2)**2/W**2
      IF(ISUB.EQ.0)THEN
       TAU=TH+(1.D0-TH)*Y(1)
       X=TAU+(1.D0-TAU)*Y(2)
       DJAC=(1.D0-TH)*(1.D0-TAU)*(1.D0-2.D0*EPS)**(N-1)
      ELSE
       YY=-DLOG(TH)*Y(1)
       TAU=DEXP(-YY)
       XY=-DLOG(TAU)*Y(2)
       X=DEXP(-XY)
       DJAC=DLOG(TH)*DLOG(TAU)*TAU*X*(1.D0-2.D0*EPS)**(N-1)
      ENDIF
      S=TAU*W**2
      WHAT = DSQRT(S)
      Q=ZSS*DSQRT(S) + ZQQ*QQ
      SIGPART=DJAC*DLUMGG(TAU,X,Q)*SIGHAT(Y(N))
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION SIG(XX)
C--FUNCTION FOR GG -> HH: BORN TERM
      PARAMETER(N=3)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DOUBLE PRECISION M1,M2,MZ
      DIMENSION XX(N),Y(N)
      COMPLEX*16 A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/SCALE/ZSS,ZPT,ZQQ,QQ,ISUB
      COMMON/ALSN/LOOP
      COMMON/CUTS/EPS,REPS
      COMMON/PARAM/S,T,U,M1,M2,MZ
      COMMON/PROC/IPROC,IAPP,ITOP,ISFAC
      COMMON/HMASS/AMA,AML,AMH,AMCH
      COMMON/SUSY/TGBET
      COMMON/FORM/A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/TOTPAR/W,GEVPB,GF,GAMZ,GAMA,GAML,GAMH
      PI=4.D0*DATAN(1.D0)
      DO 1 I=1,N
       Y(I)=EPS+(1.D0-2.D0*EPS)*XX(I)
1     CONTINUE
      TH=(M1+M2)**2/W**2
      IF(ISUB.EQ.0)THEN
       TAU=TH+(1.D0-TH)*Y(1)
       X=TAU+(1.D0-TAU)*Y(2)
       DJAC=(1.D0-TH)*(1.D0-TAU)*(1.D0-2.D0*EPS)**N
      ELSE
       YY=-DLOG(TH)*Y(1)
       TAU=DEXP(-YY)
       XY=-DLOG(TAU)*Y(2)
       X=DEXP(-XY)
       DJAC=DLOG(TH)*DLOG(TAU)*TAU*X*(1.D0-2.D0*EPS)**N
      ENDIF
      V=Y(3)/2.D0
      S=TAU*W**2
      XLAM=DSQRT((S+M1**2-M2**2)**2-4.D0*S*M1**2)
      BET=XLAM/(S+M1**2-M2**2)
      T1=-(S+M1**2-M2**2)*((1.D0+BET)/2.D0-BET*V)
      U1=-(S+M1**2-M2**2)*((1.D0-BET)/2.D0+BET*V)
      T=T1+M1**2
      U=U1+M1**2
      PT2=(T1*U1-S*M1**2)/S
      AM2=(M1+M2)**2/4.D0
      Q=ZSS*DSQRT(S) + ZPT*DSQRT(AM2+PT2) + ZQQ*QQ
      IF(IPROC.EQ.0.OR.IPROC.EQ.2.OR.IPROC.EQ.5)DJAC=DJAC/2.D0
      CALL MATRIX(DMAT)
      ALPS=ALPHAS(Q,LOOP)
      QCD=2.D0*ALPS**2/(4.D0*PI)**2
      SIG=GEVPB*QCD*DJAC*DLUMGG(TAU,X,Q)*XLAM/S**2/4096.D0/PI*DMAT
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c    .    * s
c     QCUT = 1500.D0
c     IF(S.GE.QCUT**2) SIG = 0
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION SIGRES(XX)
C--FUNCTION FOR GG -> HH: RESONANCE CONTRIBUTION
      PARAMETER(N=3)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DOUBLE PRECISION M1,M2,MZ
      DIMENSION XX(N),Y(N)
      COMPLEX*16 A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/SCALE/ZSS,ZPT,ZQQ,QQ,ISUB
      COMMON/ALSN/LOOP
      COMMON/CUTS/EPS,REPS
      COMMON/PARAM/S,T,U,M1,M2,MZ
      COMMON/PROC/IPROC,IAPP,ITOP,ISFAC
      COMMON/HMASS/AMA,AML,AMH,AMCH
      COMMON/SUSY/TGBET
      COMMON/FORM/A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/TOTPAR/W,GEVPB,GF,GAMZ,GAMA,GAML,GAMH
      PI=4.D0*DATAN(1.D0)
      IF(AMH.LE.M1+M2.OR.W.LE.AMH)THEN
       SIGRES=0
      ELSE
      DO 1 I=1,N
       Y(I)=EPS+(1.D0-2.D0*EPS)*XX(I)
1     CONTINUE
      TH=(M1+M2)**2/W**2
      IF(ISUB.EQ.0)THEN
       TAU=TH+(1.D0-TH)*Y(1)
       X=TAU+(1.D0-TAU)*Y(2)
       DJAC=(1.D0-TH)*(1.D0-TAU)*(1.D0-2.D0*EPS)**N
      ELSE
       YY=-DLOG(TH)*Y(1)
       TAU=DEXP(-YY)
       XY=-DLOG(TAU)*Y(2)
       X=DEXP(-XY)
       DJAC=DLOG(TH)*DLOG(TAU)*TAU*X*(1.D0-2.D0*EPS)**N
      ENDIF
      V=Y(3)/2.D0
      S=TAU*W**2
      XLAM=DSQRT((AMH**2-M1**2-M2**2)**2-4.D0*M1**2*M2**2)
      BET=XLAM/(AMH**2+M1**2-M2**2)
      T1=-(AMH**2+M1**2-M2**2)*((1.D0+BET)/2.D0-BET*V)
      U1=-(AMH**2+M1**2-M2**2)*((1.D0-BET)/2.D0+BET*V)
      T=T1+M1**2
      U=U1+M1**2
      PT2=(T1*U1-AMH**2*M1**2)/S
      AM2=(M1+M2)**2/4.D0
      Q=ZSS*DSQRT(AMH**2) + ZPT*DSQRT(AM2+PT2) + ZQQ*QQ
      IF(IPROC.EQ.0.OR.IPROC.EQ.2.OR.IPROC.EQ.5)DJAC=DJAC/2.D0
      CALL MATRES(DMAT)
      ALPS=ALPHAS(Q,LOOP)
      QCD=2.D0*ALPS**2/(4.D0*PI)**2
      SIGRES=GEVPB*QCD*DJAC*DLUMGG(TAU,X,Q)*XLAM/AMH**4/4096.D0/PI*DMAT
      ENDIF
      RETURN
      END
  
      SUBROUTINE MATRIX2(DMAT)
C--MATRIX ELEMENT FOR C2 TERM OF VIRTUAL GG -> HH CORRECTIONS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION M1,M2,MZ
      COMPLEX*16 A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMPLEX*16 F1,F2,PROA,PROL,PROH,PROZ,CBOX
      COMPLEX*16 F10,F20,F11,F12,F21,H10,HH10,HH20,CBOX0,CBOX1
      COMMON/PARAM/S,T,U,M1,M2,MZ
      COMMON/PROC/IPROC,IAPP,ITOP,ISFAC
      COMMON/COUP/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,
     .            GHHH,GLLL,GHLL,GLHH,GHAA,GLAA,GLVV,GHVV,
     .            GLPM,GHPM,GZT,GZB
      COMMON/HMASS/AMA,AML,AMH,AMCH
      COMMON/QMASS/AMB,AMT
      COMMON/FORM/A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/TOTPAR/W,GEVPB,GF,GAMZ,GAMA,GAML,GAMH
      COMMON/TREE/ITREE,ITRIA,IWRITE
      COMMON/TRILINEAR/FACTRIA,FACT,FACB,FACTT,FACG,FACGG,ISILH
      FACA=1
      FACL=1
      FACH=1
      FACZ=1
      IF(ISFAC.EQ.1)THEN
       FACA=S/AMA**2
       FACL=S/AML**2
       FACH=S/AMH**2
       FACZ=S/MZ**2
      ENDIF
      PROL = DCMPLX(S-AML**2,AML*GAML*FACL)
      PROH = DCMPLX(S-AMH**2,AMH*GAMH*FACH)
      PROA = DCMPLX(S-AMA**2,AMA*GAMA*FACA)
      PROZ = DCMPLX(S- MZ**2, MZ*GAMZ*FACZ)
      IF(IPROC.EQ.0)THEN
       IF(ITRIA.EQ.3) GLLL = 0
       IF(ITRIA.EQ.4) GHLL = 0
       CALL FORMFAC(AMT)
       F1 = H1*(GLT*GLLL/PROL + GHT*GHLL/PROH) + HH1*GLT**2
       F2 = H2*(GLT*GLLL/PROL + GHT*GHLL/PROH) + HH2*GLT**2
       CALL FORMFAC(AMB)
       F1 = F1 + H1*(GLB*GLLL/PROL + GHB*GHLL/PROH) + HH1*GLB**2
       F2 = F2 + H2*(GLB*GLLL/PROL + GHB*GHLL/PROH) + HH2*GLB**2
       CBOX = GLT**2
      ELSEIF(IPROC.EQ.1)THEN
       IF(ITRIA.EQ.3) GHLL = 0
       IF(ITRIA.EQ.4) GLHH = 0
       CALL FORMFAC(AMT)
       F1 = H1*(GLT*GHLL/PROL + GHT*GLHH/PROH) + HH1*GLT*GHT
       F2 = H2*(GLT*GHLL/PROL + GHT*GLHH/PROH) + HH2*GLT*GHT
       CALL FORMFAC(AMB)
       F1 = F1 + H1*(GLB*GHLL/PROL + GHB*GLHH/PROH) + HH1*GLB*GHB
       F2 = F2 + H2*(GLB*GHLL/PROL + GHB*GLHH/PROH) + HH2*GLB*GHB
       CBOX = GLT*GHT
      ELSEIF(IPROC.EQ.2)THEN
       IF(ITRIA.EQ.3) GLHH = 0
       IF(ITRIA.EQ.4) GHHH = 0
       CALL FORMFAC(AMT)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       v=1.d0/dsqrt(dsqrt(2.d0)*gf)
       h10  = h1
       hh10 = hh1
       hh20 = hh2
       hh1 = fact**2*hh1+factt*h1/amt+8*facgg*2*s/amt**2
       h1  = fact*h1+8*facg*2*s/amt
       hh2 = fact**2*hh2
       d3  = factria-1
       dt  = fact-1
       dtt = factt/amt
       dg  = facg*2*s/amt
       dgg = facgg*2*s/amt**2
       ghhh0 = 3*amh**2/v
       f10 = h10*ght*ghhh0/proh + hh10*ght*ght
       f20 = hh20*ght*ght
       f11 = ght*ghhh0/proh*(h10*(d3+dt)+8*dg)
     .     + ght*ght*(2*dt*hh10+dtt*h10+8*dgg)
       f21 = 2*dt*f20
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       F1 = H1*(GLT*GLHH/PROL + GHT*GHHH/PROH) + HH1*GHT*GHT
       F2 = H2*(GLT*GLHH/PROL + GHT*GHHH/PROH) + HH2*GHT*GHT
       CBOX = GHT*GHT
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       cbox0 = cbox
       cbox1 = (2*(dt+12*facg)-99*facgg-99*facg*ghhh0/proh*v)*cbox0
       cbox = ((fact+12*facg)**2-99*facgg-99*facg*ghhh/proh*v)*cbox
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       CALL FORMFAC(AMB)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       XHB = FACB*GHB
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       F1 = F1 + H1*(GLB*GLHH/PROL + XHB*GHHH/PROH) + HH1*XHB*XHB
       F2 = F2 + H2*(GLB*GLHH/PROL + XHB*GHHH/PROH) + HH2*XHB*XHB
      ELSEIF(IPROC.EQ.3)THEN
       IF(ITRIA.EQ.3) GZAL = 0
       IF(ITRIA.EQ.4) GLAA = 0
       CALL FORMFAC(AMT)
       F1 = A1*GAT*GLAA/PROA + Z1*GZT*GZAL/PROZ + AH1*GAT*GLT
       F2 = A2*GAT*GLAA/PROA + Z2*GZT*GZAL/PROZ + AH2*GAT*GLT
       CALL FORMFAC(AMB)
       F1 = F1 + A1*GAB*GLAA/PROA + Z1*GZB*GZAL/PROZ + AH1*GAB*GLB
       F2 = F2 + A2*GAB*GLAA/PROA + Z2*GZB*GZAL/PROZ + AH2*GAB*GLB
       CBOX = GAT*GLT
      ELSEIF(IPROC.EQ.4)THEN
       IF(ITRIA.EQ.3) GZAH = 0
       IF(ITRIA.EQ.4) GHAA = 0
       CALL FORMFAC(AMT)
       F1 = A1*GAT*GHAA/PROA + Z1*GZT*GZAH/PROZ + AH1*GAT*GHT
       F2 = A2*GAT*GHAA/PROA + Z2*GZT*GZAH/PROZ + AH2*GAT*GHT
       CALL FORMFAC(AMB)
       F1 = F1 + A1*GAB*GHAA/PROA + Z1*GZB*GZAH/PROZ + AH1*GAB*GHB
       F2 = F2 + A2*GAB*GHAA/PROA + Z2*GZB*GZAH/PROZ + AH2*GAB*GHB
       CBOX = GAT*GHT
      ELSEIF(IPROC.EQ.5)THEN
       IF(ITRIA.EQ.3) GLAA = 0
       IF(ITRIA.EQ.4) GHAA = 0
       CALL FORMFAC(AMT)
       F1 = H1*(GLT*GLAA/PROL + GHT*GHAA/PROH) + AA1*GAT*GAT
       F2 = H2*(GLT*GLAA/PROL + GHT*GHAA/PROH) + AA2*GAT*GAT
       CALL FORMFAC(AMB)
       F1 = F1 + H1*(GLB*GLAA/PROL + GHB*GHAA/PROH) + AA1*GAB*GAB
       F2 = F2 + H2*(GLB*GLAA/PROL + GHB*GHAA/PROH) + AA2*GAB*GAB
       CBOX = GAT*GAT
      ENDIF
      DMAT = 2.D0*DREAL(F1*CBOX) * 2*S/AMT**2
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      if(isilh.ne.0)then
       delta = 2*dreal((dconjg(f10)*cbox1+dconjg(f11)*cbox0))
       dmat0 = 2*dreal(f10*cbox0)
       dmat = (dmat0 + delta) * 2*s/amt**2
      endif
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      RETURN
      END

      SUBROUTINE MATRIX3(DMAT)
C--MATRIX ELEMENT FOR C3 TERM OF VIRTUAL GG -> HH CORRECTIONS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION M1,M2,MZ
      COMPLEX*16 A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMPLEX*16 F1,F2,PROA,PROL,PROH,PROZ,CBOX
      COMPLEX*16 F10,F20,F11,F12,F21,H10,HH10,HH20,CBOX0,CBOX1
      COMMON/PARAM/S,T,U,M1,M2,MZ
      COMMON/PROC/IPROC,IAPP,ITOP,ISFAC
      COMMON/COUP/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,
     .            GHHH,GLLL,GHLL,GLHH,GHAA,GLAA,GLVV,GHVV,
     .            GLPM,GHPM,GZT,GZB
      COMMON/HMASS/AMA,AML,AMH,AMCH
      COMMON/QMASS/AMB,AMT
      COMMON/FORM/A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/TOTPAR/W,GEVPB,GF,GAMZ,GAMA,GAML,GAMH
      COMMON/TREE/ITREE,ITRIA,IWRITE
      COMMON/VIRTCOF/VC1,VC2,VC3,IHA
      COMMON/TRILINEAR/FACTRIA,FACT,FACB,FACTT,FACG,FACGG,ISILH
      FACA=1
      FACL=1
      FACH=1
      FACZ=1
      IF(ISFAC.EQ.1)THEN
       FACA=S/AMA**2
       FACL=S/AML**2
       FACH=S/AMH**2
       FACZ=S/MZ**2
      ENDIF
      PROL = DCMPLX(S-AML**2,AML*GAML*FACL)
      PROH = DCMPLX(S-AMH**2,AMH*GAMH*FACH)
      PROA = DCMPLX(S-AMA**2,AMA*GAMA*FACA)
      PROZ = DCMPLX(S- MZ**2, MZ*GAMZ*FACZ)
      T1=T-M1**2
      U1=U-M1**2
      PT2=(T1*U1-S*M1**2)/S
      IF(IPROC.EQ.0)THEN
       IF(ITRIA.EQ.3) GLLL = 0
       IF(ITRIA.EQ.4) GHLL = 0
       CALL FORMFAC(AMT)
       F1 = H1*(GLT*GLLL/PROL + GHT*GHLL/PROH) + HH1*GLT**2
       F2 = H2*(GLT*GLLL/PROL + GHT*GHLL/PROH) + HH2*GLT**2
       CALL FORMFAC(AMB)
       F1 = F1 + H1*(GLB*GLLL/PROL + GHB*GHLL/PROH) + HH1*GLB**2
       F2 = F2 + H2*(GLB*GLLL/PROL + GHB*GHLL/PROH) + HH2*GLB**2
       CBOX = GLT*GLT
      ELSEIF(IPROC.EQ.1)THEN
       IF(ITRIA.EQ.3) GHLL = 0
       IF(ITRIA.EQ.4) GLHH = 0
       CALL FORMFAC(AMT)
       F1 = H1*(GLT*GHLL/PROL + GHT*GLHH/PROH) + HH1*GLT*GHT
       F2 = H2*(GLT*GHLL/PROL + GHT*GLHH/PROH) + HH2*GLT*GHT
       CALL FORMFAC(AMB)
       F1 = F1 + H1*(GLB*GHLL/PROL + GHB*GLHH/PROH) + HH1*GLB*GHB
       F2 = F2 + H2*(GLB*GHLL/PROL + GHB*GLHH/PROH) + HH2*GLB*GHB
       CBOX = GHT*GLT
      ELSEIF(IPROC.EQ.2)THEN
       IF(ITRIA.EQ.3) GLHH = 0
       IF(ITRIA.EQ.4) GHHH = 0
       CALL FORMFAC(AMT)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       v=1.d0/dsqrt(dsqrt(2.d0)*gf)
       hh10 = hh1
       hh20 = hh2
       h10 = h1
       hh1 = fact**2*hh1+factt*h1/amt+8*facgg*2*s/amt**2
       h1  = fact*h1+8*facg*2*s/amt
       hh2 = fact**2*hh2
       d3  = factria-1
       dt  = fact-1
       dtt = factt/amt
       dg  = facg*2*s/amt
       dgg = facgg*2*s/amt**2
       ghhh0 =3*amh**2/v
       f10 = h10*ght*ghhh0/proh + hh10*ght*ght
       f20 = hh20*ght*ght
       f11 = ght*ghhh0/proh*(h10*(d3+dt)+8*dg)
     .     + ght*ght*(2*dt*hh10+dtt*h10+8*dgg)
       f21 = 2*dt*hh20*ght*ght
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       F1 = H1*(GLT*GLHH/PROL + GHT*GHHH/PROH) + HH1*GHT*GHT
       F2 = H2*(GLT*GLHH/PROL + GHT*GHHH/PROH) + HH2*GHT*GHT
       CBOX = GHT*GHT
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       cbox0 = cbox
       cbox = (fact+12*facg)**2*cbox
       cbox1 = 2*(dt+12*facg)*cbox0
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c      F2 = AH2
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       CALL FORMFAC(AMB)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       XHB = FACB*GHB
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       F1 = F1 + H1*(GLB*GLHH/PROL + XHB*GHHH/PROH) + HH1*XHB*XHB
       F2 = F2 + H2*(GLB*GLHH/PROL + XHB*GHHH/PROH) + HH2*XHB*XHB
      ELSEIF(IPROC.EQ.3)THEN
       IF(ITRIA.EQ.3) GZAL = 0
       IF(ITRIA.EQ.4) GLAA = 0
       CALL FORMFAC(AMT)
       F1 = A1*GAT*GLAA/PROA + Z1*GZT*GZAL/PROZ + AH1*GAT*GLT
       F2 = A2*GAT*GLAA/PROA + Z2*GZT*GZAL/PROZ + AH2*GAT*GLT
       CALL FORMFAC(AMB)
       F1 = F1 + A1*GAB*GLAA/PROA + Z1*GZB*GZAL/PROZ + AH1*GAB*GLB
       F2 = F2 + A2*GAB*GLAA/PROA + Z2*GZB*GZAL/PROZ + AH2*GAB*GLB
       CBOX = GAT*GLT
      ELSEIF(IPROC.EQ.4)THEN
       IF(ITRIA.EQ.3) GZAH = 0
       IF(ITRIA.EQ.4) GHAA = 0
       CALL FORMFAC(AMT)
       F1 = A1*GAT*GHAA/PROA + Z1*GZT*GZAH/PROZ + AH1*GAT*GHT
       F2 = A2*GAT*GHAA/PROA + Z2*GZT*GZAH/PROZ + AH2*GAT*GHT
       CALL FORMFAC(AMB)
       F1 = F1 + A1*GAB*GHAA/PROA + Z1*GZB*GZAH/PROZ + AH1*GAB*GHB
       F2 = F2 + A2*GAB*GHAA/PROA + Z2*GZB*GZAH/PROZ + AH2*GAB*GHB
       CBOX = GAT*GHT
      ELSEIF(IPROC.EQ.5)THEN
       IF(ITRIA.EQ.3) GLAA = 0
       IF(ITRIA.EQ.4) GHAA = 0
       CALL FORMFAC(AMT)
       F1 = H1*(GLT*GLAA/PROL + GHT*GHAA/PROH) + AA1*GAT*GAT
       F2 = H2*(GLT*GLAA/PROL + GHT*GHAA/PROH) + AA2*GAT*GAT
       CALL FORMFAC(AMB)
       F1 = F1 + H1*(GLB*GLAA/PROL + GHB*GHAA/PROH) + AA1*GAB*GAB
       F2 = F2 + H2*(GLB*GLAA/PROL + GHB*GHAA/PROH) + AA2*GAB*GAB
       CBOX = GAT*GAT
      ENDIF
      DMAT = 2.D0*DREAL(F2*CBOX) * 2*S/AMT**2
     .       *PT2/2/T/U
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      if(isilh.ne.0)then
       delta = 2*dreal(dconjg(f20)*(cbox1+2*dt*cbox0))
       dmat0 = 2*dreal(f20*cbox0)
       dmat = (dmat0 + delta) * 2*s/amt**2
     .         *pt2/2/t/u
      endif
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF(IHA.EQ.0)THEN
       DMAT = DMAT*(S-M1**2-M2**2)
      ELSE
       DMAT = DMAT*(T-U)
      ENDIF
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     DMAT = 2.D0*DREAL(F2) * 2*S/AMT**2
c    .       *PT2/2/T/U
c     DMAT = DMAT*(T-U)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      RETURN
      END

      SUBROUTINE MATRIX(DMAT)
C--LO MATRIX ELEMENT FOR GG -> HH
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION M1,M2,MZ
      COMPLEX*16 A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMPLEX*16 F1,F2,PROA,PROL,PROH,PROZ
      COMPLEX*16 F10,F20,F11,F12,F21,H10,HH10,HH20
      COMMON/PARAM/S,T,U,M1,M2,MZ
      COMMON/PROC/IPROC,IAPP,ITOP,ISFAC
      COMMON/COUP/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,
     .            GHHH,GLLL,GHLL,GLHH,GHAA,GLAA,GLVV,GHVV,
     .            GLPM,GHPM,GZT,GZB
      COMMON/HMASS/AMA,AML,AMH,AMCH
      COMMON/QMASS/AMB,AMT
      COMMON/FORM/A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/TOTPAR/W,GEVPB,GF,GAMZ,GAMA,GAML,GAMH
      COMMON/TREE/ITREE,ITRIA,IWRITE
      COMMON/TRILINEAR/FACTRIA,FACT,FACB,FACTT,FACG,FACGG,ISILH
      FACA=1
      FACL=1
      FACH=1
      FACZ=1
      IF(ISFAC.EQ.1)THEN
       FACA=S/AMA**2
       FACL=S/AML**2
       FACH=S/AMH**2
       FACZ=S/MZ**2
      ENDIF
      PROL = DCMPLX(S-AML**2,AML*GAML*FACL)
      PROH = DCMPLX(S-AMH**2,AMH*GAMH*FACH)
      PROA = DCMPLX(S-AMA**2,AMA*GAMA*FACA)
      PROZ = DCMPLX(S- MZ**2, MZ*GAMZ*FACZ)
      IF(IPROC.EQ.0)THEN
       IF(ITRIA.EQ.3) GLLL = 0
       IF(ITRIA.EQ.4) GHLL = 0
       CALL FORMFAC(AMT)
       F1 = H1*(GLT*GLLL/PROL + GHT*GHLL/PROH) + HH1*GLT**2
       F2 = H2*(GLT*GLLL/PROL + GHT*GHLL/PROH) + HH2*GLT**2
       CALL FORMFAC(AMB)
       F1 = F1 + H1*(GLB*GLLL/PROL + GHB*GHLL/PROH) + HH1*GLB**2
       F2 = F2 + H2*(GLB*GLLL/PROL + GHB*GHLL/PROH) + HH2*GLB**2
      ELSEIF(IPROC.EQ.1)THEN
       IF(ITRIA.EQ.3) GHLL = 0
       IF(ITRIA.EQ.4) GLHH = 0
       CALL FORMFAC(AMT)
       F1 = H1*(GLT*GHLL/PROL + GHT*GLHH/PROH) + HH1*GLT*GHT
       F2 = H2*(GLT*GHLL/PROL + GHT*GLHH/PROH) + HH2*GLT*GHT
       CALL FORMFAC(AMB)
       F1 = F1 + H1*(GLB*GHLL/PROL + GHB*GLHH/PROH) + HH1*GLB*GHB
       F2 = F2 + H2*(GLB*GHLL/PROL + GHB*GLHH/PROH) + HH2*GLB*GHB
      ELSEIF(IPROC.EQ.2)THEN
       IF(ITRIA.EQ.3) GLHH = 0
       IF(ITRIA.EQ.4) GHHH = 0
       CALL FORMFAC(AMT)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       v=1.d0/dsqrt(dsqrt(2.d0)*gf)
       ghhh0 =3*amh**2/v
       hh10 = hh1
       hh20 = hh2
       h10 = h1
       hh1 = fact**2*hh1+factt*h1/amt+8*facgg*2*s/amt**2
       h1  = fact*h1+8*facg*2*s/amt
       hh2 = fact**2*hh2
       d3  = factria-1
       dt  = fact-1
       dtt = factt/amt
       dg  = facg*2*s/amt
       dgg = facgg*2*s/amt**2
       f10 = h10*ght*ghhh0/proh + hh10*ght*ght
       f20 = hh20*ght*ght
       f11 = ght*ghhh0/proh*(h10*(d3+dt)+8*dg)
     .     + ght*ght*(2*dt*hh10+dtt*h10+8*dgg)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       F1 = H1*(GLT*GLHH/PROL + GHT*GHHH/PROH) + HH1*GHT*GHT
       F2 = H2*(GLT*GLHH/PROL + GHT*GHHH/PROH) + HH2*GHT*GHT
       CALL FORMFAC(AMB)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       XHB = FACB*GHB
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       F1 = F1 + H1*(GLB*GLHH/PROL + XHB*GHHH/PROH) + HH1*XHB*XHB
       F2 = F2 + H2*(GLB*GLHH/PROL + XHB*GHHH/PROH) + HH2*XHB*XHB
      ELSEIF(IPROC.EQ.3)THEN
       IF(ITRIA.EQ.3) GZAL = 0
       IF(ITRIA.EQ.4) GLAA = 0
       CALL FORMFAC(AMT)
       F1 = A1*GAT*GLAA/PROA + Z1*GZT*GZAL/PROZ + AH1*GAT*GLT
       F2 = A2*GAT*GLAA/PROA + Z2*GZT*GZAL/PROZ + AH2*GAT*GLT
       CALL FORMFAC(AMB)
       F1 = F1 + A1*GAB*GLAA/PROA + Z1*GZB*GZAL/PROZ + AH1*GAB*GLB
       F2 = F2 + A2*GAB*GLAA/PROA + Z2*GZB*GZAL/PROZ + AH2*GAB*GLB
      ELSEIF(IPROC.EQ.4)THEN
       IF(ITRIA.EQ.3) GZAH = 0
       IF(ITRIA.EQ.4) GHAA = 0
       CALL FORMFAC(AMT)
       F1 = A1*GAT*GHAA/PROA + Z1*GZT*GZAH/PROZ + AH1*GAT*GHT
       F2 = A2*GAT*GHAA/PROA + Z2*GZT*GZAH/PROZ + AH2*GAT*GHT
       CALL FORMFAC(AMB)
       F1 = F1 + A1*GAB*GHAA/PROA + Z1*GZB*GZAH/PROZ + AH1*GAB*GHB
       F2 = F2 + A2*GAB*GHAA/PROA + Z2*GZB*GZAH/PROZ + AH2*GAB*GHB
      ELSEIF(IPROC.EQ.5)THEN
       IF(ITRIA.EQ.3) GLAA = 0
       IF(ITRIA.EQ.4) GHAA = 0
       CALL FORMFAC(AMT)
       F1 = H1*(GLT*GLAA/PROL + GHT*GHAA/PROH) + AA1*GAT*GAT
       F2 = H2*(GLT*GLAA/PROL + GHT*GHAA/PROH) + AA2*GAT*GAT
       CALL FORMFAC(AMB)
       F1 = F1 + H1*(GLB*GLAA/PROL + GHB*GHAA/PROH) + AA1*GAB*GAB
       F2 = F2 + H2*(GLB*GLAA/PROL + GHB*GHAA/PROH) + AA2*GAB*GAB
      ENDIF
      DMAT = 2.D0*(CDABS(F1)**2+CDABS(F2)**2)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      if(isilh.ne.0)then
       delta = 2*(2*dreal(dconjg(f10)*f11) + 4*dt*cdabs(f20)**2)
       dmat0 = 2*(cdabs(f10)**2+cdabs(f20)**2)
       dmat = dmat0 + delta
      endif
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      RETURN
      END

      SUBROUTINE MATRES(DMAT)
C--MATRIX ELEMENT FOR GG -> HH: RESONANCE CONTRIBUTION
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION M1,M2,MZ
      COMPLEX*16 A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMPLEX*16 F1,F2,PROA,PROL,PROH,PROZ
      COMMON/PARAM/S,T,U,M1,M2,MZ
      COMMON/PROC/IPROC,IAPP,ITOP,ISFAC
      COMMON/COUP/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,
     .            GHHH,GLLL,GHLL,GLHH,GHAA,GLAA,GLVV,GHVV,
     .            GLPM,GHPM,GZT,GZB
      COMMON/HMASS/AMA,AML,AMH,AMCH
      COMMON/QMASS/AMB,AMT
      COMMON/FORM/A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/TOTPAR/W,GEVPB,GF,GAMZ,GAMA,GAML,GAMH
      FACA=1
      FACL=1
      FACH=1
      FACZ=1
      IF(ISFAC.EQ.1)THEN
       FACA=S/AMA**2
       FACL=S/AML**2
       FACH=S/AMH**2
       FACZ=S/MZ**2
      ENDIF
      PROL = DCMPLX(S-AML**2,AML*GAML*FACL)
      PROH = DCMPLX(S-AMH**2,AMH*GAMH*FACH)
      PROA = DCMPLX(S-AMA**2,AMA*GAMA*FACA)
      PROZ = DCMPLX(S- MZ**2, MZ*GAMZ*FACZ)
      IF(IPROC.EQ.0)THEN
       CALL FORMRES(AMT)
       F1 = H1*GHT*GHLL/PROH
       F2 = H2*GHT*GHLL/PROH
       CALL FORMRES(AMB)
       F1 = F1 + H1*GHB*GHLL/PROH
       F2 = F2 + H2*GHB*GHLL/PROH
      ELSEIF(IPROC.EQ.5)THEN
       CALL FORMRES(AMT)
       F1 = H1*GHT*GHAA/PROH
       F2 = H2*GHT*GHAA/PROH
       CALL FORMRES(AMB)
       F1 = F1 + H1*GHB*GHAA/PROH
       F2 = F2 + H2*GHB*GHAA/PROH
      ENDIF
      DMAT = 2.D0*(CDABS(F1)**2+CDABS(F2)**2)
      RETURN
      END

      SUBROUTINE FORMFAC(AMQ)
C--FORM FACTORS FOR LO MATRIX ELEMENTS FOR GG -> HH
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION M1,M2,MZ
      COMPLEX*16 A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMPLEX*16 C0AB,C0AC,C0AD,C0BC,C0BD,C0CD,D0ABC,D0BAC,D0ACB
      COMPLEX*16 CQ2,CA5
      COMMON/PARAM/S,T,U,M1,M2,MZ
      COMMON/SCALINT/C0AB,C0AC,C0AD,C0BC,C0BD,C0CD,D0ABC,D0BAC,D0ACB
      COMMON/FORM/A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/PROC/IPROC,IAPP,ITOP,ISFAC
      COMMON/TREE/ITREE,ITRIA,IWRITE
      COMMON/CUTS/EPS,REPS
      COMMON/TRILINEAR/FACTRIA,FACT,FACB,FACTT,FACG,FACGG,ISILH
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     factor = 1.d0
c--------------------------------------------------------------------
c     rhs = 1.d0   / factor
c     rht =-1/3.d0 / factor
c     rhu =-2/3.d0 / factor
c     rhh = 0.1d0  / factor
c--------------------------------------------------------------------
c     rhs = 10.d0   / factor
c     rht =-10/3.d0 / factor
c     rhu =-20/3.d0 / factor
c     rhh = 1.d0  / factor
c--------------------------------------------------------------------
c     m1 = dsqrt(rhh) * amq
c     m2 = dsqrt(rhh) * amq
c     s = rhs * amq**2
c     t = rht * amq**2 + m1**2
c     u = rhu * amq**2 + m1**2
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      SS=S/AMQ**2
      TT=T/AMQ**2
      UU=U/AMQ**2
      R1=M1**2/AMQ**2
      R2=M2**2/AMQ**2
      RZ=MZ**2/AMQ**2
      T1=TT-R1
      U1=UU-R1
      IF(IAPP.EQ.0)THEN

      IF(ITRIA.LE.2) THEN
       CALL INISCAL(AMQ)
      ELSE
       DQ2=AMQ**2
       CQ2=AMQ**2*DCMPLX(1.D0,-REPS)
       CA5 = CDSQRT(1.D0-4.D0*CQ2/S)
       C0AB = 0.5D0*CDLOG((CA5+1.D0)/(CA5-1.D0))**2/S*DQ2
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c      C0AB = 0.5D0*CDLOG((CA5+1.D0)/(CA5-1.D0))**2/S*CQ2
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ENDIF
 
      A1 =  - 4*C0AB*SS
 
      A2 = 0
 
      H1 =  - 4*(C0AB*SS - 4*C0AB - 2)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     H1 =  - 4*(C0AB*SS*DQ2/CQ2 - 4*C0AB - 2)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 
      H2 = 0
 
      Z1 = ( - 4*(RZ - SS)*(2*C0AB + 1)*(U1 + T1 + SS))/RZ
   
      Z2 = 0
 
      AA1 = 0
      AA2 = 0
      HH1 = 0
      HH2 = 0
      AH1 = 0
      AH2 = 0

      IF(ITRIA.LE.2) THEN

      AA1 = (16*C0AB*SS + 2*C0AC*T1*(U1 + T1 + SS + 2*R1) + 2*C0AD*(
     . - U1*T1 - U1*SS - T1**2 - 2*T1*SS - 2*T1*R1 - SS**2
     . - 2*SS*R1) + 2*C0BC*U1*(U1 + T1 + SS + 2*R1) + 2*C0BD
     .*( - U1**2 - U1*T1 - 2*U1*SS - 2*U1*R1 - T1*SS - SS**2 -
     .2*SS*R1) + 4*D0BAC*SS*( - U1 - T1 - 2*R1) + 2*D0ACB*(
     . - U1**2*T1 - U1*T1**2 - U1*T1*SS - 2*U1*T1*R1 + U1*SS
     .*R1 - 2*U1*SS + T1*SS*R1 - 2*T1*SS + SS**2*R1 + 2*SS*
     .R1**2 - 4*SS*R1) + 4*D0ABC*SS*( - U1 - T1 - 2*R1) + 8*
     .SS)/SS
 
      AA2 = (2*C0AB*SS*(U1**2 + 4*U1*R1 + T1**2 + 4*T1*R1 + 2*SS*R1 +
     .4*R1**2) + 2*C0AC*T1*(U1*R1 + T1**2 + 3*T1*R1 + SS*R1
     . + 2*R1**2) + 2*C0CD*( - U1**3 - U1**2*T1 - 2*U1**2*R1
     . - U1*T1**2 + 2*U1*SS*R1 - T1**3 - 2*T1**2*R1 + 2*T1*
     .SS*R1 + 4*SS*R1**2) + 2*C0AD*( - U1**2*T1 - U1**2*SS
     . - 3*U1*T1*R1 - 3*U1*SS*R1 - T1**2*R1 - 2*T1*SS*R1 - 2
     .*T1*R1**2 - SS**2*R1 - 2*SS*R1**2) + 2*C0BC*U1*(U1**2
     . + 3*U1*R1 + T1*R1 + SS*R1 + 2*R1**2) + 2*C0BD*( - U1
     .**2*R1 - U1*T1**2 - 3*U1*T1*R1 - 2*U1*SS*R1 - 2*U1*R1
     .**2 - T1**2*SS - 3*T1*SS*R1 - SS**2*R1 - 2*SS*R1**2))
      AA2 = (AA2
     . + 2*D0BAC*( - 2*U1**2*T1 - 2*U1*T1**2 - U1*T1*SS*R1 - 4*
     .U1*T1*R1 - U1*SS*R1**2 + 2*U1*SS*R1 - T1**3*SS - 4*T1
     .**2*SS*R1 - T1*SS**2*R1 - 5*T1*SS*R1**2 + 2*T1*SS*R1
     . - SS**2*R1**2 - 2*SS*R1**3 + 4*SS*R1**2) + 4*D0ACB*(
     . - U1**2*T1 - U1*T1**2 - 2*U1*T1*R1 + U1*SS*R1 + T1*SS
     .*R1 + 2*SS*R1**2) + 2*D0ABC*( - U1**3*SS - 2*U1**2*T1
     . - 4*U1**2*SS*R1 - 2*U1*T1**2 - U1*T1*SS*R1 - 4*U1*T1*
     .R1 - U1*SS**2*R1 - 5*U1*SS*R1**2 + 2*U1*SS*R1 - T1*SS*
     .R1**2 + 2*T1*SS*R1 - SS**2*R1**2 - 2*SS*R1**3 + 4*SS*
     .R1**2))/(U1*T1 - SS*R1)
 
      HH1 = (16*C0AB*SS + 2*C0AC*T1*(U1 + T1 + SS + 2*R1 - 8) + 2*C0AD
     .*( - U1*T1 - U1*SS - T1**2 - 2*T1*SS - 2*T1*R1 + 8*T1 -
     .SS**2 - 2*SS*R1 + 8*SS) + 2*C0BC*U1*(U1 + T1 + SS + 2*
     .R1 - 8) + 2*C0BD*( - U1**2 - U1*T1 - 2*U1*SS - 2*U1*R1
     . + 8*U1 - T1*SS - SS**2 - 2*SS*R1 + 8*SS) + 4*D0BAC*SS
     .*( - U1 - T1 - 2*SS - 2*R1 + 8) + 2*D0ACB*( - U1**2*T1 -
     .U1*T1**2 - U1*T1*SS - 2*U1*T1*R1 + 8*U1*T1 + U1*SS*R1
     . - 2*U1*SS + T1*SS*R1 - 2*T1*SS + SS**2*R1 - 4*SS**2
     . + 2*SS*R1**2 - 12*SS*R1 + 16*SS) + 4*D0ABC*SS*( - U1
     . - T1 - 2*SS - 2*R1 + 8) + 8*SS)/SS
 
      HH2 = (2*C0AB*SS*(U1**2 + 4*U1*R1 - 8*U1 + T1**2 + 4*T1*R1 - 8*
     .T1 + 2*SS*R1 + 4*R1**2 - 16*R1) + 2*C0AC*T1*(U1*R1 +
     .T1**2 + 3*T1*R1 - 8*T1 + SS*R1 + 2*R1**2 - 8*R1) + 2*
     .C0CD*( - U1**3 - U1**2*T1 - 2*U1**2*R1 + 8*U1**2 - U1*T1
     .**2 + 2*U1*SS*R1 - T1**3 - 2*T1**2*R1 + 8*T1**2 + 2
     .*T1*SS*R1 + 4*SS*R1**2 - 16*SS*R1) + 2*C0AD*( - U1
     .**2*T1 - U1**2*SS - 3*U1*T1*R1 + 8*U1*T1 - 3*U1*SS*R1
     . + 8*U1*SS - T1**2*R1 - 2*T1*SS*R1 - 2*T1*R1**2 + 8*T1
     .*R1 - SS**2*R1 - 2*SS*R1**2 + 8*SS*R1) + 2*C0BC*U1*(U1
     .**2 + 3*U1*R1 - 8*U1 + T1*R1 + SS*R1 + 2*R1**2 - 8*R1))
      HH2 = (HH2
     . + 2*C0BD*( - U1**2*R1 - U1*T1**2 - 3*U1*T1*R1 + 8*U1*T1
     . - 2*U1*SS*R1 - 2*U1*R1**2 + 8*U1*R1 - T1**2*SS - 3*T1
     .*SS*R1 + 8*T1*SS - SS**2*R1 - 2*SS*R1**2 + 8*SS*R1) +
     .2*D0BAC*( - 2*U1**2*T1 - 2*U1*T1**2 - U1*T1*SS*R1 - 4*U1*
     .T1*R1 + 16*U1*T1 - U1*SS*R1**2 + 2*U1*SS*R1 - T1
     .**3*SS - 4*T1**2*SS*R1 + 8*T1**2*SS - T1*SS**2*
     .R1 - 5*T1*SS*R1**2 + 18*T1*SS*R1 - SS**2*R1**2
     . - 2*SS*R1**3 + 12*SS*R1**2 - 16*SS*R1))
      HH2 = ( HH2 + 4*
     .D0ACB*( - U1**2*T1 - U1*T1**2 - 2*U1*T1*R1 + 8*U1*T1 + U1
     .*SS*R1 + T1*SS*R1 + 2*SS*R1**2 - 8*SS*R1) + 2*
     .D0ABC*( - U1**3*SS - 2*U1**2*T1 - 4*U1**2*SS*R1 + 8*U1**2
     .*SS - 2*U1*T1**2 - U1*T1*SS*R1 - 4*U1*T1*R1 + 16*
     .U1*T1 - U1*SS**2*R1 - 5*U1*SS*R1**2 + 18*U1*SS*R1
     . - T1*SS*R1**2 + 2*T1*SS*R1 - SS**2*R1**2 - 2*SS*
     .R1**3 + 12*SS*R1**2 - 16*SS*R1))/(U1*T1 - SS*R1)
 
      AH1 = (2*C0AC*T1*(U1 + T1 + SS) + 2*C0AD*( - U1*T1 - U1*SS - T1
     .**2 - 2*T1*SS - SS**2) + 2*C0BC*U1*(U1 + T1 + SS) + 2*
     .C0BD*( - U1**2 - U1*T1 - 2*U1*SS - T1*SS - SS**2) + 4*
     .D0BAC*SS*( - U1 - T1 - 2*SS) + 2*D0ACB*( - U1**2*T1 - U1*
     .T1**2 - U1*T1*SS + U1*SS*R1 - 2*U1*SS + T1*SS*R1 - 2*
     .T1*SS + SS**2*R1 - 4*SS**2) + 4*D0ABC*SS*( - U1 - T1
     . - 2*SS))/SS
 
      AH2 = (2*C0AB*SS*( - U1**2 - 2*U1*R1 + T1**2 + 2*T1*R1) + 2*C0AC
     .*T1*( - U1*R1 + T1**2 + T1*R1 - SS*R1) + 2*C0CD*(U1**3 +
     .U1**2*T1 - U1*T1**2 - 4*U1*SS*R1 - T1**3 + 4*T1*SS*R1)
     . + 2*C0AD*(U1**2*T1 + U1**2*SS + U1*T1*R1 + U1*SS*R1 - T1
     .**2*R1 - 2*T1*SS*R1 - SS**2*R1) + 2*C0BC*U1*( - U1**2
     . - U1*R1 + T1*R1 + SS*R1) + 2*C0BD*(U1**2*R1 - U1*T1**
     .2 - U1*T1*R1 + 2*U1*SS*R1 - T1**2*SS - T1*SS*R1 + SS**
     .2*R1) + 2*D0BAC*(2*U1**2*T1 - 2*U1*T1**2 + U1*T1*SS*R1
     . + U1*SS*R1**2 - 2*U1*SS*R1 - T1**3*SS - 2*T1**2*SS*R1
     . + T1*SS**2*R1 - T1*SS*R1**2 + 2*T1*SS*R1 + SS**2*R1**2))
      AH2 = (AH2
     . + 4*D0ACB*(U1**2*T1 - U1*T1**2 - U1*SS*R1 + T1*SS*
     .R1) + 2*D0ABC*(U1**3*SS + 2*U1**2*T1 + 2*U1**2*SS*R1
     . - 2*U1*T1**2 - U1*T1*SS*R1 - U1*SS**2*R1 + U1*SS*R1**
     .2 - 2*U1*SS*R1 - T1*SS*R1**2 + 2*T1*SS*R1 - SS**2*R1**
     .2))/(U1*T1 - SS*R1)
 
      ENDIF

      ELSE

       A1 = 2*SS
       A2 = 0
       H1 = 4*SS/3.D0
       H2 = 0
       Z1 = 0
       Z2 = 0
       AA1 = 4*SS/3.D0
       AA2 = 0
       HH1 = -4*SS/3.D0
       HH2 = 0
       AH1 = -2*SS
       AH2 = 0
       IF(ITRIA.GT.2) THEN
        AA1 = 0
        AA2 = 0
        HH1 = 0
        HH2 = 0
        AH1 = 0
        AH2 = 0
       ENDIF
      ENDIF

c     write(6,*)H1/2
c     write(6,*)H2/2
c     write(6,*)HH1/2
c     write(6,*)HH2/2
c     write(6,*)D0BAC/SS
c     write(6,*)D0ACB/SS
c     write(6,*)D0ABC/SS
c     write(6,*)

C---SIGN MISTAKE IN Z CONTRIBUTIONS!!!

      A1 = A1*AMQ
      H1 = H1*AMQ
      Z1 = -Z1*AMQ**2
      A2 = A2*AMQ
      H2 = H2*AMQ
      Z2 = -Z2*AMQ**2

      RETURN
      END
  
      SUBROUTINE FORMRES(AMQ)
C--RESONANCE FORM FACTORS FOR LO MATRIX ELEMENTS FOR GG -> HH
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DOUBLE PRECISION M1,M2,MZ
      COMPLEX*16 A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMPLEX*16 C0AB
      COMMON/PARAM/S,T,U,M1,M2,MZ
      COMMON/FORM/A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/PROC/IPROC,IAPP,ITOP,ISFAC
      COMMON/HMASS/AMA,AML,AMH,AMCH
      COMMON/CUTS/EPS,REPS
      S0= S
      S = AMH**2
      SS=S/AMQ**2
      TT=T/AMQ**2
      UU=U/AMQ**2
      R1=M1**2/AMQ**2
      R2=M2**2/AMQ**2
      RZ=MZ**2/AMQ**2
      T1=TT-R1
      U1=UU-R1
      IF(IAPP.EQ.0)THEN
       DQ2=AMQ**2
       CQ2=AMQ**2*DCMPLX(1.D0,-REPS)
       CA5 = CDSQRT(1.D0-4.D0*CQ2/S)
       C0AB = 0.5D0*CDLOG((CA5+1.D0)/(CA5-1.D0))**2/S*DQ2
       H1 =  - 4*(C0AB*SS - 4*C0AB - 2)
      ELSE
       H1 = 4*SS/3.D0
      ENDIF
      H1 = H1*AMQ
      S = S0
      RETURN
      END
  
      SUBROUTINE VIRTINT
C--SCALAR INTEGRALS FOR GG -> HH AT LO
      PARAMETER(N=3)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DOUBLE PRECISION M1,M2,MZ
      DIMENSION XX(N),Y(N)
      COMPLEX*16 C0AB,C0AC,C0AD,C0BC,C0BD,C0CD,D0ABC,D0BAC,D0ACB
      COMPLEX*16 CD0AB,CD0AC,CD0AD,CD0BC,CD0BD,CD0CD,
     .           DD0ABC,DD0BAC,DD0ACB
      COMPLEX*16  A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/FORM/A1,A2,H1,H2,Z1,Z2,AA1,AA2,HH1,HH2,AH1,AH2
      COMMON/INOUT/NIN,NOUT
      COMMON/PARAM/S,T,U,M1,M2,MZ
      COMMON/SCALINT/C0AB,C0AC,C0AD,C0BC,C0BD,C0CD,D0ABC,D0BAC,D0ACB
      COMMON/SCAL0/V,TAU,AMQ
      COMMON/TOTPAR/W,GEVPB,GF,GAMZ,GAMA,GAML,GAMH
      COMMON/VEGPAR/ABSERR,IPOINT,ITERAT,NPRN
      COMMON/RESULT/RES,ERR,DCHI2,DUM
      COMMON/INT4PAR/CD1,CD2,CD3,CD12,CD13,CD23
      COMMON/INT3PAR/CC1,CC2,CC12
      COMMON/CUTS/EPS,REPS
      EXTERNAL DRINT3,DIINT3,DRINT4,DIINT4
      CQ2 = AMQ**2*DCMPLX(1.D0,-REPS)
      TH=(M1+M2)**2/W**2
      S=TAU*W**2
      XLAM=DSQRT((S+M1**2-M2**2)**2-4.D0*S*M1**2)
      BET=XLAM/(S+M1**2-M2**2)
      T1=-(S+M1**2-M2**2)*((1.D0+BET)/2.D0-BET*V)
      U1=-(S+M1**2-M2**2)*((1.D0-BET)/2.D0+BET*V)
      T=T1+M1**2
      U=U1+M1**2
      SS=S/AMQ**2
      UU1=U1/AMQ**2
      TT1=T1/AMQ**2
      R1=M1**2/AMQ**2
      R2=M2**2/AMQ**2
      WRITE(NOUT,*)'EPS   = ',EPS
      WRITE(NOUT,*)'IMAG  = ',REPS
      WRITE(NOUT,*)
      WRITE(NOUT,*)'VEGAS:'
      WRITE(NOUT,*)'======'
      WRITE(NOUT,*)'ACCURACY   = ',ABSERR
      WRITE(NOUT,*)'POINTS     = ',IPOINT
      WRITE(NOUT,*)'ITERATIONS = ',ITERAT
      WRITE(NOUT,*)
      WRITE(NOUT,*)'V     = ',V
      WRITE(NOUT,*)'TAU   = ',TAU
      WRITE(NOUT,*)'M1    = ',M1
      WRITE(NOUT,*)'M2    = ',M2
      WRITE(NOUT,*)'MQ    = ',AMQ
      WRITE(NOUT,*)'S     = ',S
      WRITE(NOUT,*)'T     = ',T
      WRITE(NOUT,*)'U     = ',U
      WRITE(NOUT,*)'T1    = ',T1
      WRITE(NOUT,*)'U1    = ',U1
      WRITE(NOUT,*)'S+T+U-M1**2-M2**2 = ',S+T+U-M1**2-M2**2
      WRITE(NOUT,*)
      C0AB = -0.5D0 - SS/24.D0
      C0AC = -0.5D0 - (R1+TT1+R1)/24.D0
      C0AD = -0.5D0 - (R2+UU1+R1)/24.D0
      C0BC = -0.5D0 - (R1+UU1+R1)/24.D0
      C0BD = -0.5D0 - (R2+TT1+R1)/24.D0
      C0CD = -0.5D0 - (R1+R2+SS)/24.D0
      D0ABC = 1.D0/6.D0 + (R1+R2+SS+UU1+R1)/60.D0
      D0BAC = 1.D0/6.D0 + (R1+R2+SS+TT1+R1)/60.D0
      D0ACB = 1.D0/6.D0 + (R1+R2+TT1+R1+UU1+R1)/60.D0
      WRITE(NOUT,*)'APPROXIMATION:'
      WRITE(NOUT,*)'=============='
      WRITE(NOUT,*)
      WRITE(NOUT,*)'C0AB  = ',DREAL(C0AB),' +I* ',DIMAG(C0AB)
      WRITE(NOUT,*)'C0CD  = ',DREAL(C0CD),' +I* ',DIMAG(C0CD)
      WRITE(NOUT,*)'C0AC  = ',DREAL(C0AC),' +I* ',DIMAG(C0AC)
      WRITE(NOUT,*)'C0AD  = ',DREAL(C0AD),' +I* ',DIMAG(C0AD)
      WRITE(NOUT,*)'C0BC  = ',DREAL(C0BC),' +I* ',DIMAG(C0BC)
      WRITE(NOUT,*)'C0BD  = ',DREAL(C0BD),' +I* ',DIMAG(C0BD)
      WRITE(NOUT,*)'D0ABC = ',DREAL(D0ABC),' +I* ',DIMAG(D0ABC)
      WRITE(NOUT,*)'D0BAC = ',DREAL(D0BAC),' +I* ',DIMAG(D0BAC)
      WRITE(NOUT,*)'D0ACB = ',DREAL(D0ACB),' +I* ',DIMAG(D0ACB)
      WRITE(NOUT,*)
      CALL INISCAL(AMQ)
      WRITE(NOUT,*)'ANALYTICALLY:'
      WRITE(NOUT,*)'============='
      WRITE(NOUT,*)
      WRITE(NOUT,*)'C0AB  = ',DREAL(C0AB),' +I* ',DIMAG(C0AB)
      WRITE(NOUT,*)'C0CD  = ',DREAL(C0CD),' +I* ',DIMAG(C0CD)
      WRITE(NOUT,*)'C0AC  = ',DREAL(C0AC),' +I* ',DIMAG(C0AC)
      WRITE(NOUT,*)'C0AD  = ',DREAL(C0AD),' +I* ',DIMAG(C0AD)
      WRITE(NOUT,*)'C0BC  = ',DREAL(C0BC),' +I* ',DIMAG(C0BC)
      WRITE(NOUT,*)'C0BD  = ',DREAL(C0BD),' +I* ',DIMAG(C0BD)
      WRITE(NOUT,*)'D0ABC = ',DREAL(D0ABC),' +I* ',DIMAG(D0ABC)
      WRITE(NOUT,*)'D0BAC = ',DREAL(D0BAC),' +I* ',DIMAG(D0BAC)
      WRITE(NOUT,*)'D0ACB = ',DREAL(D0ACB),' +I* ',DIMAG(D0ACB)
      WRITE(NOUT,*)

      IGRAPH=0
      NDIM=2
      CC1=0.D0/CQ2
      CC2=S/CQ2
      CC12=S/CQ2
      CALL VEGASN(DRINT3,ABSERR,NDIM,IPOINT,ITERAT,NPRN,IGRAPH)
      X=RES
      DX=ERR
      CALL VEGASN(DIINT3,ABSERR,NDIM,IPOINT,ITERAT,NPRN,IGRAPH)
      C0AB=DCMPLX(X,RES)
      CD0AB=DCMPLX(DX,ERR)

      CC1=M1**2/CQ2
      CC2=S/CQ2
      CC12=(S+M1**2-M2**2)/CQ2
      CALL VEGASN(DRINT3,ABSERR,NDIM,IPOINT,ITERAT,NPRN,IGRAPH)
      X=RES
      DX=ERR
      CALL VEGASN(DIINT3,ABSERR,NDIM,IPOINT,ITERAT,NPRN,IGRAPH)
      C0CD=DCMPLX(X,RES)
      CD0CD=DCMPLX(DX,ERR)

      CC1=0.D0/CQ2
      CC2=T/CQ2
      CC12=(T-M1**2)/CQ2
      CALL VEGASN(DRINT3,ABSERR,NDIM,IPOINT,ITERAT,NPRN,IGRAPH)
      X=RES
      DX=ERR
      CALL VEGASN(DIINT3,ABSERR,NDIM,IPOINT,ITERAT,NPRN,IGRAPH)
      C0AC=DCMPLX(X,RES)
      CD0AC=DCMPLX(DX,ERR)

      CC1=0.D0/CQ2
      CC2=U/CQ2
      CC12=(U-M2**2)/CQ2
      CALL VEGASN(DRINT3,ABSERR,NDIM,IPOINT,ITERAT,NPRN,IGRAPH)
      X=RES
      DX=ERR
      CALL VEGASN(DIINT3,ABSERR,NDIM,IPOINT,ITERAT,NPRN,IGRAPH)
      C0AD=DCMPLX(X,RES)
      CD0AD=DCMPLX(DX,ERR)

      CC1=0.D0/CQ2
      CC2=U/CQ2
      CC12=(U-M1**2)/CQ2
      CALL VEGASN(DRINT3,ABSERR,NDIM,IPOINT,ITERAT,NPRN,IGRAPH)
      X=RES
      DX=ERR
      CALL VEGASN(DIINT3,ABSERR,NDIM,IPOINT,ITERAT,NPRN,IGRAPH)
      C0BC=DCMPLX(X,RES)
      CD0BC=DCMPLX(DX,ERR)

      CC1=0.D0/CQ2
      CC2=T/CQ2
      CC12=(T-M2**2)/CQ2
      CALL VEGASN(DRINT3,ABSERR,NDIM,IPOINT,ITERAT,NPRN,IGRAPH)
      X=RES
      DX=ERR
      CALL VEGASN(DIINT3,ABSERR,NDIM,IPOINT,ITERAT,NPRN,IGRAPH)
      C0BD=DCMPLX(X,RES)
      CD0BD=DCMPLX(DX,ERR)

      NDIM=3
      CD1=0.D0/CQ2
      CD2=S/CQ2
      CD3=M2**2/CQ2
      CD12=S/CQ2
      CD13=-(U-M2**2)/CQ2
      CD23=(S+M2**2-M1**2)/CQ2
      CALL VEGASN(DRINT4,ABSERR,NDIM,IPOINT,ITERAT,NPRN,IGRAPH)
      X=RES
      DX=ERR
      CALL VEGASN(DIINT4,ABSERR,NDIM,IPOINT,ITERAT,NPRN,IGRAPH)
      D0ABC=DCMPLX(X,RES)
      DD0ABC=DCMPLX(DX,ERR)

      CD1=0.D0/CQ2
      CD2=S/CQ2
      CD3=M2**2/CQ2
      CD12=S/CQ2
      CD13=-(T-M2**2)/CQ2
      CD23=(S+M2**2-M1**2)/CQ2
      CALL VEGASN(DRINT4,ABSERR,NDIM,IPOINT,ITERAT,NPRN,IGRAPH)
      X=RES
      DX=ERR
      CALL VEGASN(DIINT4,ABSERR,NDIM,IPOINT,ITERAT,NPRN,IGRAPH)
      D0BAC=DCMPLX(X,RES)
      DD0BAC=DCMPLX(DX,ERR)

      CD1=0.D0/CQ2
      CD2=T/CQ2
      CD3=M2**2/CQ2
      CD12=(T-M1**2)/CQ2
      CD13=-(U-M2**2)/CQ2
      CD23=(T+M2**2)/CQ2
      CALL VEGASN(DRINT4,ABSERR,NDIM,IPOINT,ITERAT,NPRN,IGRAPH)
      X=RES
      DX=ERR
      CALL VEGASN(DIINT4,ABSERR,NDIM,IPOINT,ITERAT,NPRN,IGRAPH)
      D0ACB=DCMPLX(X,RES)
      DD0ACB=DCMPLX(DX,ERR)

      WRITE(NOUT,*)'NUMERICALLY:'
      WRITE(NOUT,*)'============'
      WRITE(NOUT,*)
      WRITE(NOUT,*)'C0AB  = ',DREAL(C0AB),' +I* ',DIMAG(C0AB)
      WRITE(NOUT,*)'C0CD  = ',DREAL(C0CD),' +I* ',DIMAG(C0CD)
      WRITE(NOUT,*)'C0AC  = ',DREAL(C0AC),' +I* ',DIMAG(C0AC)
      WRITE(NOUT,*)'C0AD  = ',DREAL(C0AD),' +I* ',DIMAG(C0AD)
      WRITE(NOUT,*)'C0BC  = ',DREAL(C0BC),' +I* ',DIMAG(C0BC)
      WRITE(NOUT,*)'C0BD  = ',DREAL(C0BD),' +I* ',DIMAG(C0BD)
      WRITE(NOUT,*)'D0ABC = ',DREAL(D0ABC),' +I* ',DIMAG(D0ABC)
      WRITE(NOUT,*)'D0BAC = ',DREAL(D0BAC),' +I* ',DIMAG(D0BAC)
      WRITE(NOUT,*)'D0ACB = ',DREAL(D0ACB),' +I* ',DIMAG(D0ACB)
      WRITE(NOUT,*)
      WRITE(NOUT,*)'ERRORS:'
      WRITE(NOUT,*)'======='
      WRITE(NOUT,*)
      WRITE(NOUT,*)'C0AB  = ',DREAL(CD0AB),' +I* ',DIMAG(CD0AB)
      WRITE(NOUT,*)'C0CD  = ',DREAL(CD0CD),' +I* ',DIMAG(CD0CD)
      WRITE(NOUT,*)'C0AC  = ',DREAL(CD0AC),' +I* ',DIMAG(CD0AC)
      WRITE(NOUT,*)'C0AD  = ',DREAL(CD0AD),' +I* ',DIMAG(CD0AD)
      WRITE(NOUT,*)'C0BC  = ',DREAL(CD0BC),' +I* ',DIMAG(CD0BC)
      WRITE(NOUT,*)'C0BD  = ',DREAL(CD0BD),' +I* ',DIMAG(CD0BD)
      WRITE(NOUT,*)'D0ABC = ',DREAL(DD0ABC),' +I* ',DIMAG(DD0ABC)
      WRITE(NOUT,*)'D0BAC = ',DREAL(DD0BAC),' +I* ',DIMAG(DD0BAC)
      WRITE(NOUT,*)'D0ACB = ',DREAL(DD0ACB),' +I* ',DIMAG(DD0ACB)
      WRITE(NOUT,*)

      CALL FORMFAC(AMQ)
      H1=H1/AMQ/SS
      H2=H2/AMQ/SS
      A1=A1/AMQ/SS
      A2=A2/AMQ/SS
      Z1=Z1/S
      Z2=Z2/S
      HH1=HH1/SS
      HH2=HH2/SS
      AA1=AA1/SS
      AA2=AA2/SS
      AH1=AH1/SS
      AH2=AH2/SS
      WRITE(NOUT,*)'FORMFACTORS:'
      WRITE(NOUT,*)'============'
      WRITE(NOUT,*)
      WRITE(NOUT,*)'H1    = ',DREAL(H1),' +I* ',DIMAG(H1)
      WRITE(NOUT,*)'H2    = ',DREAL(H2),' +I* ',DIMAG(H2)
      WRITE(NOUT,*)'A1    = ',DREAL(A1),' +I* ',DIMAG(A1)
      WRITE(NOUT,*)'A2    = ',DREAL(A2),' +I* ',DIMAG(A2)
      WRITE(NOUT,*)'Z1    = ',DREAL(Z1),' +I* ',DIMAG(Z1)
      WRITE(NOUT,*)'Z2    = ',DREAL(Z2),' +I* ',DIMAG(Z2)
      WRITE(NOUT,*)'HH1   = ',DREAL(HH1),' +I* ',DIMAG(HH1)
      WRITE(NOUT,*)'HH2   = ',DREAL(HH2),' +I* ',DIMAG(HH2)
      WRITE(NOUT,*)'AA1   = ',DREAL(AA1),' +I* ',DIMAG(AA1)
      WRITE(NOUT,*)'AA2   = ',DREAL(AA2),' +I* ',DIMAG(AA2)
      WRITE(NOUT,*)'AH1   = ',DREAL(AH1),' +I* ',DIMAG(AH1)
      WRITE(NOUT,*)'AH2   = ',DREAL(AH2),' +I* ',DIMAG(AH2)
      WRITE(NOUT,*)

      A1 = 2
      A2 = 0
      H1 = 4/3.D0
      H2 = 0
      Z1 = 0
      Z2 = 0
      AA1 = 4/3.D0
      AA2 = 0
      HH1 = -4/3.D0
      HH2 = 0
      AH1 = -2
      AH2 = 0
      WRITE(NOUT,*)'APPROXIMATION:'
      WRITE(NOUT,*)'=============='
      WRITE(NOUT,*)
      WRITE(NOUT,*)'H1    = ',DREAL(H1),' +I* ',DIMAG(H1)
      WRITE(NOUT,*)'H2    = ',DREAL(H2),' +I* ',DIMAG(H2)
      WRITE(NOUT,*)'A1    = ',DREAL(A1),' +I* ',DIMAG(A1)
      WRITE(NOUT,*)'A2    = ',DREAL(A2),' +I* ',DIMAG(A2)
      WRITE(NOUT,*)'Z1    = ',DREAL(Z1),' +I* ',DIMAG(Z1)
      WRITE(NOUT,*)'Z2    = ',DREAL(Z2),' +I* ',DIMAG(Z2)
      WRITE(NOUT,*)'HH1   = ',DREAL(HH1),' +I* ',DIMAG(HH1)
      WRITE(NOUT,*)'HH2   = ',DREAL(HH2),' +I* ',DIMAG(HH2)
      WRITE(NOUT,*)'AA1   = ',DREAL(AA1),' +I* ',DIMAG(AA1)
      WRITE(NOUT,*)'AA2   = ',DREAL(AA2),' +I* ',DIMAG(AA2)
      WRITE(NOUT,*)'AH1   = ',DREAL(AH1),' +I* ',DIMAG(AH1)
      WRITE(NOUT,*)'AH2   = ',DREAL(AH2),' +I* ',DIMAG(AH2)

      RETURN
      END
  
      DOUBLE PRECISION FUNCTION DRINT3(XX)
C--INTEGRAND FOR NUMERICAL INTEGRATION: REAL PART
      PARAMETER (N=2)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DIMENSION XX(N)
      DRINT3=DREAL(CINT3(XX))
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DIINT3(XX)
C--INTEGRAND FOR NUMERICAL INTEGRATION: IMAGINARY PART
      PARAMETER (N=2)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DIMENSION XX(N)
      DIINT3=DIMAG(CINT3(XX))
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DRINT4(XX)
C--INTEGRAND FOR NUMERICAL INTEGRATION: REAL PART
      PARAMETER (N=3)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DIMENSION XX(N)
      DRINT4=DREAL(CINT4(XX))
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DIINT4(XX)
C--INTEGRAND FOR NUMERICAL INTEGRATION: IMAGINARY PART
      PARAMETER (N=3)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DIMENSION XX(N)
      DIINT4=DIMAG(CINT4(XX))
      RETURN
      END
 
      COMPLEX*16 FUNCTION CINT3(VV)
C--COMPLEX INTEGRAND FOR NUMERICAL INTEGRATION
      PARAMETER (N=2)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DIMENSION VV(N)
      COMMON/CUTS/EPS,REPS
      COMMON/INT3PAR/C1,C2,C12
      XX = EPS + (1.D0-EPS)*VV(1)
      YY = EPS + (1.D0-EPS)*VV(2)
      X = 1.D0-XX
      Y = XX*YY
      CN = 1.D0 - C1*X*(1.D0-X) - C2*Y*(1.D0-Y) + C12*X*Y
      CINT3 = -XX/CN
      RETURN
      END
 
      COMPLEX*16 FUNCTION CINT4(VV)
C--COMPLEX INTEGRAND FOR NUMERICAL INTEGRATION
      PARAMETER (N=3)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DIMENSION VV(N)
      COMMON/CUTS/EPS,REPS
      COMMON/INT4PAR/C1,C2,C3,C12,C13,C23
      XX = EPS + (1.D0-EPS)*VV(1)
      YY = EPS + (1.D0-EPS)*VV(2)
      ZZ = EPS + (1.D0-EPS)*VV(3)
      X = 1.D0-XX
      Y = XX*(1.D0-YY)
      Z = XX*YY*ZZ
      CN = 1.D0 - C1*X*(1.D0-X) - C2*Y*(1.D0-Y)
     .          - C3*Z*(1.D0-Z) + C12*X*Y + C13*X*Z + C23*Y*Z
      CINT4 = XX**2*YY/CN**2
      RETURN
      END
 
      SUBROUTINE INISCAL(AMQ)
C--INITIALIZATION OF SCALAR INTEGRALS
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
      DOUBLE PRECISION M1,M2,MZ
      COMPLEX*16 D04
      COMPLEX*16 C03
      COMPLEX*16 LI2
      COMPLEX*16 CI,CJ
      COMPLEX*16 R,R0,X,XP,XM
      COMPLEX*16 C0AB,C0AC,C0AD,C0BC,C0BD,C0CD,D0ABC,D0BAC,D0ACB
      COMMON/PARAM/S,T,U,M1,M2,MZ
      COMMON/SCALINT/C0AB,C0AC,C0AD,C0BC,C0BD,C0CD,D0ABC,D0BAC,D0ACB
      COMMON/CUTS/EPS,REPS
      COMMON/KNIEHL/IKNIEHL
      DQ2=AMQ**2

      IF(IKNIEHL.EQ.1)THEN

      CQ2=AMQ**2*DCMPLX(1.D0,-REPS)
      CA5 = CDSQRT(1.D0-4.D0*CQ2/S)
      C0AB = 0.5D0*CDLOG((CA5+1.D0)/(CA5-1.D0))**2/S*DQ2
      S1=M1**2
      S2=M2**2
      S5=S
      CA1 = CDSQRT(1.D0-4.D0*CQ2/S1)
      CA2 = CDSQRT(1.D0-4.D0*CQ2/S2)
      CA5 = CDSQRT(1.D0-4.D0*CQ2/S5)
      XLAM=DSQRT(S1**2+S2**2+S5**2-2.D0*(S1*S2+S1*S5+S2*S5))
      CB1=1.D0/XLAM*(S1-S2-S5)
      CB2=1.D0/XLAM*(S2-S1-S5)
      CB5=1.D0/XLAM*(S5-S1-S2)
      C0CD = -(
     .       LI2((1.D0+CB1)/(CA1+CB1)) - LI2((1.D0-CB1)/(CA1-CB1))
     .     - LI2((-1.D0+CB1)/(CA1+CB1)) + LI2((-1.D0-CB1)/(CA1-CB1))
     .     + LI2((1.D0+CB2)/(CA2+CB2)) - LI2((1.D0-CB2)/(CA2-CB2))
     .     - LI2((-1.D0+CB2)/(CA2+CB2)) + LI2((-1.D0-CB2)/(CA2-CB2))
     .     + LI2((1.D0+CB5)/(CA5+CB5)) - LI2((1.D0-CB5)/(CA5-CB5))
     .     - LI2((-1.D0+CB5)/(CA5+CB5)) + LI2((-1.D0-CB5)/(CA5-CB5))
     .       )/XLAM*DQ2
      S1=0.D0
      S2=M1**2
      S5=T
      CA2 = CDSQRT(1.D0-4.D0*CQ2/S2)
      CA5 = CDSQRT(1.D0-4.D0*CQ2/S5)
      XLAM=DSQRT(S1**2+S2**2+S5**2-2.D0*(S1*S2+S1*S5+S2*S5))
      CB2=1.D0/XLAM*(S2-S1-S5)
      CB5=1.D0/XLAM*(S5-S1-S2)
      C0AC = -(
     .       LI2((1.D0+CB2)/(CA2+CB2)) - LI2((1.D0-CB2)/(CA2-CB2))
     .     - LI2((-1.D0+CB2)/(CA2+CB2)) + LI2((-1.D0-CB2)/(CA2-CB2))
     .     + LI2((1.D0+CB5)/(CA5+CB5)) - LI2((1.D0-CB5)/(CA5-CB5))
     .     - LI2((-1.D0+CB5)/(CA5+CB5)) + LI2((-1.D0-CB5)/(CA5-CB5))
     .       )/XLAM*DQ2
      S1=0.D0
      S2=M2**2
      S5=U
      CA2 = CDSQRT(1.D0-4.D0*CQ2/S2)
      CA5 = CDSQRT(1.D0-4.D0*CQ2/S5)
      XLAM=DSQRT(S1**2+S2**2+S5**2-2.D0*(S1*S2+S1*S5+S2*S5))
      CB2=1.D0/XLAM*(S2-S1-S5)
      CB5=1.D0/XLAM*(S5-S1-S2)
      C0AD = -(
     .       LI2((1.D0+CB2)/(CA2+CB2)) - LI2((1.D0-CB2)/(CA2-CB2))
     .     - LI2((-1.D0+CB2)/(CA2+CB2)) + LI2((-1.D0-CB2)/(CA2-CB2))
     .     + LI2((1.D0+CB5)/(CA5+CB5)) - LI2((1.D0-CB5)/(CA5-CB5))
     .     - LI2((-1.D0+CB5)/(CA5+CB5)) + LI2((-1.D0-CB5)/(CA5-CB5))
     .       )/XLAM*DQ2
      S1=0.D0
      S2=M1**2
      S5=U
      CA2 = CDSQRT(1.D0-4.D0*CQ2/S2)
      CA5 = CDSQRT(1.D0-4.D0*CQ2/S5)
      XLAM=DSQRT(S1**2+S2**2+S5**2-2.D0*(S1*S2+S1*S5+S2*S5))
      CB2=1.D0/XLAM*(S2-S1-S5)
      CB5=1.D0/XLAM*(S5-S1-S2)
      C0BC = -(
     .       LI2((1.D0+CB2)/(CA2+CB2)) - LI2((1.D0-CB2)/(CA2-CB2))
     .     - LI2((-1.D0+CB2)/(CA2+CB2)) + LI2((-1.D0-CB2)/(CA2-CB2))
     .     + LI2((1.D0+CB5)/(CA5+CB5)) - LI2((1.D0-CB5)/(CA5-CB5))
     .     - LI2((-1.D0+CB5)/(CA5+CB5)) + LI2((-1.D0-CB5)/(CA5-CB5))
     .       )/XLAM*DQ2
      S1=0.D0
      S2=M2**2
      S5=T
      CA2 = CDSQRT(1.D0-4.D0*CQ2/S2)
      CA5 = CDSQRT(1.D0-4.D0*CQ2/S5)
      XLAM=DSQRT(S1**2+S2**2+S5**2-2.D0*(S1*S2+S1*S5+S2*S5))
      CB2=1.D0/XLAM*(S2-S1-S5)
      CB5=1.D0/XLAM*(S5-S1-S2)
      C0BD = -(
     .       LI2((1.D0+CB2)/(CA2+CB2)) - LI2((1.D0-CB2)/(CA2-CB2))
     .     - LI2((-1.D0+CB2)/(CA2+CB2)) + LI2((-1.D0-CB2)/(CA2-CB2))
     .     + LI2((1.D0+CB5)/(CA5+CB5)) - LI2((1.D0-CB5)/(CA5-CB5))
     .     - LI2((-1.D0+CB5)/(CA5+CB5)) + LI2((-1.D0-CB5)/(CA5-CB5))
     .       )/XLAM*DQ2
      Z=M1**2
      H=M2**2
      S1=S
      S2=Z
      S5=H
      XLAM=DSQRT(S1**2+S2**2+S5**2-2.D0*(S1*S2+S1*S5+S2*S5))
      XN=T*U-Z*H
      R0=CDSQRT(1.D0+4.D0*S*CQ2/XN)
      R=(1.D0+R0)/2.D0
      RR=DREAL(R0)
      D0ACB = 2.D0/XN/RR*(CJ(Z,AMQ**2,R)+CJ(H,AMQ**2,R)
     .                   -CJ(T,AMQ**2,R)-CJ(U,AMQ**2,R))*DQ2**2
      X=CDSQRT(1.D0+4.D0*XN/S/U**2*CQ2)
      X0=DSQRT(1.D0+4.D0*XN/S/U**2*DQ2)
      XP=-U/(T-U+XLAM)*(1.D0+X)
      XM=-U/(T-U+XLAM)*(1.D0-X)
      ALP=(S+Z-H+XLAM)/2.D0/S
      BET=(U-T+XLAM)/2.D0/S
      ZERO=0.D-15
      D0ABC = 1.D0/S/U/X0*(CI(S,ZERO,XP)-CI(T-H,S,(1.D0-XM)/(1.D0-ALP))
     .      -CI(Z-U,ZERO,XP/ALP)+CI(Z-U,ZERO,XP/BET)
     .      -CI(T-H,S,XP/BET)
     .      -CJ(S,AMQ**2,XP)+CJ(H,AMQ**2,(1.D0-XM)/(1.D0-ALP))
     .      +CJ(Z,AMQ**2,XP/ALP)-CJ(U,AMQ**2,XP/BET)
     .      +CJ(H,AMQ**2,XP/BET)
     .      -(CI(S,ZERO,XM)-CI(T-H,S,(1.D0-XP)/(1.D0-ALP))
     .      -CI(Z-U,ZERO,XM/ALP)+CI(Z-U,ZERO,XM/BET)
     .      -CI(T-H,S,XM/BET)
     .      -CJ(S,AMQ**2,XM)+CJ(H,AMQ**2,(1.D0-XP)/(1.D0-ALP))
     .      +CJ(Z,AMQ**2,XM/ALP)-CJ(U,AMQ**2,XM/BET)
     .      +CJ(H,AMQ**2,XM/BET)))*DQ2**2
      X=CDSQRT(1.D0+4.D0*XN/S/T**2*CQ2)
      X0=DSQRT(1.D0+4.D0*XN/S/T**2*DQ2)
      XP=-T/(U-T+XLAM)*(1.D0+X)
      XM=-T/(U-T+XLAM)*(1.D0-X)
      ALP=(S+Z-H+XLAM)/2.D0/S
      BET=(T-U+XLAM)/2.D0/S
      D0BAC = 1.D0/S/T/X0*(CI(S,ZERO,XP)-CI(U-H,S,(1.D0-XM)/(1.D0-ALP))
     .      -CI(Z-T,ZERO,XP/ALP)+CI(Z-T,ZERO,XP/BET)
     .      -CI(U-H,S,XP/BET)
     .      -CJ(S,AMQ**2,XP)+CJ(H,AMQ**2,(1.D0-XM)/(1.D0-ALP))
     .      +CJ(Z,AMQ**2,XP/ALP)-CJ(T,AMQ**2,XP/BET)
     .      +CJ(H,AMQ**2,XP/BET)
     .      -(CI(S,ZERO,XM)-CI(U-H,S,(1.D0-XP)/(1.D0-ALP))
     .      -CI(Z-T,ZERO,XM/ALP)+CI(Z-T,ZERO,XM/BET)
     .      -CI(U-H,S,XM/BET)
     .      -CJ(S,AMQ**2,XM)+CJ(H,AMQ**2,(1.D0-XP)/(1.D0-ALP))
     .      +CJ(Z,AMQ**2,XM/ALP)-CJ(T,AMQ**2,XM/BET)
     .      +CJ(H,AMQ**2,XM/BET)))*DQ2**2
      ELSE

c     write(6,*)C0AB,C0CD,C0AC,C0AD,C0BC,C0BD

      EPM = 0.D-8
      C0AB = C03(EPM,EPM,S,AMQ,AMQ,AMQ)*DQ2
      C0AC = C03(EPM,M1**2,T,AMQ,AMQ,AMQ)*DQ2
      C0AD = C03(EPM,M2**2,U,AMQ,AMQ,AMQ)*DQ2
      C0BC = C03(EPM,M1**2,U,AMQ,AMQ,AMQ)*DQ2
      C0BD = C03(EPM,M2**2,T,AMQ,AMQ,AMQ)*DQ2
      C0CD = C03(M1**2,M2**2,S,AMQ,AMQ,AMQ)*DQ2

c     write(6,*)C0AB,C0CD,C0AC,C0AD,C0BC,C0BD
c     write(6,*)
c     write(6,*)

************************************************************************
*       FUNCTION D04(P1,P2,P3,P4,P12,P23,M1,M2,M3,M4)                  *
*  SCALAR 4-POINT FUNCTION WITH AT LEAST ONE MASS ZERO                 *
*  P1,P2,P3,P4 = SQUARED EXTERNAL MOMENTA			       *
*  P12 = (p1+p2)**2,  P23 = (p2+p3)**2				       *

      D0ABC = DQ2**2*D04(0.D0,0.D0,M1**2,M2**2,S,U,AMQ,AMQ,AMQ,AMQ)
      D0BAC = DQ2**2*D04(0.D0,0.D0,M1**2,M2**2,S,T,AMQ,AMQ,AMQ,AMQ)
      D0ACB = DQ2**2*D04(0.D0,M1**2,0.D0,M2**2,T,U,AMQ,AMQ,AMQ,AMQ)
      ENDIF

      RETURN
      END

      DOUBLE PRECISION FUNCTION DLUMGG(TAU,X,QQ)
C--GG-LUMINOSITY
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PDF(-6:6)
      X1=X
      X2=TAU/X1
      CALL STRUC(X1,QQ,PDF)
      GG=PDF(0)
      CALL STRUC(X2,QQ,PDF)
      DLUM=GG*PDF(0)
      DLUMGG=DLUM/TAU/X
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DLUMGQ(TAU,X,QQ)
C--GQ-LUMINOSITY
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PDF1(-6:6),PDF2(-6:6)
      X1=X
      X2=TAU/X1
      CALL STRUC(X1,QQ,PDF1)
      CALL STRUC(X2,QQ,PDF2)
      DLUM = 0
      DO 1 I=1,5
       DLUM = DLUM + PDF1(I)*PDF2(0) + PDF1(-I)*PDF2(0)
     .             + PDF2(I)*PDF1(0) + PDF2(-I)*PDF1(0)
1     CONTINUE
      DLUMGQ=DLUM/TAU/X
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DLUMQQ(TAU,X,QQ)
C--Q QBAR-LUMINOSITY
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PDF1(-6:6),PDF2(-6:6)
      COMMON/COLLIDER/ICOLL
      IF(ICOLL.EQ.0)THEN
       NCOLL = 1
      ELSE
       NCOLL = -1
      ENDIF
      X1=X
      X2=TAU/X1
      CALL STRUC(X1,QQ,PDF1)
      CALL STRUC(X2,QQ,PDF2)
      DLUM = 0
      DO 1 I=1,5
       DLUM=DLUM + PDF1(I)*PDF2(-NCOLL*I) + PDF1(-I)*PDF2(NCOLL*I)
1     CONTINUE
      DLUMQQ=DLUM/TAU/X
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DLUMQQZ(TAU,X,QQ)
C--Q QBAR-LUMINOSITY WITH Z-CHARGES
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PDF1(-6:6),PDF2(-6:6)
      COMMON/WEINBERG/SW2
      COMMON/COLLIDER/ICOLL
      IF(ICOLL.EQ.0)THEN
       NCOLL = 1
      ELSE
       NCOLL = -1
      ENDIF
      CU = 1.D0 + (1.D0-8.D0/3.D0*SW2)**2
      CD = 1.D0 + (1.D0-4.D0/3.D0*SW2)**2
      X1=X
      X2=TAU/X1
      CALL STRUC(X1,QQ,PDF1)
      CALL STRUC(X2,QQ,PDF2)
      DLUMU=0.D0
      DLUMD=0.D0
      DO 1 I=1,5,2
       DLUMD=DLUMD + PDF1(I)*PDF2(-NCOLL*I) + PDF1(-I)*PDF2(NCOLL*I)
       J=I+1
       IF(J.LT.6)THEN
       DLUMU=DLUMU + PDF1(J)*PDF2(-NCOLL*J) + PDF1(-J)*PDF2(NCOLL*J)
       ENDIF
1     CONTINUE
      DLUMQQZ=(CU*DLUMU+CD*DLUMD)/TAU/X
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DLUMGQZ(TAU,X,QQ)
C--G Q-LUMINOSITY WITH Z-CHARGES
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PDF1(-6:6),PDF2(-6:6)
      COMMON/WEINBERG/SW2
      CU = 1.D0 + (1.D0-8.D0/3.D0*SW2)**2
      CD = 1.D0 + (1.D0-4.D0/3.D0*SW2)**2
      X1=X
      X2=TAU/X1
      CALL STRUC(X1,QQ,PDF1)
      CALL STRUC(X2,QQ,PDF2)
      DLUMU=0.D0
      DLUMD=0.D0
      DO 1 I=1,5,2
       DLUMD=DLUMD + PDF1(0)*PDF2(I) + PDF1(I)*PDF2(0)
       DLUMD=DLUMD + PDF1(0)*PDF2(-I) + PDF1(-I)*PDF2(0)
       J=I+1
       IF(J.LT.6)THEN
       DLUMU=DLUMU + PDF1(0)*PDF2(J) + PDF1(J)*PDF2(0)
       DLUMU=DLUMU + PDF1(0)*PDF2(-J) + PDF1(-J)*PDF2(0)
       ENDIF
1     CONTINUE
      DLUMGQZ=(CU*DLUMU+CD*DLUMD)/TAU/X
      RETURN
      END
 
      COMPLEX*16 FUNCTION CI(A,B,C)
C--SCALAR INTEGRAL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 C,CC,LI2
C     X=DREAL(C)
C     CC=DCMPLX(X,0.D0)
C     IF(X.LT.0.D0.OR.X.GT.1.D0)C=CC
      CC=C
      CI = CDLOG(A*CC+B)*CDLOG(1.D0-1.D0/C)
     .    -LI2(1.D0-(A+B)/(A*C+B))+LI2(1.D0-B/(A*C+B))
      RETURN
      END
 
      COMPLEX*16 FUNCTION CJ(A,B,C)
C--SCALAR INTEGRAL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 C,CB,AP,AM,LI2
      COMMON/CUTS/EPS,REPS
      CB = B*DCMPLX(1.D0,-REPS)
      C=DCMPLX(DREAL(C),DIMAG(C)/1.D3)
      CC = C
C     CC = DREAL(C)
C     IF(CC.LT.0.D0.OR.CC.GT.1.D0)C=DCMPLX(CC,0.D0)
      AP = 0.5D0*(1.D0+CDSQRT(1.D0-4.D0*CB/A))
      AM = 0.5D0*(1.D0-CDSQRT(1.D0-4.D0*CB/A))
      CJ = CDLOG(A*CC*(1.D0-CC)-CB)*CDLOG(1.D0-1.D0/C)
     .    -LI2((CC-1.D0)/(CC-AP))+LI2(CC/(CC-AP))
     .    -LI2((CC-1.D0)/(CC-AM))+LI2(CC/(CC-AM))
      RETURN
      END
 
      COMPLEX*16 FUNCTION LI2(X)
C--COMPLEX DILOGARITHM (SPENCE-FUNCTION)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 X,Y,CLI2
      COMMON/CONST/ZETA2
      ZERO=1.D-40
      XR=DREAL(X)
      XI=DIMAG(X)
      R2=XR*XR+XI*XI
      IF(R2.LE.ZERO)THEN
        LI2=X
        RETURN
      ENDIF
      RR=XR/R2
      IF(R2.EQ.1.D0.AND.XI.EQ.0.D0)THEN
        IF(XR.EQ.1.D0)THEN
          LI2=DCMPLX(ZETA2)
        ELSE
          LI2=-DCMPLX(ZETA2/2.D0)
        ENDIF
        RETURN
      ELSEIF(R2.GT.1.D0.AND.RR.GT.0.5D0)THEN
        Y=(X-1.D0)/X
        LI2=CLI2(Y)+ZETA2-CDLOG(X)*CDLOG(1.D0-X)+0.5D0*CDLOG(X)**2
        RETURN
      ELSEIF(R2.GT.1.D0.AND.RR.LE.0.5D0)THEN
        Y=1.D0/X
        LI2=-CLI2(Y)-ZETA2-0.5D0*CDLOG(-X)**2
        RETURN
      ELSEIF(R2.LE.1.D0.AND.XR.GT.0.5D0)THEN
        Y=1.D0-X
        LI2=-CLI2(Y)+ZETA2-CDLOG(X)*CDLOG(1.D0-X)
       RETURN
      ELSEIF(R2.LE.1.D0.AND.XR.LE.0.5D0)THEN
        Y=X
        LI2=CLI2(Y)
        RETURN
      ENDIF
      END
 
      COMPLEX*16 FUNCTION CLI2(X)
C--TAYLOR-EXPANSION FOR COMPLEX DILOGARITHM (SPENCE-FUNCTION)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 X,Z
      COMMON/BERNOULLI/B(18)
      COMMON/CONST/ZETA2
      COMMON/POLY/NBER
      N=NBER-1
      Z=-CDLOG(1.D0-X)
      CLI2=B(NBER)
      DO 111 I=N,1,-1
        CLI2=Z*CLI2+B(I)
111   CONTINUE
      CLI2=Z**2*CLI2+Z
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION FAC(N)
C--MATHEMATICAL FACULTY
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      FAC=1.D0
      IF(N.EQ.0) RETURN
      DO 986 I=1,N
        FAC=FAC*DFLOAT(I)
986   CONTINUE
      RETURN
      END
 
      SUBROUTINE BERNINI(N)
C--INITIALIZATION OF COEFFICIENTS FOR POLYLOGARITHMS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/BERNOULLI/B(18)
      COMMON/CONST/ZETA2
      COMMON/POLY/NBER
 
      NBER=N
      PI=4.D0*DATAN(1.D0)
 
      B(1)=-1.D0/2.D0
      B(2)=1.D0/6.D0
      B(3)=0.D0
      B(4)=-1.D0/30.D0
      B(5)=0.D0
      B(6)=1.D0/42.D0
      B(7)=0.D0
      B(8)=-1.D0/30.D0
      B(9)=0.D0
      B(10)=5.D0/66.D0
      B(11)=0.D0
      B(12)=-691.D0/2730.D0
      B(13)=0.D0
      B(14)=7.D0/6.D0
      B(15)=0.D0
      B(16)=-3617.D0/510.D0
      B(17)=0.D0
      B(18)=43867.D0/798.D0
      ZETA2=PI**2/6.D0
      DO 995 I=1,18
        B(I)=B(I)/FAC(I+1)
995   CONTINUE
      RETURN
      END
 
      SUBROUTINE HDEC(TGBET)
C--HIGGS DECAY WIDTHS AND BRANCHING RATIOS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION LAMB
      COMPLEX*16 CF,CG,CI1,CI2,CA,CB,CTT,CTB,CTW,CLT,CLB,CLW,
     .           CAT,CAB,CAW,CTH,CLH
      COMMON/HMASS/AMA,AML,AMH,AMCH
      COMMON/MASSES/AMS,AMC,AMB,AMT,XMZ,XMW
      COMMON/WDPARAM/GF,ALPH,AMTAU,AMMUON,AMZ,AMW
      COMMON/WIDTHS/SMWDTH,HLWDTH,HHWDTH,AWDTH,HCWDTH
      COMMON/COUP/YAT,YAB,YLT,YLB,YHT,YHB,YZAH,YZAL,
     .            YHHH,YLLL,YHLL,YLHH,YHAA,YLAA,YLVV,YHVV,
     .            YLPM,YHPM,YZT,YZB
      BETA(X)=DSQRT(1.D0-4.D0*X)
      LAMB(X,Y)=DSQRT((1.D0-X-Y)**2-4.D0*X*Y)
      HVV(X,Y)= GF/(4.D0*PI*DSQRT(2.D0))*X**3/2.D0*BETA(Y)
     .            *(1.D0-4.D0*Y+12.D0*Y**2)
      AFF(X,Y)= GF/(4.*PI*DSQRT(2.D0))*X**3*Y*(BETA(Y))
      HFF(X,Y)= GF/(4.*PI*DSQRT(2.D0))*X**3*Y*(BETA(Y))**3
      CFF(Z,TB,X,Y)= GF/(4.*PI*DSQRT(2.D0))*Z**3*LAMB(X,Y)
     .              *((1.D0-X-Y)*(X*TB**2+Y/TB**2)-4.D0*X*Y)
      HV(V)=3.D0*(1.D0-8.D0*V+20.D0*V**2)/DSQRT((4.D0*V-1.D0))
     .      *DACOS((3.D0*V-1.D0)/2.D0/DSQRT(V**3))
     .      -(1.D0-V)*(47.D0/2.D0*V-13.D0/2.D0+1.D0/V)
     .      -3.D0/2.D0*(1.D0-6.D0*V+4.D0*V**2)*DLOG(V)
      QCD(X)=X-13.D0*X*DLOG(X)
      CF(CA) = -CDLOG(-(1+CDSQRT(1-CA))/(1-CDSQRT(1-CA)))**2/4
      CG(CA) = CDSQRT(1-CA)/2*CDLOG(-(1+CDSQRT(1-CA))/(1-CDSQRT(1-CA)))
      CI1(CA,CB) = CA*CB/2/(CA-CB) 
     .           + CA**2*CB**2/2/(CA-CB)**2*(CF(CA)-CF(CB))
     .           + CA**2*CB/(CA-CB)**2*(CG(CA)-CG(CB))
      CI2(CA,CB) = -CA*CB/2/(CA-CB)*(CF(CA)-CF(CB))

      PI=4D0*DATAN(1D0)
      DLT=2.D0
      SS=1.D0-(AMW/AMZ)**2
      CS=1.D0-SS

      V = 1/DSQRT(DSQRT(2D0)*GF)
      GLT=YLT/AMT*V
      GLB=YLB/AMB*V
      GHT=YHT/AMT*V
      GHB=YHB/AMB*V
      GAT=YAT/AMT*V
      GAB=YAB/AMB*V
      GLLL=YLLL/AMZ**2*V
      GHHH=YHHH/AMZ**2*V
      GLHH=YLHH/AMZ**2*V
      GHLL=YHLL/AMZ**2*V
      GLAA=YLAA/AMZ**2*V
      GHAA=YHAA/AMZ**2*V
      GZAL=YZAL/AMZ*V
      GZAH=YZAH/AMZ*V
      GLVV=YLVV
      GHVV=YHVV
      GLPM=YLPM/AMZ**2*V
      GHPM=YHPM/AMZ**2*V

C        =========================================================
C                   SM HIGGS DECAYS
C        =========================================================
C     AMXX=AMH
C     AMH=AMSM
C     =============  RUNNING MASSES ================================
      RMS = RUNM(AMH,3)
      RMC = RUNM(AMH,4)
      RMB = RUNM(AMH,5)
C     =============== PARTIAL WIDTHS ================================
      ASH=ALPHAS(AMH,2)
C  H ---> MU MU
      HMM=HFF(AMH,(AMMUON/AMH)**2)
C  H ---> TAU TAU
      HLL=HFF(AMH,(AMTAU/AMH)**2)
C  H --> SS
      HSS=3.D0*HFF(AMH,(RMS/AMH)**2)*(1.D0+1.804D0*ASH+2.953D0*ASH**2)
     .*(1.D0+4.D0/3.D0*ASH/PI*QCD(RMS**2/AMH**2))
C  H --> CC
      HCC=3.D0*HFF(AMH,(RMC/AMH)**2)*(1.D0+1.804D0*ASH+2.953D0*ASH**2)
     .*(1.D0+4.D0/3.D0*ASH/PI*QCD(RMC**2/AMH**2))
C  H --> BB :
      HBB=3.D0*HFF(AMH,(RMB/AMH)**2)*(1.D0+1.804D0*ASH+2.953D0*ASH**2)
     .*(1.D0+4.D0/3.D0*ASH/PI*QCD(RMB**2/AMH**2))
C  H ---> TT
      IF (AMH.LE.2.D0*AMT) THEN
      HTT=0.D0
      ELSE
      HTT=3.D0*HFF(AMH,(AMT/AMH)**2)
      ENDIF
C  H ---> G G
       RB=AMH**2/4.D0/AMB**2
       RT=AMH**2/4.D0/AMT**2
       IF(AMH.LE.2.D0*AMT) THEN
        FTR=DASIN(DSQRT(RT))**2
        FTI=0.D0
       ELSE
        FTR=1.D0/4*(PI**2-DLOG((1+DSQRT(1-1/RT))
     .                      /(1-DSQRT(1-1/RT)))**2)
        FTI=1.D0/2*PI*DLOG((1+DSQRT(1-1/RT))
     .                      /(1-DSQRT(1-1/RT)))
       ENDIF
       FBR=1.D0/4*(PI**2-DLOG((1+DSQRT(1-1/RB))
     .                      /(1-DSQRT(1-1/RB)))**2)
       FBI=1.D0/2*PI*DLOG((1+DSQRT(1-1/RB))
     .                      /(1-DSQRT(1-1/RB)))
       AFTR=2.D0/RT*(1.D0+(RT-1.D0)/RT*FTR)
       AFTI=2.D0/RT*((RT-1.D0)/RT*FTI)
       AFBR=2.D0/RB*(1.D0+(RB-1.D0)/RB*FBR)
       AFBI=2.D0/RB*((RB-1.D0)/RB*FBI)
       FQCD=1.D0+ASH/PI*(95.D0/4.D0-35.D0/6.D0+2.D0)
       XFAC=((AFTR+AFBR)**2+(AFTI+AFBI)**2)*FQCD
       HGG=HVV(AMH,0.D0)*(ASH/PI)**2*XFAC/8
C  H ---> GAMMA GAMMA
       RW=AMH**2/4.D0/AMW**2
       IF(AMH.LE.2.D0*AMW) THEN
        FWR=DASIN(DSQRT(RW))**2
        FWI=0.D0
       ELSE
        FWR=1.D0/4*(PI**2-DLOG((1+DSQRT(1-1/RW))
     .                      /(1-DSQRT(1-1/RW)))**2)
        FWI=1.D0/2*PI*DLOG((1+DSQRT(1-1/RW))
     .                      /(1-DSQRT(1-1/RW)))
       ENDIF
       AXWR=-1.D0/RW**2*(2.D0*RW**2+3.D0*RW+3.D0*(2.D0*RW-1.D0)*FWR)
       AXWI=-1.D0/RW**2*(3.D0*(2.D0*RW-1.D0)*FWI)
       AXTR=2.D0/RT*(1.D0+(RT-1.D0)/RT*FTR)*4/3D0
       AXTI=2.D0/RT*((RT-1.D0)/RT*FTI)*4/3D0
       AXBR=2.D0/RB*(1.D0+(RB-1.D0)/RB*FBR)/3
       AXBI=2.D0/RB*((RB-1.D0)/RB*FBI)/3
       XFAC=(AXWR+AXTR+AXBR)**2+(AXWI+AXTI+AXBI)**2
       HGA=HVV(AMH,0.D0)*(ALPH/PI)**2/16.D0*XFAC
C  H ---> Z GAMMA
      IF(AMH.LE.AMZ)THEN
       HZGA=0
      ELSE
       TS = SS/CS
       FT = -3*2D0/3*(1-4*2D0/3*SS)/DSQRT(SS*CS)
       FB = 3*1D0/3*(-1+4*1D0/3*SS)/DSQRT(SS*CS)
       EPS=1.D-8
       CTT = 4*AMT**2/AMH**2*DCMPLX(1D0,-EPS)
       CTB = 4*AMB**2/AMH**2*DCMPLX(1D0,-EPS)
       CTW = 4*AMW**2/AMH**2*DCMPLX(1D0,-EPS)
       CLT = 4*AMT**2/AMZ**2*DCMPLX(1D0,-EPS)
       CLB = 4*AMB**2/AMZ**2*DCMPLX(1D0,-EPS)
       CLW = 4*AMW**2/AMZ**2*DCMPLX(1D0,-EPS)
       CAT = FT*(CI1(CTT,CLT) - CI2(CTT,CLT))
       CAB = FB*(CI1(CTB,CLB) - CI2(CTB,CLB))
       CAW = -1/DSQRT(TS)*(4*(3-TS)*CI2(CTW,CLW)
     .     + ((1+2/CTW)*TS - (5+2/CTW))*CI1(CTW,CLW))
       XFAC = CDABS(CAT+CAB+CAW)**2
       ACOUP = DSQRT(2D0)*GF*AMZ**2*SS*CS/PI**2
       HZGA = GF/(4.D0*PI*DSQRT(2.D0))*AMH**3*(ALPH/PI)*ACOUP/16.D0
     .        *XFAC*(1-AMZ**2/AMH**2)**3
      ENDIF
C  H ---> W W
      XM1 = 2D0*AMW-DLT
      XM2 = 2D0*AMW+DLT
      IF (AMH.LE.AMW) THEN
      HWW=0
      ELSE IF (AMH.LE.XM1) THEN
      CWW=3.D0*GF**2*AMW**4/16.D0/PI**3
      HWW=HV(AMW**2/AMH**2)*CWW*AMH
      ELSE IF (AMH.LT.XM2) THEN
      CWW=3.D0*GF**2*AMW**4/16.D0/PI**3
      X1=HV(AMW**2/XM1**2)*CWW*XM1
      X2=HVV(XM2,AMW**2/XM2**2)
      HWW = X1 + (AMH-XM1)/2D0/DLT*(X2-X1)
      ELSE
      HWW=HVV(AMH,AMW**2/AMH**2)
      ENDIF
C  H ---> Z Z
      XM1 = 2D0*AMZ-DLT
      XM2 = 2D0*AMZ+DLT
      IF (AMH.LE.AMZ) THEN
      HZZ=0
      ELSE IF (AMH.LE.XM1) THEN
      CZZ=3.D0*GF**2*AMZ**4/192.D0/PI**3*(7-40/3.D0*SS+160/9.D0*SS**2)
      HZZ=HV(AMZ**2/AMH**2)*CZZ*AMH
      ELSE IF (AMH.LT.XM2) THEN
      CZZ=3.D0*GF**2*AMZ**4/192.D0/PI**3*(7-40/3.D0*SS+160/9.D0*SS**2)
      X1=HV(AMZ**2/XM1**2)*CZZ*XM1
      X2=HVV(XM2,AMZ**2/XM2**2)/2.D0
      HZZ = X1 + (AMH-XM1)/2D0/DLT*(X2-X1)
      ELSE
      HZZ=HVV(AMH,AMZ**2/AMH**2)/2.D0
      ENDIF
C    ==========  TOTAL WIDTH AND BRANCHING RATIOS ==================
      WTOT=HLL+HMM+HSS+HCC+HBB+HTT+HGG+HGA+HZGA+HWW+HZZ
      SMBRT=HTT/WTOT
      SMBRB=HBB/WTOT
      SMBRL=HLL/WTOT
      SMBRM=HMM/WTOT
      SMBRC=HCC/WTOT
      SMBRS=HSS/WTOT
      SMBRG=HGG/WTOT
      SMBRGA=HGA/WTOT
      SMBRZGA=HZGA/WTOT
      SMBRW=HWW/WTOT
      SMBRZ=HZZ/WTOT
      SMWDTH=WTOT

C     AMH=AMXX

C        =========================================================
C                   CP ODD  HIGGS DECAYS
C        =========================================================
C     =============  RUNNING MASSES ================================
      RMS = RUNM(AMA,3)
      RMC = RUNM(AMA,4)
      RMB = RUNM(AMA,5)
C     =============== PARTIAL WIDTHS ================================
      ASH=ALPHAS(AMA,2)
C  A ---> MU MU
      HMM=AFF(AMA,(AMMUON/AMA)**2)*GAB**2
C  A ---> LL
      HLL=AFF(AMA,(AMTAU/AMA)**2)*GAB**2
C  A --> SS
      HSS=3.D0*AFF(AMA,(RMS/AMA)**2)*(1.D0+1.804D0*ASH+2.953D0*ASH**2)
     .    *GAB**2
C  A --> CC
      HCC=3.D0*AFF(AMA,(RMC/AMA)**2)*(1.D0+1.804D0*ASH+2.953D0*ASH**2)
     .    *GAT**2
C  A --> BB :
      HBB=3.D0*AFF(AMA,(RMB/AMA)**2)*(1.D0+1.804D0*ASH+2.953D0*ASH**2)
     .    *GAB**2
C  A --> TT :
      IF(AMA.LE.2.D0*AMT) THEN
      HTT=0.D0
      ELSE
      HTT=3.D0*AFF(AMA,(AMT/AMA)**2)*GAT**2
      ENDIF
C  A ---> G G 
       RB=AMA**2/4.D0/AMB**2
       RT=AMA**2/4.D0/AMT**2
       IF(AMA.LE.2.D0*AMT) THEN
        FTR=DASIN(DSQRT(RT))**2
        FTI=0.D0
       ELSE
        FTR=1.D0/4*(PI**2-DLOG((1+DSQRT(1-1/RT))
     .                      /(1-DSQRT(1-1/RT)))**2)
        FTI=1.D0/2*PI*DLOG((1+DSQRT(1-1/RT))
     .                      /(1-DSQRT(1-1/RT)))
       ENDIF
       FBR=1.D0/4*(PI**2-DLOG((1+DSQRT(1-1/RB))
     .                      /(1-DSQRT(1-1/RB)))**2)
       FBI=1.D0/2*PI*DLOG((1+DSQRT(1-1/RB))
     .                      /(1-DSQRT(1-1/RB)))
       AFTR=FTR/RT*GAT
       AFTI=FTI/RT*GAT
       AFBR=FBR/RB*GAB
       AFBI=FBI/RB*GAB
       FQCD=1.D0+ASH/PI*(97.D0/4.D0-35.D0/6.D0+1.75D0)
       XFAC=((AFTR+AFBR)**2+(AFTI+AFBI)**2)*FQCD
       HGG=GF/(16.D0*PI*DSQRT(2.D0))*AMA**3*(ASH/PI)**2*XFAC
C  A ---> GAMMA GAMMA
       AXR = 4D0/3*AFTR + AFBR/3
       AXI = 4D0/3*AFTI + AFBI/3
       XFAC=AXR**2+AXI**2
       HGA=GF/(32.D0*PI*DSQRT(2.D0))*AMA**3*(ALPH/PI)**2*XFAC
C  A ---> Z GAMMA
      IF(AMA.LE.AMZ)THEN
       HZGA=0
      ELSE
       TS = SS/CS
       FT = -3*2D0/3*(1-4*2D0/3*SS)/DSQRT(SS*CS)*GAT
       FB = 3*1D0/3*(-1+4*1D0/3*SS)/DSQRT(SS*CS)*GAB
       EPS=1.D-8
       CTT = 4*AMT**2/AMA**2*DCMPLX(1D0,-EPS)
       CTB = 4*AMB**2/AMA**2*DCMPLX(1D0,-EPS)
       CLT = 4*AMT**2/AMZ**2*DCMPLX(1D0,-EPS)
       CLB = 4*AMB**2/AMZ**2*DCMPLX(1D0,-EPS)
       CAT = FT*(- CI2(CTT,CLT))
       CAB = FB*(- CI2(CTB,CLB))
       XFAC = CDABS(CAT+CAB)**2
       ACOUP = DSQRT(2D0)*GF*AMZ**2*SS*CS/PI**2
       HZGA = GF/(4.D0*PI*DSQRT(2.D0))*AMA**3*(ALPH/PI)*ACOUP/16.D0
     .        *XFAC*(1-AMZ**2/AMA**2)**3
      ENDIF
C  A ---> H Z* ---> HFF
      IF (AMA.GE.AMZ+AML) THEN
      CAZ=LAMB(AML**2/AMA**2,AMZ**2/AMA**2)
     .   *LAMB(AMA**2/AMZ**2,AML**2/AMZ**2)**2
      HAZ=GF/8D0/DSQRT(2D0)/PI*AMZ**4/AMA*GZAL**2*CAZ
      ELSE
      HAZ=0.D0
      ENDIF
C    ==========  TOTAL WIDTH AND BRANCHING RATIOS ==================
      WTOT=HLL+HMM+HSS+HCC+HBB+HGG+HGA+HZGA+HAZ+HTT
      ABRT=HTT/WTOT
      ABRB=HBB/WTOT
      ABRL=HLL/WTOT
      ABRM=HMM/WTOT
      ABRS=HSS/WTOT
      ABRC=HCC/WTOT
      ABRG=HGG/WTOT
      ABRGA=HGA/WTOT
      ABRZGA=HZGA/WTOT
      ABRZ=HAZ/WTOT
      AWDTH=WTOT

C        =========================================================
C                   HEAVY CP EVEN HIGGS DECAYS
C        =========================================================
C     =============  RUNNING MASSES ================================
      RMS = RUNM(AMH,3)
      RMC = RUNM(AMH,4)
      RMB = RUNM(AMH,5)
C     =============== PARTIAL WIDTHS ================================
      ASH=ALPHAS(AMH,2)
C  H ---> MU MU
      HMM=HFF(AMH,(AMMUON/AMH)**2)*GHB**2
C  H ---> LL
      HLL=HFF(AMH,(AMTAU/AMH)**2)*GHB**2
C  H --> SS
      HSS=3.D0*HFF(AMH,(RMS/AMH)**2)*(1.D0+1.804D0*ASH+2.953D0*ASH**2)
     .*(1.D0+4.D0/3.D0*ASH/PI*QCD(RMS**2/AMH**2))*GHB**2
C  H --> CC
      HCC=3.D0*HFF(AMH,(RMC/AMH)**2)*(1.D0+1.804D0*ASH+2.953D0*ASH**2)
     .*(1.D0+4.D0/3.D0*ASH/PI*QCD(RMC**2/AMH**2))*GHT**2
C  H --> BB :
      HBB=3.D0*HFF(AMH,(RMB/AMH)**2)*(1.D0+1.804D0*ASH+2.953D0*ASH**2)
     .*(1.D0+4.D0/3.D0*ASH/PI*QCD(RMB**2/AMH**2))*GHB**2
C  H ---> TT
      IF (AMH.LE.2.D0*AMT) THEN
      HTT=0.D0
      ELSE
      HTT=3.D0*HFF(AMH,(AMT/AMH)**2)*GHT**2
      ENDIF
C  H ---> G G
       RB=AMH**2/4.D0/AMB**2
       RT=AMH**2/4.D0/AMT**2
       IF(AMH.LE.2.D0*AMT) THEN
        FTR=DASIN(DSQRT(RT))**2
        FTI=0.D0
       ELSE
        FTR=1.D0/4*(PI**2-DLOG((1+DSQRT(1-1/RT))
     .                      /(1-DSQRT(1-1/RT)))**2)
        FTI=1.D0/2*PI*DLOG((1+DSQRT(1-1/RT))
     .                      /(1-DSQRT(1-1/RT)))
       ENDIF
       FBR=1.D0/4*(PI**2-DLOG((1+DSQRT(1-1/RB))
     .                      /(1-DSQRT(1-1/RB)))**2)
       FBI=1.D0/2*PI*DLOG((1+DSQRT(1-1/RB))
     .                      /(1-DSQRT(1-1/RB)))
       AFTR=2.D0/RT*(1.D0+(RT-1.D0)/RT*FTR)*GHT
       AFTI=2.D0/RT*((RT-1.D0)/RT*FTI)*GHT
       AFBR=2.D0/RB*(1.D0+(RB-1.D0)/RB*FBR)*GHB
       AFBI=2.D0/RB*((RB-1.D0)/RB*FBI)*GHB
       FQCD=1.D0+ASH/PI*(95.D0/4.D0-35.D0/6.D0+2.D0)
       XFAC=((AFTR+AFBR)**2+(AFTI+AFBI)**2)*FQCD
       HGG=HVV(AMH,0.D0)*(ASH/PI)**2*XFAC/8
C  H ---> GAMMA GAMMA
       RW=AMH**2/4.D0/AMW**2
       IF(AMH.LE.2.D0*AMW) THEN
        FWR=DASIN(DSQRT(RW))**2
        FWI=0.D0
       ELSE
        FWR=1.D0/4*(PI**2-DLOG((1+DSQRT(1-1/RW))
     .                      /(1-DSQRT(1-1/RW)))**2)
        FWI=1.D0/2*PI*DLOG((1+DSQRT(1-1/RW))
     .                      /(1-DSQRT(1-1/RW)))
       ENDIF
       AXWR=-1.D0/RW**2*(2.D0*RW**2+3.D0*RW+3.D0*(2.D0*RW-1.D0)*FWR)
     .     *GHVV
       AXWI=-1.D0/RW**2*(3.D0*(2.D0*RW-1.D0)*FWI)*GHVV
       AXTR=2.D0/RT*(1.D0+(RT-1.D0)/RT*FTR)*GHT*4/3D0
       AXTI=2.D0/RT*((RT-1.D0)/RT*FTI)*GHT*4/3D0
       AXBR=2.D0/RB*(1.D0+(RB-1.D0)/RB*FBR)*GHB/3
       AXBI=2.D0/RB*((RB-1.D0)/RB*FBI)*GHB/3
       XFAC=(AXWR+AXTR+AXBR)**2+(AXWI+AXTI+AXBI)**2
       HGA=HVV(AMH,0.D0)*(ALPH/PI)**2/16.D0*XFAC
C  H ---> Z GAMMA
      IF(AMH.LE.AMZ)THEN
       HZGA=0
      ELSE
       TS = SS/CS
       FT = -3*2D0/3*(1-4*2D0/3*SS)/DSQRT(SS*CS)*GHT
       FB = 3*1D0/3*(-1+4*1D0/3*SS)/DSQRT(SS*CS)*GHB
       EPS=1.D-8
       CTT = 4*AMT**2/AMH**2*DCMPLX(1D0,-EPS)
       CTB = 4*AMB**2/AMH**2*DCMPLX(1D0,-EPS)
       CTW = 4*AMW**2/AMH**2*DCMPLX(1D0,-EPS)
       CTH = 4*AMCH**2/AMH**2*DCMPLX(1D0,-EPS)
       CLT = 4*AMT**2/AMZ**2*DCMPLX(1D0,-EPS)
       CLB = 4*AMB**2/AMZ**2*DCMPLX(1D0,-EPS)
       CLW = 4*AMW**2/AMZ**2*DCMPLX(1D0,-EPS)
       CLH = 4*AMCH**2/AMZ**2*DCMPLX(1D0,-EPS)
       CAT = FT*(CI1(CTT,CLT) - CI2(CTT,CLT))
       CAB = FB*(CI1(CTB,CLB) - CI2(CTB,CLB))
       CAW = -1/DSQRT(TS)*(4*(3-TS)*CI2(CTW,CLW)
     .     + ((1+2/CTW)*TS - (5+2/CTW))*CI1(CTW,CLW))*GHVV
       CAH = (1-2*SS)/DSQRT(SS*CS)*AMW**2/AMCH**2*CI1(CTH,CLH)*GHPM
       XFAC = CDABS(CAT+CAB+CAW+CAH)**2
       ACOUP = DSQRT(2D0)*GF*AMZ**2*SS*CS/PI**2
       HZGA = GF/(4.D0*PI*DSQRT(2.D0))*AMH**3*(ALPH/PI)*ACOUP/16.D0
     .        *XFAC*(1-AMZ**2/AMH**2)**3
      ENDIF
C  H ---> W W
      XM1 = 2D0*AMW-DLT
      XM2 = 2D0*AMW+DLT
      IF (AMH.LE.AMW) THEN
      HWW=0
      ELSE IF (AMH.LE.XM1) THEN
      CWW=3.D0*GF**2*AMW**4/16.D0/PI**3
      HWW=HV(AMW**2/AMH**2)*CWW*AMH*GHVV**2
      ELSE IF (AMH.LT.XM2) THEN
      CWW=3.D0*GF**2*AMW**4/16.D0/PI**3
      X1=HV(AMW**2/XM1**2)*CWW*XM1*GHVV**2
      X2=HVV(XM2,AMW**2/XM2**2)*GHVV**2
      HWW = X1 + (AMH-XM1)/2D0/DLT*(X2-X1)
      ELSE
      HWW=HVV(AMH,AMW**2/AMH**2)*GHVV**2
      ENDIF
C  H ---> Z Z
      XM1 = 2D0*AMZ-DLT
      XM2 = 2D0*AMZ+DLT
      IF (AMH.LE.AMZ) THEN
      HZZ=0
      ELSE IF (AMH.LE.XM1) THEN
      CZZ=3.D0*GF**2*AMZ**4/192.D0/PI**3*(7-40/3.D0*SS+160/9.D0*SS**2)
      HZZ=HV(AMZ**2/AMH**2)*CZZ*AMH*GHVV**2
      ELSE IF (AMH.LT.XM2) THEN
      CZZ=3.D0*GF**2*AMZ**4/192.D0/PI**3*(7-40/3.D0*SS+160/9.D0*SS**2)
      X1=HV(AMZ**2/XM1**2)*CZZ*XM1*GHVV**2
      X2=HVV(XM2,AMZ**2/XM2**2)/2.D0*GHVV**2
      HZZ = X1 + (AMH-XM1)/2D0/DLT*(X2-X1)
      ELSE
      HZZ=HVV(AMH,AMZ**2/AMH**2)/2.D0*GHVV**2
      ENDIF
C  H ---> h h
      IF (AMH.LE.2.D0*AML) THEN
      HHH=0
      ELSE
      HHH=GF/16D0/DSQRT(2D0)/PI*AMZ**4/AMH*BETA(AML**2/AMH**2)*GHLL**2
      ENDIF
C  H ---> A A
      IF (AMH.LE.2.D0*AMA) THEN
      HAA=0
      ELSE
      HAA=GF/16.D0/DSQRT(2D0)/PI*AMZ**4/AMH*BETA(AMA**2/AMH**2)*GHAA**2
      ENDIF
C  H ---> A Z
      IF (AMH.LE.AMZ+AMA) THEN
      HAZ=0
      ELSE
      CAZ=LAMB(AMA**2/AMH**2,AMZ**2/AMH**2)
     .   *LAMB(AMH**2/AMZ**2,AMA**2/AMZ**2)**2
      HAZ=GF/8.D0/DSQRT(2D0)/PI*AMZ**4/AMH*CAZ*GZAH**2
      ENDIF
C  H ---> H+ W+ 
      IF (AMH.LE.AMW+AMCH) THEN
      HHW=0
      ELSE
      CHW=LAMB(AMCH**2/AMH**2,AMW**2/AMH**2)
     .   *LAMB(AMH**2/AMW**2,AMCH**2/AMW**2)**2
      HHW=GF/8.D0/DSQRT(2D0)/PI*AMZ**2*AMW**2/AMH*CHW*GLVV**2
      ENDIF
C    ==========  TOTAL WIDTH AND BRANCHING RATIOS ==================
      WTOT=HLL+HMM+HSS+HCC+HBB+HTT+HGG+HGA+HZGA+HWW+HZZ+HHH+HAA+HAZ+HHW
      HHBRT=HTT/WTOT
      HHBRB=HBB/WTOT
      HHBRL=HLL/WTOT
      HHBRM=HMM/WTOT
      HHBRS=HSS/WTOT
      HHBRC=HCC/WTOT
      HHBRG=HGG/WTOT
      HHBRGA=HGA/WTOT
      HHBRZGA=HZGA/WTOT
      HHBRW=HWW/WTOT
      HHBRZ=HZZ/WTOT
      HHBRH=HHH/WTOT
      HHBRA=HAA/WTOT
      HHBRAZ=HAZ/WTOT
      HHBRHW=HHW/WTOT
      HHWDTH=WTOT

C        =========================================================
C                   LIGHT CP EVEN HIGGS DECAYS
C        =========================================================
C     =============  RUNNING MASSES ================================
      RMS = RUNM(AML,3)
      RMC = RUNM(AML,4)
      RMB = RUNM(AML,5)
C     =============== PARTIAL WIDTHS ================================
      ASH=ALPHAS(AML,2)
C  H ---> MU MU
      HMM=HFF(AML,(AMMUON/AML)**2)*GLB**2
C  H ---> TAU TAU
      HLL=HFF(AML,(AMTAU/AML)**2)*GLB**2
C  H --> SS
      HSS=3.D0*HFF(AML,(RMS/AML)**2)*(1.D0+1.804D0*ASH+2.953D0*ASH**2)
     .*(1.D0+4.D0/3.D0*ASH/PI*QCD(RMS**2/AML**2))*GLB**2
C  H --> CC
      HCC=3.D0*HFF(AML,(RMC/AML)**2)*(1.D0+1.804D0*ASH+2.953D0*ASH**2)
     .*(1.D0+4.D0/3.D0*ASH/PI*QCD(RMC**2/AML**2))*GLT**2
C  H --> BB :
      HBB=3.D0*HFF(AML,(RMB/AML)**2)*(1.D0+1.804D0*ASH+2.953D0*ASH**2)
     .*(1.D0+4.D0/3.D0*ASH/PI*QCD(RMB**2/AML**2))*GLB**2
C  H ---> TT
      IF (AML.LE.2.D0*AMT) THEN
      HTT=0.D0
      ELSE
      HTT=3.D0*HFF(AML,(AMT/AML)**2)*GLT**2
      ENDIF
C  H ---> G G
       RB=AML**2/4.D0/AMB**2
       RT=AML**2/4.D0/AMT**2
       IF(AML.LE.2.D0*AMT) THEN
        FTR=DASIN(DSQRT(RT))**2
        FTI=0.D0
       ELSE
        FTR=1.D0/4*(PI**2-DLOG((1+DSQRT(1-1/RT))
     .                      /(1-DSQRT(1-1/RT)))**2)
        FTI=1.D0/2*PI*DLOG((1+DSQRT(1-1/RT))
     .                      /(1-DSQRT(1-1/RT)))
       ENDIF
       FBR=1.D0/4*(PI**2-DLOG((1+DSQRT(1-1/RB))
     .                      /(1-DSQRT(1-1/RB)))**2)
       FBI=1.D0/2*PI*DLOG((1+DSQRT(1-1/RB))
     .                      /(1-DSQRT(1-1/RB)))
       AFTR=2.D0/RT*(1.D0+(RT-1.D0)/RT*FTR)*GLT
       AFTI=2.D0/RT*((RT-1.D0)/RT*FTI)*GLT
       AFBR=2.D0/RB*(1.D0+(RB-1.D0)/RB*FBR)*GLB
       AFBI=2.D0/RB*((RB-1.D0)/RB*FBI)*GLB
       FQCD=1.D0+ASH/PI*(95.D0/4.D0-35.D0/6.D0)
       XFAC=((AFTR+AFBR)**2+(AFTI+AFBI)**2)*FQCD
       HGG=HVV(AML,0.D0)*(ASH/PI)**2*XFAC/8
C  H ---> GAMMA GAMMA
       RW=AML**2/4.D0/AMW**2
       IF(AML.LE.2.D0*AMW) THEN
        FWR=DASIN(DSQRT(RW))**2
        FWI=0.D0
       ELSE
        FWR=1.D0/4*(PI**2-DLOG((1+DSQRT(1-1/RW))
     .                      /(1-DSQRT(1-1/RW)))**2)
        FWI=1.D0/2*PI*DLOG((1+DSQRT(1-1/RW))
     .                      /(1-DSQRT(1-1/RW)))
       ENDIF
       AXWR=-1.D0/RW**2*(2.D0*RW**2+3.D0*RW+3.D0*(2.D0*RW-1.D0)*FWR)
     .     *GLVV
       AXWI=-1.D0/RW**2*(3.D0*(2.D0*RW-1.D0)*FWI)*GLVV
       AXTR=2.D0/RT*(1.D0+(RT-1.D0)/RT*FTR)*GLT*4/3D0
       AXTI=2.D0/RT*((RT-1.D0)/RT*FTI)*GLT*4/3D0
       AXBR=2.D0/RB*(1.D0+(RB-1.D0)/RB*FBR)*GLB/3
       AXBI=2.D0/RB*((RB-1.D0)/RB*FBI)*GLB/3
       XFAC=(AXWR+AXTR+AXBR)**2+(AXWI+AXTI+AXBI)**2
       HGA=HVV(AML,0.D0)*(ALPH/PI)**2/16.D0*XFAC
C  H ---> Z GAMMA
      IF(AML.LE.AMZ)THEN
       HZGA=0
      ELSE
       TS = SS/CS
       FT = -3*2D0/3*(1-4*2D0/3*SS)/DSQRT(SS*CS)*GLT
       FB = 3*1D0/3*(-1+4*1D0/3*SS)/DSQRT(SS*CS)*GLB
       EPS=1.D-8
       CTT = 4*AMT**2/AML**2*DCMPLX(1D0,-EPS)
       CTB = 4*AMB**2/AML**2*DCMPLX(1D0,-EPS)
       CTW = 4*AMW**2/AML**2*DCMPLX(1D0,-EPS)
       CTH = 4*AMCH**2/AML**2*DCMPLX(1D0,-EPS)
       CLT = 4*AMT**2/AMZ**2*DCMPLX(1D0,-EPS)
       CLB = 4*AMB**2/AMZ**2*DCMPLX(1D0,-EPS)
       CLW = 4*AMW**2/AMZ**2*DCMPLX(1D0,-EPS)
       CLH = 4*AMCH**2/AMZ**2*DCMPLX(1D0,-EPS)
       CAT = FT*(CI1(CTT,CLT) - CI2(CTT,CLT))
       CAB = FB*(CI1(CTB,CLB) - CI2(CTB,CLB))
       CAW = -1/DSQRT(TS)*(4*(3-TS)*CI2(CTW,CLW)
     .     + ((1+2/CTW)*TS - (5+2/CTW))*CI1(CTW,CLW))*GLVV
       CAH = (1-2*SS)/DSQRT(SS*CS)*AMW**2/AMCH**2*CI1(CTH,CLH)*GLPM
       XFAC = CDABS(CAT+CAB+CAW+CAH)**2
       ACOUP = DSQRT(2D0)*GF*AMZ**2*SS*CS/PI**2
       HZGA = GF/(4.D0*PI*DSQRT(2.D0))*AML**3*(ALPH/PI)*ACOUP/16.D0
     .        *XFAC*(1-AMZ**2/AML**2)**3
      ENDIF
C  H ---> W W
      XM1 = 2D0*AMW-DLT
      XM2 = 2D0*AMW+DLT
      IF (AML.LE.AMW) THEN
      HWW=0
      ELSE IF (AML.LE.XM1) THEN
      CWW=3.D0*GF**2*AMW**4/16.D0/PI**3
      HWW=HV(AMW**2/AML**2)*CWW*AML*GLVV**2
      ELSE IF (AML.LT.XM2) THEN
      CWW=3.D0*GF**2*AMW**4/16.D0/PI**3
      X1=HV(AMW**2/XM1**2)*CWW*XM1*GLVV**2
      X2=HVV(XM2,AMW**2/XM2**2)*GLVV**2
      HWW = X1 + (AML-XM1)/2D0/DLT*(X2-X1)
      ELSE
      HWW=HVV(AML,AMW**2/AML**2)*GLVV**2
      ENDIF
C  H ---> Z Z
      XM1 = 2D0*AMZ-DLT
      XM2 = 2D0*AMZ+DLT
      IF (AML.LE.AMZ) THEN
      HZZ=0
      ELSE IF (AML.LE.XM1) THEN
      CZZ=3.D0*GF**2*AMZ**4/192.D0/PI**3*(7-40/3.D0*SS+160/9.D0*SS**2)
      HZZ=HV(AMZ**2/AML**2)*CZZ*AML*GLVV**2
      ELSE IF (AML.LT.XM2) THEN
      CZZ=3.D0*GF**2*AMZ**4/192.D0/PI**3*(7-40/3.D0*SS+160/9.D0*SS**2)
      X1=HV(AMZ**2/XM1**2)*CZZ*XM1*GLVV**2
      X2=HVV(XM2,AMZ**2/XM2**2)/2.D0*GLVV**2
      HZZ = X1 + (AML-XM1)/2D0/DLT*(X2-X1)
      ELSE
      HZZ=HVV(AML,AMZ**2/AML**2)/2.D0*GLVV**2
      ENDIF
C  H ---> A A
      IF (AML.LE.2.D0*AMA) THEN
      HAA=0
      ELSE
      HAA=GF/16.D0/DSQRT(2D0)/PI*AMZ**4/AML*BETA(AMA**2/AML**2)*GLAA**2
      ENDIF
C  H ---> A Z
      IF (AML.LE.AMZ+AMA) THEN
      HAZ=0
      ELSE
      CAZ=LAMB(AMA**2/AML**2,AMZ**2/AML**2)
     .   *LAMB(AML**2/AMZ**2,AMA**2/AMZ**2)**2
      HAZ=GF/8.D0/DSQRT(2D0)/PI*AMZ**4/AML*CAZ*GZAL**2
      ENDIF
C  H ---> H+ W+ 
      IF (AML.LE.AMW+AMCH) THEN
      HHW=0
      ELSE
      CHW=LAMB(AMCH**2/AML**2,AMW**2/AML**2)
     .   *LAMB(AML**2/AMW**2,AMCH**2/AMW**2)**2
      HHW=GF/8.D0/DSQRT(2D0)/PI*AMZ**2*AMW**2/AML*CHW*GHVV**2
      ENDIF
C    ==========  TOTAL WIDTH AND BRANCHING RATIOS ==================
      WTOT=HLL+HMM+HSS+HCC+HBB+HTT+HGG+HGA+HZGA+HWW+HZZ+HAA+HAZ+HHW
      HLBRT=HTT/WTOT
      HLBRB=HBB/WTOT
      HLBRL=HLL/WTOT
      HLBRM=HMM/WTOT
      HLBRS=HSS/WTOT
      HLBRC=HCC/WTOT
      HLBRG=HGG/WTOT
      HLBRGA=HGA/WTOT
      HLBRZGA=HZGA/WTOT
      HLBRW=HWW/WTOT
      HLBRZ=HZZ/WTOT
      HLBRA=HAA/WTOT
      HLBRAZ=HAZ/WTOT
      HLBRHW=HHW/WTOT
      HLWDTH=WTOT

C        =========================================================
C                CHARGED  HIGGS DECAYS
C        =========================================================
      TB=TGBET
C     =============  RUNNING MASSES ================================
      RMS = RUNM(AMCH,3)
      RMC = RUNM(AMCH,4)
      RMB = RUNM(AMCH,5)
C     =============== PARTIAL WIDTHS ================================
      ASH=ALPHAS(AMCH,2)
C  H+ ---> MU NMU
      HMN=CFF(AMCH,TB,(AMMUON/AMCH)**2,0.D0)
C  H+ ---> TAU NTAU
      HLN=CFF(AMCH,TB,(AMTAU/AMCH)**2,0.D0)
C  H+ --> SU
      VUS=0.2205D0
      HSU=3.D0*VUS**2*CFF(AMCH,TB,(RMS/AMCH)**2,0D0)
C  H+ --> CS
      HSC=3.D0*CFF(AMCH,TB,(RMS/AMCH)**2,(RMC/AMCH)**2)
C  H+ --> CB
      VCB=0.04D0
      HBC=3.D0*VCB**2*CFF(AMCH,TB,(RMB/AMCH)**2,(RMC/AMCH)**2)
C  H+ --> TB :
      IF (AMCH.LE.AMT+AMB) THEN
      HBT=0
      ELSE
      HBT=3.D0*CFF(AMCH,TB,(RMB/AMCH)**2,(AMT/AMCH)**2)
      ENDIF
C  H+ ---> W H
      IF (AMCH.LE.AMW+AML) THEN
      HWH=0
      ELSE
      CWH=LAMB(AML**2/AMCH**2,AMW**2/AMCH**2)
     .   *LAMB(AMCH**2/AMW**2,AML**2/AMW**2)**2
      HWH=GF/8.D0/DSQRT(2D0)/PI*AMZ**2*AMW**2/AMCH*GHVV**2*CWH
      ENDIF
C    ==========  TOTAL WIDTH AND BRANCHING RATIOS ==================
      WTOT=HLN+HMN+HSU+HSC+HBC+HBT+HWH
      HCBRL=HLN/WTOT
      HCBRM=HMN/WTOT
      HCBRS=HSU/WTOT
      HCBRC=HSC/WTOT
      HCBRB=HBC/WTOT
      HCBRT=HBT/WTOT
      HCBRW=HWH/WTOT
      HCWDTH=WTOT

      RETURN
      END

      DOUBLE PRECISION FUNCTION RUNM(Q,NF)
C--RUNNING QUARK MASSES: Q = SCALE, NF = 3,..,6 QUARK TYPE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NN=6)
      PARAMETER (ZETA3 = 1.202056903159594D0)
      DIMENSION AM(NN),YMSB(NN)
      COMMON/ALS/XLAMBDA,AMCA,AMBA,AMTA,N0A
      COMMON/MASSES/AMS,AMC,AMB,AMT,AMZ,AMW
      COMMON/STRANGE/AMSB
      SAVE ISTRANGE
      B0(NF)=(33.D0-2.D0*NF)/12D0
      B1(NF) = (102D0-38D0/3D0*NF)/16D0
      B2(NF) = (2857D0/2D0-5033D0/18D0*NF+325D0/54D0*NF**2)/64D0
      G0(NF) = 1D0
      G1(NF) = (202D0/3D0-20D0/9D0*NF)/16D0
      G2(NF) = (1249D0-(2216D0/27D0+160D0/3D0*ZETA3)*NF
     .       - 140D0/81D0*NF**2)/64D0
      C1(NF) = G1(NF)/B0(NF) - B1(NF)*G0(NF)/B0(NF)**2
      C2(NF) = ((G1(NF)/B0(NF) - B1(NF)*G0(NF)/B0(NF)**2)**2
     .       + G2(NF)/B0(NF) + B1(NF)**2*G0(NF)/B0(NF)**3
     .       - B1(NF)*G1(NF)/B0(NF)**2 - B2(NF)*G0(NF)/B0(NF)**2)/2D0
      TRAN(X,XK)=1D0+4D0/3D0*ALPHAS(X,2)/PI+XK*(ALPHAS(X,2)/PI)**2
      CQ(X,NF)=(2D0*B0(NF)*X)**(G0(NF)/B0(NF))
     .            *(1D0+C1(NF)*X+C2(NF)*X**2)
      DATA ISTRANGE/0/
      PI=4D0*DATAN(1D0)
      ACC = 1.D-8
      AM(1) = 0
      AM(2) = 0
C--------------------------------------------
      IMSBAR = 0
      NNLO = 0
      IF(IMSBAR.EQ.1)THEN
       IF(ISTRANGE.EQ.0)THEN
C--STRANGE POLE MASS FROM MSBAR-MASS AT 1 GEV
        AMSD = XLAMBDA
        AMSU = 1.D8
123     AMS  = (AMSU+AMSD)/2
        AM(3) = AMS
        XMSB = AMS/CQ(ALPHAS(AMS,2)/PI,3)
     .            *CQ(ALPHAS(1.D0,2)/PI,3)/TRAN(AMS,0D0)
        DD = (XMSB-AMSB)/AMSB
        IF(DABS(DD).GE.ACC)THEN
         IF(DD.LE.0.D0)THEN
          AMSD = AM(3)
         ELSE
          AMSU = AM(3)
         ENDIF
         GOTO 123
        ENDIF
        ISTRANGE=1
       ENDIF
       AM(3) = AMSB
      ELSE
       AMS=AMSB
       AM(3) = AMS
      ENDIF
C--------------------------------------------
      AM(3) = AMSB
      AM(4) = AMC
      AM(5) = AMB
      AM(6) = AMT
      XK = 16.11D0
      DO 1 I=1,NF-1
       XK = XK - 1.04D0*(1.D0-AM(I)/AM(NF))
1     CONTINUE
      IF(NF.GE.4)THEN
       XMSB = AM(NF)/TRAN(AM(NF),0D0)
       XMHAT = XMSB/CQ(ALPHAS(AM(NF),2)/PI,NF)
      ELSE
       XMSB = 0
       XMHAT = 0
      ENDIF
      YMSB(3) = AMSB
      IF(NF.EQ.3)THEN
       YMSB(4) = YMSB(3)*CQ(ALPHAS(AM(4),2)/PI,3)/
     .                   CQ(ALPHAS(1.D0,2)/PI,3)
       YMSB(5) = YMSB(4)*CQ(ALPHAS(AM(5),2)/PI,4)/
     .                   CQ(ALPHAS(AM(4),2)/PI,4)
       YMSB(6) = YMSB(5)*CQ(ALPHAS(AM(6),2)/PI,5)/
     .                   CQ(ALPHAS(AM(5),2)/PI,5)
      ELSEIF(NF.EQ.4)THEN
       YMSB(4) = XMSB
       YMSB(5) = YMSB(4)*CQ(ALPHAS(AM(5),2)/PI,4)/
     .                   CQ(ALPHAS(AM(4),2)/PI,4)
       YMSB(6) = YMSB(5)*CQ(ALPHAS(AM(6),2)/PI,5)/
     .                   CQ(ALPHAS(AM(5),2)/PI,5)
      ELSEIF(NF.EQ.5)THEN
       YMSB(5) = XMSB
       YMSB(4) = YMSB(5)*CQ(ALPHAS(AM(4),2)/PI,4)/
     .                   CQ(ALPHAS(AM(5),2)/PI,4)
       YMSB(6) = YMSB(5)*CQ(ALPHAS(AM(6),2)/PI,5)/
     .                   CQ(ALPHAS(AM(5),2)/PI,5)
      ELSEIF(NF.EQ.6)THEN
       YMSB(6) = XMSB
       YMSB(5) = YMSB(6)*CQ(ALPHAS(AM(5),2)/PI,5)/
     .                   CQ(ALPHAS(AM(6),2)/PI,5)
       YMSB(4) = YMSB(5)*CQ(ALPHAS(AM(4),2)/PI,4)/
     .                   CQ(ALPHAS(AM(5),2)/PI,4)
      ENDIF
      IF(Q.LT.AMC)THEN
       N0=3
       Q0 = 1.D0
      ELSEIF(Q.LE.AMB)THEN
       N0=4
       Q0 = AMC
      ELSEIF(Q.LE.AMT)THEN
       N0=5
       Q0 = AMB
      ELSE
       N0=6
       Q0 = AMT
      ENDIF
      IF(NNLO.EQ.1.AND.NF.GT.3)THEN
       XKFAC = TRAN(AM(NF),0D0)/TRAN(AM(NF),XK)
      ELSE
       XKFAC = 1D0
      ENDIF
      RUNM = YMSB(N0)*CQ(ALPHAS(Q,2)/PI,N0)/
     .               CQ(ALPHAS(Q0,2)/PI,N0)
     .       * XKFAC
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     
C
C     THIS PROGRAM COMPUTES THE RENORMALIZATION GROUP IMPROVED 
C     VALUES OF HIGGS MASSES AND COUPLINGS IN THE MSSM. 
C
C     INPUT: MA,TANB = TAN(BETA),MQ,MUR,MTOP,AU,AD,MU.
C 
C     ALL MASSES IN GEV UNITS. MA IS THE CP-ODD HIGGS MASS,
C     MTOP IS THE PHYSICAL TOP MASS, MQ AND MUR ARE THE SOFT
C     SUPERSYMMETRY BREAKING MASS PARAMETERS OF LEFT HANDED
C     AND RIGHT HANDED STOPS RESPECTIVELY, AU AND AD ARE THE
C     STOP AND SBOTTOM TRILINEAR SOFT BREAKING TERMS,  
C     RESPECTIVELY,  AND MU IS THE SUPERSYMMETRIC
C     HIGGS MASS PARAMETER. WE USE THE  CONVENTIONS FROM 
C     THE PHYSICS REPORT OF HABER AND KANE: LEFT RIGHT
C     STOP MIXING TERM PROPORTIONAL TO (AU - MU/TANB).
C
C     WE USE AS INPUT TANB DEFINED AT THE SCALE MTOP.
C
C     OUTPUT: MH,HM,MHCH, SA = SIN(ALPHA), CA= COS(ALPHA), TANBA
C     
C
C     WHERE MH AND HM ARE THE LIGHTEST AND HEAVIEST CP-EVEN 
C     HIGGS MASSES, MHCH IS THE CHARGED HIGGS MASS AND 
C     ALPHA IS THE HIGGS MIXING ANGLE. 
C
C     TANBA IS THE ANGLE TANB AT THE CP-ODD HIGGS MASS SCALE.
C
C     RANGE OF VALIDITY: 
C
C    (STOP1**2 - STOP2**2)/(STOP2**2 + STOP1**2) < 0.5
C    (SBOT1**2 - SBOT2**2)/(SBOT2**2 + SBOT2**2) < 0.5
C
C     WHERE STOP1, STOP2, SBOT1 AND SBOT2 ARE THE STOP AND
C     ARE THE SBOTTOM  MASS EIGENVALUES, RESPECTIVELY. THIS
C     RANGE AUTOMATICALLY EXCLUDES THE EXISTENCE OF TACHYONS. 
C
C     
C     FOR THE CHARGED HIGGS MASS COMPUTATION, THE METHOD IS 
C     VALID IF
C 
C     2 * |MB * AD* TANB|  < M_SUSY**2,  2 * |MTOP * AU| < M_SUSY**2
C
C     2 * |MB * MU * TANB| < M_SUSY**2,  2 * |MTOP * MU| < M_SUSY**2     
C
C     WHERE M_SUSY**2 IS THE AVERAGE OF THE SQUARED STOP MASS
C     EIGENVALUES, M_SUSY**2 = (STOP1**2 + STOP2**2)/2. THE SBOTTOM
C     MASSES HAVE BEEN ASSUMED TO BE OF ORDER OF THE STOP ONES.
C
C     M_SUSY**2 = (MQ**2 + MUR**2)*0.5 + MTOP**2
C
C     PROGRAM BASED ON THE WORK BY M. CARENA, J.R. ESPINOSA,
C     M. QUIROS AND C.E.M. WAGNER, CERN-PREPRINT CERN-TH/95-45.
C     
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE SUBH (MA,TANB,MQ,MUR,MTOP,AU,AD,MU,MH,HM,
     * MHCH,SA,CA,TANBA)
      IMPLICIT REAL*8(A-H,L,M,O-Z)
      COMMON/HSELF/LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,LAMBDA5,
     .             LAMBDA6,LAMBDA7
      COMMON/TREE/ITREE,ITRIA,IWRITE
      MZ = 91.18
      ALPHA1 = 0.0101
      ALPHA2 = 0.0337
      ALPHA3Z = 0.12    
      V = 174.1
      PI = 3.14159
      TANBA = TANB
      TANBT = TANB

C     MBOTTOM(MTOP) = 3. GEV
      MB = 3.
      ALPHA3 = ALPHA3Z/(1. +(11. - 10./3.)/4./PI*ALPHA3Z*
     *LOG(MTOP**2/MZ**2))

C     RMTOP= RUNNING TOP QUARK MASS     
      RMTOP = MTOP/(1.+4.*ALPHA3/3./PI)
      MS = ((MQ**2 + MUR**2)/2. + MTOP**2)**.5
      T = LOG(MS**2/MTOP**2)
      SINB = TANB/((1. + TANB**2)**.5)
      COSB = SINB/TANB
C      IF(MA.LE.MTOP) TANBA = TANBT
      IF(MA.GT.MTOP)
     *TANBA = TANBT*(1.-3./32./PI**2*
     *(RMTOP**2/V**2/SINB**2-MB**2/V**2/COSB**2)*
     *LOG(MA**2/MTOP**2))
     
      SINBT = TANBT/((1. + TANBT**2)**.5)
      COSBT = 1./((1. + TANBT**2)**.5)
      COS2BT = (TANBT**2 - 1.)/(TANBT**2 + 1.)
      G1 = (ALPHA1*4.*PI)**.5 
      G2 = (ALPHA2*4.*PI)**.5 
      G3 = (ALPHA3*4.*PI)**.5 
      HU = RMTOP/V/SINBT
      HD =  MB/V/COSBT

C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF(ITREE.EQ.1)THEN
       HU=0
       HD=0
       G3=0
      ENDIF
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      XAU = (2.*AU**2/MS**2)*(1. - AU**2/12./MS**2)
      XAD = (2.*AD**2/MS**2)*(1. - AD**2/12./MS**2)
      AUD = (-6.*MU**2/MS**2 - ( MU**2- AD*AU)**2/MS**4.
     *+ 3.*(AU + AD)**2/MS**2)/6.
      LAMBDA1 = ((G1**2 + G2**2)/4.)*(1.-3.*HD**2*T/8./PI**2)
     *+(3.*HD**4./8./PI**2) * (T + XAD/2. + (3.*HD**2/2. + HU**2/2.       
     *- 8.*G3**2) * (XAD*T + T**2)/16./PI**2) 
     *-(3.*HU**4.* MU**4./96./PI**2/MS**4.) * (1+ (9.*HU**2 -5.* HD**2
     *-  16.*G3**2) *T/16./PI**2)
      LAMBDA2 = ((G1**2 + G2**2)/4.)*(1.-3.*HU**2*T/8./PI**2)
     *+(3.*HU**4./8./PI**2) * (T + XAU/2. + (3.*HU**2/2. + HD**2/2.       
     *- 8.*G3**2) * (XAU*T + T**2)/16./PI**2) 
     *-(3.*HD**4.* MU**4./96./PI**2/MS**4.) * (1+ (9.*HD**2 -5.* HU**2
     *-  16.*G3**2) *T/16./PI**2)
      LAMBDA3 = ((G2**2 - G1**2)/4.)*(1.-3.*
     *(HU**2 + HD**2)*T/16./PI**2)
     *+(6.*HU**2*HD**2/16./PI**2) * (T + AUD/2. + (HU**2 + HD**2       
     *- 8.*G3**2) * (AUD*T + T**2)/16./PI**2) 
     *+(3.*HU**4./96./PI**2) * (3.*MU**2/MS**2 - MU**2*AU**2/
     *MS**4.)* (1.+ (6.*HU**2 -2.* HD**2/2.
     *-  16.*G3**2) *T/16./PI**2)
     *+(3.*HD**4./96./PI**2) * (3.*MU**2/MS**2 - MU**2*AD**2/
     *MS**4.)*(1.+ (6.*HD**2 -2.* HU**2
     *-  16.*G3**2) *T/16./PI**2)
      LAMBDA4 = (- G2**2/2.)*(1.-3.*(HU**2 + HD**2)*T/16./PI**2)
     *-(6.*HU**2*HD**2/16./PI**2) * (T + AUD/2. + (HU**2 + HD**2       
     *- 8.*G3**2) * (AUD*T + T**2)/16./PI**2) 
     *+(3.*HU**4./96./PI**2) * (3.*MU**2/MS**2 - MU**2*AU**2/
     *MS**4.)*
     *(1+ (6.*HU**2 -2.* HD**2
     *-  16.*G3**2) *T/16./PI**2)
     *+(3.*HD**4./96./PI**2) * (3.*MU**2/MS**2 - MU**2*AD**2/
     *MS**4.)*
     *(1+ (6.*HD**2 -2.* HU**2/2.
     *-  16.*G3**2) *T/16./PI**2)
      LAMBDA5 = -(3.*HU**4.* MU**2*AU**2/96./PI**2/MS**4.) *
     * (1- (2.*HD**2 -6.* HU**2 + 16.*G3**2) *T/16./PI**2)
     *-(3.*HD**4.* MU**2*AD**2/96./PI**2/MS**4.) *
     * (1- (2.*HU**2 -6.* HD**2 + 16.*G3**2) *T/16./PI**2)
      LAMBDA6 = (3.*HU**4.* MU**3.*AU/96./PI**2/MS**4.) *
     * (1- (7.*HD**2/2. -15.* HU**2/2. + 16.*G3**2) *T/16./PI**2)
     *+(3.*HD**4.* MU *(AD**3./MS**3. - 6.*AD/MS )/96./PI**2/MS) *
     * (1- (HU**2/2. -9.* HD**2/2. + 16.*G3**2) *T/16./PI**2)
      LAMBDA7 = (3.*HD**4.* MU**3.*AD/96./PI**2/MS**4.) *
     * (1- (7.*HU**2/2. -15.* HD**2/2. + 16.*G3**2) *T/16./PI**2)
     *+(3.*HU**4.* MU *(AU**3./MS**3. - 6.*AU/MS )/96./PI**2/MS) *
     * (1- (HD**2/2. -9.* HU**2/2. + 16.*G3**2) *T/16./PI**2)
      TRM2 = MA**2 + 2.*V**2* (LAMBDA1* COSBT**2 + 
     *2.* LAMBDA6*SINBT*COSBT
     *+ LAMBDA5*SINBT**2 + LAMBDA2* SINBT**2 + 2.* LAMBDA7*SINBT*COSBT
     *+ LAMBDA5*COSBT**2) 
      DETM2 = 4.*V**4.*(-(SINBT*COSBT*(LAMBDA3 + LAMBDA4) + 
     *LAMBDA6*COSBT**2
     *+ LAMBDA7* SINBT**2)**2 + (LAMBDA1* COSBT**2 +
     *2.* LAMBDA6* COSBT*SINBT 
     *+ LAMBDA5*SINBT**2)*(LAMBDA2* SINBT**2 +2.* LAMBDA7* COSBT*SINBT 
     *+ LAMBDA5*COSBT**2)) + MA**2*2.*V**2 * 
     *((LAMBDA1* COSBT**2 +2.* 
     *LAMBDA6* COSBT*SINBT + LAMBDA5*SINBT**2)*COSBT**2 + 
     *(LAMBDA2* SINBT**2 +2.* LAMBDA7* COSBT*SINBT + LAMBDA5*COSBT**2)
     **SINBT**2
     * +2.*SINBT*COSBT* (SINBT*COSBT*(LAMBDA3 
     * + LAMBDA4) + LAMBDA6*COSBT**2
     *+ LAMBDA7* SINBT**2))  

      MH2 = (TRM2 - (TRM2**2 - 4.* DETM2)**.5)/2.
      HM2 = (TRM2 + (TRM2**2 - 4.* DETM2)**.5)/2.
      HM = HM2**.5 
      MH = MH2**.5 
      MHCH2 = MA**2 + (LAMBDA5 - LAMBDA4)* V**2
      MHCH = MHCH2**.5
      MHCH = MHCH2**.5

      SINALPHA = SQRT(((TRM2**2 - 4.* DETM2)**.5) -
     * ((2.*V**2*(LAMBDA1* COSBT**2 + 2.* 
     *LAMBDA6* COSBT*SINBT 
     *+ LAMBDA5*SINBT**2) + MA**2*SINBT**2)
     *- (2.*V**2*(LAMBDA2* SINBT**2 +2.* LAMBDA7* COSBT*SINBT 
     *+ LAMBDA5*COSBT**2) + MA**2*COSBT**2)))/
     *SQRT(((TRM2**2 - 4.* DETM2)**.5))/2.**.5

      COSALPHA = (2.*(2.*V**2*(SINBT*COSBT*(LAMBDA3 + LAMBDA4) + 
     *LAMBDA6*COSBT**2 + LAMBDA7* SINBT**2) -
     *MA**2*SINBT*COSBT))/2.**.5/
     *SQRT(((TRM2**2 - 4.* DETM2)**.5)* 
     *(((TRM2**2 - 4.* DETM2)**.5) -
     * ((2.*V**2*(LAMBDA1* COSBT**2 + 2.* 
     *LAMBDA6* COSBT*SINBT 
     *+ LAMBDA5*SINBT**2) + MA**2*SINBT**2)
     *- (2.*V**2*(LAMBDA2* SINBT**2 +2.* LAMBDA7* COSBT*SINBT  
     *+ LAMBDA5*COSBT**2) + MA**2*COSBT**2))))

      SA = -SINALPHA
      CA = -COSALPHA
 
 2242 RETURN
      END

C *****************************************************************
C ************* SUBROUTINE FOR THE SUSY COUPLINGS *****************
C *****************************************************************
      SUBROUTINE SUSYCP(TGBET)
C--CALCULATION OF SUSY-COUPLINGS FROM PSEUDOSCALAR HIGGS MASS
C--AND TAN(BETA)
C--AML: LIGHT SCALAR HIGGS MASS
C--AMH: HEAVY SCALAR HIGGS MASS
C--AMA: PSEUDOSCALAR HIGGS MASS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION LA1,LA2,LA3,LA4,LA5,LA6,LA7,LA3T
      COMMON/SUSYPAR/GF,ALPHA,AMZ,AMSQ
      COMMON/QMASS/AMB,AMT
      COMMON/HMASS/AMA,AML,AMH,AMCH
      COMMON/HSELF/LA1,LA2,LA3,LA4,LA5,LA6,LA7
      COMMON/BREAK/AMUR,AU,AD,AMU
      COMMON/COUP/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,
     .            GHHH,GLLL,GHLL,GLHH,GHAA,GLAA,GLVV,GHVV,
     .            GLPM,GHPM,GZT,GZB
      COMMON/PROC/IPROC,IAPP,ITOP,ISFAC
      COMMON/TREE/ITREE,ITRIA,IWRITE
      V=1.D0/DSQRT(DSQRT(2.D0)*GF)
      BET=DATAN(TGBET)
      SB = DSIN(BET)
      CB = DCOS(BET)
      CALL SUBH(AMA,TGBET,AMSQ,AMUR,AMT,AU,AD,AMU,AML,AMH,AMCH,SA,CA,
     .          TBA)
      LA3T=LA3+LA4+LA5
      AMA2=AMA**2
      AML2=AML**2
      AMH2=AMH**2
      AMP2=AMCH**2
      SBMA = SB*CA-CB*SA
      CBMA = CB*CA+SB*SA
      SBPA = SB*CA+CB*SA
      CBPA = CB*CA-SB*SA
      S2A = 2*SA*CA
      C2A = CA**2-SA**2
      S2B = 2*SB*CB
      C2B = CB**2-SB**2
      GLZZ = 1/V/2*AML2*SBMA
      GHZZ = 1/V/2*AMH2*CBMA
      GLWW = 2*GLZZ
      GHWW = 2*GHZZ
      GLAZ = 1/V*(AML2-AMA2)*CBMA
      GHAZ = -1/V*(AMH2-AMA2)*SBMA
      GLPW = -1/V*(AMP2-AML2)*CBMA
      GLMW = GLPW
      GHPW = 1/V*(AMP2-AMH2)*SBMA
      GHMW = GHPW
      GAPW = 1/V*(AMP2-AMA2)
      GAMW = -GAPW
      GHHH = V/2*(LA1*CA**3*CB + LA2*SA**3*SB + LA3T*SA*CA*SBPA
     .     + LA6*CA**2*(3*SA*CB+CA*SB) + LA7*SA**2*(3*CA*SB+SA*CB))
      GLLL = -V/2*(LA1*SA**3*CB - LA2*CA**3*SB + LA3T*SA*CA*CBPA
     .     - LA6*SA**2*(3*CA*CB-SA*SB) + LA7*CA**2*(3*SA*SB-CA*CB))
      GLHH = -3*V/2*(LA1*CA**2*CB*SA - LA2*SA**2*SB*CA
     .     + LA3T*(SA**3*CB-CA**3*SB+2*SBMA/3)
     .     - LA6*CA*(CB*C2A-SA*SBPA) - LA7*SA*(C2A*SB+CA*SBPA))
      GHLL = 3*V/2*(LA1*SA**2*CB*CA + LA2*CA**2*SB*SA
     .     + LA3T*(SA**3*SB+CA**3*CB-2*CBMA/3)
     .     - LA6*SA*(CB*C2A+CA*CBPA) + LA7*CA*(C2A*SB+SA*CBPA))
      GLAA = -V/2*(LA1*SB**2*CB*SA - LA2*CB**2*SB*CA
     .     - LA3T*(SB**3*CA-CB**3*SA) + 2*LA5*SBMA
     .     - LA6*SB*(CB*SBPA+SA*C2B) - LA7*CB*(C2B*CA-SB*SBPA))
      GHAA = V/2*(LA1*SB**2*CB*CA + LA2*CB**2*SB*SA
     .     + LA3T*(SB**3*SA+CB**3*CA) - 2*LA5*CBMA
     .     - LA6*SB*(CB*CBPA+CA*C2B) + LA7*CB*(SB*CBPA+SA*C2B))
      GLPM = 2*GLAA + V*(LA5 - LA4)*SBMA
      GHPM = 2*GHAA + V*(LA5 - LA4)*CBMA

      GLZZ = 2*GLZZ
      GHZZ = 2*GHZZ
      GLLL = 6*GLLL
      GHHH = 6*GHHH
      GLHH = 2*GLHH
      GHLL = 2*GHLL
      GLAA = 2*GLAA
      GHAA = 2*GHAA
      XNORM = AMZ**2/V
      GLLL = GLLL/XNORM
      GHLL = GHLL/XNORM
      GLHH = GLHH/XNORM
      GHHH = GHHH/XNORM
      GHAA = GHAA/XNORM
      GLAA = GLAA/XNORM
      GLPM = GLPM/XNORM
      GHPM = GHPM/XNORM

      GAT=1.D0/TGBET
      GAB=TGBET
      GLT=CA/SB
      GLB=-SA/CB
      GHT=SA/SB
      GHB=CA/CB
      GZAH=SBMA
      GZAL=-CBMA
      GHVV=CBMA
      GLVV=SBMA
      IF(IWRITE.EQ.1)THEN
       WRITE(6,*)
       WRITE(6,*)'TGB  = ',TGBET,'     MA   = ',AMA
       WRITE(6,*)'ML   = ',AML ,'     MH   = ',AMH
       WRITE(6,*)'SA   = ',SA  ,'     CA   = ',CA,SA**2+CA**2
       WRITE(6,*)'SB   = ',SB  ,'     CB   = ',CB,SB**2+CB**2
       WRITE(6,*)'SBPA = ',SBPA,'     CBPA = ',CBPA,SBPA**2+CBPA**2
       WRITE(6,*)'SBMA = ',SBMA,'     CBMA = ',CBMA,SBMA**2+CBMA**2
       WRITE(6,*)'S2A  = ',S2A ,'     C2A  = ',C2A ,S2A**2+C2A**2
       WRITE(6,*)'S2B  = ',S2B ,'     C2B  = ',C2B ,S2B**2+C2B**2
       WRITE(6,*)'GLT  = ',GLT ,'     GLB  = ',GLB
       WRITE(6,*)'GHT  = ',GHT ,'     GHB  = ',GHB
       WRITE(6,*)'GAT  = ',GAT ,'     GAB  = ',GAB
       WRITE(6,*)'GLLL = ',GLLL,'     GHLL = ',GHLL
       WRITE(6,*)'GLHH = ',GLHH,'     GHHH = ',GHHH
       WRITE(6,*)'GLAA = ',GLAA,'     GHAA = ',GHAA
       WRITE(6,*)'GZAL = ',GZAL,'     GZAH = ',GZAH
      ENDIF
      GLT=GLT*AMT/V
      GLB=GLB*AMB/V
      GHT=GHT*AMT/V
      GHB=GHB*AMB/V
      GAT=GAT*AMT/V
      GAB=GAB*AMB/V
      GLLL=GLLL*AMZ**2/V
      GHHH=GHHH*AMZ**2/V
      GLHH=GLHH*AMZ**2/V
      GHLL=GHLL*AMZ**2/V
      GLAA=GLAA*AMZ**2/V
      GHAA=GHAA*AMZ**2/V
      GLPM=GLPM*AMZ**2/V
      GHPM=GHPM*AMZ**2/V
      GZAL=GZAL*AMZ/V
      GZAH=GZAH*AMZ/V
      ZFAC=DSQRT(GF*AMZ**2/2.D0/DSQRT(2.D0))
      GZT=ZFAC
      GZB=-ZFAC
      IF(IWRITE.EQ.1)THEN
       WRITE(6,*)'GZT  = ',GZT ,'     GZB  = ',GZB
      ENDIF
      IF(ITOP.EQ.1)THEN
       GAB=0.D0
       GLB=0.D0
       GHB=0.D0
       GZB=0.D0
      ENDIF
      IF(ITOP.EQ.2)THEN
       GAT=0.D0
       GLT=0.D0
       GHT=0.D0
       GZT=0.D0
      ENDIF
      IF(ITRIA.EQ.1)THEN
       GLLL=0
       GLHH=0
       GHLL=0
       GHHH=0
       GLAA=0
       GHAA=0
       GZAH=0
       GZAL=0
      ELSEIF(ITRIA.EQ.2)THEN
       GZAH=0
       GZAL=0
      ENDIF
      RETURN
      END

      SUBROUTINE SUSYCP0(TGBET)
C--CALCULATION OF SUSY-COUPLINGS FROM PSEUDOSCALAR HIGGS MASS
C--AND TAN(BETA)
C--AML: LIGHT SCALAR HIGGS MASS
C--AMH: HEAVY SCALAR HIGGS MASS
C--AMA: PSEUDOSCALAR HIGGS MASS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION M1,M2,MZ
      COMMON/PROC/IPROC,IAPP,ITOP,ISFAC
      COMMON/PARAM/S,T,U,M1,M2,MZ
      COMMON/TOTPAR/W,GEVPB,GF,GAMZ,GAMA,GAML,GAMH
      COMMON/QMASS/AMB,AMT
      COMMON/HMASS/AMA,AML,AMH,AMCH
      COMMON/WDPARAM/WGF,WALPH,WMTAU,WMMUON,WMZ,AMW
      COMMON/COUP/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,
     .            GHHH,GLLL,GHLL,GLHH,GHAA,GLAA,GLVV,GHVV,
     .            GLPM,GHPM,GZT,GZB
      COMMON/TREE/ITREE,ITRIA,IWRITE
      V=1.D0/DSQRT(DSQRT(2.D0)*GF)
      ZFAC=DSQRT(GF*MZ**2/2.D0/DSQRT(2.D0))
      AMZ=MZ
      PI=4.D0*DATAN(1.D0)
      BET=DATAN(TGBET)
C     CW2=1.D0-SW2
C     SUSYEPS=3.D0/2.D0*ALPHA/PI/SW2/CW2/DSIN(BET)**2*AMQ**4/AMZ**2
C    .             *DLOG(1.D0+AMS**2/AMQ**2)
      SUSYEPS=0.D0
      DIS=(AMA**2+AMZ**2+SUSYEPS)**2    
     .    -4.D0*AMA**2*AMZ**2*DCOS(2.D0*BET)**2
     .    -4.D0*SUSYEPS*(AMA**2*DSIN(BET)**2+AMZ**2*DCOS(BET)**2)
      IF(DIS.GE.0.D0)THEN
        DELTA=DSQRT(DIS)
        AML=DSQRT((AMA**2+AMZ**2+SUSYEPS-DELTA)/2.D0)
        AMH=DSQRT((AMA**2+AMZ**2+SUSYEPS+DELTA)/2.D0)
        AMCH=DSQRT(AMA**2+AMW**2)
        ALP=DATAN(DTAN(2.D0*BET)*(AMA**2+AMZ**2)/(AMA**2-AMZ**2+
     .      SUSYEPS/DCOS(2.D0*BET)))/2.D0
        IF(ALP.GT.0.D0)ALP=ALP-PI/2.D0
        GAT=1.D0/TGBET
        GAB=TGBET
        GLT=DCOS(ALP)/DSIN(BET)
        GLB=-DSIN(ALP)/DCOS(BET)
        GHT=DSIN(ALP)/DSIN(BET)
        GHB=DCOS(ALP)/DCOS(BET)
        GZAH=DSIN(BET-ALP)
        GZAL=-DCOS(BET-ALP)
        GHVV=DCOS(BET-ALP)
        GLVV=DSIN(BET-ALP)
        GHHH=3.D0*DCOS(2.D0*ALP)*DCOS(ALP+BET)
        GLLL=3.D0*DCOS(2.D0*ALP)*DSIN(ALP+BET)
        GHLL=2.D0*DSIN(2.D0*ALP)*DSIN(ALP+BET)
     .       -DCOS(2.D0*ALP)*DCOS(ALP+BET)
        GLHH=-2.D0*DSIN(2.D0*ALP)*DCOS(ALP+BET)
     .       -DCOS(2.D0*ALP)*DSIN(ALP+BET)
        GHAA=-DCOS(2.D0*BET)*DCOS(ALP+BET)
        GLAA=DCOS(2.D0*BET)*DSIN(ALP+BET)
        GLPM=2*AMW**2/AMZ**2*DSIN(BET-ALP)+DCOS(2.D0*BET)*DSIN(ALP+BET)
        GHPM=2*AMW**2/AMZ**2*DCOS(BET-ALP)-DCOS(2.D0*BET)*DCOS(ALP+BET)
      ELSE
        GAT=0.D0
        GAB=0.D0
        GLT=0.D0
        GLB=0.D0
        GHT=0.D0
        GHB=0.D0
        GZAH=0.D0
        GZAL=0.D0
        GHHH=0.D0
        GLLL=0.D0
        GHLL=0.D0
        GLHH=0.D0
        GHAA=0.D0
        GLAA=0.D0
      ENDIF
      IF(IWRITE.EQ.1)THEN
       SA=DSIN(ALP)
       CA=DCOS(ALP)
       SB=DSIN(BET)
       CB=DCOS(BET)
       WRITE(6,*)'TGB  = ',TGBET,'     MA   = ',AMA
       WRITE(6,*)'ML   = ',AML ,'     MH   = ',AMH
       WRITE(6,*)'SA   = ',SA  ,'     CA   = ',CA,SA**2+CA**2
       WRITE(6,*)'SB   = ',SB  ,'     CB   = ',CB,SB**2+CB**2
       WRITE(6,*)'GLT  = ',GLT ,'     GLB  = ',GLB
       WRITE(6,*)'GHT  = ',GHT ,'     GHB  = ',GHB
       WRITE(6,*)'GAT  = ',GAT ,'     GAB  = ',GAB
       WRITE(6,*)'GLLL = ',GLLL,'     GHLL = ',GHLL
       WRITE(6,*)'GLHH = ',GLHH,'     GHHH = ',GHHH
       WRITE(6,*)'GLAA = ',GLAA,'     GHAA = ',GHAA
       WRITE(6,*)'GZAL = ',GZAL,'     GZAH = ',GZAH
      ENDIF
      GLT=GLT*AMT/V
      GLB=GLB*AMB/V
      GHT=GHT*AMT/V
      GHB=GHB*AMB/V
      GAT=GAT*AMT/V
      GAB=GAB*AMB/V
      GZT=ZFAC
      GZB=-ZFAC
      GLLL=GLLL*MZ**2/V
      GHHH=GHHH*MZ**2/V
      GLHH=GLHH*MZ**2/V
      GHLL=GHLL*MZ**2/V
      GLAA=GLAA*MZ**2/V
      GHAA=GHAA*MZ**2/V
      GLPM=GLPM*MZ**2/V
      GHPM=GHPM*MZ**2/V
      GZAL=GZAL*MZ/V
      GZAH=GZAH*MZ/V
      IF(ITOP.EQ.1)THEN
       GAB=0.D0
       GLB=0.D0
       GHB=0.D0
       GZB=0.D0
      ENDIF
      IF(ITOP.EQ.2)THEN
       GAT=0.D0
       GLT=0.D0
       GHT=0.D0
       GZT=0.D0
      ENDIF
      IF(ITRIA.EQ.1)THEN
       GLLL=0
       GLHH=0
       GHLL=0
       GHHH=0
       GLAA=0
       GHAA=0
       GZAH=0
       GZAL=0
      ELSEIF(ITRIA.EQ.2)THEN
       GZAH=0
       GZAL=0
      ENDIF
      RETURN
      END
 
      SUBROUTINE SMCP
C--CALCULATION OF SM-COUPLINGS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION M1,M2,MZ
      COMMON/PROC/IPROC,IAPP,ITOP,ISFAC
      COMMON/PARAM/S,T,U,M1,M2,MZ
      COMMON/TOTPAR/W,GEVPB,GF,GAMZ,GAMA,GAML,GAMH
      COMMON/QMASS/AMB,AMT
      COMMON/HMASS/AMA,AML,AMH,AMCH
      COMMON/TREE/ITREE,ITRIA,IWRITE
      COMMON/TRILINEAR/FACTRIA,FACT,FACB,FACTT,FACG,FACGG,ISILH
      COMMON/COUP/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,
     .            GHHH,GLLL,GHLL,GLHH,GHAA,GLAA,GLVV,GHVV,
     .            GLPM,GHPM,GZT,GZB
      AML=AMA
      AMH=AMA
      AMCH=AMA
      AMZ=MZ
      V=1.D0/DSQRT(DSQRT(2.D0)*GF)
      ZFAC=DSQRT(GF*MZ**2/2.D0/DSQRT(2.D0))
      HFAC=3.D0*AMH**2/V
      GAT=0.D0
      GAB=0.D0
      GLT=0.D0
      GLB=0.D0
      GZAH=0.D0
      GZAL=0.D0
      GLLL=0.D0
      GHLL=0.D0
      GLHH=0.D0
      GHAA=0.D0
      GLAA=0.D0
      GLPM=0.D0
      GHPM=0.D0

      GHT=AMT/V
      GHB=AMB/V
      GZT=ZFAC
      GZB=-ZFAC
      GHHH=HFAC * FACTRIA
      IF(ITOP.EQ.1)THEN
       GHB=0.D0
       GZB=0.D0
      ENDIF
      IF(ITOP.EQ.2)THEN
       GHT=0.D0
       GZT=0.D0
      ENDIF
      IF(ITRIA.EQ.1)THEN
       GHHH=0
      ENDIF
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION ALPHAS(Q,N)
C--ALPHA_S: Q = SCALE, N = 1 -> LO, N = 2 -> NLO
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XLB(6)
      COMMON/ALSLAM/XLB1(6),XLB2(6)
      COMMON/ALS/XLAMBDA,AMC,AMB,AMT,N0
      B0(NF)=33.D0-2.D0*NF
      B1(NF)=6.D0*(153.D0-19.D0*NF)/B0(NF)**2
      ALS1(NF,X)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB(NF)**2))
      ALS2(NF,X)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB(NF)**2))
     .          *(1.D0-B1(NF)*DLOG(DLOG(X**2/XLB(NF)**2))
     .           /DLOG(X**2/XLB(NF)**2))
      PI=4.D0*DATAN(1.D0)
      IF(N.EQ.1)THEN
       DO 1 I=1,6
        XLB(I)=XLB1(I)
1      CONTINUE
      ELSE
       DO 2 I=1,6
        XLB(I)=XLB2(I)
2      CONTINUE
      ENDIF

      IF(Q.LT.AMC)THEN
       NF=3
      ELSEIF(Q.LE.AMB)THEN
       NF=4
      ELSEIF(Q.LE.AMT)THEN
       NF=5
      ELSE
       NF=6
      ENDIF
      IF(N.EQ.1)THEN
        ALPHAS=ALS1(NF,Q)
      ELSE
        ALPHAS=ALS2(NF,Q)
      ENDIF
      RETURN
      END
 
      SUBROUTINE ALSINI(ACC)
C--ALPHA_S: INITIALIZATION OF LAMBDA_NF, ACC = ACCURAY AT THRESHOLDS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XLB(6)
      COMMON/ALSLAM/XLB1(6),XLB2(6)
      COMMON/ALS/XLAMBDA,AMC,AMB,AMT,N0
      B0(NF)=33.D0-2.D0*NF
      B1(NF)=6.D0*(153.D0-19.D0*NF)/B0(NF)**2
      ALS2(NF,X)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB(NF)**2))
     .          *(1.D0-B1(NF)*DLOG(DLOG(X**2/XLB(NF)**2))
     .           /DLOG(X**2/XLB(NF)**2))
      PI=4.D0*DATAN(1.D0)
      XLB1(1)=0D0
      XLB1(2)=0D0
      XLB2(1)=0D0
      XLB2(2)=0D0
      IF(N0.EQ.4)THEN
       XLB(4)=XLAMBDA
       XLB(5)=XLB(4)*(XLB(4)/AMB)**(2.D0/23.D0)
      ELSEIF(N0.EQ.5)THEN
       XLB(5)=XLAMBDA
       XLB(4)=XLB(5)*(XLB(5)/AMB)**(-2.D0/25.D0)
      ENDIF
      XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
      XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
      DO 1 I=1,6
       XLB1(I)=XLB(I)
1     CONTINUE
      IF(N0.EQ.4)THEN
       XLB(4)=XLAMBDA
       XLB(5)=XLB(4)*(XLB(4)/AMB)**(2.D0/23.D0)
     .             *(2.D0*DLOG(AMB/XLB(4)))**(-963.D0/13225.D0)
       XLB(5)=XITER(AMB,XLB(4),4,XLB(5),5,ACC)
      ELSEIF(N0.EQ.5)THEN
       XLB(5)=XLAMBDA
       XLB(4)=XLB(5)*(XLB(5)/AMB)**(-2.D0/25.D0)
     .             *(2.D0*DLOG(AMB/XLB(5)))**(963.D0/14375.D0)
       XLB(4)=XITER(AMB,XLB(5),5,XLB(4),4,ACC)
      ENDIF
      XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
     .            *(2.D0*DLOG(AMC/XLB(4)))**(107.D0/2025.D0)
      XLB(3)=XITER(AMC,XLB(4),4,XLB(3),3,ACC)
      XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
     .           *(2.D0*DLOG(AMT/XLB(5)))**(-321.D0/3381.D0)
      XLB(6)=XITER(AMT,XLB(5),5,XLB(6),6,ACC)
      DO 2 I=1,6
       XLB2(I)=XLB(I)
2     CONTINUE
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION XITER(Q,XLB1,NF1,XLB,NF2,ACC)
C--ALPHA_S: ITERATION FOR ALPHA_S(M_Z) -> LAMBDA_5
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      B0(NF)=33.D0-2.D0*NF
      B1(NF)=6.D0*(153.D0-19.D0*NF)/B0(NF)**2
      ALS2(NF,X,XLB)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB**2))
     .              *(1.D0-B1(NF)*DLOG(DLOG(X**2/XLB**2))
     .              /DLOG(X**2/XLB**2))
      AA(NF)=12D0*PI/B0(NF)
      BB(NF)=B1(NF)/AA(NF)
      XIT(A,B,X)=A/2.D0*(1D0+DSQRT(1D0-4D0*B*DLOG(X)))
      PI=4.D0*DATAN(1.D0)
      XLB2=XLB
      IF(ACC.GE.1.D0) GOTO 111
      II=0
1     II=II+1
      X=DLOG(Q**2/XLB2**2)
      ALP=ALS2(NF1,Q,XLB1)
      A=AA(NF2)/ALP
      B=BB(NF2)*ALP
      XX=XIT(A,B,X)
      XLB2=Q*DEXP(-XX/2.D0)
      Y1=ALS2(NF1,Q,XLB1)
      Y2=ALS2(NF2,Q,XLB2)
      DY=DABS(Y2-Y1)/Y1
      IF(DY.GE.ACC) GOTO 1
111   XITER=XLB2
      RETURN
      END

      DOUBLE PRECISION FUNCTION XITLA(NO,ALP,ACC)
C--ITERATION ROUTINE TO DETERMINE IMPROVED LAMBDAS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/SUSYPAR/GF,ALPHA,AMZ,AMSQ
      B0(NF)=33.D0-2.D0*NF
      B1(NF)=6.D0*(153.D0-19.D0*NF)/B0(NF)**2
      B2(NF)=27/2.D0*(2857-5033/9.D0*NF+325/27.D0*NF**2)/B0(NF)**3
      B3(NF)= 81*(149753/6.d0+3564*zeta3-(1078361/162.d0+6508*zeta3/27)
     .        *nf+(50065/162.d0+6472*zeta3/81)*nf**2+1093/729.d0*nf**3)
     .      / B0(NF)**4
      ALS2(NF,X,XLB)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB**2))
     .              *(1.D0-B1(NF)*DLOG(DLOG(X**2/XLB**2))
     .              /DLOG(X**2/XLB**2))
      ALS3(NF,X,XLB)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB**2))
     .          *(1.D0-B1(NF)*DLOG(DLOG(X**2/XLB**2))
     .           /DLOG(X**2/XLB**2)
     .           +(B1(NF)**2*(DLOG(DLOG(X**2/XLB**2))**2
     .                      -DLOG(DLOG(X**2/XLB**2))-1)+B2(NF))
     .           /DLOG(X**2/XLB**2)**2)
      ALS4(NF,X,XLB)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB**2))
     .          *(1.D0-B1(NF)*DLOG(DLOG(X**2/XLB**2))
     .           /DLOG(X**2/XLB**2)
     .           +(B1(NF)**2*(DLOG(DLOG(X**2/XLB**2))**2
     .                      -DLOG(DLOG(X**2/XLB**2))-1)+B2(NF))
     .           /DLOG(X**2/XLB**2)**2
     .           -(B1(NF)**3*(DLOG(DLOG(X**2/XLB**2))**3
     .                      -5*DLOG(DLOG(X**2/XLB**2))**2/2
     .                      -2*DLOG(DLOG(X**2/XLB**2))+1/2.d0)
     .            +3*B1(NF)*B2(NF)*DLOG(DLOG(X**2/XLB**2))
     .            -B3(NF)/2)/DLOG(X**2/XLB**2)**3)
      AA(NF)=12D0*PI/B0(NF)
      BB(NF)=B1(NF)/AA(NF)
      CC(NF)=B2(NF)/AA(NF)
      DD(NF)=B3(NF)/AA(NF)
      XIT(A,B,X)=A/2.D0*(1D0+DSQRT(1D0-4D0*B*DLOG(X)))
      XIT3(A,B,C,X)=A/2.D0*(1D0+DSQRT(1D0-4D0*B*DLOG(X)
     .          *(1-(A*B*(DLOG(X)**2-DLOG(X)-1)+C/B)/X/DLOG(X))))
      XIT4(A,B,C,D,X)=A/2.D0*(1D0+DSQRT(1D0-4D0*B*DLOG(X)
     .          *(1-(A*B*(DLOG(X)**2-DLOG(X)-1)+C/B)/X/DLOG(X)
     .          +(A**2*B**2*(DLOG(X)**3-5*DLOG(X)**2/2-2*DLOG(X)+1/2.D0)
     .           +3*A*C*DLOG(X)-D/B/2)/X**2/DLOG(X))))
      PI=4.D0*DATAN(1.D0)
      N3LO = 0
      ZETA2 = PI**2/6
      ZETA3 = 1.2020569031595942853997381D0
      NF=5
      Q=AMZ
      XLB=Q*DEXP(-AA(NF)/ALP/2.D0)
      IF(NO.EQ.1)GOTO 111
      II=0
1     II=II+1
      X=DLOG(Q**2/XLB**2)
      A=AA(NF)/ALP
      B=BB(NF)*ALP
      C=CC(NF)*ALP
      D=DD(NF)*ALP
      IF(NO.EQ.2)THEN
       XX=XIT(A,B,X)
      ELSEIF(NO.EQ.3)THEN
       XX=XIT3(A,B,C,X)
      ELSE
       XX=XIT4(A,B,C,D,X)
      ENDIF
      IF(N3LO.NE.0) XX=XIT4(A,B,C,D,X)
      XLB=Q*DEXP(-XX/2.D0)
      Y1=ALP
      IF(NO.EQ.2)THEN
       Y2=ALS2(NF,Q,XLB)
      ELSEIF(NO.EQ.3)THEN
       Y2=ALS3(NF,Q,XLB)
      ELSE
       Y2=ALS4(NF,Q,XLB)
      ENDIF
      IF(N3LO.NE.0) Y2=ALS4(NF,Q,XLB)
      DY=DABS(Y2-Y1)/Y1
      IF(DY.GE.ACC) GOTO 1
111   XITLA=XLB
      RETURN
      END

      SUBROUTINE VEGASN(FXN,ACC,NDIM,NPOINT,NITT,NPRN,IGRAPH)
C--VEGAS: WARMUP AND MAIN INTEGRATION
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL FXN
      NPOINT0=NPOINT/10
      NITT0=5
      CALL VEGAS(FXN,ACC,NDIM,NPOINT0,NITT0,NPRN,IGRAPH)
      CALL VEGAS1(FXN,ACC,NDIM,NPOINT,NITT,NPRN,IGRAPH)
      RETURN
      END

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        FUNCTION CSPEN(Z)                                              
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SPENCE-FUNKTION KOMPLEX, FREI NACH HOLLIK                     C
C---------------------------------------------------------------------C
C       20.07.83    LAST CHANGED 10.05.89        ANSGAR DENNER        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        COMPLEX*16 CSPEN,W,SUM,Z,U                                     
        REAL*8 RZ,AZ,A1                                                
        REAL*8 B(9)/                                                   
     1   0.1666666666666666666666666667D0,                             
     2  -0.0333333333333333333333333333D0,                             
     3   0.0238095238095238095238095238D0,                             
     4  -0.0333333333333333333333333333D0,                             
     5   0.0757575757575757575757575758D0,                             
     6  -0.2531135531135531135531135531D0,                             
     7   1.1666666666666666666666666667D0,                             
     8  -7.09215686274509804D0         ,                               
     9  54.97117794486215539D0         /                               
C     BEACHTE:                 B(N)=B2N                                
C     B(1)=1./6.                                                       
C     B(2)=-1./30.                                                     
C     B(3)=1./42.                                                      
C     B(4)=-1./30.                                                     
C     B(5)=5./66.                                                      
C     B(6)=-691./2730.                                                 
C     B(7)=7./6.                                                       
C     B(8)=-3617./510.                                                 
C     B(9)=43867./798.                                                 
C     B(10)=-174611./330.                                              
C     B(11)=854513./138.                                               
C     PI=3.1415926535897932384                                         
C     PI*PI/6.=1.6449..., PI*PI/3=3.28986...                           
C                                                                      
      Z =Z*DCMPLX(1D0)                                                 
      RZ=DREAL(Z)                                                      
      AZ=CDABS(Z)                                                      
      A1=CDABS(1D0-Z)                                                  
C     IF((SNGL(RZ) .EQ. 0.0) .AND. (SNGL(DIMAG(Z)) .EQ. 0.0)) THEN     
C ---> CHANGED  10.5.89                                                
      IF(AZ .LT. 1D-20) THEN                                           
        CSPEN=-CDLOG(1D0-Z)                                            
        RETURN                                                         
      END IF                                                           
c      IF((SNGL(RZ) .EQ. 1.0) .AND. (SNGL(DIMAG(Z)) .EQ. 0.0)) THEN     
c ---> changed 5.7.94
       IF( (ABS(RZ-1D0).LT.1D-18) .AND. (ABS(DIMAG(Z)).LT.1D-18) ) THEN     
        CSPEN=1.64493406684822643D0                                    
        RETURN                                                         
      END IF                                                           
      IF(RZ.GT.5D-1) GOTO 20                                           
      IF(AZ.GT.1D0) GOTO 10                                            
      W=-CDLOG(1D0-Z)                                                  
      SUM=W-0.25D0*W*W                                                 
      U=W                                                              
      IF(CDABS(U).LT.1D-10) GOTO 2                                     
      DO 1 K=1,9                                                       
      U=U*W*W/DFLOAT(2*K*(2*K+1))                                      
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 2                            
      SUM=SUM+U*B(K)                                                   
 1    CONTINUE                                                         
 2    CSPEN=SUM                                                        
      RETURN                                                           
10    W=-CDLOG(1D0-1D0/Z)                                              
      SUM=W-0.25D0*W*W                                                 
      U=W                                                              
      IF(CDABS(U).LT.1D-10) GOTO 12                                    
                                                                       
      DO 11 K=1,9                                                      
      U=U*W*W/DFLOAT(2*K*(2*K+1))                                      
      IF(CDABS(B(K)*U/SUM).LT.1D-20) GOTO 12                           
      SUM=SUM+U*B(K)                                                   
11    CONTINUE                                                         
12    CSPEN=-SUM-1.64493406684822643D0-.5D0*CDLOG(-Z)**2               
      RETURN                                                           
20    IF(A1.GT.1D0) GOTO 30                                            
      W=-CDLOG(Z)                                                      
      SUM=W-0.25D0*W*W                                                 
      U=W                                                              
      IF(CDABS(U).LT.1D-10) GOTO 22                                    
      DO 21 K=1,9                                                      
      U=U*W*W/DFLOAT(2*K*(2*K+1))                                      
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 22                           
      SUM=SUM+U*B(K)                                                   
21    CONTINUE                                                         
22    CSPEN=-SUM+1.64493406684822643D0-CDLOG(Z)*CDLOG(1D0-Z)           
      RETURN                                                           
30    W=CDLOG(1D0-1D0/Z)                                               
      SUM=W-0.25D0*W*W                                                 
      U=W                                                              
      IF(CDABS(U).LT.1D-10) GOTO 32                                    
      DO 31 K=1,9                                                      
      U=U*W*W/DFLOAT(2*K*(2*K+1))                                      
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 32                           
      SUM=SUM+U*B(K)                                                   
31    CONTINUE                                                         
32    CSPEN=SUM+3.28986813369645287D0                                  
     *               +.5D0*CDLOG(Z-1D0)**2-CDLOG(Z)*CDLOG(1D0-Z)       
50    CONTINUE                                                         
        END                                                            

***********************************************************************
        FUNCTION ETA(C1,C2)                                            
***********************************************************************
*       COMPLEX ETA-FUNKTION                                           
*---------------------------------------------------------------------*
*       8.06.90    ANSGAR DENNER                                       
***********************************************************************
        IMPLICIT   LOGICAL(A-Z)                                        
        COMPLEX*16 ETA,C1,C2                                           
        REAL*8     PI,IM1,IM2,IM12                                     
                                                                       
        PI     = 4D0*DATAN(1D0)                                        
        IM1    = DIMAG(C1)                                             
        IM2    = DIMAG(C2)                                             
        IM12   = DIMAG(C1*C2)                                          
                                                                       
        IF(IM1.LT.0D0.AND.IM2.LT.0D0.AND.IM12.GT.0D0) THEN             
            ETA = DCMPLX(0D0,2D0*PI)                                   
        ELSE IF (IM1.GT.0D0.AND.IM2.GT.0D0.AND.IM12.LT.0D0) THEN       
            ETA = DCMPLX(0D0,-2D0*PI)                                  
        ELSE                                                           
            ETA = DCMPLX(0D0)                                          
        END IF                                                         
        END                                                            

***********************************************************************
        FUNCTION ETAS(Y,R,RS)                                            
***********************************************************************
*       MODIFIED ETA-FUNKTION                                           
*---------------------------------------------------------------------*
*       18.1.94   SD                                       
***********************************************************************
        IMPLICIT   LOGICAL(A-Z)                                        
        COMPLEX*16 ETA,ETAS,Y,R,RS
        REAL*8     PI,IMY,IMRS
                                                                       
        PI     = 4D0*DATAN(1D0)                                        

	IF( DIMAG(R).NE.0D0 ) THEN
	    ETAS = ETA(Y,R)
	ELSE	    
	    IF( DREAL(R).GT.0D0 ) THEN
		ETAS = DCMPLX(0D0,0D0)
	    ELSE
	 	IMY  = DIMAG(Y)
		IMRS = DIMAG(RS)
		ETAS = 2D0*DCMPLX(0D0,PI)*(
     *			(1D0+SIGN(1D0,-IMY))*(1D0+SIGN(1D0,-IMRS))-
     *			(1D0+SIGN(1D0, IMY))*(1D0+SIGN(1D0, IMRS))
     *					  )/4D0
	    ENDIF
	ENDIF
        END                                                            

***********************************************************************
        FUNCTION SQE(A,B,C)                                            
***********************************************************************
*       SOLUTION OF QUADRATIC EQUATION				      *
*---------------------------------------------------------------------*
*       13.1.92  SD						      *
***********************************************************************
        IMPLICIT REAL*8 (A-Z)                                        
        COMPLEX*16 A,B,C,SQE,X1,X2

	X1=(-B+SQRT(B**2-4D0*A*C))/2D0/A
	X2=(-B-SQRT(B**2-4D0*A*C))/2D0/A

	IF (ABS(X1).GT.ABS(X2)) THEN
	   SQE=X1
	ELSE
	   SQE=X2
	ENDIF

        END                                                            

************************************************************************
        FUNCTION D04(P1,P2,P3,P4,P12,P23,M1,M2,M3,M4)
************************************************************************
*  SCALAR 4-POINT FUNCTION WITH AT LEAST ONE MASS ZERO                 *
*  P1,P2,P3,P4 = SQUARED EXTERNAL MOMENTA			       *
*  P12 = (p1+p2)**2,  P23 = (p2+p3)**2				       *
*----------------------------------------------------------------------*
*  2.1.92  SD	         					       *
************************************************************************
        IMPLICIT REAL*8 (A-Z)
	REAL*8 M(4),P(4,4),K(4,4)
	COMPLEX*16 A1,A2,A3,A4,SWAP
	COMPLEX*16 SS(4), XX(2), X(2,4),RS(4,4)
	COMPLEX*16 S0(4),XX0(2),X0(2,4), R(4,4),G(2)
        COMPLEX*16 C04,D04,CSPEN,ETA,SQE,ETAS
	COMPLEX*16 AA,BB,CC,DD,IEPS,H,HH,L1,L2,L3,L4
        COMPLEX*16 Z2,B,SC,TC,WP,WM,BS,XS
	INTEGER GEN,I,J

        MM1=M1
        MM2=M2
        MM3=M3
        MM4=M4
        M12=M1*M1
        M22=M2*M2
        M32=M3*M3
        M42=M4*M4
        Q1=P1
        Q2=P2
        Q3=P3
	Q4=P4
        Q12=P12
        Q23=P23

C	IS AT LEAST ONE MASS ZERO ???
	IF (MM1*MM2*MM3*MM4.NE.0D0) GOTO 130

C	PERMUTATE UNTIL MM3=0D0
	GOTO 20
10	CONTINUE
	MM0=MM1
	MM1=MM2
	MM2=MM3
	MM3=MM4
	MM4=MM0
	M02=M12
	M12=M22
	M22=M32
	M32=M42
	M42=M02
	Q00=Q12
	Q12=Q23
	Q23=Q00
	Q0=Q1
	Q1=Q2
	Q2=Q3
	Q3=Q4
	Q4=Q0
20	IF (MM3.NE.0D0) GOTO 10
C	ONLY MM3 IS ZERO
	IF (MM1*MM2*MM4.NE.0D0) GOTO 30
C	ONLY MM3 AND MM4 ARE ZERO ==> 3->2, 4->3...
	IF ((MM1*MM2.NE.0D0).AND.(MM4.EQ.0D0)) GOTO 10
C	ONLY MM2 AND MM3 ARE ZERO
	IF ((MM1*MM4.NE.0D0).AND.(MM2.EQ.0D0)) GOTO 40
	WRITE(*,*)'CASE OF THIS SPECIAL D0-FUNCTION NOT IMPLEMENTED!'
	STOP

C	****** NO MASS EQUAL TO ZERO ******
130	CONTINUE
	EPS=1D-18
	IEPS=DCMPLX(0D0,EPS)

	IF( ABS((MM1**2+MM3**2-Q12)/MM1/MM3).LT.2D0 ) THEN
C	R13 WOULD BE NOT REAL. -> PERMUTATION! -> R(2,4) IS NOT REAL.
	   M(1)=MM2
	   M(2)=MM3
	   M(3)=MM4
	   M(4)=MM1
	   P(1,2)=Q2
	   P(1,3)=Q23
	   P(1,4)=Q1
	   P(2,3)=Q3
	   P(2,4)=Q12
	   P(3,4)=Q4
	ELSE
C	R(1,3) IS REAL.
	   M(1)=MM1
	   M(2)=MM2
	   M(3)=MM3
	   M(4)=MM4
	   P(1,2)=Q1
	   P(1,3)=Q12
	   P(1,4)=Q4
	   P(2,3)=Q2
	   P(2,4)=Q23
	   P(3,4)=Q3
	ENDIF

	DO 11 J=2,4
	DO 11 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
11	CONTINUE

	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	S0(1)=R(1,2)
	S0(2)=R(2,3)
	S0(3)=R(3,4)
	S0(4)=R(1,4)
	AA=K(3,4)/R(2,4)+R(1,3)*K(1,2)-K(1,4)*R(1,3)/R(2,4)-K(2,3)
	BB=(R(2,4)-1D0/R(2,4))*(R(1,3)-1D0/R(1,3))
     *		+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)/R(1,3)+R(2,4)*K(3,4)-K(1,4)*R(2,4)/R(1,3)-K(2,3)
	DD=K(2,3)-R(1,3)*K(1,2)-R(2,4)*K(3,4)+R(1,3)*R(2,4)*K(1,4)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	XX0(1)=SQE(AA,BB,CC)
	XX0(2)=CC/AA/XX0(1)
c	IF (ABS(DREAL(XX0(1)-XX(2))).LT.ABS(DREAL(XX0(1)-XX(1)))) THEN
	IF (ABS(XX0(1)-XX(2)).LT.ABS(XX0(1)-XX(1))) THEN
	  SWAP  =XX0(1)
	  XX0(1)=XX0(2)
	  XX0(2)=SWAP
	ENDIF

	DO 12 I=1,2
	G(I)  =SIGN( 1D0,DREAL(AA*(XX(I)-XX(3-I))) )
	 X(I,1)= XX(I)/R(2,4)
	X0(I,1)=XX0(I)/R(2,4)
	 X(I,2)= XX(I)/R(2,4)*R(1,3)
	X0(I,2)=XX0(I)/R(2,4)*R(1,3)
	 X(I,3)= XX(I)*R(1,3)
	X0(I,3)=XX0(I)*R(1,3)
	 X(I,4)= XX(I)
	X0(I,4)=XX0(I)
12	CONTINUE

	D04 = DCMPLX(0D0,0D0)
	DO 13 I=1,2
	DO 13 J=1,4
	A1 = 1D0+X0(I,J)*S0(J) + ABS(1D0+X0(I,J)*S0(J))*IEPS*
     *				  SIGN(1D0,DIMAG(X(I,J)*SS(J)))
	A2 = 1D0+X0(I,J)/S0(J) + ABS(1D0+X0(I,J)/S0(J))*IEPS*
     *				  SIGN(1D0,DIMAG(X(I,J)/SS(J)))
	D04 = D04 + (-1D0)**(I+J)*(
     *		CSPEN(A1)+ETA(-X(I,J),SS(J))*LOG(A1)
     *	       +CSPEN(A2)+ETA(-X(I,J),1D0/SS(J))*LOG(A2)     )
13	CONTINUE

	IF( DIMAG(R(1,3)).EQ.0D0 ) THEN
	DO 14 I=1,2
	   A1 = (K(1,3)-2D0*R(1,3))/XX0(I)
     *		      -R(1,3)*K(1,4)+K(3,4)
     	   A2 = ((K(2,4)-2D0*R(2,4))*R(1,3)*XX0(I)
     *		      -R(2,4)*K(3,4)+K(2,3))/DD
	   A3 = (K(1,3)-2D0*R(1,3))*R(2,4)/XX0(I)
     *		      -R(1,3)*K(1,2)+K(2,3)
	   A4 = ((K(2,4)-2D0*R(2,4))*XX0(I)
     *		      -R(2,4)*K(1,4)+K(1,2))/DD
	   L1 = LOG( A1-ABS(A1)*IEPS )
     	   L2 = LOG( A2+ABS(A2)*IEPS*G(I)*SIGN(1D0,DREAL(R(1,3))
     *				        	  *DIMAG(RS(2,4))) )
	   L3 = LOG( A3-ABS(A3)*IEPS )
	   L4 = LOG( A4+ABS(A4)*IEPS*G(I)*SIGN(1D0,DIMAG(RS(2,4))) )

	   D04 = D04 + (3D0-2D0*I)*(
     *		 ETAS(-XX(I),R(1,3),RS(1,3))
     *		   *( LOG(R(1,3)*XX(I)) + L1 + L2 )
     *		+ETAS(-XX(I),1D0/R(2,4),1D0/RS(2,4))
     *		   *( LOG(XX(I)/R(2,4)) + L3 + L4 )
     *		-( ETAS(-XX(I),R(1,3)/R(2,4),RS(1,3)/RS(2,4))
     *		  +ETA(RS(1,3),1D0/RS(2,4)) )
     *		   *( LOG(XX(I)*R(1,3)/R(2,4)) + L3 + L2 )
     *	  	+ETA(RS(1,3),1D0/RS(2,4))
     *		   *ETAS(-XX(I),-R(1,3)/R(2,4),-RS(1,3)/RS(2,4))   )
14	CONTINUE
	ELSE
	DO 15 I=1,2
	   L1 = LOG( R(2,4)/XX0(I)+XX0(I)/R(2,4)+K(1,2)
     *		     -XX0(I)/R(2,4)*EPS*BB*G(I) )
	   L2 = LOG( R(1,3)*XX0(I)+1D0/XX0(I)/R(1,3)+K(3,4)
     *		     -XX0(I)*R(1,3)*EPS*BB*G(I) )
	   L3 = LOG( R(1,3)/R(2,4)*XX0(I)+R(2,4)/XX0(I)/R(1,3)+K(2,3)
     *		     -XX0(I)*R(1,3)/R(2,4)*EPS*BB*G(I) )

	   D04 = D04 + (3D0-2D0*I)*(
     *		+ETA(-XX(I),1D0/R(2,4))
     *		   *( LOG(XX(I)/R(2,4)) + L1 )
     *		+ETA(-XX(I),R(1,3))
     *		   *( LOG(R(1,3)*XX(I)) + L2 )
     *		-( ETA(-XX(I),R(1,3)/R(2,4))
     *		  +ETA(R(1,3),1D0/R(2,4)) )
     *		   *( LOG(XX(I)*R(1,3)/R(2,4)) + L3 )
     *	  	+ETA(R(1,3),1D0/R(2,4))
     *		   *ETA(-XX(I),-R(1,3)/R(2,4))
     *		   *(1D0-G(I)*SIGN(1D0,DREAL(BB)))	    )
15	CONTINUE
	ENDIF

	D04 = D04/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
	RETURN


C--->	***************** SPEZIELL ( --> T.SACK-PROMOTION )
C	D1=Q12-M12
C	D2=Q2 -M22
C	D3=Q3 -M42
C	IF ((D1*D2.LE.0D0).OR.(D2*D3.LE.0D0)) THEN
C	   WRITE(*,*) 'THE CASE OF DIFFERENT SIGNS OF THE D1,D2,D3'
C	   WRITE(*,*) 'IN D04(...) IS NOT IMPLEMENTED !!!'
C	   STOP
C	ENDIF
C	NM1=ABS(MM1/D1)
C	NM2=ABS(MM2/D2)
C	NM3=ABS(MM4/D3)
C	NP1=Q2/D2**2+Q12/D1**2+(Q1-Q2-Q12)/D1/D2
C	NP2=Q2/D2**2+ Q3/D3**2+(Q23-Q2-Q3)/D2/D3
C	NP3=Q3/D3**2+Q12/D1**2+(Q4-Q3-Q12)/D1/D3
C	D04=C04(NP1,NP2,NP3,NM1,NM2,NM3)/D1/D2/D3

C	*************** ALLGEMEIN


C	****** ONLY MM3 IS ZERO ******
30	CONTINUE
	EPS=1D-17
	IEPS=DCMPLX(0D0,EPS)
	M(1)=MM1
	M(2)=MM2
	M(3)=10D0
	M(4)=MM4
	P(1,2)=Q1
	P(1,3)=Q12
	P(1,4)=Q4
	P(2,3)=Q2
	P(2,4)=Q23
	P(3,4)=Q3
	DO 1 J=2,4
	DO 1 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	IF (I.EQ.3) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.3) K(I,J)=K(I,J)-M(J)/M(I)
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
1	CONTINUE
	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	AA=K(3,4)/R(2,4)-K(2,3)
	BB=K(1,3)*(1D0/R(2,4)-R(2,4))+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)*K(1,3)-K(1,3)*K(1,4)*R(2,4)+R(2,4)*K(3,4)-K(2,3)
	DD=K(2,3)-R(2,4)*K(3,4)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	DO 2 I=1,2
	X(I,1)=XX(I)/R(2,4)
	X(I,2)=XX(I)/R(2,4)*R(1,3)
	X(I,3)=XX(I)*R(1,3)
	X(I,4)=XX(I)
2	CONTINUE
	D04 = DCMPLX(0D0,0D0)
	DO 3 I=1,2
	D04 = D04 + (2D0*I-3D0)*(
     *		CSPEN(1D0+SS(4)*X(I,4))
     *	       -CSPEN(1D0+SS(1)*X(I,1))
     *	       +CSPEN(1D0+X(I,4)/SS(4))
     *	       -CSPEN(1D0+X(I,1)/SS(1))
     *	       +ETA(-X(I,4),SS(4))*LOG(1D0+SS(4)*X(I,4))
     *	       -ETA(-X(I,1),SS(1))*LOG(1D0+SS(1)*X(I,1))
     *	       +ETA(-X(I,4),1D0/SS(4))*LOG(1D0+X(I,4)/SS(4))
     *	       -ETA(-X(I,1),1D0/SS(1))*LOG(1D0+X(I,1)/SS(1))
     *	       -CSPEN(1D0+X(I,4)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       +CSPEN(1D0+X(I,1)*(K(2,3)-IEPS)/(K(1,3)-IEPS))
     *	       -ETA(-X(I,4),(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+X(I,4)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       +ETA(-X(I,1),(K(2,3)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+X(I,1)*(K(2,3)-IEPS)/(K(1,3)-IEPS))   )
	IF (DIMAG(R(2,4)).NE.0D0) THEN
	   H=ETA(-1D0/XX(I),R(2,4))
	ELSE
	   H=DCMPLX(0D0,0D0)
	   IF (DREAL(R(2,4)).LT.0D0) THEN
	      HH=-1D0/XX(I)
	      IM1=DIMAG(HH)
	      IM2=DIMAG(RS(2,4))
	      IF ((IM1.GT.0D0).AND.(IM2.GT.0D0)) THEN
	         H=-DCMPLX(0D0,2D0*PI)
	      ENDIF
	      IF ((IM1.LT.0D0).AND.(IM2.LT.0D0)) THEN
	         H=+DCMPLX(0D0,2D0*PI)
	      ENDIF
	   ENDIF
	ENDIF
	D04 = D04 + (2D0*I-3D0)*
     *	          H*( LOG( (K(1,2)-R(2,4)*K(1,4)
     *			  +XX(I)*(1D0/R(2,4)-R(2,4)))/DD )
     *		     +LOG(K(1,3)-IEPS) )
3	CONTINUE
	D04 = D04/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
	RETURN

C	****** ONLY MM2 AND MM3 ARE ZERO ******
40	CONTINUE
	EPS=1D-17
	IEPS=DCMPLX(0D0,EPS)

	M(1)=MM1
	M(2)=10D0
	M(3)=10D0
	M(4)=MM4
	P(1,2)=Q1
	P(1,3)=Q12
	P(1,4)=Q4
	P(2,3)=Q2
	P(2,4)=Q23
	P(3,4)=Q3
	DO 4 J=2,4
	DO 4 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	IF (I.EQ.2) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.2) K(I,J)=K(I,J)-M(J)/M(I)
	IF (I.EQ.3) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.3) K(I,J)=K(I,J)-M(J)/M(I)
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
4	CONTINUE
	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	AA=K(2,4)*K(3,4)-K(2,3)
	BB=K(1,3)*K(2,4)+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)*K(1,3)-K(2,3)
	DD=K(2,3)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	DO 5 I=1,2
	X(I,1)=XX(I)/R(2,4)
	X(I,2)=XX(I)/R(2,4)*R(1,3)
	X(I,3)=XX(I)*R(1,3)
	X(I,4)=XX(I)
5	CONTINUE
	D04 = DCMPLX(0D0,0D0)
	DO 6 I=1,2
	D04 = D04 + (2D0*I-3D0)*(
     *		CSPEN(1D0+SS(4)*X(I,4))
     *	       +CSPEN(1D0+X(I,4)/SS(4))
     *	       +ETA(-X(I,4),SS(4))*LOG(1D0+SS(4)*X(I,4))
     *	       +ETA(-X(I,4),1D0/SS(4))*LOG(1D0+X(I,4)/SS(4))
     *	       -CSPEN(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       -CSPEN(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	       -ETA(-XX(I),(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       -ETA(-XX(I),(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	           *LOG(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	       +LOG(-XX(I))*( LOG(K(1,2)-IEPS)
     *			     +LOG(K(1,3)-IEPS)-LOG(K(2,3)-IEPS) ) )
6	CONTINUE
	D04 = D04/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))

	RETURN

	END

************************************************************************
        FUNCTION C03(P1,P2,P3,M1,M2,M3)
************************************************************************
*  SCALAR 3-POINT FUNCTION                                             *
*  P1,P2,P3 = SQUARED EXTERNAL MOMENTA  			       *
*----------------------------------------------------------------------*
*  5.12.96  M. SPIRA    					       *
************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 M1,M2,M3
      REAL*8 R(0:2)
      COMPLEX*16 C03,CSPEN,ETA,IEPS,IM
      COMPLEX*16 ALP(0:2),X(0:2,2),Y0(0:2),Y(0:2,2)
      COMPLEX*16 CDUM
C     REAL*8 KAPPA
      COMPLEX*16 KAPPA
C     KAPPA(A,B,C) = DSQRT(A**2+B**2+C**2-2*(A*B+A*C+B*C))
C     KAPPA(A,B,C) = DSQRT(DABS(A**2+B**2+C**2-2*(A*B+A*C+B*C)))
      KAPPA(A,B,C) = CDSQRT(DCMPLX(A**2+B**2+C**2-2*(A*B+A*C+B*C)))
      EPS = 1.D-8
      IM = DCMPLX(0.D0,1.D0)
      IEPS = DCMPLX(0.D0,1.D-17)
      PI = 4*DATAN(1.D0)
      XX = 0.D0
C     IF(P1.LT.0.D0.OR.P2.LT.0.D0.OR.P3.LT.0.D0) XX=1.D0
      IF(P1.NE.0.D0.OR.XX.NE.0.D0)THEN
       Q10 = P1
      ELSE
       Q10 = EPS
      ENDIF
      IF(P3.NE.0.D0.OR.XX.NE.0.D0)THEN
       Q20 = P3
      ELSE
       Q20 = EPS
      ENDIF
      IF(P2.NE.0.D0.OR.XX.NE.0.D0)THEN
       Q21 = P2
      ELSE
       Q21 = EPS
      ENDIF
      R(0) = P2
      R(1) = P3
      R(2) = P1
      SM0 = M1**2
      SM1 = M2**2
      SM2 = M3**2
      ALPHA = KAPPA(Q10,Q21,Q20)
      ALP(0) = KAPPA(Q21,SM1,SM2)*(1+IEPS*Q21)
      ALP(1) = KAPPA(Q20,SM2,SM0)*(1+IEPS*Q20)
      ALP(2) = KAPPA(Q10,SM0,SM1)*(1+IEPS*Q10)
      X(0,1) = (Q21 - SM1 + SM2 + ALP(0))/2/Q21
      X(0,2) = (Q21 - SM1 + SM2 - ALP(0))/2/Q21
      X(1,1) = (Q20 - SM2 + SM0 + ALP(1))/2/Q20
      X(1,2) = (Q20 - SM2 + SM0 - ALP(1))/2/Q20
      X(2,1) = (Q10 - SM0 + SM1 + ALP(2))/2/Q10
      X(2,2) = (Q10 - SM0 + SM1 - ALP(2))/2/Q10
      Y0(0) = (Q21*(Q21-Q20-Q10+2*SM0-SM1-SM2) - (Q20-Q10)*(SM1-SM2)
     .      + ALPHA*(Q21-SM1+SM2))/2/ALPHA/Q21
      Y0(1) = (Q20*(Q20-Q10-Q21+2*SM1-SM2-SM0) - (Q10-Q21)*(SM2-SM0)
     .      + ALPHA*(Q20-SM2+SM0))/2/ALPHA/Q20
      Y0(2) = (Q10*(Q10-Q21-Q20+2*SM2-SM0-SM1) - (Q21-Q20)*(SM0-SM1)
     .      + ALPHA*(Q10-SM0+SM1))/2/ALPHA/Q10
      Y(0,1) = Y0(0) - X(0,1)
      Y(0,2) = Y0(0) - X(0,2)
      Y(1,1) = Y0(1) - X(1,1)
      Y(1,2) = Y0(1) - X(1,2)
      Y(2,1) = Y0(2) - X(2,1)
      Y(2,2) = Y0(2) - X(2,2)
      CDUM=0.D0
      DO I=0,2
       DO J=1,2
        CDUM = CDUM + CSPEN((Y0(I)-1)/Y(I,J)) - CSPEN(Y0(I)/Y(I,J))
        CX = ETA(1-X(I,J),1/Y(I,J))
        IF(CX.NE.0.D0)THEN
         CDUM = CDUM + CX*CDLOG((Y0(I)-1)/Y(I,J))
        ENDIF
        CY = ETA(-X(I,J),1/Y(I,J))
        IF(CY.NE.0.D0)THEN
         CDUM = CDUM - CY*CDLOG(Y0(I)/Y(I,J))
        ENDIF
       ENDDO
       CX = ETA(-X(I,1),-X(I,2))
       IF(CX.NE.0.D0)THEN
        CDUM = CDUM - CX*CDLOG((1-Y0(I))/(-Y0(I)))
       ENDIF
       CY = ETA(Y(I,1),Y(I,2))
       IF(CY.NE.0.D0)THEN
        CDUM = CDUM + CY*CDLOG((1-Y0(I))/(-Y0(I)))
       ENDIF
       A = -R(I)
       B = -DIMAG(Y(I,1)*Y(I,2))
       IF(A.GT.0.D0.AND.B.GT.0.D0) THEN
        CDUM = CDUM + 2*PI*IM*CDLOG((1-Y0(I))/(-Y0(I)))
       ENDIF
      ENDDO
      C03 = CDUM/ALPHA
      RETURN
      END

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE VEGAS(FXN,BCC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
C--VEGAS SUBROUTINE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/RESU/RES
      COMMON/VEGOUT/NV
      COMMON/BVEG2/NDO,IT,SI,SI2,SWGT,SCHI,XI(50,10),SCALLS
     1,D(50,10),DI(50,10),NXI(50,10)
      DIMENSION XIN(50),R(50),DX(10),IA(10),KG(10),DT(10)
      DIMENSION XL(10),XU(10),QRAN(10),X(10)
      COMMON/RESULT/S1,S2,S3,S4
       EXTERNAL FXN
      DATA XL,XU/10*0.D0,10*1.D0/
      DATA NDMX/50/,ALPH/1.5D0/,ONE/1.D0/,MDS/1/
      IPR=1
      IF(NPRN.GT.0)IPR=0
      NDO=1
      DO 1 J=1,NDIM
1     XI(1,J)=ONE
      ENTRY VEGAS1(FXN,BCC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
      NOW=IGRAPH
CS    IF(IGRAPH.GT.0)CALL INPLOT(NOW,F1,W)
      IT=0
      SI=0.D0
      SI2=SI
      SWGT=SI
      SCHI=SI
      SCALLS=SI
      ENTRY VEGAS2(FXN,BCC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
      ND=NDMX
      NG=1
      IF(MDS.EQ.0) GO TO 2
      NG=(NCALL*0.5)**(1./NDIM)
      MDS=1
      IF((2*NG-NDMX).LT.0) GO TO 2
      MDS=-1
      NPG=NG/NDMX+1
      ND=NG/NPG
      NG=NPG*ND
2     K=NG**NDIM
      NPG=NCALL/K
      IF(NPG.LT.2)NPG=2
      CALLS=NPG*K
      DXG=ONE/NG
      DV2G=DXG**(2*NDIM)/NPG/NPG/(NPG-ONE)
      XND=ND
      NDM=ND-1
      DXG=DXG*XND
      XJAC=ONE
      DO 3 J=1,NDIM
      DX(J)=XU(J)-XL(J)
3     XJAC=XJAC*DX(J)
      IF(ND.EQ.NDO) GO TO 8
      RC=NDO/XND
      DO 7 J=1,NDIM
      K=0
      XN=0.D0
      DR=XN
      I=K
4     K=K+1
      DR=DR+ONE
      XO=XN
      XN=XI(K,J)
5     IF(RC.GT.DR) GO TO 4
      I=I+1
      DR=DR-RC
      XIN(I)=XN-(XN-XO)*DR
      IF(I.LT.NDM) GO TO 5
      DO 6  I=1,NDM
6     XI(I,J)=XIN(I)
7     XI(ND,J)=ONE
      NDO=ND
      ACC=BCC
      IF(NPRN.NE.0.AND.NPRN.NE.10)WRITE(NV,200)NDIM,CALLS,IT,ITMX
     1,ACC,MDS,ND
8     CONTINUE
      IF(NPRN.EQ.10)WRITE(NV,290)NDIM,CALLS,ITMX,ACC,MDS,ND
      ENTRY VEGAS3(FXN,BCC,NDIM,NCALL,ITMX,NPRN,IGRAPH)
9     IT=IT+1
      TI=0.D0
      TSI=TI
CS    IF(IGRAPH.GT.0)CALL REPLOT(NOW,F1,W)
      DO 10 J=1,NDIM
      KG(J)=1
      DO 10 I=1,ND
      NXI(I,J)=0
      D(I,J)=TI
10    DI(I,J)=TI
11    FB=0.D0
      F2B=FB
      K=0
12    K=K+1
      DO 121 J=1,NDIM
121   QRAN(J)=RANDM(0)
      WGT=XJAC
      DO 15 J=1,NDIM
      XN=(KG(J)-QRAN(J))*DXG+ONE
      IA(J)=XN
      IAJ=IA(J)
      IAJ1=IAJ-1
      IF(IAJ.GT.1) GO TO 13
      XO=XI(IAJ,J)
      RC=(XN-IAJ)*XO
      GO TO 14
13    XO=XI(IAJ,J)-XI(IAJ1,J)
      RC=XI(IAJ1,J)+(XN-IAJ)*XO
14    X(J)=XL(J)+RC*DX(J)
15    WGT=WGT*XO*XND
      F=FXN(X)*WGT
      F1=F/CALLS
      W=WGT/CALLS
CS    IF(IGRAPH.GT.0)CALL XPLOT(NOW,F1,W)
      F2=F**2
      FB=FB+F
      F2B=F2B+F2
      DO 16 J=1,NDIM
      IAJ=IA(J)
      NXI(IAJ,J)=NXI(IAJ,J)+1
      DI(IAJ,J)=DI(IAJ,J)+F/CALLS
16    IF(MDS.GE.0)  D(IAJ,J)=D(IAJ,J)+F2
      IF(K.LT.NPG) GO TO 12
      F2B=F2B*NPG
      F2B=SQRT(F2B)
      F2B=(F2B-FB)*(F2B+FB)
      TI=TI+FB
      TSI=TSI+F2B
      IF(MDS.GE.0) GO TO 18
      DO 17 J=1,NDIM
      IAJ=IA(J)
17    D(IAJ,J)=D(IAJ,J)+F2B
18    K=NDIM
19    KG(K)=MOD(KG(K),NG)+1
      IF(KG(K).NE.1) GO TO 11
      K=K-1
      IF(K.GT.0) GO TO 19
      TI=TI/CALLS
      TSI=TSI*DV2G
      TI2=TI*TI
      WGT=TI2/TSI
      SI=SI+TI*WGT
      SI2=SI2+TI2
      SWGT=SWGT+WGT
      SCHI=SCHI+TI2*WGT
      SCALLS=SCALLS+CALLS
      AVGI=SI/SWGT
      SD=SWGT*IT/SI2
      CHI2A=0.D0
      IF(IT.GT.1)CHI2A=SD*(SCHI/SWGT-AVGI*AVGI)/(IT-1)
      SD=ONE/SD
      SD=SQRT(SD)
      IF(NPRN.EQ.0) GO TO 21
      TSI=SQRT(TSI)
      IF(NPRN.NE.10)WRITE(NV,201)IPR,IT,TI,TSI,AVGI,SD,CHI2A
      RES=AVGI
      IF(NPRN.EQ.10)WRITE(NV,203)IT,TI,TSI,AVGI,SD,CHI2A
      IF(NPRN.GE.0) GO TO 21
      DO 20 J=1,NDIM
      WRITE(NV,202)J
20    WRITE(NV,204)(XI(I,J),DI(I,J),D(I,J),I=1,ND)
21    IF(ABS(SD/AVGI).LE.ABS(ACC).OR.IT.GE.ITMX)NOW=2
      S1=AVGI
      S2=SD
      S3=TI
      S4=TSI
CS    IF(IGRAPH.GT.0)CALL PLOTIT(NOW,F1,W)
C      DO 23 J=1,NDIM
C      XO=D(1,J)
C      XN=D(2,J)
C      D(1,J)=(XO+XN)*0.5D0
C      DT(J)=D(1,J)
C      DO 22 I=2,NDM
C      D(I,J)=XO+XN
C      XO=XN
C      XN=D(I+1,J)
C      D(I,J)=(D(I,J)+XN)/3.D0
C22    DT(J)=DT(J)+D(I,J)
C      D(ND,J)=(XN+XO)*0.5D0
C23    DT(J)=DT(J)+D(ND,J)
C-----THIS PART OF THE VEGAS-ALGORITHM IS UNSTABLE
C-----IT SHOULD BE REPLACED BY
      DO 23 J=1,NDIM
      DT(J)=0.D0
      DO 23 I=1,ND
      IF(NXI(I,J).GT.0)D(I,J)=D(I,J)/NXI(I,J)
23    DT(J)=DT(J)+D(I,J)
      DO 28 J=1,NDIM
      RC=0.D0
      DO 24 I=1,ND
      R(I)=0.D0
      IF(D(I,J).LE.0.D0)GO TO 24
      XO=DT(J)/D(I,J)
      R(I)=((XO-ONE)/XO/LOG(XO))**ALPH
24    RC=RC+R(I)
      RC=RC/XND
      K=0
      XN=0.D0
      DR=XN
      I=K
25    K=K+1
      DR=DR+R(K)
      XO=XN
      XN=XI(K,J)
26    IF(RC.GT.DR) GO TO 25
      I=I+1
      DR=DR-RC
      XIN(I)=XN-(XN-XO)*DR/R(K)
      IF(I.LT.NDM) GO TO 26
      DO 27 I=1,NDM
27    XI(I,J)=XIN(I)
28    XI(ND,J)=ONE
      IF(IT.LT.ITMX.AND.ABS(ACC).LT.ABS(SD/AVGI))GO TO 9
200   FORMAT(35H0INPUT PARAMETERS FOR VEGAS   NDIM=,I3
     1,8H  NCALL=,F8.0/28X,5H  IT=,I5,8H  ITMX =,I5/28X
     2,6H  ACC=,G9.3/28X,6H  MDS=,I3,6H   ND=,I4//)
290   FORMAT(13H0VEGAS  NDIM=,I3,8H  NCALL=,F8.0,8H  ITMX =,I5
     1,6H  ACC=,G9.3,6H  MDS=,I3,6H   ND=,I4)
201   FORMAT(/I1,20HINTEGRATION BY VEGAS/13H0ITERATION NO,I3,
     114H.   INTEGRAL =,G14.8/20X,10HSTD DEV  =,G10.4/
     234H ACCUMULATED RESULTS.   INTEGRAL =,G14.8/
     324X,10HSTD DEV  =,G10.4 / 24X,18HCHI**2 PER ITN   =,G10.4)
202   FORMAT(14H0DATA FOR AXIS,I2 / 7X,1HX,7X,10H  DELT I  ,
     12X,11H CONVCE    ,11X,1HX,7X,10H  DELT I  ,2X,11H CONVCE
     2,11X,1HX,7X,10H  DELT I  ,2X,11H CONVCE     /)
204   FORMAT(1X,3G12.4,5X,3G12.4,5X,3G12.4)
203   FORMAT(1H ,I3,G20.8,G12.4,G20.8,G12.4,G12.4)
      S1=AVGI
      S2=SD
      S3=CHI2A
      RETURN
      END

C----------------------------------------------------------------------
C  A UNIVERSAL RANDOM NUMBER GENERATOR
 
        DOUBLE PRECISION FUNCTION RANDM(IDMY)
C--RANDOM NUMBER GENERATOR
        IMPLICIT REAL*8(A-H,O-Z)
        REAL*4 UNIV
        RANDM=DBLE(UNIV(1))
        RETURN 
        END

C ---------------------------------------------------------------------

        FUNCTION UNIV(IDUM)
C--FUNCTION FOR RANDOM NUMBER GENERATOR
        REAL U(97)
        COMMON /SET1/ U,C,CD,CM,I,J
        UNIV=U(I)-U(J)
        IF(UNIV.LT.0.) UNIV=UNIV+1.
        U(I)=UNIV
        I=I-1
        IF(I.EQ.0) I=97
        J=J-1
        IF(J.EQ.0) J=97
        C=C-CD
        IF(C.LT.0.) C=C+CM
        UNIV=UNIV-C
        IF(UNIV.LT.0.) UNIV=UNIV+1
        RETURN
        END
 
C----------------------------------------------------------------------
C INITIALIZING THE RANDOM NUMBER GENERATOR
C TO INITIALIZE CALL RSTART(12,34,56,78)


        SUBROUTINE RSTART(I,J,K,L)
C--INITIALIZATION ROUTINE FOR RANDOM NUMBER GENERATOR
        REAL U(97)
        COMMON /SET1/ U,C,CD,CM,ISTART,JSTART
        IF ((I.LT.0).OR.(I.GT.178 )) STOP 'FIRST SEED .LT.0 OR .GT.178'
        IF ((J.LT.0).OR.(J.GT.178 )) STOP 'SECOND SEED .LT.0 OR .GT.178'
        IF ((K.LT.0).OR.(K.GT.178 )) STOP 'THIRD SEED .LT.0 OR .GT.178'
        IF ((L.LT.0).OR.(L.GT.168 )) STOP 'FOURTH SEED .LT.0 OR .GT.168'
        IF ( (I.EQ.1).AND.(J.EQ.1).AND.(K.EQ.1) ) STOP
     &     'FIRST, SECOND AND THIRD SEEDS ARE ALL EQUAL TO 1'
        ISTART=97
        JSTART=33
        IDUM=I
        JDUM=J
        KDUM=K
        LDUM=L
        DO 2 II=1,97
        S=0.
        T=.5
        DO 3 JJ=1,24
          M=MOD(MOD(IDUM*JDUM,179)*K,179)
          IDUM=JDUM
          JDUM=KDUM
          KDUM=M
          LDUM=MOD(53*LDUM+1,169)
          IF(MOD(LDUM*M,64).GE.32) S=S+T
3         T=.5*T
2         U(II)=S
        C=362436./16777216.
        CD=7654321./16777216.
        CM=16777213./16777216.
        RETURN
        END

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C*******************************************************************
C* PAKPDF 1.1 /21 August 1991/                                     *
C*----------------------------                                     *
C*                                                                 *
C* authors: D.Duke & J.Owens (DO)                                  *
C*          /Phys. Rev. D30 (1984) 49/                             *
C*          J.Owens (set 1.1)                                      *
C*          /preprint FSU-HEP-910606/                              *
C*                                                                 *
C* prepared by: K.Charchula, DESY                                  *
C*              bitnet: F1PCHA@DHHDESY3                            *
C*              decnet: 13313::CHARCHULA                           *
C*                                                                 *
C*part of the code adapted from PYTHIA program (author:T.Sjostrand)*
C*******************************************************************
 
      SUBROUTINE  PDDO(ISET,X,Q2,XPDF)
 
C...ISET = 1 - DO(1)  , Lambda=0.200 GeV, Nfl=4
C...       2 - DO(2)  ,       =0.400 GeV,    =4
C...       3 - DO(1.1),       =0.177 GeV,    =4
C...X          - Bjorken x
C...Q2         - square of the momentum scale  (in GeV**2)
C...XPDF(-6:6) - matrix containing  x*p(x,Q2)
C...     IPDF = -6 ,  -5 ,  -4 ,  -3 ,  -2 ,  -1 ,0 ,1,2,3,4,5,6
C...          t_bar,b_bar,c_bar,s_bar,u_bar,d_bar,gl,d,u,s,c,b,t
C...range of validity:
C...     D-04 < X  < 1
C...      4   < Q2 < D4  GeV^2
C...REAL*8 version
      REAL*8 X,Q2,XPDF(-6:6)
      DIMENSION XQ(6),TS(6)
      DIMENSION CDO(4,6,5,3)
 
C...expansion coefficients for (up+down) valence quark distribution
      DATA ((CDO(IP,IS,1,1),IS=1,6),IP=1,4)/
     1 4.190E-01, 3.460E+00, 4.400E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     2 4.000E-03, 7.240E-01,-4.860E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     3-7.000E-03,-6.600E-02, 1.330E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     4 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA ((CDO(IP,IS,1,2),IS=1,6),IP=1,4)/
     1 3.740E-01, 3.330E+00, 6.030E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     2 1.400E-02, 7.530E-01,-6.220E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     3 0.000E+00,-7.600E-02, 1.560E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     4 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA ((CDO(IP,IS,1,3),IS=1,6),IP=1,4)/
     1 6.6500E-1, 3.6140E+0, 8.6730E-1, 0.0000E+0, 0.0000E+0, 0.0000E+0,
     2-1.0970E-1, 8.3950E-1,-1.6637E+0, 1.1049E+0, 0.0000E+0, 0.0000E+0,
     3-2.4420E-3,-2.1860E-2, 3.4200E-1,-2.3690E-1, 0.0000E+0, 0.0000E+0,
     4 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00/
C...expansion coefficients for down valence quark distribution
      DATA ((CDO(IP,IS,2,1),IS=1,6),IP=1,4)/
     1 7.630E-01, 4.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     2-2.370E-01, 6.270E-01,-4.210E-01, 0.000E+00, 0.000E+00, 0.000E+00,
     3 2.600E-02,-1.900E-02, 3.300E-02, 0.000E+00, 0.000E+00, 0.000E+00,
     4 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA ((CDO(IP,IS,2,2),IS=1,6),IP=1,4)/
     1 7.610E-01, 3.830E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     2-2.320E-01, 6.270E-01,-4.180E-01, 0.000E+00, 0.000E+00, 0.000E+00,
     3 2.300E-02,-1.900E-02, 3.600E-02, 0.000E+00, 0.000E+00, 0.000E+00,
     4 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA ((CDO(IP,IS,2,3),IS=1,6),IP=1,4)/
     1 8.3880E-1, 4.6670E+0, 0.0000E+0, 0.0000E+0, 0.0000E+0, 0.0000E+0,
     2-2.0920E-1, 7.9510E-1,-1.0232E+0, 8.6160E-1, 0.0000E+0, 0.0000E+0,
     3 2.6570E-2, 1.0810E-1, 5.7990E-2, 1.5300E-1, 0.0000E+0, 0.0000E+0,
     4 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00/
C...expansion coefficients for (up+down+strange) sea quark distribution
      DATA ((CDO(IP,IS,3,1),IS=1,6),IP=1,4)/
     1 1.265E+00, 0.000E+00, 8.050E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     2-1.132E+00,-3.720E-01, 1.590E+00, 6.310E+00,-1.050E+01, 1.470E+01,
     3 2.930E-01,-2.900E-02,-1.530E-01,-2.730E-01,-3.170E+00, 9.800E+00,
     4 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA ((CDO(IP,IS,3,2),IS=1,6),IP=1,4)/
     1 1.670E+00, 0.000E+00, 9.150E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     2-1.920E+00,-2.730E-01, 5.300E-01, 1.570E+01,-1.010E+02, 2.230E+02,
     3 5.820E-01,-1.640E-01,-7.630E-01,-2.830E+00, 4.470E+01,-1.170E+02,
     4 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA ((CDO(IP,IS,3,3),IS=1,6),IP=1,4)/
     1 9.0900E-1, 0.0000E+0, 7.2780E+0, 0.0000E+0, 0.0000E+0, 0.0000E+0,
     2-4.0230E-1,-3.8230E-1,-7.9042E-1,-1.6629E+0,-1.3330E-2, 1.2110E-1,
     3 6.3050E-3, 2.7660E-2, 8.1080E-1, 5.7190E-1, 5.2990E-1,-1.7390E-1,
     4 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00/
c...expansion coefficients for charm sea quark distribution
      DATA ((CDO(IP,IS,4,1),IS=1,6),IP=1,4)/
     1 0.000E+00,-3.600E-02, 6.350E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     2 1.350E-01,-2.220E-01, 3.260E+00,-3.030E+00, 1.740E+01,-1.790E+01,
     3-7.500E-02,-5.800E-02,-9.090E-01, 1.500E+00,-1.130E+01, 1.560E+01,
     4 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00/
       DATA ((CDO(IP,IS,4,2),IS=1,6),IP=1,4)/
     1 0.000E+00,-1.200E-01, 3.510E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     2 6.700E-02,-2.330E-01, 3.660E+00,-4.740E-01, 9.500E+00,-1.660E+01,
     3-3.100E-02,-2.300E-02,-4.530E-01, 3.580E-01,-5.430E+00, 1.550E+01,
     4 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA ((CDO(IP,IS,4,3),IS=1,6),IP=1,4)/
     1 0.0000E+0,-1.4470E-1, 6.7599E+0, 0.0000E+0, 0.0000E+0, 0.0000E+0,
     2 9.4690E-2,-4.0200E-1, 1.6596E+0,-4.4559E+0, 7.8620E+0,-2.4720E-1,
     3-7.0660E-2, 1.5330E-1, 6.7980E-1, 3.3756E+0,-3.6591E+0,-7.5100E-1,
     4 1.2360E-2,-6.4790E-2,-8.5250E-1,-9.4680E-1, 3.6720E-2, 4.8700E-2/
C...expansion coefficients for gluon distribution
      DATA ((CDO(IP,IS,5,1),IS=1,6),IP=1,4)/
     1 1.560E+00, 0.000E+00, 6.000E+00, 9.000E+00, 0.000E+00, 0.000E+00,
     2-1.710E+00,-9.490E-01, 1.440E+00,-7.190E+00,-1.650E+01, 1.530E+01,
     3 6.380E-01, 3.250E-01,-1.050E+00, 2.550E-01, 1.090E+01,-1.010E+01,
     4 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA ((CDO(IP,IS,5,2),IS=1,6),IP=1,4)/
     1 8.790E-01, 0.000E+00, 4.000E+00, 9.000E+00, 0.000E+00, 0.000E+00,
     2-9.710E-01,-1.160E+00, 1.230E+00,-5.640E+00,-7.540E+00,-5.960E-01,
     3 4.340E-01, 4.760E-01,-2.540E-01,-8.170E-01, 5.500E+00, 1.260E-01,
     4 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA ((CDO(IP,IS,5,3),IS=1,6),IP=1,4)/
     1 3.0170E+0, 0.0000E+0, 5.3040E+0, 0.0000E+0, 0.0000E+0, 0.0000E+0,
     2-4.7347E+0,-9.3420E-1, 1.4654E+0,-3.9141E+0, 9.0176E+0,-5.9602E+0,
     3 3.3594E+0, 5.4540E-1,-1.4292E+0, 2.8445E+0,-1.0426E+1, 7.5150E+0,
     4-9.4430E-1,-1.6680E-1, 7.5690E-1,-8.4110E-1, 4.0983E+0,-2.7329E+0/
C...euler's beta function, requires ordinary gamma function
      PDEBE(Z,Y)=PDGAMM(Z)*PDGAMM(Y)/PDGAMM(Z+Y)
 
      XR=X
      Q2R=Q2
C...reset structure functions, check x and hadron flavour
      ALAM=0.
      DO 100 IFL=-6,6
  100 XPDF(IFL)=0.
 
C...determine set, lambda and s expansion parameter
        IF(ISET.EQ.1) ALAM=0.2
        IF(ISET.EQ.2) ALAM=0.4
        IF(ISET.EQ.3) ALAM=0.177
       SD=LOG(LOG(Q2R/ALAM**2)/LOG(4./ALAM**2))
 
C...calculate structure functions
      DO 11  IFL=1,5
       IF (IFL.EQ.1) FNO=3.
       IF (IFL.EQ.2) FNO=1.
       DO 15  IS=1,6
        TS(IS)=CDO(1,IS,IFL,ISET)+CDO(2,IS,IFL,ISET)*SD
     &         + CDO(3,IS,IFL,ISET)*SD**2 + CDO(4,IS,IFL,ISET)*SD**3
 15    CONTINUE
C...valence
       IF(IFL.LE.2) THEN
        IF(ISET.LE.2) THEN
          AN=PDEBE(TS(1),TS(2)+1.)*(1.+TS(3)*TS(1)/(TS(1)+TS(2)+1.))
        ELSE
          AN=PDEBE(TS(1),TS(2)+1.) + TS(3)*PDEBE(TS(1)+1.,TS(2)+1.)
     &          + TS(4)*PDEBE(TS(1)+2.,TS(2)+1.)
        ENDIF
        ANO=FNO/AN
        XQ(IFL)=ANO*XR**TS(1)*(1.-XR)**TS(2)*(1.+TS(3)*XR+TS(4)*XR*XR)
       ELSE
C...sea, charm and gluons
        XQ(IFL)=TS(1)*XR**TS(2)*(1.-XR)**TS(3)*(1.+TS(4)*XR+TS(5)*XR**2+
     &  TS(6)*XR**3)
       ENDIF
  11  CONTINUE
C...change of u <--> d  code
        XPDF(0)=XQ(5)
        XPDF(2)=XQ(1)-XQ(2)+XQ(3)/6.
        XPDF(1)=XQ(2)+XQ(3)/6.
        XPDF(3)=XQ(3)/6.
        XPDF(4)=XQ(4)
        XPDF(-2)=XQ(3)/6.
        XPDF(-1)=XQ(3)/6.
        XPDF(-3)=XQ(3)/6.
        XPDF(-4)=XQ(4)
 
      RETURN
      END
C---------------------------------------------
      FUNCTION PDGAMM(X)
      DIMENSION C(13)
      DATA C
     1/ 0.00053 96989 58808, 0.00261 93072 82746, 0.02044 96308 23590,
     2  0.07309 48364 14370, 0.27964 36915 78538, 0.55338 76923 85769,
     3  0.99999 99999 99998,-0.00083 27247 08684, 0.00469 86580 79622,
     4  0.02252 38347 47260,-0.17044 79328 74746,-0.05681 03350 86194,
     5  1.13060 33572 86556/
      Z=X
      IF(X .GT. 0.0) GO TO 1
      IF(X .EQ. AINT(X)) GO TO 5
      Z=1.0-Z
    1 F=1.0/Z
      IF(Z .LE. 1.0) GO TO 4
      F=1.0
    2 IF(Z .LT. 2.0) GO TO 3
      Z=Z-1.0
      F=F*Z
      GO TO 2
    3 Z=Z-1.0
    4 PDGAMM=
     1 F*((((((C(1)*Z+C(2))*Z+C(3))*Z+C(4))*Z+C(5))*Z+C(6))*Z+C(7))/
     2   ((((((C(8)*Z+C(9))*Z+C(10))*Z+C(11))*Z+C(12))*Z+C(13))*Z+1.0)
      IF(X .GT. 0.0) RETURN
      PDGAMM=3.141592653589793/(SIN(3.141592653589793*X)*PDGAMM)
      RETURN
    5 PDGAMM=0.
C     WRITE(*,10)X
      RETURN
   10 FORMAT(1X,45HGAMMA ... ARGUMENT IS NON-POSITIVE INTEGER = ,E20.5)
      END

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C******************************************************************
C* last change: 17/09/92                                          *
C*----------------------                                          *
C*                                                                *
C* authors: M.Glueck, E.Reya & A.Vogt (GRV)                       *
C*          /Z.Phys. C53 (1992) 127/                              *
C*                                                                *
C* prepared by: K.Charchula, DESY                                 *
C*              bitnet: F1PCHA@DHHDESY3                           *
C*              decnet: 13313::CHARCHULA                          *
C*                                                                *
C* code adapted from the original program provided by the authors *
C******************************************************************
 
      SUBROUTINE PDGRV(ISET,X,Q2,XPDF)
 
C...ISET = 1 - LO,          Lambda_4=0.20 GeV, N_f=5
C...       2 - NLO, MS_bar, Lambda_4=0.20 GeV, N_f=5
C...X          - Bjorken x
C...Q2         - square of the momentum scale  (in GeV**2)
C...XPDF(-6:6) - matrix containing  x*p(x,Q2)
C...     IPDF = -6 ,  -5 ,  -4 ,  -3 ,  -2 ,  -1 ,0 ,1,2,3,4,5,6
C...          t_bar,b_bar,c_bar,s_bar,u_bar,d_bar,gl,d,u,s,c,b,t
C...range of validity:
C...     D-5  < X  < 1
C...      0.3 < Q2 < D8  GeV^2
C...REAL*8 version
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XPDF(-6:6)
 
C...LO PARAMETRIZATIONS :
      IF (ISET.EQ.1) THEN
        AMU2  = 0.25
        ALAM2 = 0.232 * 0.232
        S  = LOG (LOG(Q2/ALAM2)/LOG(AMU2/ALAM2))
        S2 = S * S
        S3 = S2 * S
 
C...X * (UV + DV) :
        DNUD = 0.663 + 0.191 * S - 0.041 * S2 + 0.031 * S3
        AKUD = 0.326
        AGUD = -1.97 +  6.74 * S -  1.96 * S2
        BUD  =  24.4 -  20.7 * S +  4.08 * S2
        DUD  =  2.86 +  0.70 * S -  0.02 * S2
        UDV  = PDFV (X, DNUD, AKUD, AGUD, BUD, DUD)
 
C...X * DV :
        DND  = 0.579 + 0.283 * S + 0.047 * S2
        AKD = 0.523 - 0.015 * S
        AGD =  2.22 -  0.59 * S -  0.27 * S2
        BD  =  5.95 -  6.19 * S +  1.55 * S2
        DD  =  3.57 +  0.94 * S -  0.16 * S2
        DV  = PDFV (X, DND, AKD, AGD, BD, DD)
 
C...X * G :
        ALG =  0.558
        BEG =  1.218
        AKG =   1.00 -  0.17 * S
        BKG =   0.0
        AGG =   0.0  + 4.879 * S - 1.383 * S2
        BGG =  25.92 - 28.97 * S + 5.596 * S2
        CG  = -25.69 + 23.68 * S - 1.975 * S2
        DG  =  2.537 + 1.718 * S + 0.353 * S2
        EG  =  0.595 + 2.138 * S
        ESG =  4.066
        GL =PDFW (X, S, ALG, BEG, AKG, BKG, AGG, BGG, CG, DG, EG, ESG)
 
C...X * UBAR = X * DBAR :
        ALU =  1.396
        BEU =  1.331
        AKU =  0.412 - 0.171 * S
        BKU =  0.566 - 0.496 * S
        AGU =  0.363
        BGU = -1.196
        CU  =  1.029 + 1.785 * S - 0.459 * S2
        DU  =  4.696 + 2.109 * S
        EU  =  3.838 + 1.944 * S
        ESU =  2.845
        UDB=PDFW (X, S, ALU, BEU, AKU, BKU, AGU, BGU, CU, DU, EU, ESU)
 
C...X * SBAR = X * S :
        SS  =   0.0
        ALS =  0.803
        BES =  0.563
        AKS =  2.082 - 0.577 * S
        AGS = -3.055 + 1.024 * S **  0.67
        BS  =   27.4 -  20.0 * S ** 0.154
        DS  =   6.22
        EST =   4.33 + 1.408 * S
        ESS =   8.27 - 0.437 * S
        SB =PDFWS (X, S, SS, ALS, BES, AKS, AGS, BS, DS, EST, ESS)
 
C...X * CBAR = X * C :
        SC  =  0.888
        ALC =   1.01
        BEC =   0.37
        AKC =   0.0
        AGC =   0.0
        BC  =   4.24 - 0.804 * S
        DC  =   3.46 + 1.076 * S
        EC  =   4.61 + 1.490 * S
        ESC =  2.555 + 1.961 * S
        CB =PDFWS (X, S, SC, ALC, BEC, AKC, AGC, BC, DC, EC, ESC)
 
C...X * BBAR = X * B :
        SBO =  1.351
        ALB =   1.00
        BEB =   0.51
        AKB =   0.0
        AGB =   0.0
        BBO =  1.848
        DB  =  2.929 + 1.396 * S
        EB  =   4.71 + 1.514 * S
        ESB =   4.02 + 1.239 * S
        BB =PDFWS (X, S, SBO, ALB, BEB, AKB, AGB, BBO, DB, EB, ESB)
 
C...HO parametrization:
      ELSEIF(ISET.EQ.2) THEN
        AMU2  = 0.3
        ALAM2 = 0.248 * 0.248
        S  = LOG (LOG(Q2/ALAM2)/LOG(AMU2/ALAM2))
        DS = SQRT (S)
        S2 = S * S
        S3 = S2 * S
 
C...X * (UV + DV) :
        DNUD  = 0.330 + 0.151 * S - 0.059 * S2 + 0.027 * S3
        AKUD = 0.285
        AGUD = -2.28 + 15.73 * S -  4.58 * S2
        BUD  =  56.7 -  53.6 * S + 11.21 * S2
        DUD  =  3.17 +  1.17 * S -  0.47 * S2 +  0.09 * S3
        UDV  = PDFV (X, DNUD, AKUD, AGUD, BUD, DUD)
 
C...X * DV :
        DND  = 0.459 + 0.315 * DS + 0.515 * S
        AKD = 0.624              - 0.031 * S
        AGD =  8.13 -  6.77 * DS +  0.46 * S
        BD  =  6.59 - 12.83 * DS +  5.65 * S
        DD  =  3.98              +  1.04 * S  -  0.34 * S2
        DV  = PDFV (X, DND, AKD, AGD, BD, DD)
 
C...X * G :
        ALG =  1.128
        BEG =  1.575
        AKG =  0.323 + 1.653 * S
        BKG =  0.811 + 2.044 * S
        AGG =   0.0  + 1.963 * S - 0.519 * S2
        BGG =  0.078 +  6.24 * S
        CG  =  30.77 - 24.19 * S
        DG  =  3.188 + 0.720 * S
        EG  = -0.881 + 2.687 * S
        ESG =  2.466
        GL =PDFW (X, S, ALG, BEG, AKG, BKG, AGG, BGG, CG, DG, EG, ESG)
 
C...X * UBAR = X * DBAR :
        ALU =  0.594
        BEU =  0.614
        AKU =  0.636 - 0.084 * S
        BKU =   0.0
        AGU =  1.121 - 0.193 * S
        BGU =  0.751 - 0.785 * S
        CU  =   8.57 - 1.763 * S
        DU  =  10.22 + 0.668 * S
        EU  =  3.784 + 1.280 * S
        ESU =  1.808 + 0.980 * S
        UDB=PDFW (X, S, ALU, BEU, AKU, BKU, AGU, BGU, CU, DU, EU, ESU)
 
C...X * SBAR = X * S :
        SS  =   0.0
        ALS =  0.756
        BES =  0.101
        AKS =  2.942 - 1.016 * S
        AGS =  -4.60 + 1.167 * S
        BS  =   9.31 - 1.324 * S
        DS  =  11.49 - 1.198 * S + 0.053 * S2
        EST =  2.630 + 1.729 * S
        ESS =   8.12
        SB =PDFWS (X, S, SS, ALS, BES, AKS, AGS, BS, DS, EST, ESS)
 
C...X * CBAR = X * C :
        SC  =  0.820
        ALC =   0.98
        BEC =   0.0
        AKC = -0.625 - 0.523 * S
        AGC =   0.0
        BC  =  1.896 + 1.616 * S
        DC  =   4.12 + 0.683 * S
        EC  =   4.36 + 1.328 * S
        ESC =  0.677 + 0.679 * S
        CB =PDFWS (X, S, SC, ALC, BEC, AKC, AGC, BC, DC, EC, ESC)
 
C...X * BBAR = X * B :
        SBO =  1.297
        ALB =   0.99
        BEB =   0.0
        AKB =   0.0  - 0.193 * S
        AGB =   0.0
        BBO =   0.0
        DB  =  3.447 + 0.927 * S
        EB  =   4.68 + 1.259 * S
        ESB =  1.892 + 2.199 * S
        BB =PDFWS (X, S, SBO, ALB, BEB, AKB, AGB, BBO, DB, EB, ESB)
      ELSE
       WRITE(*,*) ' error in PDGRV: wrong ISET value'
      ENDIF
 
C...final results
      XPDF(0)=GL
      XPDF(1)=DV+UDB
      XPDF(2)=UDV-DV+UDB
      XPDF(3)=SB
      XPDF(4)=CB
      XPDF(5)=BB
      XPDF(6)=0.
      XPDF(-1)=UDB
      XPDF(-2)=UDB
      XPDF(-3)=SB
      XPDF(-4)=CB
      XPDF(-5)=BB
      XPDF(-6)=0.
 
      RETURN
      END
C-----------------------------------------------------
      FUNCTION PDFV (X, DN, AK, AG, B, D)
 
C...functional forms for ho and lo parametrizations :
      IMPLICIT REAL*8 (A-H,O-Z)
       DX = SQRT (X)
       PDFV = DN * X**AK * (1.+ AG*DX + B*X) * (1.- X)**D
      RETURN
      END
C
      FUNCTION PDFW (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
      IMPLICIT REAL*8 (A-H,O-Z)
       ALX = LOG (1./X)
       PDFW = (X**AK * (AG + X * (BG + X*C)) *ALX**BK + S**AL
     1      * EXP (-E + SQRT (ES * S**BE *ALX))) * (1.- X)**D
      RETURN
      END
C-----------------------------------------------------
      FUNCTION PDFWS (X, S, ST, AL, BE, AK, AG, B, D, E, ES)
      IMPLICIT REAL*8 (A-H,O-Z)
       DX = SQRT (X)
       ALX = LOG (1./X)
       IF (S .LE. ST) THEN
         FWS = 0.0
       ELSE
         FWS = (S-ST)**AL / ALX**AK * (1.+ AG*DX + B*X) * (1.- X)**D
     1          * EXP (-E + SQRT (ES * S**BE *ALX))
       ENDIF
       PDFWS=FWS
      RETURN
      END
 
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c     SUBROUTINE STRUC(X,Q,PDF)
C--PARTON DENSITIES
c     IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
c     DIMENSION PDF(-6:6), VALUE(20)
c     CHARACTER*20 PARM(20)
c     COMMON/PDFLIB/NGROUP,NSET,SCALFAC
c     QQ=SCALFAC*Q
c     PARM(1)='NPTYPE'
c     PARM(2)='NGROUP'
c     PARM(3)='NSET'
c     VALUE(1)=1.D0
c     VALUE(2)=DFLOAT(NGROUP)
c     VALUE(3)=DFLOAT(NSET)
c     IF(NGROUP.GE.0)THEN
c      CALL PDFSET(PARM,VALUE)
c      CALL PFTOPDG(X,QQ,PDF)
c     ELSEIF(NGROUP.EQ.-1)THEN
c      I=NSET
c      CALL PDDO(I,X,QQ**2,PDF)
c     ELSEIF(NGROUP.EQ.-2)THEN
c      I=NSET
c      CALL PDGRV(I,X,QQ**2,PDF)
c     ENDIF
c     RETURN
c     END

c     subroutine struc(x,q,pdf)
c     implicit double precision (a-h,o-z)
c     dimension pdf(-6:6), value(20)
c     character*20 parm(20)
c     character prefix*50
c     common/pdflib/ipdflib,ngroup,nset,iseterr
c     common/pdflib/ngroup,nset,scalfac
c     common/pdflib0/iseterr
c     ipdflib=ngroup
c     icase=ipdflib
c     ngroup = ipdflib
c     iseterr=0
c     if(icase.eq.1)then
c      q2=q**2
c      call pdgrv(iset,x,q2,pdf)
c      return
c     endif
c     if(ngroup.gt.0)then
c      parm(1)='nptype'
c      parm(2)='ngroup'
c      parm(3)='nset'
c      value(1)=1.d0
c      value(2)=dfloat(ngroup)
c      value(3)=dfloat(nset)
c      call pdfset(parm,value)
c      call pftopdg(x,q,pdf)
c     elseif(ngroup.eq.-1)then
c      call SetCtq6(nset)
c      pdf(6)  = 0
c      pdf(-6) = 0
c      do i=-5,5
c       j = i
c       if(i.eq.1)j=2
c       if(i.eq.2)j=1
c       if(i.eq.-1)j=-2
c       if(i.eq.-2)j=-1
c       pdf(j) = x*Ctq6Pdf(i,x,q)
c      enddo
c     elseif(ngroup.eq.-2)then
c      mode = nset
c      call mrst2001(x,q,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
c      pdf(-6) = 0
c      pdf(-5) = bot
c      pdf(-4) = chm
c      pdf(-3) = str
c      pdf(-2) = usea
c      pdf(-1) = dsea
c      pdf(0) = glu
c      pdf(1) = dnv + dsea
c      pdf(2) = upv + usea
c      pdf(3) = str
c      pdf(4) = chm
c      pdf(5) = bot
c      pdf(6) = 0
c     else
c      if(iseterr.eq.0)then
c       iset = 0
c       if(nset.eq.0)then
c        prefix = "Grids/mstw2008lo"
c       elseif(nset.eq.1)then
c        prefix = "Grids/mstw2008nlo"
c       else
c        prefix = "Grids/mstw2008nnlo"
c       endif
c      elseif(iseterr.gt.0)then
c       iset = iseterr
c       if(nset.eq.0)then
c        prefix = "Grids/mstw2008lo.90cl"
c       else
c        prefix = "Grids/mstw2008nlo.90cl"
c       endif
c      else
c       iset =-iseterr
c       if(nset.eq.0)then
c        prefix = "Grids/mstw2008lo.68cl"
c       else
c        prefix = "Grids/mstw2008nlo.68cl"
c       endif
c      endif
C--   First the traditional MRST-like interface
C--   (but note the "sbar", "cbar", "bbar" and "phot").
c      CALL GetAllPDFs(prefix,iset,x,q,upv,dnv,usea,dsea,str,sbar,
c    &        chm,cbar,bot,bbar,glu,phot)
c      pdf(-6) = 0
c      pdf(-5) = bbar
c      pdf(-4) = cbar
c      pdf(-3) = sbar
c      pdf(-2) = usea
c      pdf(-1) = dsea
c      pdf(0) = glu
c      pdf(1) = dnv + dsea
c      pdf(2) = upv + usea
c      pdf(3) = str
c      pdf(4) = chm
c      pdf(5) = bot
c      pdf(6) = 0
c     endif
c     pdf( 6) = 0
c     pdf(-6) = 0
c     return
c     end

c     subroutine struc(x,q,pdf)
c     implicit double precision (a-h,o-z)
c     dimension pdf(-6:6), value(20)
c     character*20 parm(20)
c     character*50 prefix
c     common/pdflib/ngroup,nset,scalfac
c     common/pdflib0/iseterr,ials
c     common/nloop/loop
c     if(ngroup.gt.0)then
c      parm(1)='nptype'
c      parm(2)='ngroup'
c      parm(3)='nset'
c      value(1)=1.d0
c      value(2)=dfloat(ngroup)
c      value(3)=dfloat(nset)
c      call pdfset(parm,value)
c      call pftopdg(x,q,pdf)
c     elseif(ngroup.eq.-1)then
c      call SetCtq6(nset)
c      pdf(6)  = 0
c      pdf(-6) = 0
c      do i=-5,5
c       j = i
c       if(i.eq.1)j=2
c       if(i.eq.2)j=1
c       if(i.eq.-1)j=-2
c       if(i.eq.-2)j=-1
c       pdf(j) = x*Ctq6Pdf(i,x,q)
c      enddo
c     elseif(ngroup.eq.-2)then
c      mode = nset
c      call   mrst2001(x,q,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
c      call mrst2004f4(x,q,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
c      pdf(-6) = 0
c      pdf(-5) = bot
c      pdf(-4) = chm
c      pdf(-3) = str
c      pdf(-2) = usea
c      pdf(-1) = dsea
c      pdf(0) = glu
c      pdf(1) = dnv + dsea
c      pdf(2) = upv + usea
c      pdf(3) = str
c      pdf(4) = chm
c      pdf(5) = bot
c      pdf(6) = 0
c     else
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c      if(ials.eq.0)then
c       if(iseterr.eq.0)then
c        iset = 0
c        if(nset.eq.0)then
c         prefix = "Grids/mstw2008lo"
c        else
c         prefix = "Grids/mstw2008nlo"
c        endif
c       elseif(iseterr.gt.0)then
c        iset = iseterr
c        if(nset.eq.0)then
c         prefix = "Grids/mstw2008lo.90cl"
c        else
c         prefix = "Grids/mstw2008nlo.90cl"
c        endif
c       else
c        iset =-iseterr
c        if(nset.eq.0)then
c         prefix = "Grids/mstw2008lo.68cl"
c        else
c         prefix = "Grids/mstw2008nlo.68cl"
c        endif
c       endif
c      else
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c       if(iseterr.gt.0)then
c        iset = iseterr
c        if(ials.eq.1)then
c         prefix = "Grids/mstw2008nlo_asmz-90cl.90cl"
c        elseif(ials.eq.2)then
c         prefix = "Grids/mstw2008nlo_asmz-90clhalf.90cl"
c        elseif(ials.eq.3)then
c         prefix = "Grids/mstw2008nlo_asmz+90clhalf.90cl"
c        else
c         prefix = "Grids/mstw2008nlo_asmz+90cl.90cl"
c        endif
c       else
c        iset =-iseterr
c        if(ials.eq.1)then
c         prefix = "Grids/mstw2008nlo_asmz-68cl.68cl"
c        elseif(ials.eq.2)then
c         prefix = "Grids/mstw2008nlo_asmz-68clhalf.68cl"
c        elseif(ials.eq.3)then
c         prefix = "Grids/mstw2008nlo_asmz+68clhalf.68cl"
c        else
c         prefix = "Grids/mstw2008nlo_asmz+68cl.68cl"
c        endif
c       endif
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c       if(iseterr.gt.0)then
c        iset = 0
c        if(ials.eq.1)then
c         prefix = "Grids/mstw2008nlo_asmz-90cl"
c        elseif(ials.eq.2)then
c         prefix = "Grids/mstw2008nlo_asmz-90clhalf"
c        elseif(ials.eq.3)then
c         prefix = "Grids/mstw2008nlo_asmz+90clhalf"
c        else
c         prefix = "Grids/mstw2008nlo_asmz+90cl"
c        endif
c       else
c        if(ials.eq.1)then
c         prefix = "Grids/mstw2008nlo_asmz-68cl"
c        elseif(ials.eq.2)then
c         prefix = "Grids/mstw2008nlo_asmz-68clhalf"
c        elseif(ials.eq.3)then
c         prefix = "Grids/mstw2008nlo_asmz+68clhalf"
c        else
c         prefix = "Grids/mstw2008nlo_asmz+68cl"
c        endif
c       endif
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c      endif
C--   First the traditional MRST-like interface
C--   (but note the "sbar", "cbar", "bbar" and "phot").
c      CALL GetAllPDFs(prefix,iset,x,q,upv,dnv,usea,dsea,str,sbar,
c    &        chm,cbar,bot,bbar,glu,phot)
c      pdf(-6) = 0
c      pdf(-5) = bbar
c      pdf(-4) = cbar
c      pdf(-3) = sbar
c      pdf(-2) = usea
c      pdf(-1) = dsea
c      pdf(0) = glu
c      pdf(1) = dnv + dsea
c      pdf(2) = upv + usea
c      pdf(3) = str
c      pdf(4) = chm
c      pdf(5) = bot
c      pdf(6) = 0
c     endif
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     pdf( 6) = 0
c     pdf(-6) = 0
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     return
c     end

      subroutine struc(x,q,pdf)
      implicit double precision (a-h,o-z)
      dimension pdf(-6:6), value(20)
      common/pdflib/ngroup,nset,scalfac
      ipdflib = ngroup
      if(ipdflib.eq.1)then
       q2=q**2
       call pdgrv(nset,x,q2,pdf)
      elseif(ipdflib.eq.2)then
       call SetCtq6(nset)
       pdf(6)  = 0
       pdf(-6) = 0
       do i=-5,5
        j = i
        if(i.eq.1)j=2
        if(i.eq.2)j=1
        if(i.eq.-1)j=-2
        if(i.eq.-2)j=-1
        pdf(j) = x*Ctq6Pdf(i,x,q)
       enddo
      else
       call evolvePDF(x,q,pdf)
      endif
      pdf( 6) = 0
      pdf(-6) = 0
      return
      end

      subroutine pdfset(pathname,pdfname)
      implicit double precision (a-h,o-z)
      character*100 pdfname, pathname
      common/pdflib/ngroup,nset,scalfac
      common/pdflib0/iseterr,ials
      if(ngroup.eq.0)then
       call SetPDFpath(pathname)
       call InitPDFsetByName(pdfname)
       call InitPDF(nset)
      endif
      return
      end
