
MODULE iffi_types
  USE kinds,        ONLY: real_8
  
  IMPLICIT NONE
 
  !FIXME make sizes dynamic and arrays ALLOCATABLE
  INTEGER, PARAMETER :: maxind          = 250          !max number of QM atoms
  INTEGER, PARAMETER :: maxnear         = 10000       !max length of merged nearlist (=C_-3 'NEAR' and C_-2 'OLE' atoms from iffi), MXNO_MERGED_NEAR_ATOMS in iffi  SET ALSO IN epot.inc!!!!
  INTEGER, PARAMETER :: maxnearpervo    = 8000  !max number of near atoms (C_-3) per voxel
  INTEGER, PARAMETER :: maxolepervo     = 8000   !max number of ole atoms (C_-2) per voxel
  
  !.....Coefficients of the local expansion
!.....for readable access to LOCEXP entries
   INTEGER, PARAMETER :: KPHI       = 1
   INTEGER, PARAMETER :: KX         = 2
   INTEGER, PARAMETER :: KY         = 3
   INTEGER, PARAMETER :: KZ         = 4
   INTEGER, PARAMETER :: KXX        = 5
   INTEGER, PARAMETER :: KYY        = 6
   INTEGER, PARAMETER :: KXY        = 7
   INTEGER, PARAMETER :: KXZ        = 8
   INTEGER, PARAMETER :: KYZ        = 9
   INTEGER, PARAMETER :: KXXX       = 10
   INTEGER, PARAMETER :: KYYY       = 11
   INTEGER, PARAMETER :: KXYY       = 12
   INTEGER, PARAMETER :: KYXX       = 13
   INTEGER, PARAMETER :: KZXX       = 14
   INTEGER, PARAMETER :: KZYY       = 15
   INTEGER, PARAMETER :: KXYZ       = 16
   INTEGER, PARAMETER :: KXXYY      = 17
   INTEGER, PARAMETER :: KXXZZ      = 18
   INTEGER, PARAMETER :: KYYZZ      = 19
   INTEGER, PARAMETER :: KXXXY      = 20
   INTEGER, PARAMETER :: KXXXZ      = 21
   INTEGER, PARAMETER :: KYYYX      = 22
   INTEGER, PARAMETER :: KYYYZ      = 23
   INTEGER, PARAMETER :: KZZZX      = 24
   INTEGER, PARAMETER :: KZZZY      = 25
   INTEGER, PARAMETER :: dimlocexp  = KZZZY

   INTEGER, PARAMETER :: LAD        = 1
   INTEGER, PARAMETER :: PX         = 2
   INTEGER, PARAMETER :: PY         = 3
   INTEGER, PARAMETER :: PZ         = 4
   INTEGER, PARAMETER :: QXX        = 5
   INTEGER, PARAMETER :: QYY        = 6
   INTEGER, PARAMETER :: QXY        = 7
   INTEGER, PARAMETER :: QXZ        = 8
   INTEGER, PARAMETER :: QYZ        = 9
   INTEGER, PARAMETER :: OXYZ       = 10
   INTEGER, PARAMETER :: OXXY       = 11
   INTEGER, PARAMETER :: OXXZ       = 12
   INTEGER, PARAMETER :: OYYX       = 13
   INTEGER, PARAMETER :: OYYZ       = 14
   INTEGER, PARAMETER :: OZZX       = 15
   INTEGER, PARAMETER :: OZZY       = 16
   INTEGER, PARAMETER :: HXXYY      = 17
   INTEGER, PARAMETER :: HXXZZ      = 18
   INTEGER, PARAMETER :: HYYZZ      = 19
   INTEGER, PARAMETER :: HXXXY      = 20
   INTEGER, PARAMETER :: HYYYX      = 21
   INTEGER, PARAMETER :: HXXXZ      = 22
   INTEGER, PARAMETER :: HZZZX      = 23
   INTEGER, PARAMETER :: HYYYZ      = 24
   INTEGER, PARAMETER :: HZZZY      = 25
   INTEGER, PARAMETER :: dimmultexp = HZZZY
   
   INTEGER, PARAMETER :: N_KX       = 1
   INTEGER, PARAMETER :: N_KY       = 2
   INTEGER, PARAMETER :: N_KZ       = 3
   INTEGER, PARAMETER :: PCHARGE    = 4
   INTEGER, PARAMETER :: SIGCHARGE  = 5      
   INTEGER, PARAMETER :: SIGDIPOLE  = 6
   INTEGER, PARAMETER :: CORECHARGE = 7
   INTEGER, PARAMETER :: CO_WIDTH   = 8
   INTEGER, PARAMETER :: CO_VALUE   = 9
   INTEGER, PARAMETER :: sizeneardata = CO_VALUE

   INTEGER, PARAMETER :: N_PX       = 1
   INTEGER, PARAMETER :: N_PY       = 2
   INTEGER, PARAMETER :: N_PZ       = 3
   INTEGER, PARAMETER :: sizeneardata_var = N_PZ
   
!.....for readable access to MMPOTEXP entries
! ....first the fixed quantities..
   INTEGER, PARAMETER :: MPOT   = 1
   INTEGER, PARAMETER :: MEX    = 2
   INTEGER, PARAMETER :: MEY    = 3
   INTEGER, PARAMETER :: MEZ    = 4
   INTEGER, PARAMETER :: MEXX   = 5
   INTEGER, PARAMETER :: MEYY   = 6
   INTEGER, PARAMETER :: MEZZ   = 7
   INTEGER, PARAMETER :: MEXY   = 8
   INTEGER, PARAMETER :: MEXZ   = 9
   INTEGER, PARAMETER :: MEYZ   = 10
   INTEGER, PARAMETER :: dimpotexp  = MEYZ

! ....then the variable quantities..
   INTEGER, PARAMETER :: MEXPFF = 1
   INTEGER, PARAMETER :: MEYPFF = 2
   INTEGER, PARAMETER :: MEZPFF = 3
   INTEGER, PARAMETER :: dimpotexp_var  = MEZPFF

   INTEGER, PARAMETER :: RGZERO = 1
   INTEGER, PARAMETER :: RGX    = 2
   INTEGER, PARAMETER :: RGY    = 3
   INTEGER, PARAMETER :: RGZ    = 4
   INTEGER, PARAMETER :: RGXX   = 5
   INTEGER, PARAMETER :: RGYY   = 6
   INTEGER, PARAMETER :: RGZZ   = 7
   INTEGER, PARAMETER :: RGXY   = 8
   INTEGER, PARAMETER :: RGXZ   = 9
   INTEGER, PARAMETER :: RGYZ   = 10
   INTEGER, PARAMETER :: NGAMMA = 11
   INTEGER, PARAMETER :: dimrgyrexp = NGAMMA
  
!.....fixed quantities (=fixed during pff iteration)............................................................
!     interaction lists (changes each listupdate step)
   TYPE :: epot1_list_t
     INTEGER :: nnear            !number of near atoms
     INTEGER :: nloc             !number of local atomic expansions (=qm atoms)
     INTEGER :: locid(maxind)    !data for the local expansions
     INTEGER :: nincl(maxind)    !number of near atoms
     INTEGER :: idnear(maxnear)  !
     INTEGER :: lincl(maxind,maxnear) !list of near atoms to add to each local expansion
   END TYPE  epot1_list_t
   TYPE(epot1_list_t) epot1_list

   TYPE :: epot1_fix_t
     REAL(real_8) :: koloc(3,maxind)   ! Position of the Local Taylor Expansion (= qm atoms)
     REAL(real_8) :: NEARDATA    (SIZENEARDATA    ,MAXNEAR)  ! data of near atoms (charges,dipoles etc)   
     REAL(real_8) :: boxoffs(3)        ! vector from iffi to cpmd coordinate frame (needed for output of potential/rhoe)
     REAL(real_8) :: boxtrans(3)       ! translation vector when box has been recentered by IPHIGENIE
   END TYPE  epot1_fix_t
   TYPE(epot1_fix_t) epot1_fix
   
!.....variable quantities (=change during pff iteration)............................................................
   TYPE :: epot1_var_t
     LOGICAL :: TWRITERESTARTFILE
     LOGICAL :: TGRIDPART
     LOGICAL :: TDOEXTRAP
     LOGICAL :: TCALCONLYPOTFROMPFF
     LOGICAL :: TNOPOTANDWFCTUPD
     LOGICAL :: TCALCONLYPFFFIELD
     INTEGER :: CALCESPCHARGES
     LOGICAL :: TUPDATEQMINTLISTS
     LOGICAL :: TUPDATERGYR
     LOGICAL :: UPDATEMEANFIELD
     INTEGER :: NMEAN
   
     REAL(real_8) :: NEARDATA_VAR    (SIZENEARDATA_VAR    ,MAXNEAR) 
     REAL(real_8) :: LOCEXP(DIMLOCEXP,MAXIND)
   END TYPE  epot1_var_t
   TYPE(epot1_var_t) epot1_var
  
  !.....quantities that do not change during the calculation
  TYPE :: runinfo_t
    INTEGER :: nrqmatoms
    INTEGER :: sammtype
    LOGICAL :: pff
    INTEGER :: voxnum(3)
    LOGICAL :: normalmodeanalysis
    LOGICAL :: meanfieldmode
    LOGICAL :: tverbose
    
    REAL(real_8) :: boxqm(6)          !           
    REAL(real_8) :: boxdum(6)         ! box edges 
    REAL(real_8) :: OLERADIUS
    REAL(real_8) :: OLEEPS
    REAL(real_8) :: OLENSIG
    REAL(real_8) :: MYRAG(100)
  END TYPE runinfo_t
  TYPE(runinfo_t) :: runinfo
  
      REAL(real_8), ALLOCATABLE :: EXTFMM(:)
      REAL(real_8), ALLOCATABLE :: VOXEXP(:,:)
      REAL(real_8), ALLOCATABLE :: VOXEXP_MM(:,:)
      INTEGER     , ALLOCATABLE :: NEARLIST(:,:)
      INTEGER     , ALLOCATABLE :: NEARLISTLEN(:)
      INTEGER     , ALLOCATABLE :: OLELIST(:,:)
      INTEGER     , ALLOCATABLE :: OLELISTLEN(:)
      REAL(real_8), ALLOCATABLE :: QMPOTANDFLD(:,:)
      REAL(real_8), ALLOCATABLE :: ECHRG(:)
      REAL(real_8), ALLOCATABLE :: MMPOTEXP(:,:,:)
      REAL(real_8), ALLOCATABLE :: MMPOTEXP_VAR(:,:,:)
      REAL(real_8), ALLOCATABLE :: VOMULTEXP(:,:,:)
      REAL(real_8), ALLOCATABLE :: ATMULTEXP(:,:,:)
      REAL(real_8), ALLOCATABLE :: VOXRGYREXP(:,:,:)
      REAL(real_8), ALLOCATABLE :: ATRGYREXP(:,:,:)
      REAL(real_8), ALLOCATABLE :: EXTFMEAN(:)
      REAL(real_8), ALLOCATABLE :: EXTFLAST(:)
      REAL(real_8), ALLOCATABLE :: MMPOTEXPGLOB(:,:)
      REAL(real_8), ALLOCATABLE :: MMPOTEXPGLOB_VAR(:,:)
      REAL(real_8), ALLOCATABLE :: ATMULTEXPGLOB(:,:)
      REAL(real_8), ALLOCATABLE :: ATRGYREXPGLOB(:,:)
        
!....variables to handle voxels
      INTEGER :: MYNVOX
!     voxel coordinates r_lambda      
      REAL(real_8), ALLOCATABLE :: KOVOX(:,:)
!     voxel partitioning of grid
      INTEGER     , ALLOCATABLE :: VOXPART(:)
!     mapping voxels->qmatoms
      INTEGER     , ALLOCATABLE :: VO2IQ(:)
 
  END MODULE iffi_types
