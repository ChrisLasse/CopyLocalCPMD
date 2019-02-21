#include "cpmd_global.h"

MODULE elpa_utils
#ifdef _HAS_LIBELPA
  USE error_handling,                  ONLY: stopgm
  USE mp_interface,                    ONLY: mp_min,&
                                            mp_group
  USE parac,                           ONLY: parai, paral
  use elpa


  IMPLICIT NONE

  PRIVATE
  !     ==--------------------------------------------------------------==
  INTEGER, SAVE   ::      elpa_ctxt,elpagrp,elpa_nprows
  INTEGER, SAVE   ::      elpa_npcols,elpa_my_prow,elpa_my_pcol,num_elpa
  INTEGER, PARAMETER    :: elpa_max_proc=16, elpa_apiversion=20170403
  PRIVATE :: scalapack_init
  PUBLIC  :: elpa_ob_create
 

CONTAINS
  ! ==================================================================
  SUBROUTINE elpa_ob_create(nst,el,idx,na_rows,na_cols)
    IMPLICIT NONE
    INTEGER, INTENT(IN)               :: nst  ! dimension of the problem
    INTEGER, INTENT(OUT),ALLOCATABLE  :: idx(:,:,:)  ! global/local index array
    INTEGER, INTENT(OUT)              :: na_rows,na_cols  ! rows and cols
    class(elpa_t),INTENT(OUT),pointer :: el

    !local variables
    INTEGER                          :: blocksize,colmin,rowmin,success,&
                                        i,j,k
    LOGICAL                          :: success1
    LOGICAL,SAVE                     :: initialized=.false.
    ! externals
    INTEGER, external                  :: numroc, indxl2g

    ! set up the communicator and process grid
    if(.not.initialized) then
      call scalapack_init(nst)
      initialized=.true.
    endif

    IF (parai%cp_me .LT. num_elpa ) THEN
      blocksize=64
      na_rows= numroc(nst,blocksize,ELPA_my_prow,0,ELPA_nprows)
      na_cols= numroc(nst,blocksize,ELPA_my_pcol,0,ELPA_npcols)
      colmin=na_cols
      rowmin=na_rows
      call mp_min(colmin,elpagrp)
      call mp_min(rowmin,elpagrp)
      if (colmin .eq. 0 .or. rowmin .eq.0 ) then
        if (paral%io_parent)write(6,*)&
          'WARNING: RESETTING ELPA BLOCKSIZE FROM 64 TO 32 CONSIDER USING FEWER MPI PROCS'
        blocksize=32
        na_rows= numroc(nst,blocksize,ELPA_my_prow,0,ELPA_nprows)
        na_cols= numroc(nst,blocksize,ELPA_my_pcol,0,ELPA_npcols)
      end if
      colmin=na_cols
      rowmin=na_rows
      call mp_min(colmin,elpagrp)
      call mp_min(rowmin,elpagrp)
      if (colmin .eq. 0 .or. rowmin .eq.0 ) then
        if (paral%io_parent)write(6,*) &
          'LAST WARNING: RESETTING ELPA BLOCKSIZE FROM 64 TO 16 USE FEWER MPI PROCS'
        blocksize=16
        na_rows= numroc(nst,blocksize,ELPA_my_prow,0,ELPA_nprows)
        na_cols= numroc(nst,blocksize,ELPA_my_pcol,0,ELPA_npcols)
      end if
      colmin=na_cols
      rowmin=na_rows
      call mp_min(colmin,elpagrp)
      call mp_min(rowmin,elpagrp)
      if (colmin .eq. 0 .or. rowmin .eq.0 ) &
        call stopgm('crotwf','use fewer mpi procs', &
            __LINE__,__FILE__)
      
      success1=elpa_init(elpa_apiversion)
      el => elpa_allocate()
      call el%set("na", nst, success)
      call el%set("nev", nst, success)
      call el%set("local_nrows", na_rows, success)
      call el%set("local_ncols", na_cols, success)
      call el%set("nblk", blocksize, success)
      call el%set("mpi_comm_parent", elpagrp, success)
      call el%set("process_row", ELPA_my_prow, success)
      call el%set("process_col", ELPA_my_pcol, success)
      success = el%setup()
      if (nst .ge. 1000 )then
        call el%set("solver", elpa_solver_2stage, success)
        call el%set("qr", 1, success)
        if (success .ne. 1 )call el%set("qr", 0, success)
      end if
      
      ! set up index array
      ALLOCATE(idx(na_rows,na_cols,2))
      !$OMP parallel do private(i,j,k)
      DO i=1,na_rows
        k=indxl2g(i,blocksize,ELPA_my_prow,0,ELPA_nprows)
        do j=1,na_cols
          idx(i,j,1)=k
          idx(i,j,2)=indxl2g(j,blocksize,ELPA_my_pcol,0,ELPA_NPCOLS)
        end do
      end DO
    else      
      na_rows=0
      na_cols=0
      ALLOCATE(idx(na_rows,na_cols,0))
    end if
  END SUBROUTINE

!     ==================================================================
  SUBROUTINE SCALAPACK_INIT(nst)
!     ==--------------------------------------------------------------==
    IMPLICIT NONE
    INTEGER,INTENT(IN)        :: nst
    INTEGER                   :: i
    INTEGER,ALLOCATABLE       :: elpagrouplist(:)

    if (nst .lt. 100 .or. parai%cp_nproc .lt. 4) then
       elpa_nprows=1
    elseif (nst .lt. 400) then
       elpa_nprows=2
    else
       elpa_nprows=min(elpa_max_proc,int(sqrt(dble(parai%cp_nproc))))
    end if
    elpa_npcols=elpa_nprows
    num_elpa=elpa_nprows*elpa_npcols
    ALLOCATE(elpagrouplist(num_elpa))
    
    DO i=1,num_elpa
      elpagrouplist(i)=i-1
    ENDDO
    CALL mp_group(num_elpa,elpagrouplist,elpagrp,parai%cp_grp)
    IF (parai%cp_me .LT. NUM_ELPA) THEN
      elpa_ctxt = elpagrp
      CALL blacs_gridinit(elpa_ctxt,'C',elpa_nprows,elpa_npcols)
      CALL blacs_gridinfo(elpa_ctxt,elpa_nprows,elpa_npcols,elpa_my_prow,elpa_my_pcol)
    ENDIF
    DEALLOCATE(elpagrouplist)
!       ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE
!     ==================================================================
#endif
END MODULE elpa_utils
