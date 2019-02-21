MODULE nort
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================


  ! ==================================================================
  TYPE :: nort_com_t
     REAL(real_8) :: slimit = HUGE(0.0_real_8)
     REAL(real_8) :: scond = HUGE(0.0_real_8)
     LOGICAL :: lcon
     LOGICAL :: ncon
     LOGICAL :: fcon
  END TYPE nort_com_t
  TYPE(nort_com_t), SAVE :: nort_com
  real(real_8),allocatable, save :: nort_ovlap(:,:)
END MODULE nort
