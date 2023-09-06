!=----------------------------------------------------------------------------=!
MODULE fftpw_param
!=----------------------------------------------------------------------------=!

!  USE mpi_f08

  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  INTEGER, PARAMETER :: ndims = 10
  REAL(SELECTED_REAL_KIND ( 14, 200 )), PARAMETER :: scal = 1.0

END MODULE fftpw_param
