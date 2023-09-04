!=----------------------------------------------------------------------------=!
MODULE fftpw_param
!=----------------------------------------------------------------------------=!

!  USE mpi_f08

  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  INTEGER, PARAMETER :: ndims = 10

END MODULE fftpw_param
