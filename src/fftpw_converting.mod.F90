!=----------------------------------------------------------------------------=!
MODULE fftpw_converting
!=----------------------------------------------------------------------------=!

  USE fftpw_base,                    ONLY: smap
  USE fftpw_param
  USE fftpw_stick_base,              ONLY: sticks_map
  USE fftpw_types,                   ONLY: fft_type_init,&
                                           PW_fft_type_descriptor,&
                                           create_mpi_communicators
  USE gvec,                          ONLY: gvec_com
  USE parac,                         ONLY: parai
  USE system,                        ONLY: spar 

  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: Create_PwFFT_datastructure

CONTAINS

  SUBROUTINE Create_PwFFT_datastructure( dfft )
    IMPLICIT NONE
  
    TYPE( PW_fft_type_descriptor ), INTENT(INOUT) :: dfft

    TYPE( sticks_map ) :: smap
    INTEGER :: i
    REAL(DP) :: gcutw, gcutp
  
    dfft%nr1 = spar%nr1s
    dfft%nr2 = spar%nr2s
    dfft%nr3 = spar%nr3s
    dfft%nnr = dfft%nr1 * dfft%nr2 * dfft%nr3

    dfft%nproc = parai%nproc
    dfft%mype  = parai%me

    dfft%comm            = MPI_COMM_NULL
    dfft%node_comm       = MPI_COMM_NULL
    dfft%inter_node_comm = MPI_COMM_NULL

!    dfft%comm = parai%cp_grp

    gcutw = gvec_com%gcutw
    gcutp = gvec_com%gcut

    DO i = 1, 3
       dfft%bg(i,1) = gvec_com%b1(i)
       dfft%bg(i,2) = gvec_com%b2(i)
       dfft%bg(i,3) = gvec_com%b3(i)
    ENDDO

    CALL fft_type_init( dfft, smap, .true., parai%cp_grp, dfft%bg, gcutw, gcutp )
  
    CALL create_mpi_communicators( dfft )
  
  END SUBROUTINE Create_PwFFT_datastructure


END MODULE fftpw_converting
