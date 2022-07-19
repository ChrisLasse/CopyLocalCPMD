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

!    IF( dfft%rsactive ) THEN
!       DO i = 1, dfft%nr1
!          DO j = 1, dfft%nr2 * dfft%nr3
!             
!          ENDDO
!       ENDDO
  
  END SUBROUTINE Create_PwFFT_datastructure

  SUBROUTINE ConvertFFTCoeffs( dfft, c0, c0_pw )
    IMPLICIT NONE

    TYPE( PW_fft_type_descriptor ), INTENT(INOUT) :: dfft
    COMPLEX(DP), INTENT(IN) :: c0
    COMPLEX(DP), INTENT(OUT) :: c0_pw
    
    
    inyh_pw = inyh

    ind = 0

    DO i = 1, dfft%ngw
       
       IF( dfft%isind .eq. 0 .or. dfft%isind .gt. dfft%ngw ) CYCLE

       y = position bestimmen anhand von position array (maybe in smap)
       z = position bestimmen anhand von position array (maybe in smap)

       DO j = 1, dfft%... 506
       
          IF( g(2,j) .eq. y .and. g(3,j) .eq. z ) THEN
          
             ind = ind + 1
             g_pw(1,ind) = g(1,j)
             g_pw(2,ind) = g(2,j)
             g_pw(3,ind) = g(3,j)

          END IF

       ENDDO

    ENDDO

  END SUBROUTINE ConvertFFTCoeffs

END MODULE fftpw_converting
