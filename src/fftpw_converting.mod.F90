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
  USE mp_interface,                  ONLY: mp_sum
  USE parac,                         ONLY: parai
  USE system,                        ONLY: spar 

  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: Create_PwFFT_datastructure
  PUBLIC :: ConvertFFT_array
  PUBLIC :: ConvertFFT_Coeffs
  PUBLIC :: ConvertFFT_v

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

    ALLOCATE( dfft%time_adding( 100 ) )
    ALLOCATE( dfft%averaged_times( 100 ) )
    ALLOCATE( dfft%nnr_all( dfft%nproc ) )
    dfft%nnr_all = 0
    dfft%nnr_all( dfft%mype+1 ) = dfft%nnr
    CALL mp_sum( dfft%nnr_all, dfft%nproc, dfft%comm )
    dfft%nnr_offset = SUM( dfft%nnr_all( 1:dfft%mype ) )
    dfft%nnr_total  = SUM( dfft%nnr_all )
  
    CALL create_mpi_communicators( dfft )

!    IF( dfft%rsactive ) THEN
!       DO i = 1, dfft%nr1
!          DO j = 1, dfft%nr2 * dfft%nr3
!             
!          ENDDO
!       ENDDO
  
  END SUBROUTINE Create_PwFFT_datastructure

  SUBROUTINE ConvertFFT_array( dfft, g_pw, g_cpmd, ngw_cpmd )
    IMPLICIT NONE

    TYPE( PW_fft_type_descriptor ), INTENT(INOUT) :: dfft
    
    REAL(DP), INTENT(IN) :: g_pw(:,:), g_cpmd(:,:)
    INTEGER, INTENT(IN)  :: ngw_cpmd

    INTEGER :: c0_total( 3, dfft%ngw_total )
    REAL(DP), PARAMETER :: eps8=1.0E-8_DP
    INTEGER :: i,j, offset
    LOGICAL :: found 
   
    iloop: DO i = 1, dfft%ngw_total

       found = .false.

       jloop: DO j = 1, ngw_cpmd

          IF( abs(g_pw(1,i) - g_cpmd(1,j)) .le. eps8 .and. &
              abs(g_pw(2,i) - g_cpmd(2,j)) .le. eps8 .and. &
              abs(g_pw(3,i) - g_cpmd(3,j)) .le. eps8 ) THEN
             dfft%conv_inv( i ) = j
             dfft%conv_fw(  j ) = i
!             CYCLE iloop
             IF( found ) write(6,*) "same gvector twice"
             found = .true.
          END IF

       ENDDO jloop

    ENDDO iloop

    DO i = 1, dfft%ngw_total
       IF( dfft%conv_inv( i ) .gt. 0 ) THEN
          c0_total( :, i ) = g_cpmd( :, dfft%conv_inv( i ) )
       ELSE
          c0_total( :, i ) = (0.0_DP,0.0_DP)
       END IF
    ENDDO
   
    Call mp_sum( c0_total, 3*dfft%ngw_total, dfft%comm )

    offset = SUM( dfft%ngw_all( 1:dfft%mype ) )
    DO i = 1, dfft%ngw
       DO j = 1, 3
          IF( c0_total( j, i + offset ) .ne. g_pw( j, i + offset ) ) write(6,*) "warning missing gvec"
       ENDDO
    ENDDO

!    i = dfft%conv_inv( 1 )
!    dfft%conv_inv( 1 ) = dfft%conv_inv( 2 )
!    dfft%conv_inv( 2 ) = i
!
!    i = dfft%conv_fw( 1 )
!    dfft%conv_fw( 1 ) = dfft%conv_fw( 4 )
!    dfft%conv_fw( 4 ) = i

  END SUBROUTINE ConvertFFT_array

  SUBROUTINE ConvertFFT_Coeffs( dfft, isign, c0_in, c0_out, ngw )
    IMPLICIT NONE

    TYPE( PW_fft_type_descriptor ), INTENT(INOUT) :: dfft
    COMPLEX(DP), INTENT(IN)  :: c0_in(:)
    COMPLEX(DP), INTENT(OUT) :: c0_out(:)
    INTEGER, INTENT(IN) :: ngw, isign

    COMPLEX(DP) :: c0_total( dfft%ngw_total )
    INTEGER     :: i, offset
    
    IF( isign .eq. -1 ) THEN ! CPMD in ; PW out

       DO i = 1, dfft%ngw_total
          IF( dfft%conv_inv( i ) .gt. 0 ) THEN
             c0_total( i ) = c0_in( dfft%conv_inv( i ) )
          ELSE
             c0_total( i ) = (0.0_DP,0.0_DP)
          END IF
       ENDDO
   
       Call mp_sum( c0_total, dfft%ngw_total, dfft%comm )
   
       offset = SUM( dfft%ngw_all( 1:dfft%mype ) )
       DO i = 1, dfft%ngw
          c0_out( i ) = c0_total( i + offset )
       ENDDO

    ELSE ! PW in ; CPMD out

       c0_total = (0.0_DP,0.0_DP)       

       offset = SUM( dfft%ngw_all( 1:dfft%mype ) )
       DO i = 1, dfft%ngw
          c0_total( i + offset ) = c0_in( i )
       ENDDO

       Call mp_sum( c0_total, dfft%ngw_total, dfft%comm )

       DO i = 1, ngw
          c0_out( i ) = c0_total( dfft%conv_fw( i ) )
       ENDDO

    END IF

  END SUBROUTINE ConvertFFT_Coeffs

  SUBROUTINE ConvertFFT_v( dfft, v_in, v_out )
    IMPLICIT NONE

    TYPE( PW_fft_type_descriptor ), INTENT(INOUT) :: dfft
    REAL(DP), INTENT(IN)  :: v_in (:)
    REAL(DP), INTENT(OUT) :: v_out(:)

    REAL(DP) :: v_total(dfft%nnr_total)
    INTEGER :: i, j

    v_total = 0.0d0
    DO i = 1, dfft%nr1*dfft%nr2
       DO j = 1, dfft%my_nr3p
          v_total( j + (i-1)*dfft%nr3 + dfft%mype*dfft%my_nr3p ) = v_in( j + (i-1)*dfft%my_nr3p )
       ENDDO
    ENDDO

    CALL mp_sum( v_total, dfft%nnr_total, dfft%comm )

    DO i = 1, dfft%nnr
       v_out( i ) = v_total( i + dfft%nnr_offset )
    ENDDO

  END SUBROUTINE ConvertFFT_v

END MODULE fftpw_converting
