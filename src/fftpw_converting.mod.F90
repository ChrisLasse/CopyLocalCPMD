!=----------------------------------------------------------------------------=!
MODULE fftpw_converting
!=----------------------------------------------------------------------------=!

  USE elct,                          ONLY: crge
  USE error_handling,                ONLY: stopgm
  USE fft,                           ONLY: plac
  USE fftpw_base,                    ONLY: smap
  USE fftpw_ggen,                    ONLY: ggen_pw,&
                                           fft_set_nl
  USE fftpw_stick_base,              ONLY: sticks_map
  USE fftpw_types,                   ONLY: fft_type_init,&
                                           PW_fft_type_descriptor,&
                                           create_mpi_communicators
  USE gvec,                          ONLY: gvec_com
  USE mp_interface,                  ONLY: mp_sum,&
                                           mp_bcast
  USE parac,                         ONLY: parai
  USE rswfmod,                       ONLY: rsactive
  USE system,                        ONLY: cntl,&
                                           spar,&
                                           ncpw 
  USE zeroing_utils,                 ONLY: zeroing

  USE mpi_f08

  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: Create_PwFFT_datastructure
!  PUBLIC :: ConvertFFT_array
!  PUBLIC :: ConvertFFT_Coeffs
!  PUBLIC :: ConvertFFT_v
  PUBLIC :: Prep_pwFFT_Wave
  PUBLIC :: Prep_pwFFT_Rho

  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
      
CONTAINS

  SUBROUTINE Create_PwFFT_datastructure( dfft, which )
    IMPLICIT NONE
  
    TYPE( PW_fft_type_descriptor ), INTENT(INOUT) :: dfft
    CHARACTER(LEN=*), INTENT(IN) :: which
    CHARACTER(*), PARAMETER                  :: procedureN = 'create_pwfft'

    TYPE( sticks_map ) :: smap
    INTEGER :: i, gstart, offset, offset2, ierr
    INTEGER, ALLOCATABLE  :: ig_l2g_pw(:)
    REAL(DP), ALLOCATABLE :: g_pw(:,:)
    REAL(DP), ALLOCATABLE :: g_pw_all(:,:)
    REAL(DP) :: gcutw, gcutp
  
    dfft%nr1 = spar%nr1s
    dfft%nr2 = spar%nr2s
    dfft%nr3 = spar%nr3s
    dfft%nnr = dfft%nr1 * dfft%nr2 * dfft%nr3

    dfft%nstate = crge%n
    IF( cntl%krwfn ) dfft%rsactive = .true.
!    dfft%rsactive = .true.

    dfft%nproc = parai%nproc
    dfft%mype  = parai%me
    dfft%nthreads = parai%ncpus

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

    CALL fft_type_init( dfft, which, smap, .true., parai%cp_grp, dfft%bg, gcutw, gcutp )

!    IF( dfft%exact_batchsize ) THEN
!       dfft%batch_size_save = dfft%given_batchsize
!       dfft%max_batch_size = dfft%given_batchsize
!    ELSE
!       dfft%batch_size_save = dfft%max_batch_size
!    END IF
!    IF( dfft%nstate/2 .lt. dfft%max_batch_size ) THEN
!       dfft%max_batch_size = dfft%nstate/2
!       dfft%batch_size_save = dfft%max_batch_size
!    END IF
    dfft%max_batch_size = 50
    dfft%batch_size_save = 1
    dfft%buffer_size_save = 1
    dfft%x_group_autosize = dfft%max_batch_size
    dfft%y_group_autosize = dfft%max_batch_size
    dfft%z_group_autosize = dfft%max_batch_size

    dfft%overlapp = cntl%overlapp_comm_comp
    dfft%fft_tuning = cntl%fft_tune_batchsize

    dfft%nthreads = parai%ncpus

    dfft%z_set_size_save = 1
    dfft%y_set_size_save = 1
    dfft%scatter_set_size_save = 1
    dfft%x_set_size_save = 1
    dfft%apply_set_size_save = 1

    IF( dfft%rsactive ) dfft%max_buffer_size = 2

    ALLOCATE( dfft%z_loop_size( 2 ) )
    ALLOCATE( dfft%y_loop_size( 2 ) )
    ALLOCATE( dfft%x_loop_size( 4 ) )
    ALLOCATE( dfft%z_groups( dfft%max_batch_size, 2 ) )
    ALLOCATE( dfft%y_groups( dfft%max_batch_size, 2 ) )
    ALLOCATE( dfft%x_groups( dfft%max_batch_size, 4 ) )
       
    ALLOCATE( dfft%time_adding( 100 ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE( dfft%averaged_times( 100 ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE( dfft%auto_timings( 4 ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    dfft%auto_timings = 0
    ALLOCATE( dfft%optimal_groups( dfft%max_batch_size, dfft%max_buffer_size, 3 ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    dfft%optimal_groups = 0
    ALLOCATE( dfft%auto_4Stimings( 3 ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    dfft%auto_4Stimings = 0
    ALLOCATE( dfft%nnr_all( dfft%nproc ), STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    dfft%nnr_all = 0
    dfft%nnr_all( dfft%mype+1 ) = dfft%nnr
    CALL mp_sum( dfft%nnr_all, dfft%nproc, dfft%comm )
    dfft%nnr_offset = SUM( dfft%nnr_all( 1:dfft%mype ) )
    dfft%nnr_total  = SUM( dfft%nnr_all )
  
    CALL create_mpi_communicators( dfft )

    dfft%ngm_l  = ( dfft%ngl( dfft%mype + 1 ) + 1 ) / 2

    dfft%nr3px = MAXVAL ( dfft%nr3p )
    IF( dfft%what .eq. 1 ) THEN
       dfft%my_nr1p = count ( dfft%ir1p > 0 )
       dfft%small_chunks = dfft%nr3px * MAXVAL( dfft%nsp )
    ELSE
       dfft%my_nr1p = count ( dfft%ir1w > 0 )
       dfft%small_chunks = dfft%nr3px * MAXVAL( dfft%nsw )
    END IF
    dfft%big_chunks = dfft%small_chunks * dfft%node_task_size * dfft%node_task_size
    dfft%tscale = 1.0d0 / dble( dfft%nr1 * dfft%nr2 * dfft%nr3 )

    ALLOCATE( dfft%ngm_all( dfft%nproc ),STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing( dfft%ngm_all )
    dfft%ngm_all( dfft%mype + 1 ) = dfft%ngm_l
    CALL mp_sum( dfft%ngm_all, dfft%nproc, dfft%comm )

    dfft%ngm_gl = SUM( dfft%ngm_all )
    dfft%ngm_max = MAXVAL( dfft%ngm_all )

    ALLOCATE(g_pw(3,dfft%ngm_gl),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(dfft%g_pw(3,dfft%ngm_gl),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(dfft%gg_pw(dfft%ngm_gl),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(ig_l2g_pw(dfft%ngm_gl),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(g_pw)
    CALL zeroing(dfft%gg_pw)
    CALL zeroing(ig_l2g_pw)
    
    CALL ggen_pw( dfft, dfft%bg, dfft%bg, gvec_com%gcut, dfft%ngm_gl, dfft%ngm_l, & 
                  g_pw, dfft%gg_pw, ig_l2g_pw, gstart )

    dfft%g_pw = NINT( g_pw + ( ( spar%nr1s / 2 ) + 1 ) )

!    ALLOCATE( dfft%ng_all( dfft%nproc ) ,STAT=ierr)
!    CALL zeroing( dfft%ng_all )
!    IF( which == 'wave' ) THEN
!       dfft%ng_all( dfft%mype + 1 ) = dfft%ngw
!    ELSE IF( which == 'rho' ) THEN
!       dfft%ng_all( dfft%mype + 1 ) = dfft%ngm
!    END IF
!    CALL mp_sum( dfft%ng_all, dfft%nproc, dfft%comm )
!    dfft%ng_total = SUM( dfft%ng_all )
!    dfft%ng_max   = MAXVAL( dfft%ng_all )
!
!    ALLOCATE( g_pw_all( 3, dfft%ng_total ) ,STAT=ierr)
!    CALL zeroing( g_pw_all )
!    offset  = SUM( dfft%ng_all( 1:dfft%mype ) )
!    offset2 = SUM( dfft%ng_all( 1:dfft%mype+1 ) )
!    IF( which == 'wave' ) THEN
!       DO i = 1, dfft%ngw
!          g_pw_all( 1:3, i + offset ) = NINT(g_pw( 1:3, i ))
!       ENDDO
!    ELSE IF( which == 'rho' ) THEN
!       DO i = 1, dfft%ngm
!          g_pw_all( 1:3, i + offset ) = NINT(g_pw( 1:3, i ))
!       ENDDO
!    END IF
!    CALL mp_sum( g_pw_all, dfft%ng_total*3, dfft%comm )
!
!    ALLOCATE( dfft%conv_inv( dfft%ng_total ),STAT=ierr )
!    ALLOCATE( dfft%conv_fw(  dfft%ng_total ),STAT=ierr )
!    CALL zeroing( dfft%conv_inv )
!    CALL zeroing( dfft%conv_fw  )
!
!    IF( which == 'wave' ) THEN
!       CALL ConvertFFT_array( dfft, g_pw_all, dfft%g_cpmd, ncpw%ngw, dfft%ngw ) 
!    ELSE IF( which == 'rho' ) THEN
!       CALL ConvertFFT_array( dfft, g_pw_all, dfft%g_cpmd, ncpw%nhg, dfft%ngm ) 
!    END IF

    DEALLOCATE( g_pw,STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
!    DEALLOCATE( dfft%gg_pw )
    DEALLOCATE( ig_l2g_pw,STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE( dfft%ngm_all,STAT=ierr )
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
!    DEALLOCATE( g_pw_all )

  END SUBROUTINE Create_PwFFT_datastructure

  SUBROUTINE Prep_pwFFT_Wave( dfft )
    IMPLICIT NONE
  
    TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  
    CHARACTER(*), PARAMETER :: procedureN = 'Prep_pwFFT_Maps'

  
  END SUBROUTINE Prep_pwFFT_Wave

  SUBROUTINE Prep_pwFFT_Rho( dfft )
    IMPLICIT NONE
  
    TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  
    CHARACTER(*), PARAMETER :: procedureN = 'Prep_pwFFT_Maps'

  
  END SUBROUTINE Prep_pwFFT_Rho





!  SUBROUTINE ConvertFFT_array( dfft, g_pw, g_cpmd, ng_cpmd, ng_pw )
!    IMPLICIT NONE
!
!    TYPE( PW_fft_type_descriptor ), INTENT(INOUT) :: dfft
!    
!    REAL(DP), INTENT(IN) :: g_pw(:,:), g_cpmd(:,:)
!    INTEGER, INTENT(IN)  :: ng_cpmd, ng_pw
!
!    INTEGER :: c0_total( 3, dfft%ng_total )
!    REAL(DP), PARAMETER :: eps8=1.0E-8_DP
!    INTEGER :: i,j, offset
!    LOGICAL :: found 
!   
!    iloop: DO i = 1, dfft%ng_total
!
!       found = .false.
!
!       jloop: DO j = 1, ng_cpmd
!
!          IF( abs(g_pw(1,i) - g_cpmd(1,j)) .le. eps8 .and. &
!              abs(g_pw(2,i) - g_cpmd(2,j)) .le. eps8 .and. &
!              abs(g_pw(3,i) - g_cpmd(3,j)) .le. eps8 ) THEN
!             dfft%conv_inv( i ) = j
!             dfft%conv_fw(  j ) = i
!!             CYCLE iloop
!             IF( found ) write(6,*) "same gvector twice"
!             found = .true.
!          END IF
!
!       ENDDO jloop
!
!    ENDDO iloop
!
!    DO i = 1, dfft%ng_total
!       IF( dfft%conv_inv( i ) .gt. 0 ) THEN
!          c0_total( :, i ) = g_cpmd( :, dfft%conv_inv( i ) )
!       ELSE
!          c0_total( :, i ) = (0.0_DP,0.0_DP)
!       END IF
!    ENDDO
!   
!    Call mp_sum( c0_total, 3*dfft%ng_total, dfft%comm )
!
!    offset = SUM( dfft%ng_all( 1:dfft%mype ) )
!    DO i = 1, ng_pw
!       DO j = 1, 3
!          IF( c0_total( j, i + offset ) .ne. g_pw( j, i + offset ) ) write(6,*) "warning missing gvec"
!       ENDDO
!    ENDDO
!
!!    i = dfft%conv_inv( 1 )
!!    dfft%conv_inv( 1 ) = dfft%conv_inv( 2 )
!!    dfft%conv_inv( 2 ) = i
!!
!!    i = dfft%conv_fw( 1 )
!!    dfft%conv_fw( 1 ) = dfft%conv_fw( 4 )
!!    dfft%conv_fw( 4 ) = i
!
!  END SUBROUTINE ConvertFFT_array
!
!  SUBROUTINE ConvertFFT_Coeffs( dfft, isign, c0_in, c0_out, ng_cpmd, ng_pw )
!    IMPLICIT NONE
!
!    TYPE( PW_fft_type_descriptor ), INTENT(INOUT) :: dfft
!    COMPLEX(DP), INTENT(IN)  :: c0_in(:)
!    COMPLEX(DP), INTENT(OUT) :: c0_out(:)
!    INTEGER, INTENT(IN) :: ng_cpmd, ng_pw, isign
!
!    COMPLEX(DP) :: c0_total( dfft%ng_total )
!    INTEGER     :: i, offset
!    
!    IF( isign .eq. -1 ) THEN ! CPMD in ; PW out
!
!       DO i = 1, dfft%ng_total
!          IF( dfft%conv_inv( i ) .gt. 0 ) THEN
!             c0_total( i ) = c0_in( dfft%conv_inv( i ) )
!          ELSE
!             c0_total( i ) = (0.0_DP,0.0_DP)
!          END IF
!       ENDDO
!   
!       Call mp_sum( c0_total, dfft%ng_total, dfft%comm )
!   
!       offset = SUM( dfft%ng_all( 1:dfft%mype ) )
!       DO i = 1, ng_pw
!          c0_out( i ) = c0_total( i + offset )
!       ENDDO
!
!    ELSE ! PW in ; CPMD out
!
!       c0_total = (0.0_DP,0.0_DP)       
!
!       offset = SUM( dfft%ng_all( 1:dfft%mype ) )
!       DO i = 1, ng_pw
!          c0_total( i + offset ) = c0_in( i )
!       ENDDO
!
!       Call mp_sum( c0_total, dfft%ng_total, dfft%comm )
!
!       DO i = 1, ng_cpmd
!          c0_out( i ) = c0_total( dfft%conv_fw( i ) )
!       ENDDO
!
!    END IF
!
!  END SUBROUTINE ConvertFFT_Coeffs
!
!  SUBROUTINE ConvertFFT_v( dfft, v_in, v_out )
!    IMPLICIT NONE
!
!    TYPE( PW_fft_type_descriptor ), INTENT(INOUT) :: dfft
!    REAL(DP), INTENT(IN)  :: v_in (:)
!    REAL(DP), INTENT(OUT) :: v_out(:)
!
!    REAL(DP) :: v_total(dfft%nnr_total)
!    INTEGER :: i, j
!
!    v_total = 0.0d0
!    DO i = 1, dfft%nr1*dfft%nr2
!       DO j = 1, dfft%my_nr3p
!          v_total( j + (i-1)*dfft%nr3 + dfft%mype*dfft%my_nr3p ) = v_in( j + (i-1)*dfft%my_nr3p )
!       ENDDO
!    ENDDO
!
!    CALL mp_sum( v_total, dfft%nnr_total, dfft%comm )
!
!    DO i = 1, dfft%nnr
!       v_out( i ) = v_total( i + dfft%nnr_offset )
!    ENDDO
!
!  END SUBROUTINE ConvertFFT_v

END MODULE fftpw_converting
