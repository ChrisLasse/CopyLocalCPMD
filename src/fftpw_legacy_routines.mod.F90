!=----------------------------------------------------------------------=!
   MODULE fftpw_legacy_routines
!=----------------------------------------------------------------------=!
!! iso_c_binding provides C_PTR, C_NULL_PTR, C_ASSOCIATED
       USE, intrinsic :: iso_c_binding
       USE error_handling,                ONLY: stopgm
       USE fftpw_param
       USE fftpw_types,                         ONLY: PW_fft_type_descriptor
       IMPLICIT NONE
       SAVE
       PRIVATE


       PUBLIC :: fft_1D

!=----------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------=!

!=----------------------------------------------------------------------=!
!
   SUBROUTINE fft_1D( c, howmany, length, isign )

!     driver routine for "howmany" 1d complex fft's of length "length"
!     input/output  :  c(howmany * length)   (complex)
!     NOTA BENE: transform is in-place!
!     isign > 0 : forward (f(G)=>f(R)), isign <0 backward (f(R) => f(G))
!     Up to "ndims" initializations (for different combinations of input
!     parameters: howmany, length) are stored and re-used if available

     INTEGER, INTENT(IN) :: isign
     INTEGER, INTENT(IN) :: howmany, length

     COMPLEX (DP), INTENT(INOUT) :: c( length , howmany )

     INTEGER    :: i, err, idir, ip
     INTEGER, SAVE :: dims_1d( 2, ndims ) = -1
     INTEGER, SAVE :: icurrent = 1
     LOGICAL :: done
 
     !   Pointers to the "C" structures containing FFT factors ( PLAN )

     TYPE(C_PTR), SAVE :: fw_plan_1d( ndims ) = C_NULL_PTR
     TYPE(C_PTR), SAVE :: bw_plan_1d( ndims ) = C_NULL_PTR
     !
     !   Here initialize table only if necessary
     !

     CALL lookup()

     IF( .NOT. done ) THEN

       !   no table exist for these parameters
       !   initialize a new one

       CALL init_plan()

     END IF

     !
     !   Now perform the FFTs using machine specific drivers
     !
     IF (isign < 0) THEN
        CALL dfftw_execute_dft( fw_plan_1d( ip ), c, c)
     ELSE IF (isign > 0) THEN
        CALL dfftw_execute_dft( bw_plan_1d( ip ), c, c)
     END IF

     RETURN

   CONTAINS

     SUBROUTINE lookup()
        ! lookup for stored plan
        DO ip = 1, ndims
           !   first check if there is already a table initialized
           !   for this combination of parameters
           done = ( howmany == dims_1d(1,ip) ) .AND. ( length == dims_1d(2,ip) )
           IF (done) EXIT
        END DO
     END SUBROUTINE lookup

     SUBROUTINE init_plan()
       logical :: FFTW_IN_PLACE
       logical :: FFTW_ESTIMATE
       INTEGER :: void
       INTEGER, EXTERNAL :: omp_get_max_threads
       !
       CALL dfftw_cleanup_threads()
       CALL dfftw_init_threads( void )
       CALL dfftw_plan_with_nthreads(omp_get_max_threads())

       IF( C_ASSOCIATED(fw_plan_1d( icurrent)) ) CALL dfftw_destroy_plan( fw_plan_1d( icurrent) )
       IF( C_ASSOCIATED(bw_plan_1d( icurrent)) ) CALL dfftw_destroy_plan( bw_plan_1d( icurrent) )
       idir = -1
       CALL dfftw_plan_many_dft( fw_plan_1d( icurrent), 1, length, howmany, c, &
            (/SIZE(c)/), 1, length, c, (/SIZE(c)/), 1, length, idir, FFTW_ESTIMATE, FFTW_IN_PLACE)
       idir = 1
       CALL dfftw_plan_many_dft( bw_plan_1d( icurrent), 1, length, howmany, c, &
            (/SIZE(c)/), 1, length, c, (/SIZE(c)/), 1, length, idir, FFTW_ESTIMATE, FFTW_IN_PLACE)

       dims_1d(1,icurrent) = howmany; dims_1d(2,icurrent) = length;
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1
     END SUBROUTINE init_plan

   END SUBROUTINE fft_1D

END MODULE
