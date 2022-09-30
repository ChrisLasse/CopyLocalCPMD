!=----------------------------------------------------------------------=!
   MODULE fftpw_legacy_routines
!=----------------------------------------------------------------------=!
!! iso_c_binding provides C_PTR, C_NULL_PTR, C_ASSOCIATED
       USE, intrinsic :: iso_c_binding
       USE fftpw_param
       USE fftpw_types,                         ONLY: PW_fft_type_descriptor
       IMPLICIT NONE
       SAVE
       PRIVATE


       PUBLIC :: fft_1D
       PUBLIC :: fft_scatter_xy

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

     COMPLEX (DP), INTENT(INOUT), CONTIGUOUS :: c(:)

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

   !-----------------------------------------------------------------------
   SUBROUTINE fft_scatter_xy ( desc, f_in, f_aux, nxx_, isgn )
     !-----------------------------------------------------------------------
     !
     ! Transpose of the fft xy planes across the comm communicator.
     ! If the optional comm is not provided as input, the transpose is made
     ! across desc%comm2 communicator.
     !
     ! a) From Y-oriented columns to X-oriented partial slices (isgn > 0)
     !    Active columns along the Y direction corresponding to a subset of the
     !    active X values and a range of Z values (in this order) are stored
     !    consecutively for each processor and are such that the subgroup owns
     !    all data for a range of Z values.
     !
     !    The Y pencil -> X-oriented partial slices transposition is performed
     !    in the subgroup of processors (desc%comm2) owning this range of Z values.
     !
     !    The transpose takes place in two steps:
     !    1) on each processor the columns are sliced into sections along Y
     !       that are stored one after the other. On each processor, slices for
     !       processor "iproc2" are desc%nr2p(iproc2)*desc%nr1p(me2)*desc%my_nr3p big.
     !    2) all processors communicate to exchange slices (all sectin of columns with
     !       Y in the slice belonging to "me" must be received, all the others
     !       must be sent to "iproc2")
     !
     !    Finally one gets the "partial slice" representation: each processor has
     !    all the X values of desc%my_nr2p Y and desc%my_nr3p Z values.
     !    Data are organized with the X index running fastest, then Y, then Z.
     !
     !    f_in  contains the input Y columns, is destroyed on output
     !    f_aux contains the output X-oriented partial slices.
     !
     !  b) From planes to columns (isgn < 0)
     !
     !    Quite the same in the opposite direction
     !    f_aux contains the input X-oriented partial slices, is destroyed on output
     !    f_in  contains the output Y columns.
     !
     IMPLICIT NONE
   
     TYPE (PW_fft_type_descriptor), INTENT(in) :: desc
     INTEGER, INTENT(in)           :: nxx_, isgn
     COMPLEX (DP), INTENT(inout)   :: f_in (nxx_), f_aux (nxx_)
     INTEGER :: nr1_temp(1)
   
     !
     if ( .not. desc%vpsi ) then          ! It's a potential FFT
        CALL impl_xy( MAXVAL ( desc%nr2p ), desc%nproc2, desc%my_nr2p, desc%nr1p, desc%indp, desc%iplp, 1 )
     else     ! It's a wavefunction FFT
        CALL impl_xy( MAXVAL ( desc%nr2p ), desc%nproc2, desc%my_nr2p, desc%nr1w, desc%indw, desc%iplw, 2 )
     end if
     !
   
     RETURN
   
     CONTAINS
   
     SUBROUTINE impl_xy(nr2px, nproc2, my_nr2p, nr1p_, indx, iplx, which )
     IMPLICIT NONE
     !
     INTEGER, INTENT(in):: nr2px, nproc2, my_nr2p, which
     INTEGER, INTENT(in):: nr1p_(nproc2), indx(desc%nr1,nproc2), iplx(desc%nr1)
     !
     INTEGER :: ierr, me2, iproc2, ncpx
     INTEGER :: i, it, j, k, kfrom, kdest, mc, m1, m3, i1, icompact, sendsize, l
     Logical, allocatable, save :: map( : , : )
     Integer, allocatable, save :: map_source( : , : )
     Integer, save :: zero_start( 2 ), zero_end( 2 )
     Logical :: first
     LOGICAL, save :: already( 2 ) = .false.
     !
     me2    = desc%mype2 + 1
     ncpx = MAXVAL(nr1p_) * desc%my_nr3p       ! maximum number of Y columns to be disributed
     ! calculate the message size
     sendsize = ncpx * nr2px       ! dimension of the scattered chunks (safe value)
   
     ierr = 0
     IF (isgn.gt.0) THEN
   
        IF( .true. ) THEN !desc%use_maps ) THEN
        
           IF( .not. already( which ) ) THEN 

              IF( .not. allocated( map ) ) THEN
                 Allocate( map( nxx_, 2 ) )
                 Allocate( map_source( nxx_, 2 ) )
              END IF

              already(which) = .true.              
           
              do i = 1, nxx_
                 map(i,which) = .false.
              end do
              !$omp parallel private(it,m3,i1,m1,icompact)
              DO iproc2 = 1, nproc2
                 !$omp do 
                 DO i = 0, ncpx-1
                    IF(i>=nr1p_(iproc2)*desc%my_nr3p) CYCLE ! control i from 0 to nr1p_(iproc2)*desc%my_nr3p-1
                    it = ( iproc2 - 1 ) * sendsize + nr2px * i
                    m3 = i/nr1p_(iproc2)+1
                    i1 = mod(i,nr1p_(iproc2))+1
                    m1 = indx(i1,iproc2)
                    icompact = m1 + (m3-1)*desc%nr1*my_nr2p
                    DO j = 1, my_nr2p
                       map( icompact, which ) = .true.
                       map_source( icompact, which ) = j + it
                       icompact = icompact + desc%nr1
                    ENDDO
                 ENDDO
                 !$omp end do nowait
              ENDDO
              !$omp end parallel
              
              first = .true.
              
              do l = 1, desc%nr1
              
                 IF( map( l, which ) .eqv. .true. ) THEN
                    IF( first .eqv. .false. ) THEN
                       zero_end( which ) = l-1
                       EXIT
                    END IF
                 ELSE
                    IF( first .eqv. .true. ) THEN
                       zero_start( which ) = l
                       first = .false.
                    END IF
                 END IF
              
              end do
           
           END IF

           !$omp parallel do private( i, j, k )
           DO i = 1, desc%my_nr3p
              DO j = 1, my_nr2p
                 DO k = 1, zero_start(which)-1
                    f_aux( (i-1)*my_nr2p*desc%nr1 + (j-1)*desc%nr1 + k ) = f_in( map_source( (i-1)*my_nr2p*desc%nr1 + (j-1)*desc%nr1 + k, which ) )
                 END DO
                 DO k = zero_start(which), zero_end(which)
                    f_aux( (i-1)*my_nr2p*desc%nr1 + (j-1)*desc%nr1 + k ) = (0.0_DP, 0.0_DP)
                 END DO
                 DO k = zero_end(which)+1, desc%nr1
                    f_aux( (i-1)*my_nr2p*desc%nr1 + (j-1)*desc%nr1 + k ) = f_in( map_source( (i-1)*my_nr2p*desc%nr1 + (j-1)*desc%nr1 + k, which ) )
                 END DO
              END DO
           END DO
           !$omp end parallel do
        
        ELSE
        
           !$omp parallel
           !$omp do collapse(2) private(it,m3,i1,m1,icompact)
                DO iproc2 = 1, nproc2
                   DO i = 0, ncpx-1
                      IF(i>=nr1p_(iproc2)*desc%my_nr3p) CYCLE ! control i from 0 to nr1p_(iproc2)*desc%my_nr3p-1
                      it = ( iproc2 - 1 ) * sendsize + nr2px * i
                      m3 = i/nr1p_(iproc2)+1
                      i1 = mod(i,nr1p_(iproc2))+1
                      m1 = indx(i1,iproc2)
                      icompact = m1 + (m3-1)*desc%nr1*my_nr2p
                      DO j = 1, my_nr2p
                         !f_aux( m1 + (j-1)*desc%nr1x + (m3-1)*desc%nr1x*my_nr2p ) = f_in( j + it )
                         f_aux( icompact ) = f_in( j + it )
                         icompact = icompact + desc%nr1
                      ENDDO
                   ENDDO
                ENDDO
           !$omp end do nowait
           !$omp end parallel
!                !
!                ! clean extra array elements in each stick
!                !
!           !$omp do
!                DO k = 1, my_nr2p*desc%my_nr3p
!                   DO i1 = 1, desc%nr1x
!                      IF(iplx(i1)==0) f_aux(desc%nr1x*(k-1)+i1) = (0.0_DP, 0.0_DP)
!                   ENDDO
!                ENDDO
!           !$omp end do nowait
!           !$omp end parallel
        
        END IF
   
        !
     ELSE
              !
              !  "backward" scatter from planes to columns
              !
         !$omp parallel do collapse(2) private(it,m3,i1,m1,icompact)
              DO iproc2 = 1, nproc2
                 DO i = 0, ncpx-1
                    IF(i>=nr1p_(iproc2)*desc%my_nr3p) CYCLE ! control i from 0 to nr1p_(iproc2)*desc%my_nr3p-1
                    it = ( iproc2 - 1 ) * sendsize + nr2px * i
                    m3 = i/nr1p_(iproc2)+1
                    i1 = mod(i,nr1p_(iproc2))+1
                    m1 = indx(i1,iproc2)
                    icompact = m1 + (m3-1)*desc%nr1*my_nr2p
                    DO j = 1, my_nr2p
                       !f_in( j + it ) = f_aux( m1 + (j-1)*desc%nr1x + (m3-1)*desc%nr1x*my_nr2p )
                       f_in( j + it ) = f_aux( icompact )
                       icompact = icompact + desc%nr1
                    ENDDO
                 ENDDO
              ENDDO
         !$omp end parallel do
              !
!         !$omp parallel do collapse(2) private(kdest,kfrom)
!              DO iproc2 = 1, nproc2
!                 DO k = 0, nr1p_(me2)*desc%my_nr3p-1
!                    kdest = ( iproc2 - 1 ) * sendsize + nr2px * k
!                    kfrom = desc%nr2p_offset(iproc2) + desc%nr2 * k
!                    DO i = 1, desc%nr2p( iproc2 )
!                       f_in ( kfrom + i ) = f_aux ( kdest + i )
!                    ENDDO
!                 ENDDO
!              ENDDO
!         !$omp end parallel do
!              !
!              ! clean extra array elements in each stick
!              !
!              IF( desc%nr2x /= desc%nr2 ) THEN
!                 DO k = 1, nr1p_(me2)*desc%my_nr3p
!                    f_in(desc%nr2x*(k-1)+desc%nr2+1:desc%nr2x*k) = (0.0_DP, 0.0_DP)
!                 ENDDO
!              ENDIF
   
     ENDIF
   
     END SUBROUTINE impl_xy
   
   END SUBROUTINE fft_scatter_xy

END MODULE
