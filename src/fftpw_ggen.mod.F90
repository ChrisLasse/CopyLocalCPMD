!
! Copyright (C) 2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!=----------------------------------------------------------------------=
MODULE fftpw_ggen
!=----------------------------------------------------------------------=

  !  ... subroutines generating variables nl* needed to map G-vector
  !  ... components onto the FFT grid(s) in reciprocal space
  !
   USE error_handling,                  ONLY: stopgm
   USE fftpw_param
   USE fftpw_types,  ONLY : PW_fft_type_descriptor
   PRIVATE
   SAVE

   PUBLIC :: ggen_pw, fft_set_nl

!=----------------------------------------------------------------------=
CONTAINS
!=----------------------------------------------------------------------=
  SUBROUTINE ggen_pw ( dfft, at, bg,  gcutm, ngm_g, ngm, &
       g, gg, ig_l2g, gstart )
    !----------------------------------------------------------------------
    !! This routine generates all the reciprocal lattice vectors
    !! contained in the sphere of radius gcutm. Furthermore it
    !! computes the indices nl which give the correspondence
    !! between the fft mesh points and the array of g vectors.
    !
!    USE mp, ONLY: mp_rank, mp_size, mp_sum
    !
    IMPLICIT NONE
    !
    TYPE(PW_fft_type_descriptor),INTENT(INOUT) :: dfft
    REAL(DP), INTENT(IN) :: at(3,3), bg(3,3), gcutm
    INTEGER, INTENT(IN) :: ngm_g
    INTEGER, INTENT(INOUT) :: ngm
    REAL(DP), INTENT(OUT) :: g(:,:), gg(:)
    INTEGER, INTENT(OUT) :: ig_l2g(:), gstart
    !  if no_global_sort is present (and it is true) G vectors are sorted only
    !  locally and not globally. In this case no global array needs to be
    !  allocated and sorted: saves memory and a lot of time for large systems.
    !
    !     here a few local variables
    !
    REAL(DP), PARAMETER :: eps8=1.0E-8_DP
    REAL(DP) :: tx(3), ty(3), t(3)
    REAL(DP), ALLOCATABLE :: tt(:)
    INTEGER :: ngm_save, n1, n2, n3, ngm_offset, ngm_max, ngm_local
    !
    REAL(DP), ALLOCATABLE :: g2sort_g(:)
    ! array containing only g vectors for the current processor
    INTEGER, ALLOCATABLE :: mill_unsorted(:,:)
    ! array containing all g vectors generators, on all processors
    ! (replicated data). When no_global_sort is present and .true.,
    ! only g-vectors for the current processor are stored
    INTEGER, ALLOCATABLE :: igsrt(:), g2l(:)
    !
    INTEGER :: ni, nj, nk, i, j, k, ipol, ng, igl, indsw
    INTEGER :: istart, jstart, kstart
    INTEGER :: mype, npe
    LOGICAL :: global_sort, is_local
    INTEGER, ALLOCATABLE :: ngmpe(:)
    !
    ngm_max = ngm_g
    !
    ! save current value of ngm
    !
    ngm_save  = ngm
    !
    ngm = 0
    ngm_local = 0
    !
    !    set the total number of fft mesh points and and initial value of gg
    !    The choice of gcutm is due to the fact that we have to order the
    !    vectors after computing them.
    !
    gg(:) = gcutm + 1.d0
    !
    !    and computes all the g vectors inside a sphere
    !
    ALLOCATE( mill_unsorted( 3, ngm_save ) )
    ALLOCATE( igsrt( ngm_max ) )
    ALLOCATE( g2l( ngm_max ) )
    ALLOCATE( g2sort_g( ngm_max ) )
    !
    g2sort_g(:) = 1.0d20
    !
    ! allocate temporal array
    !
    ALLOCATE( tt( dfft%nr3 ) )
    !
    ! max miller indices (same convention as in module stick_set)
    !
    ni = (dfft%nr1-1)/2
    nj = (dfft%nr2-1)/2
    nk = (dfft%nr3-1)/2
    !
    ! gamma-only: exclude space with x < 0
    !
    istart = 0
    !
    iloop: DO i = istart, ni
       !
       ! gamma-only: exclude plane with x = 0, y < 0
       !
       IF ( i == 0 ) THEN
          jstart = 0
       ELSE
          jstart = -nj
       ENDIF
       !
       tx(1:3) = i * bg(1:3,1)
       !
       jloop: DO j = jstart, nj
          !
          IF ( dfft%lpara .AND. fft_stick_index( dfft, i, j ) == 0) THEN
             is_local = .FALSE.
          ELSE
             is_local = .TRUE.
          END IF
          !
          ! gamma-only: exclude line with x = 0, y = 0, z < 0
          !
          IF ( i == 0 .and. j == 0 ) THEN
             kstart = 0
          ELSE
             kstart = -nk
          ENDIF
          !
          ty(1:3) = tx(1:3) + j * bg(1:3,2)
          !
          !  compute all the norm square
          !
          DO k = kstart, nk
             !
             t(1) = ty(1) + k * bg(1,3)
             t(2) = ty(2) + k * bg(2,3)
             t(3) = ty(3) + k * bg(3,3)
             tt(k-kstart+1) = t(1)**2 + t(2)**2 + t(3)**2
          ENDDO
          !
          !  save all the norm square within cutoff
          !
          DO k = kstart, nk
             IF (tt(k-kstart+1) <= gcutm) THEN
                ngm = ngm + 1
                IF (ngm > ngm_max) write(6,*) "GVEC ERROR"
                IF ( tt(k-kstart+1) > eps8 ) THEN
                   g2sort_g(ngm) = tt(k-kstart+1)
                ELSE
                   g2sort_g(ngm) = 0.d0
                ENDIF
                IF (is_local) THEN
                  ngm_local = ngm_local + 1
                  mill_unsorted( :, ngm_local ) = (/ i,j,k /)
                  g2l(ngm) = ngm_local
                ELSE
                  g2l(ngm) = 0
                ENDIF
             ENDIF
          ENDDO
       ENDDO jloop
    ENDDO iloop
    IF (ngm  /= ngm_max) write(6,*) "GVEC ERROR"
    !
    igsrt(1) = 0
    CALL hpsort_eps( ngm_g, g2sort_g, igsrt, eps8 )
    DEALLOCATE( g2sort_g, tt )
    
    ngm = 0
    !
    ngloop: DO ng = 1, ngm_max
       !
       IF (g2l(igsrt(ng))>0) THEN
          ! fetch the indices
          i = mill_unsorted(1, g2l(igsrt(ng)))
          j = mill_unsorted(2, g2l(igsrt(ng)))
          k = mill_unsorted(3, g2l(igsrt(ng)))
          !
          ngm = ngm + 1
          !
          !  Here map local and global g index !!! N.B: :
          !  the global G vectors arrangement depends on the number of processors
          !
          ig_l2g( ngm ) = ng
       
          g(1:3, ngm) = i * bg (:, 1) + j * bg (:, 2) + k * bg (:, 3)
          gg(ngm) = sum(g(1:3, ngm)**2)
       ENDIF
    ENDDO ngloop

    DEALLOCATE( igsrt, g2l )

    IF (ngm /= ngm_save) write(6,*) "GVEC ERROR"
    !
    !     determine first nonzero g vector
    !
    IF (gg(1).le.eps8) THEN
       gstart=2
    ELSE
       gstart=1
    ENDIF
    !
    !     Now set nl and nls with the correct fft correspondence
    !
    CALL fft_set_nl( dfft, at, g )
    !
  END SUBROUTINE ggen_pw
!
  PURE FUNCTION fft_stick_index( desc, i, j )
     IMPLICIT NONE
     TYPE(PW_fft_type_descriptor), INTENT(IN) :: desc
     INTEGER :: fft_stick_index
     INTEGER, INTENT(IN) :: i, j
     INTEGER :: mc, m1, m2
     m1 = mod (i, desc%nr1) + 1
     IF (m1 < 1) m1 = m1 + desc%nr1
     m2 = mod (j, desc%nr2) + 1
     IF (m2 < 1) m2 = m2 + desc%nr2
     mc = m1 + (m2 - 1) * desc%nr1
     fft_stick_index = desc%isind ( mc )
  END FUNCTION
!-----------------------------------------------------------------------
   SUBROUTINE fft_set_nl ( dfft, at, g )
!----------------------------------------------------------------------
   !
   ! Input:  FFT descriptor dfft, lattice vectors at, list of G-vectors g
   ! Output: indices nl such that G_fft(nl(i)) = G(i)
   !         indices nlm such that G_fft(nlm(i)) = -G(i) only if lgamma=.true.
   !         optionally, Miller indices: if bg = reciprocal lattice vectors,
   ! G(:,i) = mill(1,i)*bg(:,1) + mill(2,i)*bg(:,2) + mill(3,i)*bg(:,3)
   !  
   !
   IMPLICIT NONE
   CHARACTER(*), PARAMETER                  :: procedureN = 'fft_set_nl'
   !
   TYPE (PW_fft_type_descriptor), INTENT(inout) :: dfft
   REAL(DP), INTENT(IN) :: g(:,:)
   REAL(DP), INTENT(IN) :: at(:,:)
   INTEGER :: ng, n1, n2, n3, ierr, ngm_max
   Integer :: l
   !
   ngm_max = MAXVAL(dfft%ngl)
   l=0
   !
   IF( ALLOCATED( dfft%nlnew ) ) DEALLOCATE( dfft%nlnew )
   ALLOCATE( dfft%nlnew( ngm_max, dfft%nproc3 ) )
   dfft%nlnew = 0

   IF( ALLOCATED( dfft%nl ) ) DEALLOCATE( dfft%nl )
   ALLOCATE( dfft%nl( dfft%ngm ),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
   IF( ALLOCATED( dfft%nl_r ) ) DEALLOCATE( dfft%nl_r )
   ALLOCATE( dfft%nl_r( dfft%nnr ) )
   dfft%nl_r = 0
   IF( ALLOCATED( dfft%nlm ) ) DEALLOCATE( dfft%nlm )
   ALLOCATE( dfft%nlm( dfft%ngm ) )
   IF( ALLOCATED( dfft%nlm_r ) ) DEALLOCATE( dfft%nlm_r )
   ALLOCATE( dfft%nlm_r( dfft%nnr ) )
   !
   DO ng = 1, dfft%ngm
      n1 = nint (sum(g (:, ng) * at (:, 1)))
      IF (n1<0) n1 = n1 + dfft%nr1

      n2 = nint (sum(g (:, ng) * at (:, 2)))
      IF (n2<0) n2 = n2 + dfft%nr2

      n3 = nint (sum(g (:, ng) * at (:, 3)))
      IF (n3<0) n3 = n3 + dfft%nr3

      IF ( dfft%lpara ) THEN
         dfft%nl (ng) = 1 + n3 + ( dfft%isind ( 1 + n1 + n2*dfft%nr1) - 1) * dfft%nr3
         dfft%nl_r( dfft%nl (ng) ) = ng
         dfft%nlnew (ng, dfft%mype3+1) = dfft%nl (ng) 
      ELSE
         dfft%nl (ng) = 1 + n1 + n2 * dfft%nr1 + n3 * dfft%nr1 * dfft%nr2
         dfft%nl_r( dfft%nl (ng) ) = ng
      ENDIF

      n1 = - n1 ; IF (n1<0) n1 = n1 + dfft%nr1
      n2 = - n2 ; IF (n2<0) n2 = n2 + dfft%nr2
      n3 = - n3 ; IF (n3<0) n3 = n3 + dfft%nr3

      IF ( dfft%lpara ) THEN
         dfft%nlm(ng) = 1 + n3 + ( dfft%isind ( 1 + n1 + n2*dfft%nr1) - 1) * dfft%nr3
         dfft%nlm_r( dfft%nlm( ng ) ) = ng
      ELSE
         dfft%nlm(ng) = 1 + n1 + n2 * dfft%nr1 + n3 * dfft%nr1 * dfft%nr2
         dfft%nlm_r( dfft%nlm( ng ) ) = ng
      ENDIF

   ENDDO
   !
   CALL MPI_ALLREDUCE(MPI_IN_PLACE, dfft%nlnew, ngm_max*dfft%nproc3, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

   END SUBROUTINE fft_set_nl 
   !
   !---------------------------------------------------------------------
   SUBROUTINE hpsort_eps (n, ra, ind, eps)
     !---------------------------------------------------------------------
     !! Sort an array ra(1:n) into ascending order using heapsort algorithm,
     !! and considering two elements being equal if their values differ
     !! for less than "eps".  
     !! \(\text{n}\) is input, \(\text{ra}\) is replaced on output by its 
     !! sorted rearrangement.  
     !! Create an index table (ind) by making an exchange in the index array
     !! whenever an exchange is made on the sorted data array (\(\text{ra}\)).  
     !! In case of equal values in the data array (\(\text{ra}\)) the values
     !! in the index array (ind) are used to order the entries.  
     !! If on input ind(1) = 0 then indices are initialized in the routine,
     !! if on input ind(1) != 0 then indices are assumed to have been
     !! initialized before entering the routine and these indices are carried
     !! around during the sorting process.
     !
     ! no work space needed !
     ! free us from machine-dependent sorting-routines !
     !
     ! adapted from Numerical Recipes pg. 329 (new edition)
     !
     implicit none  
     !-input/output variables
     integer, intent(in) :: n  
     integer, intent(inout) :: ind (*)  
     real(DP), intent(inout) :: ra (*)
     real(DP), intent(in) :: eps
     !-local variables
     integer :: i, ir, j, l, iind  
     real(DP) :: rra  
     ! initialize index array
     if (ind (1) .eq.0) then  
        do i = 1, n  
           ind (i) = i  
        enddo
     endif
     ! nothing to order
     if (n.lt.2) return  
     ! initialize indices for hiring and retirement-promotion phase
     l = n / 2 + 1  
   
     ir = n  
   
     sorting: do 
     
       ! still in hiring phase
       if ( l .gt. 1 ) then  
          l    = l - 1  
          rra  = ra (l)  
          iind = ind (l)  
          ! in retirement-promotion phase.
       else  
          ! clear a space at the end of the array
          rra  = ra (ir)  
          !
          iind = ind (ir)  
          ! retire the top of the heap into it
          ra (ir) = ra (1)  
          !
          ind (ir) = ind (1)  
          ! decrease the size of the corporation
          ir = ir - 1  
          ! done with the last promotion
          if ( ir .eq. 1 ) then  
             ! the least competent worker at all !
             ra (1)  = rra  
             !
             ind (1) = iind  
             exit sorting  
          endif
       endif
       ! wheter in hiring or promotion phase, we
       i = l  
       ! set up to place rra in its proper level
       j = l + l  
       !
       do while ( j .le. ir )  
          if ( j .lt. ir ) then  
             ! compare to better underling
             if ( abs(ra(j)-ra(j+1)).ge.eps ) then  
                if (ra(j).lt.ra(j+1)) j = j + 1
             else
                ! this means ra(j) == ra(j+1) within tolerance
                if (ind (j) .lt.ind (j + 1) ) j = j + 1
             endif
          endif
          ! demote rra
          if ( abs(rra - ra(j)).ge.eps ) then  
             if (rra.lt.ra(j)) then
                ra (i) = ra (j)  
                ind (i) = ind (j)  
                i = j  
                j = j + j  
             else
                ! set j to terminate do-while loop
                j = ir + 1  
             end if
          else
             !this means rra == ra(j) within tolerance
             ! demote rra
             if (iind.lt.ind (j) ) then
                ra (i) = ra (j)
                ind (i) = ind (j)
                i = j
                j = j + j
             else
                ! set j to terminate do-while loop
                j = ir + 1
             endif
          end if
       enddo
       ra (i) = rra  
       ind (i) = iind  
   
     end do sorting    
     !
   END SUBROUTINE hpsort_eps


   !
!=----------------------------------------------------------------------=
   END MODULE fftpw_ggen
!=----------------------------------------------------------------------=
