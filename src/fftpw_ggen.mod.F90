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
   USE fftpw_param
   USE fftpw_types,  ONLY : PW_fft_type_descriptor
   PRIVATE
   SAVE

   PUBLIC :: fft_set_nl

!=----------------------------------------------------------------------=
CONTAINS
!=----------------------------------------------------------------------=
!
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
   ALLOCATE( dfft%nlnew( ngm_max, dfft%nproc ) )
   dfft%nlnew = 0

   IF( ALLOCATED( dfft%nl ) ) DEALLOCATE( dfft%nl )
   ALLOCATE( dfft%nl( dfft%ngm ) )
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
         dfft%nlnew (ng, dfft%mype+1) = dfft%nl (ng) 
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
   CALL MPI_ALLREDUCE(MPI_IN_PLACE, dfft%nlnew, ngm_max*dfft%nproc, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

   END SUBROUTINE fft_set_nl 
   !
!=----------------------------------------------------------------------=
   END MODULE fftpw_ggen
!=----------------------------------------------------------------------=
