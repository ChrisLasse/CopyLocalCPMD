!=----------------------------------------------------------------------=
MODULE fftpw_single_calls
!=----------------------------------------------------------------------=

! This module was mainly brought into existance in order to transfer the 
! FFT code into CPMD and therefore a comparison with the singular routines
! was thought necessary. Things here are mostly just a repeat or remixing 
! of things out of fft_batching.

  USE fftpw_batching
  USE fftpw_param
  USE fftpw_types,                              ONLY: PW_fft_type_descriptor

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Prepare_Psi_single
  PUBLIC :: invfft_single

CONTAINS

SUBROUTINE Prepare_Psi_single( dfft, psi )
  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN)  :: psi( : , : )
  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft

  INTEGER :: j, i, howmany
  INTEGER :: offset

  howmany = 2

  IF( howmany .eq. 2 ) THEN
     !$omp parallel do private( i, j, offset )
     DO i = 1, dfft%nsw( dfft%mype+1 )
        offset = (i-1)*dfft%nr3
        !$omp simd
        DO j = 1, dfft%zero_prep_start(i,1)-1
           dfft%aux( offset + j ) = psi( dfft%nl_r( offset + j ), 1 ) + (0.0d0,1.0d0) * psi( dfft%nl_r( offset + j ), 2 )
        ENDDO
        !$omp simd
        DO j = dfft%zero_prep_start(i,2), dfft%zero_prep_end(i,2)
           dfft%aux( offset + j ) = (0.d0, 0.d0)
        ENDDO
        !$omp simd
        DO j = dfft%zero_prep_end(i,1)+1, dfft%nr3
           dfft%aux( offset + j ) = psi( dfft%nl_r( offset + j ), 1 ) + (0.0d0,1.0d0) * psi( dfft%nl_r( offset + j ), 2 )
        ENDDO
     ENDDO
     !$omp end parallel do
  
     !$omp parallel do private( i, j, offset )
     DO i = 1, dfft%nsw( dfft%mype+1 )
        offset = (i-1)*dfft%nr3
        !$omp simd
        DO j = 1, dfft%zero_prep_start(i,3)-1
           dfft%aux( offset + j ) = conjg( psi( dfft%nlm_r( offset + j ), 1 ) - (0.0d0,1.0d0) * psi( dfft%nlm_r( offset + j ), 2 ) )
        ENDDO
        !$omp simd
        DO j = dfft%zero_prep_end(i,3)+1, dfft%nr3
           dfft%aux( offset + j ) = conjg( psi( dfft%nlm_r( offset + j ), 1 ) - (0.0d0,1.0d0) * psi( dfft%nlm_r( offset + j ), 2 ) )
        ENDDO
     ENDDO
     !$omp end parallel do
  ELSE
     !$omp parallel do private( i, j, offset )
     DO i = 1, dfft%nsw( dfft%mype+1 )
        offset = (i-1)*dfft%nr3
        !$omp simd
        DO j = 1, dfft%zero_prep_start(i,1)-1
           dfft%aux( offset + j ) = psi( dfft%nl_r( offset + j ), 1 )
        ENDDO
        !$omp simd
        DO j = dfft%zero_prep_start(i,2), dfft%zero_prep_end(i,2)
           dfft%aux( offset + j ) = (0.d0, 0.d0)
        ENDDO
        !$omp simd
        DO j = dfft%zero_prep_end(i,1)+1, dfft%nr3
           dfft%aux( offset + j ) = psi( dfft%nl_r( offset + j ), 1 )
        ENDDO
     ENDDO
     !$omp end parallel do
  
     !$omp parallel do private( i, j, offset )
     DO i = 1, dfft%nsw( dfft%mype+1 )
        offset = (i-1)*dfft%nr3
        !$omp simd
        DO j = 1, dfft%zero_prep_start(i,3)-1
           dfft%aux( offset + j ) = conjg( psi( dfft%nlm_r( offset + j ), 1 ) )
        ENDDO
        !$omp simd
        DO j = dfft%zero_prep_end(i,3)+1, dfft%nr3
           dfft%aux( offset + j ) = conjg( psi( dfft%nlm_r( offset + j ), 1 ) )
        ENDDO
     ENDDO
     !$omp end parallel do
  END IF

END SUBROUTINE Prepare_Psi_single

SUBROUTINE invfft_single( dfft, psi, comm_send, comm_recv, sendsize )
  IMPLICIT NONE


  TYPE(PW_fft_type_descriptor), INTENT(INOUT) :: dfft
  COMPLEX(DP), INTENT(OUT) :: psi( dfft%nnr )
  INTEGER, INTENT(IN)    :: sendsize
  COMPLEX(DP), INTENT(INOUT) :: comm_send( : ), comm_recv( : )

  

  CALL invfft_pre_com( dfft, comm_send, comm_recv, 1, 1 )
!  CALL fft_com( dfft, comm_send, comm_recv, sendsize, dfft%my_node_rank, &
!                dfft%inter_node_comm, dfft%nodes_numb, dfft%my_inter_node_rank, dfft%non_blocking )
  CALL invfft_after_com( dfft, psi, comm_recv, 1 )

END SUBROUTINE invfft_single

END MODULE fftpw_single_calls
