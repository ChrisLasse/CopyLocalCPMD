#include "cpmd_global.h"
#if defined(__FFT_HASNT_THREADED_COPIES)
#define HASNT_THREADED_COPIES .TRUE.
#else
#define HASNT_THREADED_COPIES .FALSE.
#endif

#if defined(__FFT_HAS_LOW_LEVEL_TIMERS)
#define HAS_LOW_LEVEL_TIMERS .TRUE.
#else
#define HAS_LOW_LEVEL_TIMERS .FALSE.
#endif

#if defined(__FFT_HAS_SPECIAL_COPY)
#define HAS_SPECIAL_COPY __FFT_HAS_SPECIAL_COPY
#else
#define HAS_SPECIAL_COPY 0
#endif

#if defined(__FFT_HAS_OMP_COLLAPSE)
#define __COLLAPSE2 collapse(2)
#else
#define __COLLAPSE2
#endif


MODULE fftutil_utils
  USE cppt,                            ONLY: nzh_r,&
                                             indz_r
  USE fft,                             ONLY: FFT_TYPE_DESCRIPTOR
  USE fft_maxfft,                      ONLY: maxfft
  USE kinds,                           ONLY: real_4,&
                                             real_8
  USE mltfft_utils,                    ONLY: mltfft_fftw_t
  USE mp_interface,                    ONLY: mp_all2all,&
                                             mp_startall,&
                                             mp_waitall
  USE parac,                           ONLY: parai
  USE reshaper,                        ONLY: type_cast,&
                                             reshape_inplace
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             parap,&
                                             spar,&
                                             parm,&
                                             cntl
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zgthr_no_omp,&
                                             zsctr_no_omp
  USE zeroing_utils,                   ONLY: zeroing
  USE iso_fortran_env

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: phase
  PUBLIC :: putz
  PUBLIC :: getz
  PUBLIC :: unpack_y2x
  PUBLIC :: pack_y2x
  PUBLIC :: fft_comm
  PUBLIC :: pack_x2y
  PUBLIC :: unpack_x2y
  PUBLIC :: phasen
!TK special routines for batched fft
  PUBLIC :: getz_n
  PUBLIC :: putz_n
  PUBLIC :: pack_x2y_n
  PUBLIC :: unpack_x2y_n
  PUBLIC :: pack_y2x_n
  PUBLIC :: unpack_y2x_n
!TK
!CR special routines for improved FFT
  LOGICAL, ALLOCATABLE, SAVE :: locks_omp(:,:,:)
  PUBLIC :: locks_omp
  PUBLIC :: Prepare_Psi
  PUBLIC :: fft_com
  PUBLIC :: invfft_z_section
  PUBLIC :: invfft_y_section
  PUBLIC :: invfft_x_section
  PUBLIC :: fwfft_x_section
  PUBLIC :: fwfft_y_section
  PUBLIC :: fwfft_z_section
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  REAL(DP), PARAMETER :: scal = 1.0
!CR
CONTAINS


  ! ==================================================================
  SUBROUTINE phase(f)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8) :: f(fpar%kr1,fpar%kr2s,fpar%kr3s)

    INTEGER                                  :: ii, ijk1, ijk2, isub, j, k

! ==--------------------------------------------------------------==

    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset('     PHASE',isub)
    ijk2=parap%nrxpl(parai%mepos,2)-parap%nrxpl(parai%mepos,1)+1 ! jh-mb
    !$omp parallel do private (J,K,II,IJK1) shared(IJK2) __COLLAPSE2
    DO k=1,spar%nr3s
       DO j=1,spar%nr2s
          ii=k+j+parap%nrxpl(parai%mepos,1)
          ijk1=MOD(ii,2)+1
          DO ii=ijk1,ijk2,2
             f(ii,j,k)=-f(ii,j,k)
          ENDDO
       ENDDO
    ENDDO
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt('     PHASE',isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE phase
  ! ==================================================================
  SUBROUTINE putz(a,b,krmin,krmax,kr,m)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: krmin, krmax
    COMPLEX(real_8)                          :: a(krmax-krmin+1,*)
    INTEGER                                  :: kr, m
    COMPLEX(real_8)                          :: b(kr,m)

    CHARACTER(*), PARAMETER                  :: procedureN = 'putz'

    INTEGER                                  :: isub, n

    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset(procedureN,isub)
    n=krmax-krmin+1
    CALL zeroing(b)!,m*kr)
    CALL matmov(n,m,a,n,b(krmin,1),kr)
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE putz
  ! ==================================================================
  SUBROUTINE getz(a,b,krmin,krmax,kr,m)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: krmin, krmax
    COMPLEX(real_8)                          :: b(krmax-krmin+1,*)
    INTEGER                                  :: kr
    COMPLEX(real_8)                          :: a(kr,*)
    INTEGER                                  :: m

    CHARACTER(*), PARAMETER                  :: procedureN = 'getz'

    INTEGER                                  :: isub, n

! ==--------------------------------------------------------------==

    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset(procedureN,isub)
    n=krmax-krmin+1
    CALL matmov(n,m,a(krmin,1),kr,b,n)
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE getz
  ! ==================================================================
  ! ==================================================================
  ! CODE FOR NEW FFT ROUTINES
  ! ==================================================================


  SUBROUTINE unpack_y2x(xf,yf,m,nrays,lda,jrxpl,sp5,maxfft,mproc,&
       tr4a2a)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: xf(*), yf(*)
    INTEGER                                  :: m, nrays, lda, maxfft, mproc, &
                                                sp5(0:mproc-1), &
                                                jrxpl(0:mproc-1)
    LOGICAL                                  :: tr4a2a

    COMPLEX(real_4), POINTER                 :: yf4(:)
    INTEGER                                  :: ip, ipp, isub1, k, nrs, nrx

    !$    INTEGER   max_threads
    !$    INTEGER, EXTERNAL :: omp_get_max_threads
    CHARACTER(*),PARAMETER :: procedureN='UNPACK_Y2X'
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    ! ..Pack the data for sending
    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset(procedureN,isub1)
    !$    IF(HASNT_THREADED_COPIES) THEN
    !$       max_threads = omp_get_max_threads()
    !$       call omp_set_num_threads(1)
    !$    ENDIF
    IF (tr4a2a) THEN
       CALL type_cast(yf, maxfft, yf4)
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays + 1
             ipp = ip*lda + 1
             CALL dcopy_s(2*nrx,yf4(ipp),1,xf(nrs),1)
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP5,NRAYS,JRXPL,LDA) &
          !$omp             private(NRX,NRS,IPP,K)
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays
             ipp = ip*lda
             DO k=1,nrx
                xf(nrs+k)=yf4(ipp+k)
             ENDDO
          ENDDO
       END SELECT
    ELSE
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays + 1
             ipp = ip*lda + 1
             CALL dcopy(2*nrx,yf(ipp),1,xf(nrs),1)
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP5,NRAYS,JRXPL,LDA) &
          !$omp             private(NRX,NRS,IPP,K)
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays
             ipp = ip*lda
             DO k=1,nrx
                xf(nrs+k)=yf(ipp+k)
             ENDDO
          ENDDO
       END SELECT
    ENDIF
    !$    IF(HASNT_THREADED_COPIES) call omp_set_num_threads(max_threads)
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt(procedureN,isub1)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE unpack_y2x
  ! ==================================================================
  SUBROUTINE pack_y2x(xf,yf,m,lr1,lda,msp,lmsp,sp8,maxfft,mproc,&
       tr4a2a)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: xf(*), yf(*)
    INTEGER                                  :: m, lr1, lda, lmsp, &
                                                msp(lmsp,*), maxfft, mproc, &
                                                sp8(0:mproc-1)
    LOGICAL                                  :: tr4a2a

    COMPLEX(real_4), POINTER                 :: xf4(:)
    INTEGER                                  :: i, ii, ip, isub1, jj, k, &
                                                mxrp, nrx

    !$    INTEGER   max_threads
    !$    INTEGER, EXTERNAL :: omp_get_max_threads
    CHARACTER(*),PARAMETER :: procedureN='PACK_Y2X'
    ! ==--------------------------------------------------------------==
    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset(procedureN,isub1)
    ! ..Pack the data for sending
    ! IF(HAS_LOW_LEVEL_TIMERS) CALL TISET(procedureN,ISUB1)
    !$    IF(HASNT_THREADED_COPIES) THEN
    !$       max_threads = omp_get_max_threads()
    !$       call omp_set_num_threads(1)
    !$    ENDIF
    IF (tr4a2a) THEN
       CALL type_cast(xf, maxfft, xf4)
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp + 1
                jj = (i-1)*m + 1
                CALL cgthr_z(mxrp,yf(jj),xf4(ii),msp(1,ip+1))
             ENDDO
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP8,LR1,M,LDA) &
          !$omp             private(MXRP,I,II,JJ,K)
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp
                jj = (i-1)*m
                DO k=1,mxrp
                   xf4(ii+k) = yf(jj+msp(k,ip+1))
                ENDDO
             ENDDO
          ENDDO
       END SELECT
    ELSE
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default

          !$omp parallel do shared(SP8,NRX,M,LDA) &
          !$omp             private(MXRP,I,II,JJ)
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp + 1
                jj = (i-1)*m + 1
                CALL zgthr_no_omp(mxrp,yf(jj),xf(ii),msp(1,ip+1))
             ENDDO
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP8,LR1,M,LDA) &
          !$omp             private(MXRP,I,II,JJ,K)
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp
                jj = (i-1)*m
                DO k=1,mxrp
                   xf(ii+k) = yf(jj+msp(k,ip+1))
                ENDDO
             ENDDO
          ENDDO
       END SELECT
    ENDIF
    !$    IF(HASNT_THREADED_COPIES) call omp_set_num_threads(max_threads)
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt(procedureN,isub1)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE pack_y2x
  ! ==================================================================
  SUBROUTINE fft_comm(xf,yf,lda,tr4a2a, comm )
    ! ==--------------------------------------------------------------==
#ifdef __PARALLEL
    USE mpi_f08
#endif
    COMPLEX(real_8), TARGET                  :: xf(*), yf(*)
    INTEGER, INTENT(IN)                      :: lda
    LOGICAL, INTENT(IN)                      :: tr4a2a
#ifdef __PARALLEL
    type(MPI_COMM), INTENT(IN)                      :: comm
#else
    INTEGER, INTENT(IN)                      :: comm
#endif

    CHARACTER(*), PARAMETER                  :: procedureN = 'fft_comm'

    COMPLEX(real_4), POINTER                 :: xf4(:), yf4(:)
    INTEGER                                  :: isub1, isub2

! Variables
! ==--------------------------------------------------------------==
! ..All to all communication

    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tiset(procedureN//'_tuning',isub1)
       ELSE
          CALL tiset(procedureN,isub2)
       END IF
    END IF
    IF (tr4a2a) THEN
       ! is that needed on P?
       CALL type_cast(xf, maxfft, xf4)
       CALL type_cast(yf, maxfft, yf4)       
       CALL mp_all2all( xf4, yf4, lda, comm )
    ELSE
       CALL mp_all2all( xf, yf, lda, comm )
    ENDIF
    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tihalt(procedureN//'_tuning',isub1)
       ELSE
          CALL tihalt(procedureN,isub2)
       END IF
    END IF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE fft_comm
  ! ==================================================================
  SUBROUTINE pack_x2y(xf,yf,nrays,lda,jrxpl,sp5,&
       maxfft,mproc,tr4a2a)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: xf(*), yf(*)
    INTEGER                                  :: nrays, lda, maxfft, mproc, &
                                                sp5(0:mproc-1), &
                                                jrxpl(0:mproc-1)
    LOGICAL                                  :: tr4a2a

    COMPLEX(real_4), POINTER                 :: yf4(:)
    INTEGER                                  :: ip, ipp, isub1, k, nrs, nrx

    !$    INTEGER   max_threads
    !$    INTEGER, EXTERNAL :: omp_get_max_threads
    CHARACTER(*),PARAMETER :: procedureN='PACK_X2Y'
    ! ==--------------------------------------------------------------==
    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset(procedureN,isub1)
    ! ..Prepare data for sending
    !$    IF(HASNT_THREADED_COPIES) THEN
    !$       max_threads = omp_get_max_threads()
    !$       call omp_set_num_threads(1)
    !$    ENDIF
    IF (tr4a2a) THEN
       CALL type_cast(yf, maxfft, yf4)
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays + 1
             ipp = ip*lda + 1
             CALL scopy_d(2*nrx,xf(nrs),1,yf4(ipp),1)
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP5,NRAYS,JRXPL,LDA) &
          !$omp             private(NRX,NRS,IPP,K)
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays
             ipp = ip*lda
             DO k=1,nrx
                yf4(ipp+k)=xf(nrs+k)
             ENDDO
          ENDDO
       END SELECT
    ELSE
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays + 1
             ipp = ip*lda + 1
             CALL dcopy(2*nrx,xf(nrs),1,yf(ipp),1)
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP5,NRAYS,JRXPL,LDA) &
          !$omp             private(NRX,NRS,IPP,K)
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays
             ipp = ip*lda
             DO k=1,nrx
                yf(ipp+k)=xf(nrs+k)
             ENDDO
          ENDDO
       END SELECT
    ENDIF
    !$    IF(HASNT_THREADED_COPIES) call omp_set_num_threads(max_threads)
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt(procedureN,isub1)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE pack_x2y
  ! ==================================================================
  SUBROUTINE unpack_x2y(xf,yf,m,lr1,lda,msp,lmsp,sp8,maxfft,mproc,tr4a2a)
    ! ==--------------------------------------------------------------==
    ! include 'parac.inc'
    COMPLEX(real_8)                          :: xf(*)
    INTEGER                                  :: m, lr1, lda, lmsp, &
                                                msp(lmsp,*), maxfft
    COMPLEX(real_8)                          :: yf(maxfft)
    INTEGER                                  :: mproc, sp8(0:mproc-1)
    LOGICAL                                  :: tr4a2a

    COMPLEX(real_4), POINTER                 :: xf4(:)
    INTEGER                                  :: i, ii, ip, isub1, jj, k, mxrp

    !$    INTEGER   max_threads
    !$    INTEGER, EXTERNAL :: omp_get_max_threads
    CHARACTER(*),PARAMETER :: procedureN='UNPACK_X2Y'
    ! ==--------------------------------------------------------------==
    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset(procedureN,isub1)
    ! ..Unpacking the data
    CALL zeroing(yf)!,maxfft)
    !$    IF(HASNT_THREADED_COPIES) THEN
    !$       max_threads = omp_get_max_threads()
    !$       call omp_set_num_threads(1)
    !$    ENDIF
    IF (tr4a2a) THEN
       CALL type_cast(xf, maxfft, xf4)
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp + 1
                jj = (i-1)*m + 1
                CALL zsctr_c(mxrp,xf4(ii),msp(1,ip+1),yf(jj))
             ENDDO
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP8,LR1,M,LDA) &
          !$omp             private(MXRP,I,II,JJ,K)
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp
                jj = (i-1)*m
                DO k=1,mxrp
                   yf(jj+msp(k,ip+1)) = xf4(ii+k)
                ENDDO
             ENDDO
          ENDDO
       END SELECT
    ELSE
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default
          !$omp parallel do shared(SP8,LR1,M,LDA) &
          !$omp             private(MXRP,I,II,JJ)
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp + 1
                jj = (i-1)*m + 1
                CALL zsctr_no_omp(mxrp,xf(ii),msp(1,ip+1),yf(jj))
             ENDDO
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP8,LR1,M,LDA) &
          !$omp             private(MXRP,I,II,JJ,K)
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp
                jj = (i-1)*m
                DO k=1,mxrp
                   yf(jj+msp(k,ip+1)) = xf(ii+k)
                ENDDO
             ENDDO
          ENDDO
       END SELECT
    ENDIF
    !$    IF(HASNT_THREADED_COPIES) call omp_set_num_threads(max_threads)
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt(procedureN,isub1)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE unpack_x2y
  ! ==================================================================
  SUBROUTINE phasen(f,kr1,kr2s,kr3s,n1u,n1o,nr2s,nr3s)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: kr1, kr2s, kr3s
    COMPLEX(real_8)                          :: f(kr1,kr2s,kr3s)
    INTEGER                                  :: n1u, n1o, nr2s, nr3s

    INTEGER                                  :: i, ii, ijk, isub, j, k
    REAL(real_8), DIMENSION(2)               :: pf = (/1._real_8,-1._real_8/)

    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset('     PHASE',isub)
    !$omp parallel do default(none) __COLLAPSE2 &
    !$omp             private(K,J,I,II,IJK) &
    !$omp             shared(F,PF,NR3S,NR2S,N1U,N1O)
    DO k=1,nr3s
       DO j=1,nr2s
          DO i=n1u,n1o
             ii=i-n1u+1
             ijk=MOD(k+j+i+1,2)+1
             f(ii,j,k)=f(ii,j,k)*pf(ijk)
          ENDDO
       ENDDO
    ENDDO
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt('     PHASE',isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE phasen
  ! ==================================================================

  !TK
  ! Rewritten routines taking care of batches of states
  ! no support for A2A in single precision
  ! Author:
  ! Tobias Kloeffel, CCC,FAU Erlangen-Nuernberg tobias.kloeffel@fau.de
  ! Gerald Mathias, LRZ, Garching Gerald.Mathias@lrz.de
  ! Bernd Meyer, CCC, FAU Erlangen-Nuernberg bernd.meyer@fau.de

  SUBROUTINE putz_n(a,b,krmin,krmax,kr,kr1,kr2s,nperbatch)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: krmin, krmax, kr, kr1, kr2s, nperbatch
    COMPLEX(real_8),INTENT(IN)               :: a(krmax-krmin+1,kr1,nperbatch,*)
    COMPLEX(real_8),INTENT(OUT)              :: b(kr,kr1,kr2s,*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'putz_n'
    REAL(real_8), POINTER __CONTIGUOUS       :: b_r(:,:,:,:)
    REAL(real_8), POINTER __CONTIGUOUS       :: a_r(:,:,:,:)
    INTEGER                                  :: isub1,isub2

    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tiset(procedureN//'_tuning',isub1)
       ELSE
          CALL tiset(procedureN,isub2)
       END IF
    END IF
    
    CALL reshape_inplace(a,(/(krmax-krmin+1)*2,kr1,nperbatch,kr2s/),a_r)
    CALL reshape_inplace(b,(/kr*2,kr1,kr2s,nperbatch/),b_r)
    CALL putz_n_r(a_r,b_r,krmin,krmax,kr,kr1,kr2s,nperbatch)

    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tihalt(procedureN//'_tuning',isub1)
       ELSE
          CALL tihalt(procedureN,isub2)
       END IF
    END IF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE putz_n

  SUBROUTINE getz_n(a,b,krmin,krmax,kr,kr1,kr2s,nperbatch)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: krmin, krmax, kr, kr1, kr2s, nperbatch
    COMPLEX(real_8),INTENT(OUT)              :: b(krmax-krmin+1,kr1,nperbatch,*)
    COMPLEX(real_8),INTENT(IN)               :: a(kr,kr1,kr2s,*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'getz_n'

    INTEGER                                  :: isub1,isub2,n,is,i,j,k,n1

! ==--------------------------------------------------------------==
    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tiset(procedureN//'_tuning',isub1)
       ELSE
          CALL tiset(procedureN,isub2)
       END IF
    END IF

    n=krmax-krmin+1
    n1=krmin-1
    !$omp parallel private (i,is,k,j) proc_bind(close)
    DO is=1,nperbatch
       !$omp do
       DO i=1,kr2s
          DO k=1,kr1
             !$omp simd
             DO j=1,n
                b(j,k,is,i)=a(n1+j,k,i,is)
             END DO
          END DO
       END DO
       !$omp end do nowait
    END DO
    !$omp end parallel

    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tihalt(procedureN//'_tuning',isub1)
       ELSE
          CALL tihalt(procedureN,isub2)
       END IF
    END IF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE getz_n

  SUBROUTINE pack_y2x_n(xf,yf,m,lr1,lda,msp,lmsp,sp8,maxfft,mproc,&
       tr4a2a,nperbatch,offset_in)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8),INTENT(OUT)              :: xf(*)
    COMPLEX(real_8),INTENT(IN)               :: yf(*)
    INTEGER,INTENT(IN)                       :: m, lr1, lda, lmsp, &
                                                msp(lmsp,*), mproc, &
                                                sp8(0:mproc-1), maxfft
    LOGICAL,INTENT(IN)                       :: tr4a2a
    INTEGER,INTENT(IN),OPTIONAL              :: nperbatch,offset_in

    INTEGER                                  :: i, ii, ip, isub1, isub2, &
                                                jj, k, mxrp, is, nstate, offset

    CHARACTER(*),PARAMETER :: procedureN='PACK_Y2X_n'
    ! ==--------------------------------------------------------------==
    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tiset(procedureN//'_tuning',isub1)
       ELSE
          CALL tiset(procedureN,isub2)
       END IF
    END IF

    IF(PRESENT(nperbatch))THEN
       nstate=nperbatch
    ELSE
       nstate=1
    END IF
    IF(PRESENT(offset_in))THEN
       offset=offset_in
    ELSE
       offset=0
    END IF

    !$omp parallel private(is,ip,mxrp,i,ii,jj,k) proc_bind(close)
    DO is=1,nstate
       !$omp do
       DO ip=0,mproc-1
          mxrp = sp8(ip)
          DO i=1,lr1
             ii = ip*lda*nstate + (i-1)*mxrp + (is-1)*lda
             jj = (i-1)*m + (is-1)*offset
             DO k=1,mxrp
                xf(ii+k) = yf(jj+msp(k,ip+1))
             ENDDO
          ENDDO
       ENDDO
       !$omp end do nowait
    end do
    !$omp end parallel

    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tihalt(procedureN//'_tuning',isub1)
       ELSE
          CALL tihalt(procedureN,isub2)
       END IF
    END IF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE pack_y2x_n

  SUBROUTINE unpack_y2x_n(xf,yf,m,nrays,lda,jrxpl,sp5,maxfft,mproc,tr4a2a,count)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8),INTENT(IN)               :: yf(*)
    COMPLEX(real_8),INTENT(OUT)              :: xf(*)
    LOGICAL,INTENT(IN)                       :: tr4a2a
    INTEGER,INTENT(IN)                       :: m, nrays, lda, mproc, &
                                                sp5(0:mproc-1), &
                                                jrxpl(0:mproc-1), maxfft
    INTEGER,INTENT(IN),OPTIONAL              :: count

    INTEGER                                  :: ip, ipp, isub1, isub2, k, nrs, &
                                                nrx, i,is,nstate

    CHARACTER(*),PARAMETER :: procedureN='UNPACK_Y2X_n'
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    ! ..Pack the data for sending
    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tiset(procedureN//'_tuning',isub1)
       ELSE
          CALL tiset(procedureN,isub2)
       END IF
    END IF

    If(PRESENT(count))THEN
       nstate=count
    ELSE
       nstate=1
    END IF

    !$omp parallel do private(ip,i,is,nrs,ipp,k) proc_bind(close)
    DO ip=0,mproc-1
       DO i=1,sp5(ip)
          DO is=1,nstate
             nrs = (jrxpl(ip)-1)*nrays*nstate + (i-1)*nrays*nstate +(is-1)*nrays
             ipp = ip*lda*nstate + (i-1)*nrays + (is-1)*lda
             !$omp simd
             DO k=1,nrays
                xf(nrs+k)=yf(ipp+k)
             END DO
          END DO
       END DO
    END DO

    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tihalt(procedureN//'_tuning',isub1)
       ELSE
          CALL tihalt(procedureN,isub2)
       END IF
    END IF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE unpack_y2x_n
  ! ==================================================================

  SUBROUTINE pack_x2y_n(xf,yf,nrays,lda,jrxpl,sp5,maxfft,mproc,tr4a2a,count)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8),INTENT(IN)               :: xf(*)
    COMPLEX(real_8),INTENT(OUT)              :: yf(*)
    INTEGER,INTENT(IN)                       :: nrays, lda, maxfft, mproc, &
                                                sp5(0:mproc-1), &
                                                jrxpl(0:mproc-1)
    LOGICAL,INTENT(IN)                       :: tr4a2a
    INTEGER,INTENT(IN),OPTIONAL              :: count

    INTEGER                                  :: ip, ipp, isub1, isub2, &
                                                k, nrs, nrx, i,is,nstate
    CHARACTER(*),PARAMETER :: procedureN='PACK_X2Y_n'
    ! ==--------------------------------------------------------------==
    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tiset(procedureN//'_tuning',isub1)
       ELSE
          CALL tiset(procedureN,isub2)
       END IF
    END IF

    IF(PRESENT(count))THEN
       nstate=count
    ELSE
       nstate=1
    END IF

    !$omp parallel do private(ip,i,is,nrs,ipp,k) proc_bind(close)
    DO ip=0,mproc-1
       DO i=1,sp5(ip)
          DO is=1,nstate
             nrs = (jrxpl(ip)-1)*nrays*nstate + (i-1)*nrays*nstate +(is-1)*nrays
             ipp = ip*lda*nstate + (i-1)*nrays + (is-1)*lda
             !$omp simd
             DO k=1,nrays
                yf(ipp+k)=xf(nrs+k)
             END DO
          END DO
       END DO
    END DO

    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tihalt(procedureN//'_tuning',isub1)
       ELSE
          CALL tihalt(procedureN,isub2)
       END IF
    END IF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE pack_x2y_n

  ! ==================================================================
  SUBROUTINE unpack_x2y_n(xf,yf,m,lr1,lda,msp,lmsp,sp8,maxfft,mproc,tr4a2a,&
       count,offset_in)
    ! ==--------------------------------------------------------------==
    ! include 'parac.inc'
    COMPLEX(real_8),INTENT(IN)               :: xf(*)
    COMPLEX(real_8),INTENT(OUT)              :: yf(*)
    INTEGER, INTENT(IN)                      :: m, lr1, lda, lmsp, &
                                                msp(lmsp,*), maxfft,&
                                                mproc, sp8(0:mproc-1)
    LOGICAL, INTENT(IN)                      :: tr4a2a
    INTEGER, INTENT(IN), OPTIONAL            :: count,offset_in
    INTEGER                                  :: isub1, isub2, nstate, offset
    REAL(real_8), POINTER __CONTIGUOUS       :: yf_r(:),xf_r(:)
    
    !$    INTEGER   max_threads
    !$    INTEGER, EXTERNAL :: omp_get_max_threads
    CHARACTER(*),PARAMETER :: procedureN='UNPACK_X2Y_n'
    ! ==--------------------------------------------------------------==
    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tiset(procedureN//'_tuning',isub1)
       ELSE
          CALL tiset(procedureN,isub2)
       END IF
    END IF
    ! ..Unpacking the data
    IF(PRESENT(count))THEN
       nstate=count
    ELSE
       nstate=1
    ENDIF
    IF(PRESENT(offset_in))THEN
       offset=offset_in
    ELSE
       offset=0
    END IF
    CALL reshape_inplace(yf,(/nstate*offset*2/),yf_r)
    CALL reshape_inplace(xf,(/nstate*offset*2/),xf_r)
    CALL unpack_x2y_n_r(xf_r,yf_r,m,lr1,lda,msp,lmsp,sp8,maxfft,mproc,&
       offset,nstate)
    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tihalt(procedureN//'_tuning',isub1)
       ELSE
          CALL tihalt(procedureN,isub2)
       END IF
    END IF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE unpack_x2y_n
  ! ==================================================================
  SUBROUTINE unpack_x2y_n_r(xf_r,yf_r,m,lr1,lda,msp,lmsp,sp8,maxfft,mproc,&
       offset,nstate)
    ! ==--------------------------------------------------------------==
    ! include 'parac.inc'
    real(real_8),INTENT(IN)               :: xf_r(*)
    real(real_8),INTENT(OUT)              :: yf_r(*)
    INTEGER, INTENT(IN)                      :: m, lr1, lda, lmsp, &
                                                msp(lmsp,*), maxfft,&
                                                mproc, sp8(0:mproc-1)
    INTEGER, INTENT(IN)                    :: offset,nstate
    INTEGER                                  :: i, ii, ip, isub1, jj, k, mxrp,is,kk

    !$omp parallel private(is,ip,mxrp,i,ii,jj,k,kk) proc_bind(close)
    !    call zero(yf_r,nstate*offset*2)
    !$omp do
    do is=1,nstate*offset*2
       yf_r(is)=0._real_8
    end do
    DO is=1,nstate
       !$omp do schedule (static)
       DO ip=0,mproc-1
          mxrp = sp8(ip)
          DO i=1,lr1
             ii = ip*lda*nstate + (i-1)*mxrp + (is-1)*lda
             jj = (i-1)*m + (is-1)*offset
             DO k=1,mxrp
                do kk=-1,0
                   yf_r((jj+msp(k,ip+1))*2+kk) = xf_r((ii+k)*2+kk)
                end do
             END DO
          END DO
       END DO
       !$omp end do nowait
    END DO
    !$omp end parallel

  END SUBROUTINE unpack_x2y_n_r
  SUBROUTINE putz_n_r(a_r,b_r,krmin,krmax,kr,kr1,kr2s,nperbatch)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: krmin, krmax, kr, kr1, kr2s, nperbatch
    REAL(real_8),INTENT(IN)               :: a_r((krmax-krmin+1)*2,kr1,nperbatch,*)
    REAL(real_8),INTENT(OUT)              :: b_r(kr*2,kr1,kr2s,*)

    INTEGER                                  :: isub,n,n1,n2,n3,n4,is,i,j,k,krmin_loc,kr_loc

    
    n=krmax-krmin+1
    n1=(krmin-1)*2
    n2=(krmin-1+n)*2
    n3=(krmin+n)*2-1
    n4=(krmin-1)*2
    kr_loc=kr*2
    krmin_loc=krmin*2-1

    !$omp parallel private (i,is,k,j) proc_bind(close)
    DO is=1,nperbatch
       !$omp do
       DO i=1,kr2s
          DO k=1,kr1
             !$omp simd
             DO j=1,n1
                b_r(j,k,i,is)=0.0_real_8
             END DO
             !$omp simd
             DO j=krmin_loc,n2
                b_r(j,k,i,is)=a_r(j-n4,k,is,i)
             END DO
             !$omp simd
             DO j=n3,kr_loc
                b_r(j,k,i,is)=0.0_real_8
             END DO
          END DO
       END DO
       !$omp end do nowait
    END DO
    !$omp end parallel

    ! ==--------------------------------------------------------------==
  END SUBROUTINE putz_n_r

  subroutine zero(in,len)
    real(real_8), intent(out) __CONTIGUOUS :: in(:)
    integer, intent(in) :: len
    integer :: i
    !$omp do simd
    do i=1,len
       in(i)=0.0_real_8
    end do
  end subroutine zero
  subroutine zero_noomp(in,len)
    real(real_8), intent(out) __CONTIGUOUS :: in(:)
    integer, intent(in) :: len
    integer :: i
    !$omp simd
    do i=1,len
       in(i)=0.0_real_8
    end do
  end subroutine zero_noomp
  !TK
  !CR
  SUBROUTINE Prepare_Psi( tfft, psi, aux, remswitch, mythread )
    IMPLICIT NONE
  
    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) ::tfft 
    INTEGER, INTENT(IN) :: remswitch, mythread
    COMPLEX(DP), INTENT(IN)  :: psi ( : , : )
    COMPLEX(DP), INTENT(OUT)  :: aux ( tfft%nr3 , * ) !ns(parai%me+1)*z_group_size )
  
    INTEGER :: j, i, iter
    INTEGER :: offset, offset2
  !  INTEGER :: isub, isub4
  !  CHARACTER(*), PARAMETER :: procedureN = 'Prepare_Psi'
  
    INTEGER(INT64) :: time(2), cr
  
  !  IF( cntl%fft_tune_batchsize ) THEN
  !     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN//'_tuning',isub4)
  !  ELSE
  !     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN,isub)
  !  END IF
  
    IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(1) )
  !------------------------------------------------------
  !----------Prepare_Psi Start---------------------------
  
    ! parai%me+1 in tfft%thread_z_start eigentlich (parai%my_node-1)*parai%node_nproc+parai%node_me+1
  
    DO i = tfft%thread_z_start( mythread+1, remswitch, parai%me+1, tfft%which ), tfft%thread_z_end( mythread+1, remswitch, parai%me+1, tfft%which )
       iter = mod( i-1, tfft%nsw(parai%me+1) ) + 1
       offset  = ( iter - 1 ) * tfft%nr3
       offset2 = 2 * ( ( (i-1) / tfft%nsw(parai%me+1) ) + 1 )
  
       DO j = 1, tfft%prep_map(1,iter)-1
          aux( j, i ) = conjg( psi( indz_r( offset + j ), offset2 - 1 ) - (0.0d0,1.0d0) * psi( indz_r( offset + j ), offset2 ) )
       ENDDO
       DO j = 1, tfft%prep_map(2,iter)-1
          aux( j, i ) = psi( nzh_r( offset + j ), offset2 - 1 ) + (0.0d0,1.0d0) * psi( nzh_r( offset + j ), offset2 )
       ENDDO
  
       DO j = tfft%prep_map(3,iter), tfft%prep_map(4,iter)
          aux( j, i ) = (0.d0, 0.d0)
       ENDDO
  
       DO j = tfft%prep_map(5,iter)+1, tfft%nr3
          aux( j, i ) = psi( nzh_r( offset + j ), offset2 - 1 ) + (0.0d0,1.0d0) * psi( nzh_r( offset + j ), offset2 )
       ENDDO
       DO j = tfft%prep_map(6,iter)+1, tfft%nr3
          aux( j, i ) = conjg( psi( indz_r( offset + j ), offset2 - 1 ) - (0.0d0,1.0d0) * psi( indz_r( offset + j ), offset2 ) )
       ENDDO
  
    ENDDO
  
  !----------Prepare_Psi End-----------------------------
  !------------------------------------------------------
    IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(2) )
    IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) tfft%time_adding( 1 ) = tfft%time_adding( 1 ) + ( time(2) - time(1) )
  
  !  CALL SYSTEM_CLOCK( count_rate = cr )
  !  write(6,*) "ACTUAL", REAL( tfft%time_adding( 1 ) / REAL( cr ) )
  
  !  IF( cntl%fft_tune_batchsize ) THEN
  !     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN//'_tuning',isub4)
  !  ELSE
  !     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN,isub)
  !  END IF
  
  END SUBROUTINE Prepare_Psi
  
  SUBROUTINE fft_com( tfft, remswitch, work_buffer, which )
    IMPLICIT NONE
  
    INTEGER, INTENT(IN)                         :: remswitch, work_buffer, which
    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT)    :: tfft 
  
    CHARACTER(*), PARAMETER :: procedureN = 'fft_com'
  
    INTEGER :: ierr, isub, isub4
  
    IF( cntl%fft_tune_batchsize ) THEN
       CALL tiset(procedureN//'_tuning',isub4)
    ELSE
       CALL tiset(procedureN,isub)
    END IF
  
    !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, tfft%mpi_window( 1 ), ierr )
    !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, tfft%mpi_window( 2 ), ierr )
    
    CALL MP_STARTALL( tfft%comm_sendrecv(1,tfft%which), parai%send_handle(:,work_buffer,remswitch,which) )
    CALL MP_STARTALL( tfft%comm_sendrecv(2,tfft%which), parai%recv_handle(:,work_buffer,remswitch,which) )
    
    CALL MP_WAITALL( tfft%comm_sendrecv(1,tfft%which), parai%send_handle(:,work_buffer,remswitch,which) )
    CALL MP_WAITALL( tfft%comm_sendrecv(2,tfft%which), parai%recv_handle(:,work_buffer,remswitch,which) )
  
    !CALL mpi_win_unlock_all( tfft%mpi_window( 2 ), ierr )
    !CALL mpi_win_unlock_all( tfft%mpi_window( 1 ), ierr )
  
    IF( cntl%fft_tune_batchsize ) THEN
       CALL tihalt(procedureN//'_tuning',isub4)
    ELSE
       CALL tihalt(procedureN,isub)
    END IF
  
  END SUBROUTINE fft_com
  
  SUBROUTINE invfft_z_section( tfft, aux, comm_mem_send, comm_mem_recv, batch_size, remswitch, mythread, nss )
    IMPLICIT NONE
  
    INTEGER, INTENT(IN) :: batch_size, remswitch, mythread
    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) :: tfft 
    COMPLEX(DP), INTENT(INOUT) :: comm_mem_send( * ), comm_mem_recv( * )
    COMPLEX(DP), INTENT(INOUT)  :: aux ( tfft%nr3 , * ) !ns(parai%me+1)*z_group_size )
    INTEGER, INTENT(IN) :: nss(*)
  
    INTEGER :: l, m, j, k, i
    INTEGER :: offset, kdest, ierr
  !  INTEGER :: isub, isub4
  !  CHARACTER(*), PARAMETER :: procedureN = 'invfft_z_section'
  
    INTEGER(INT64) :: time(3)
  
  !  IF( cntl%fft_tune_batchsize ) THEN
  !     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN//'_tuning',isub4)
  !  ELSE
  !     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN,isub)
  !  END IF
  
    IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(1) )
  !------------------------------------------------------
  !------------z-FFT Start-------------------------------
  
    CALL mltfft_fftw_t('n','n',aux( : , tfft%thread_z_start( mythread+1, remswitch, parai%me+1, tfft%which ) : tfft%thread_z_end( mythread+1, remswitch, parai%me+1, tfft%which ) ), &
                     tfft%nr3, tfft%thread_z_sticks(mythread+1,remswitch,parai%me+1,tfft%which), &
                     aux( : , tfft%thread_z_start( mythread+1, remswitch, parai%me+1, tfft%which ) : tfft%thread_z_end( mythread+1, remswitch, parai%me+1, tfft%which ) ), &
                     tfft%nr3, tfft%thread_z_sticks(mythread+1,remswitch,parai%me+1,tfft%which), &
                     tfft%nr3, tfft%thread_z_sticks(mythread+1,remswitch,parai%me+1,tfft%which),-1,scal,.FALSE.,mythread,parai%ncpus)
  
  !-------------z-FFT End--------------------------------
  !------------------------------------------------------
    IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(2) )
  !------------------------------------------------------
  !---------Pre-Com-Copy Start---------------------------
  
    IF( parai%nnode .ne. 1 ) THEN
  
       !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, tfft%mpi_window( 1 ), ierr )
  
       DO l = 1, parai%nnode
          IF( l .eq. parai%my_node+1 ) CYCLE 
          DO m = 1, parai%node_nproc
             j = (l-1)*parai%node_nproc + m
             offset = ( parai%node_me + (j-1)*parai%node_nproc ) * tfft%small_chunks(tfft%which) + (l-1)*(batch_size-1) * tfft%big_chunks(tfft%which)
             DO k = tfft%thread_z_start( mythread+1, remswitch, parai%me+1, tfft%which ), tfft%thread_z_end( mythread+1, remswitch, parai%me+1, tfft%which )
                kdest = offset + tfft%nr3px * mod( (k-1), nss(parai%me+1) ) + ( (k-1) / nss(parai%me+1) ) * tfft%big_chunks(tfft%which)
                DO i = 1, tfft%nr3p( j )
                   comm_mem_send( kdest + i ) = aux( i + tfft%nr3p_offset( j ), k )
                ENDDO
             ENDDO
          ENDDO
       ENDDO
  
       !CALL mpi_win_unlock_all( tfft%mpi_window( 1 ), ierr )
  
    END IF
  
    !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, tfft%mpi_window( 2 ), ierr )
  
    DO m = 1, parai%node_nproc
       j = parai%my_node*parai%node_nproc + m
       offset = ( parai%node_me + (j-1)*parai%node_nproc ) * tfft%small_chunks(tfft%which) + parai%my_node*(batch_size-1) * tfft%big_chunks(tfft%which)
       DO k = tfft%thread_z_start( mythread+1, remswitch, parai%me+1, tfft%which ), tfft%thread_z_end( mythread+1, remswitch, parai%me+1, tfft%which )
          kdest = offset + tfft%nr3px * mod( k-1, nss(parai%me+1) ) + ( (k-1) / nss(parai%me+1) ) * tfft%big_chunks(tfft%which)
          DO i = 1, tfft%nr3p( j )
             comm_mem_recv( kdest + i ) = aux( i + tfft%nr3p_offset( j ), k )
          ENDDO
       ENDDO
    ENDDO
  
    !CALL mpi_win_unlock_all( tfft%mpi_window( 2 ), ierr )
  
  !----------Pre-Com-Copy End----------------------------
  !------------------------------------------------------
    IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(3) )
    IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) tfft%time_adding( 2 ) = tfft%time_adding( 2 ) + ( time(2) - time(1) )
    IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) tfft%time_adding( 3 ) = tfft%time_adding( 3 ) + ( time(3) - time(2) )
  
  !  IF( cntl%fft_tune_batchsize ) THEN
  !     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN//'_tuning',isub4)
  !  ELSE
  !     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN,isub)
  !  END IF
  
  END SUBROUTINE invfft_z_section 
  
  SUBROUTINE invfft_y_section( tfft, aux, comm_mem_recv, aux2_r, map_acinv, map_acinv_rem, counter, remswitch, mythread, my_nr1s )
    !$ USE omp_lib
    IMPLICIT NONE
  
    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) ::tfft 
    INTEGER, INTENT(IN) :: remswitch, mythread, counter, my_nr1s
    COMPLEX(DP), INTENT(IN)  :: comm_mem_recv( * )
    COMPLEX(DP), INTENT(INOUT) :: aux ( tfft%my_nr3p * tfft%nr2 * tfft%nr1 , * ) !y_group_size
    COMPLEX(DP), INTENT(INOUT) :: aux2_r( : , : )
    INTEGER, INTENT(IN) :: map_acinv( * ), map_acinv_rem( * )
  
    INTEGER :: i, k, offset, iter
  !  INTEGER :: isub, isub4
  !  CHARACTER(*), PARAMETER :: procedureN = 'invfft_y_section'
  
    INTEGER(INT64) :: time(5)
  
  !  IF( cntl%fft_tune_batchsize ) THEN
  !     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN//'_tuning',isub4)
  !  ELSE
  !     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN,isub)
  !  END IF
  
    IF( remswitch .eq. 1 ) THEN 
       Call First_Part_y_section( aux2_r, map_acinv )
    ELSE
       Call First_Part_y_section( aux2_r, map_acinv_rem )
    END IF
  
    !$  locks_omp( mythread+1, counter, 7 ) = .false.
    !$omp flush( locks_omp )
    !$  DO WHILE( ANY( locks_omp( :, counter, 7 ) ) )
    !$omp flush( locks_omp )
    !$  END DO
  
    Call Second_Part_y_section( aux2_r )
  
  !  IF( cntl%fft_tune_batchsize ) THEN
  !     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN//'_tuning',isub4)
  !  ELSE
  !     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN,isub)
  !  END IF
  
    CONTAINS
  
      SUBROUTINE First_Part_y_section( aux2, map )
        
        Implicit NONE
        COMPLEX(DP), INTENT(INOUT) :: aux2( tfft%nr2 , * ) !nr1s(parai%me+1) * tfft%my_nr3p *  y_group_size )
        INTEGER, INTENT(IN)        :: map( * )
  
  
          IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(1) )
        !------------------------------------------------------
        !--------After-Com-Copy Start--------------------------
        
          !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, tfft%mpi_window( 2 ), ierr )
          
          DO i = tfft%thread_y_start( mythread+1, remswitch, tfft%which ), tfft%thread_y_end( mythread+1, remswitch, tfft%which )
             iter = mod( i-1, my_nr1s ) + 1
             offset = ( mod( i-1, my_nr1s * tfft%my_nr3p ) + ( (i-1) / ( my_nr1s * tfft%my_nr3p ) ) * tfft%my_nr3p * my_nr1s ) * tfft%nr2
             DO k = 1, tfft%zero_acinv_start( iter, tfft%which ) - 1
                aux2( k, i ) = comm_mem_recv( map( offset + k ) )
             END DO
             DO k = tfft%zero_acinv_start( iter, tfft%which ), tfft%zero_acinv_end( iter, tfft%which )
                aux2( k, i ) = (0.0_DP,0.0_DP)
             END DO
             DO k = tfft%zero_acinv_end( iter, tfft%which ) + 1, tfft%nr2
                aux2( k, i ) = comm_mem_recv( map( offset + k ) )
             END DO
          END DO
        
          !CALL mpi_win_unlock_all( tfft%mpi_window( 2 ), ierr )
        
        !---------After-Com-Copy End---------------------------
        !------------------------------------------------------
          IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(2) )
        !------------------------------------------------------
        !------------y-FFT Start-------------------------------
  
          CALL mltfft_fftw_t('n','n',aux2( : , tfft%thread_y_start( mythread+1, remswitch, tfft%which ) : tfft%thread_y_end( mythread+1, remswitch, tfft%which ) ), &
                           tfft%nr2, tfft%thread_y_sticks(mythread+1,remswitch,tfft%which), &
                           aux2( : , tfft%thread_y_start( mythread+1, remswitch, tfft%which ) : tfft%thread_y_end( mythread+1, remswitch, tfft%which ) ), &
                           tfft%nr2, tfft%thread_y_sticks(mythread+1,remswitch,tfft%which), &
                           tfft%nr2, tfft%thread_y_sticks(mythread+1,remswitch,tfft%which),-1,scal,.FALSE.,mythread,parai%ncpus)
        
        !-------------y-FFT End--------------------------------
        !------------------------------------------------------
          IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(3) )
          IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) tfft%time_adding( 4 ) = tfft%time_adding( 4 ) + ( time(2) - time(1) )
          IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) tfft%time_adding( 5 ) = tfft%time_adding( 5 ) + ( time(3) - time(2) )
  
      END SUBROUTINE First_Part_y_section
  
      SUBROUTINE Second_Part_y_section( aux2 )
  
        IMPLICIT NONE
        COMPLEX(DP), INTENT(INOUT) :: aux2  ( my_nr1s * tfft%my_nr3p * tfft%nr2, * ) !y_group_size
        INTEGER :: ibatch
  
          IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(4) )
        !------------------------------------------------------
        !-------------yx-scatter-------------------------------
        
          DO i = tfft%thread_x_start( mythread+1, remswitch, tfft%which ), tfft%thread_x_end( mythread+1, remswitch, tfft%which )
             offset = mod( i-1, tfft%my_nr3p * tfft%nr2 ) * tfft%nr1
             ibatch = ( ( (i-1) / ( tfft%my_nr3p * tfft%nr2 ) ) + 1 )
             DO k = 1, tfft%zero_scatter_start( tfft%which ) - 1
                aux( offset + k, ibatch ) = aux2( tfft%map_scatter_inv( offset + k, tfft%which ), ibatch )
             END DO
             DO k = tfft%zero_scatter_start( tfft%which ), tfft%zero_scatter_end( tfft%which )
                aux( offset + k, ibatch ) = (0.0_DP, 0.0_DP)
             END DO
             DO k = tfft%zero_scatter_end( tfft%which ) + 1, tfft%nr1
                aux( offset + k, ibatch ) = aux2( tfft%map_scatter_inv( offset + k, tfft%which ), ibatch )
             END DO
          END DO
        
        !-------------yx-scatter-------------------------------
        !------------------------------------------------------
          IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(5) )
          IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) tfft%time_adding( 6 ) = tfft%time_adding( 6 ) + ( time(5) - time(4) )
  
      END SUBROUTINE Second_Part_y_section
  
  END SUBROUTINE invfft_y_section
  
  SUBROUTINE invfft_x_section( tfft, aux, remswitch, mythread )
    IMPLICIT NONE
  
    INTEGER, INTENT(IN) :: remswitch, mythread
    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) :: tfft
    COMPLEX(DP), INTENT(INOUT) :: aux( tfft%nr1, * ) !tfft%my_nr3p * tfft%nr2 * x_group_size )
  
  !  INTEGER :: isub, isub4
  !  CHARACTER(*), PARAMETER :: procedureN = 'invfft_x_section'
  
    INTEGER(INT64) :: time(2)
  
  !  IF( cntl%fft_tune_batchsize ) THEN
  !     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN//'_tuning',isub4)
  !  ELSE
  !     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN,isub)
  !  END IF
  
    IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(1) )
  !------------------------------------------------------
  !------------x-FFT Start-------------------------------
  
    CALL mltfft_fftw_t('n','n',aux( : , tfft%thread_x_start( mythread+1, remswitch, tfft%which ) : tfft%thread_x_end( mythread+1, remswitch, tfft%which ) ), &
                     tfft%nr1, tfft%thread_x_sticks(mythread+1,remswitch, tfft%which), &
                     aux( : , tfft%thread_x_start( mythread+1, remswitch, tfft%which ) : tfft%thread_x_end( mythread+1, remswitch, tfft%which ) ), &
                     tfft%nr1, tfft%thread_x_sticks(mythread+1,remswitch, tfft%which), &
                     tfft%nr1, tfft%thread_x_sticks(mythread+1,remswitch, tfft%which),-1,scal,.FALSE.,mythread,parai%ncpus)
  
  !-------------x-FFT End--------------------------------
  !------------------------------------------------------
    IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(2) )
    IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) tfft%time_adding( 7 ) = tfft%time_adding( 7 ) + ( time(2) - time(1) )
  
  !  IF( cntl%fft_tune_batchsize ) THEN
  !     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN//'_tuning',isub4)
  !  ELSE
  !     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN,isub)
  !  END IF
  
  END SUBROUTINE invfft_x_section
  
  SUBROUTINE fwfft_x_section( tfft, aux_r, aux2, counter, remswitch, mythread, my_nr1s )
    IMPLICIT NONE
  
    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) :: tfft
    INTEGER, INTENT(IN) :: counter, remswitch, mythread, my_nr1s
    COMPLEX(DP), INTENT(INOUT)  :: aux_r( : ) !tfft%my_nr3p * tfft%nr2 * tfft%nr1 , * ) !x_group_size )
    COMPLEX(DP), INTENT(INOUT) :: aux2( my_nr1s * tfft%my_nr3p * tfft%nr2 , * ) !x_group_size )
  
  !  INTEGER :: isub, isub4
  !  CHARACTER(*), PARAMETER :: procedureN = 'fwfft_x_section'
  
    INTEGER(INT64) :: time(4)
  
  !  IF( cntl%fft_tune_batchsize ) THEN
  !     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN//'_tuning',isub4)
  !  ELSE
  !     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN,isub)
  !  END IF
  
    Call First_Part_x_section( aux_r )
  
    !$  locks_omp( mythread+1, counter, 8 ) = .false.
    !$omp flush( locks_omp )
    !$  DO WHILE( ANY( locks_omp( :, counter, 8 ) ) )
    !$omp flush( locks_omp )
    !$  END DO
   
    Call Second_Part_x_section( aux_r )
  
  !  IF( cntl%fft_tune_batchsize ) THEN
  !     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN//'_tuning',isub4)
  !  ELSE
  !     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN,isub)
  !  END IF
  
    CONTAINS
  
      SUBROUTINE First_Part_x_section( aux )
  
        IMPLICIT NONE
        COMPLEX(DP), INTENT(INOUT) :: aux( tfft%nr1 , * ) !tfft%my_nr3p * tfft%nr2 * x_group_size )
  
        !------------------------------------------------------
        !------------x-FFT Start-------------------------------
        
          !$  locks_omp( mythread+1, counter, 9 ) = .false.
          !$omp flush( locks_omp )
          !$  DO WHILE( ANY( locks_omp( :, counter, 9 ) ) )
          !$omp flush( locks_omp )
          !$  END DO
  
          IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(1) )
  
          CALL mltfft_fftw_t('n','n',aux( : , tfft%thread_x_start( mythread+1, remswitch, tfft%which ) : tfft%thread_x_end( mythread+1, remswitch, tfft%which ) ), &
                           tfft%nr1, tfft%thread_x_sticks(mythread+1,remswitch, tfft%which), &
                           aux( : , tfft%thread_x_start( mythread+1, remswitch, tfft%which ) : tfft%thread_x_end( mythread+1, remswitch, tfft%which ) ), &
                           tfft%nr1, tfft%thread_x_sticks(mythread+1,remswitch, tfft%which), &
                           tfft%nr1, tfft%thread_x_sticks(mythread+1,remswitch, tfft%which),1,scal,.FALSE.,mythread,parai%ncpus)
        
        !-------------x-FFT End--------------------------------
        !------------------------------------------------------
          IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(2) )
          IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) tfft%time_adding( 9 ) = tfft%time_adding( 9 ) + ( time(2) - time(1) )
  
      END SUBROUTINE First_Part_x_section
  
      SUBROUTINE Second_Part_x_section( aux )
  
        IMPLICIT NONE
        COMPLEX(DP), INTENT(INOUT)  :: aux( tfft%my_nr3p * tfft%nr2 * tfft%nr1 , * ) !x_group_size )
        INTEGER :: i, j, offset, ibatch
  
          IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(3) )
        !------------------------------------------------------
        !------Forward xy-scatter Start------------------------
        
          DO i = tfft%thread_y_start( mythread+1, remswitch, tfft%which ), tfft%thread_y_end( mythread+1, remswitch, tfft%which )
             ibatch = ( ( (i-1) / ( tfft%my_nr3p * my_nr1s ) ) + 1 )
             offset = tfft%nr2 * mod( i-1, my_nr1s * tfft%my_nr3p )
             DO j = 1, tfft%nr2
                aux2( j + offset , ibatch ) = aux( tfft%map_scatter_fw( j + offset, tfft%which ), ibatch )
             ENDDO
          ENDDO
        
        !-------Forward xy-scatter End-------------------------
        !------------------------------------------------------
          IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(4) )
          IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) tfft%time_adding( 10 ) = tfft%time_adding( 10 ) + ( time(4) - time(3) )
  
      END SUBROUTINE Second_Part_x_section
  
  END SUBROUTINE fwfft_x_section
  
  SUBROUTINE fwfft_y_section( tfft, aux, comm_mem_send, comm_mem_recv, map_pcfw, batch_size, counter, remswitch, mythread, my_nr1s, nss )
    IMPLICIT NONE
  
    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) :: tfft
    INTEGER, INTENT(IN) :: counter, batch_size, remswitch, mythread, my_nr1s
    COMPLEX(DP), INTENT(INOUT)  :: comm_mem_send( * )
    COMPLEX(DP), INTENT(INOUT)  :: comm_mem_recv( * )
    COMPLEX(DP), INTENT(INOUT) :: aux( : , : )
    INTEGER, INTENT(IN) :: map_pcfw( * ), nss( * )
  
    INTEGER :: l, m, i, offset, j, k, ibatch, jter, offset2
  !  INTEGER :: isub, isub4
    CHARACTER(*), PARAMETER :: procedureN = 'fwfft_y_section'
  
    INTEGER(INT64) :: time(4)
  
  !  IF( cntl%fft_tune_batchsize ) THEN
  !     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN//'_tuning',isub4)
  !  ELSE
  !     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN,isub)
  !  END IF
  
    Call First_Part_y_section( aux )
  
    !$  locks_omp( mythread+1, counter, 10 ) = .false.
    !$omp flush( locks_omp )
    !$  DO WHILE( ANY( locks_omp( :, counter, 10 ) ) )
    !$omp flush( locks_omp )
    !$  END DO
   
    Call Second_Part_y_section( aux )
  
  !  IF( cntl%fft_tune_batchsize ) THEN
  !     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN//'_tuning',isub4)
  !  ELSE
  !     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN,isub)
  !  END IF
  
    CONTAINS
  
      SUBROUTINE First_Part_y_section( aux2 )
        
        Implicit NONE
        COMPLEX(DP), INTENT(INOUT) :: aux2( tfft%nr2 , * ) !nr1s(parai%me+1) * tfft%my_nr3p *  y_group_size )
  
          IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(1) )
        !------------------------------------------------------
        !------------y-FFT Start-------------------------------
  
          CALL mltfft_fftw_t('n','n',aux2( : , tfft%thread_y_start( mythread+1, remswitch, tfft%which ) : tfft%thread_y_end( mythread+1, remswitch, tfft%which ) ), &
                           tfft%nr2, tfft%thread_y_sticks(mythread+1,remswitch,tfft%which), &
                           aux2( : , tfft%thread_y_start( mythread+1, remswitch, tfft%which ) : tfft%thread_y_end( mythread+1, remswitch, tfft%which ) ), &
                           tfft%nr2, tfft%thread_y_sticks(mythread+1,remswitch,tfft%which), &
                           tfft%nr2, tfft%thread_y_sticks(mythread+1,remswitch,tfft%which),1,scal,.FALSE.,mythread,parai%ncpus)
        
        !-------------y-FFT End--------------------------------
        !------------------------------------------------------
          IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(2) )
          IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) tfft%time_adding( 11 ) = tfft%time_adding( 11 ) + ( time(2) - time(1) )
  
      END SUBROUTINE First_Part_y_section
  
      SUBROUTINE Second_Part_y_section( aux2 )
  
        IMPLICIT NONE
        COMPLEX(DP), INTENT(INOUT) :: aux2  ( my_nr1s * tfft%my_nr3p * tfft%nr2, * ) !y_group_size
  
          IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(3) )
  
        !------------------------------------------------------
        !---------Pre-Com-Copy Start---------------------------
        
          IF( parai%nnode .ne. 1 ) THEN
        
             !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, tfft%mpi_window( 1 ), ierr )
        
             DO l = 1, parai%nnode
                IF( l .eq. parai%my_node+1 ) CYCLE
                DO m = 1, parai%node_nproc
                   i = (l-1)*parai%node_nproc + m
                   offset = ( parai%node_me + (i-1)*parai%node_nproc ) * tfft%small_chunks(tfft%which) + (l-1)*(batch_size-1) * tfft%big_chunks(tfft%which)
                   DO j = tfft%thread_z_start( mythread+1, remswitch, i, tfft%which ), tfft%thread_z_end( mythread+1, remswitch, i, tfft%which )
                      jter = mod( j-1, nss( i ) ) 
                      ibatch = ( (j-1) / nss( i ) )
                      offset2 =  ibatch * tfft%big_chunks(tfft%which)
                      DO k = 1, tfft%my_nr3p
                         comm_mem_send( offset + offset2 + jter*tfft%nr3px + k ) = &
                         aux2( map_pcfw( (i-1)*tfft%small_chunks(tfft%which) + jter*tfft%nr3px + k ), ibatch+1 )
                      END DO
                   END DO
                END DO
             END DO
        
             !CALL mpi_win_unlock_all( tfft%mpi_window( 1 ), ierr )
        
          END IF
        
          !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, tfft%mpi_window( 2 ), ierr )
        
          DO m = 1, parai%node_nproc
             i = parai%my_node*parai%node_nproc + m
             offset = ( parai%node_me + (i-1)*parai%node_nproc ) * tfft%small_chunks(tfft%which) + parai%my_node*(batch_size-1) * tfft%big_chunks(tfft%which)
             DO j = tfft%thread_z_start( mythread+1, remswitch, i, tfft%which ), tfft%thread_z_end( mythread+1, remswitch, i, tfft%which )
                jter = mod( j-1, nss( i ) ) 
                ibatch = ( (j-1) / nss( i ) )
                offset2 =  ibatch * tfft%big_chunks(tfft%which)
                DO k = 1, tfft%my_nr3p
                   comm_mem_recv( offset + offset2 + jter*tfft%nr3px + k ) = &
                   aux2( map_pcfw( (i-1)*tfft%small_chunks(tfft%which) + jter*tfft%nr3px + k ), ibatch+1 )
                END DO
             END DO
          END DO
        
          !CALL mpi_win_unlock_all( tfft%mpi_window( 2 ), ierr )
        
        !----------Pre-Com-Copy End----------------------------
        !------------------------------------------------------
          IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(4) )
          IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) tfft%time_adding( 12 ) = tfft%time_adding( 12 ) + ( time(4) - time(3) )
  
      END SUBROUTINE Second_Part_y_section
        
  END SUBROUTINE fwfft_y_section
  
  SUBROUTINE fwfft_z_section( tfft, comm_mem_recv, aux, counter, batch_size, remswitch, mythread, nss, factor_in )
    IMPLICIT NONE
  
    INTEGER, INTENT(IN) :: counter, batch_size, remswitch, mythread
    TYPE(FFT_TYPE_DESCRIPTOR), INTENT(INOUT) :: tfft
    COMPLEX(DP), INTENT(IN)  :: comm_mem_recv( * )
    COMPLEX(DP), INTENT(INOUT)  :: aux ( tfft%nr3 , * ) !ns(parai%me+1)*z_group_size )
    INTEGER, INTENT(IN) :: nss( * )
    DOUBLE PRECISION, OPTIONAL, INTENT(IN) :: factor_in
  
    DOUBLE PRECISION :: factor
    INTEGER :: j, l, k, i
    INTEGER :: offset, kfrom, ierr
  !  INTEGER :: isub, isub4
  !  CHARACTER(*), PARAMETER :: procedureN = 'fwfft_z_section'
  
    INTEGER(INT64) :: time(4)
  
  !  IF( cntl%fft_tune_batchsize ) THEN
  !     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN//'_tuning',isub4)
  !  ELSE
  !     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tiset(procedureN,isub)
  !  END IF
  
    IF( present( factor_in ) ) THEN
       factor = factor_in
    ELSE
       factor = 1
    END IF
  
    IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(1) )
  !------------------------------------------------------
  !--------After-Com-Copy Start--------------------------
  
    !CALL mpi_win_lock_all( MPI_MODE_NOCHECK, tfft%mpi_window( 2 ), ierr )
   
    DO j = 1, parai%nnode
       DO l = 1, parai%node_nproc
          offset = ( parai%node_me*parai%node_nproc + (l-1) ) * tfft%small_chunks(tfft%which) + (j-1)*batch_size * tfft%big_chunks(tfft%which)
          DO k = tfft%thread_z_start( mythread+1, remswitch, parai%me+1, tfft%which ), tfft%thread_z_end( mythread+1, remswitch, parai%me+1, tfft%which )
             kfrom = offset + tfft%nr3px * mod( k-1, nss(parai%me+1) ) + ( (k-1) / nss(parai%me+1) ) * tfft%big_chunks(tfft%which)
             DO i = 1, tfft%nr3p( (j-1)*parai%node_nproc + l )
                aux( tfft%nr3p_offset( (j-1)*parai%node_nproc + l ) + i, k ) = comm_mem_recv( kfrom + i ) * factor
             ENDDO
          ENDDO
       ENDDO
    ENDDO
  
    !CALL mpi_win_unlock_all( tfft%mpi_window( 2 ), ierr )
  
  !---------After-Com-Copy End---------------------------
  !------------------------------------------------------
    IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(2) )
  !------------------------------------------------------
  !------------z-FFT Start-------------------------------
  
    !$  locks_omp( mythread+1, counter, 11 ) = .false.
    !$omp flush( locks_omp )
    !$  DO WHILE( ANY( locks_omp( :, counter, 11 ) ) )
    !$omp flush( locks_omp )
    !$  END DO
  
    IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(3) )
  
  
    CALL mltfft_fftw_t('n','n',aux( : , tfft%thread_z_start( mythread+1, remswitch, parai%me+1, tfft%which ) : tfft%thread_z_end( mythread+1, remswitch, parai%me+1, tfft%which ) ), &
                     tfft%nr3, tfft%thread_z_sticks(mythread+1,remswitch,parai%me+1,tfft%which), &
                     aux( : , tfft%thread_z_start( mythread+1, remswitch, parai%me+1, tfft%which ) : tfft%thread_z_end( mythread+1, remswitch, parai%me+1, tfft%which ) ), &
                     tfft%nr3, tfft%thread_z_sticks(mythread+1,remswitch,parai%me+1,tfft%which), &
                     tfft%nr3, tfft%thread_z_sticks(mythread+1,remswitch,parai%me+1,tfft%which),1,scal,.FALSE.,mythread,parai%ncpus)
  
  !-------------z-FFT End--------------------------------
  !------------------------------------------------------
    IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) CALL SYSTEM_CLOCK( time(4) )
    IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) tfft%time_adding( 13 ) = tfft%time_adding( 13 ) + ( time(2) - time(1) )
    IF( mythread .eq. 1 .or. parai%ncpus .eq. 1 ) tfft%time_adding( 14 ) = tfft%time_adding( 14 ) + ( time(4) - time(3) )
  
  !  IF( cntl%fft_tune_batchsize ) THEN
  !     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN//'_tuning',isub4)
  !  ELSE
  !     IF( parai%ncpus .eq. 1 .or. mythread .eq. 1 ) CALL tihalt(procedureN,isub)
  !  END IF
  
  END SUBROUTINE fwfft_z_section
!CR
END MODULE fftutil_utils
