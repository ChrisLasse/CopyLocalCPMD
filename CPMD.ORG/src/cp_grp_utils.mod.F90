MODULE cp_grp_utils
  USE distribution_utils,              ONLY: dist_size
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE kpts,                            ONLY: tkpts
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: imagp,&
                                             nlps_com
  USE part_1d,                         ONLY: part_1d_get_blk_bounds
  USE parac,                           ONLY: parai
  USE system,                          ONLY: ncpw,&
                                             nkpt,&
                                             maxsys,&
                                             parap,&
                                             norbpe
  USE sfac,                            ONLY: dfnl,&
                                             fnl,&
                                             dfnla,&
                                             fnla
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: cp_grp_get_sizes
  PUBLIC :: cp_grp_redist
  PUBLIC :: cp_grp_redist_array
  PUBLIC :: cp_grp_split_states
  PUBLIC :: cp_grp_split_atoms
  PUBLIC :: cp_grp_redist_dfnl_fnl

  INTERFACE cp_grp_redist
     MODULE PROCEDURE cp_grp_redist_d
     MODULE PROCEDURE cp_grp_redist_z
  END INTERFACE cp_grp_redist
  INTERFACE cp_grp_redist_array
     MODULE PROCEDURE cp_grp_redist_array_r1_d
     MODULE PROCEDURE cp_grp_redist_array_r1_z
     MODULE PROCEDURE cp_grp_redist_array_r2_d
     MODULE PROCEDURE cp_grp_redist_array_r2_z
  END INTERFACE cp_grp_redist_array
      
  INTEGER,DIMENSION(:,:),ALLOCATABLE,PUBLIC :: cp_grp_get_cp_rank


CONTAINS
  ! ==================================================================
  SUBROUTINE cp_grp_get_sizes(ngwk_l,ngw_l,&
       geq0_l,igeq0_l,firstk_g,lastk_g,first_g,last_g,i_g0_l)
    ! ==--------------------------------------------------------------==
    INTEGER, INTENT(out), OPTIONAL           :: ngwk_l, ngw_l
    LOGICAL, INTENT(out), OPTIONAL           :: geq0_l
    INTEGER, INTENT(out), OPTIONAL           :: igeq0_l, firstk_g, lastk_g, &
                                                first_g, last_g, i_g0_l

    INTEGER                                  :: first, firstk, i_g0, last, &
                                                lastk, n, nk

    nk = CEILING(REAL(nkpt%ngwk,kind=real_8)/REAL(parai%cp_nogrp,kind=real_8))
    firstk = parai%cp_inter_me*nk + 1
    lastk = MIN((parai%cp_inter_me+1)*nk,nkpt%ngwk)
    ! protect for OOB
    IF (firstk>nkpt%ngwk) THEN
       firstk = 1
       lastk = 0
    ENDIF
    IF (PRESENT(firstk_g)) firstk_g = firstk
    IF (PRESENT(lastk_g)) lastk_g = lastk
    IF (PRESENT(ngwk_l)) ngwk_l = lastk - firstk + 1
    ! 
    ! index of the G=0 element
    i_g0 = 0
    IF (tkpts%tkpnt) THEN
       IF (GEq0) i_g0 = 1 + ncpw%ngw
    ELSE
       IF (GEq0) i_g0 = 1
    ENDIF
    ! 
    ! 
    n = CEILING(REAL(ncpw%ngw,kind=real_8)/REAL(parai%cp_nogrp,kind=real_8))
    first = parai%cp_inter_me*n + 1
    last = MIN((parai%cp_inter_me+1)*n,ncpw%ngw)
    ! protect for OOB
    IF (first>ncpw%ngw) THEN
       first = 1
       last = 0
    ENDIF
    IF (PRESENT(first_g)) first_g = first
    IF (PRESENT(last_g)) last_g = last
    IF (PRESENT(ngw_l)) ngw_l = last - first + 1
    ! 
    ! 
    IF (PRESENT(geq0_l)) THEN
       geq0_l = i_g0.GE.firstk.AND.i_g0.LE.lastk
    ENDIF

    IF (PRESENT(igeq0_l)) THEN
       igeq0_l = 0
       IF (i_g0.GE.firstk.AND.i_g0.LE.lastk) igeq0_l = parai%cp_me
       CALL mp_sum(igeq0_l,parai%cp_grp)
    ENDIF

    IF (PRESENT(i_g0_l)) THEN
       i_g0_l = -HUGE(0)
       IF (i_g0.GE.firstk.AND.i_g0.LE.lastk) i_g0_l = i_g0 - firstk + &
            1
    ENDIF

       ! ==--------------------------------------------------------------==
    RETURN
    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_grp_get_sizes
  ! ==================================================================
  SUBROUTINE cp_grp_redist_z(DATA,ld,n)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: DATA(:,:)
    INTEGER                                  :: ld, n

    IF (parai%cp_nogrp.EQ.1) RETURN
    CALL mp_sum(DATA,ld*n,parai%cp_inter_grp)
    ! ==--------------------------------------------------------------==
    RETURN
    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_grp_redist_z
  ! ==================================================================
  SUBROUTINE cp_grp_redist_d(DATA,ld,n)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: DATA(:,:)
    INTEGER                                  :: ld, n

    IF (parai%cp_nogrp.EQ.1) RETURN
    CALL mp_sum(DATA,ld*n,parai%cp_inter_grp)
    ! ==--------------------------------------------------------------==
    RETURN
    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_grp_redist_d
  ! ==--------------------------------------------------------------==
  SUBROUTINE cp_grp_redist_array_r1_z(data,n)
    ! ==--------------------------------------------------------------==
    INTEGER :: n
    COMPLEX*16 :: data(*)
    INTEGER :: n_states_proc(3,parai%cp_nogrp),revcnt(parai%cp_nogrp),displ(parai%cp_nogrp),i
    !  ==--------------------------------------------------------------==
    IF(parai%CP_NOGRP.EQ.1) RETURN 
    do i = 1,parai%cp_nogrp
       CALL part_1d_get_blk_bounds(n,i-1,parai%cp_nogrp,n_states_proc(1,i),n_states_proc(2,i))
       n_states_proc(3,i)=n_states_proc(2,i)- n_states_proc(1,i)+1
    END DO
    do i = 1,parai%cp_nogrp
       revcnt(i)=n_states_proc(3,i)*2
    END DO
    displ=0
    do i= 2,parai%cp_nogrp
       displ(i)=displ(i-1)+n_states_proc(3,i-1)*2
    END DO
    call my_concatv_inplace(data,revcnt,displ,parai%cp_inter_grp)
    ! ==--------------------------------------------------------------==
    RETURN
    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_grp_redist_array_r1_z
  ! ==================================================================
  SUBROUTINE cp_grp_redist_array_r1_d(data,n)
    ! ==--------------------------------------------------------------==
    INTEGER :: n
    REAL*8 :: data(*)
    INTEGER :: n_states_proc(3,parai%cp_nogrp),revcnt(parai%cp_nogrp),displ(parai%cp_nogrp),i
    ! ==--------------------------------------------------------------==
    IF(parai%CP_NOGRP.EQ.1) RETURN 
    do i = 1,parai%cp_nogrp
       CALL part_1d_get_blk_bounds(n,i-1,parai%cp_nogrp,n_states_proc(1,i),n_states_proc(2,i))
       n_states_proc(3,i)=n_states_proc(2,i)- n_states_proc(1,i)+1
    END DO
    do i = 1,parai%cp_nogrp
       revcnt(i)=n_states_proc(3,i)
    END DO
    displ=0
    do i= 2,parai%cp_nogrp
       displ(i)=displ(i-1)+n_states_proc(3,i-1)
    END DO
    call my_concatv_inplace(data,revcnt,displ,parai%cp_inter_grp)
    ! ==--------------------------------------------------------------==
    RETURN
    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_grp_redist_array_r1_d
  ! ==================================================================
  SUBROUTINE cp_grp_redist_array_r2_z(data,ld,n)
    ! ==--------------------------------------------------------------==
    USE part_1d, ONLY: part_1d_get_blk_bounds
    IMPLICIT NONE
    INTEGER :: ld,n
    COMPLEX*16 :: data(ld,*)
    INTEGER :: n_states_proc(3,parai%cp_nogrp),revcnt(parai%cp_nogrp),displ(parai%cp_nogrp),i
    !  ==--------------------------------------------------------------==
    IF(parai%CP_NOGRP.EQ.1) RETURN 
    do i = 1,parai%cp_nogrp
       CALL part_1d_get_blk_bounds(n,i-1,parai%cp_nogrp,n_states_proc(1,i),n_states_proc(2,i))
       n_states_proc(3,i)=n_states_proc(2,i)- n_states_proc(1,i)+1
    END DO
    do i = 1,parai%cp_nogrp
       revcnt(i)=n_states_proc(3,i)*ld*2
    END DO
    displ=0
    do i= 2,parai%cp_nogrp
       displ(i)=displ(i-1)+n_states_proc(3,i-1)*ld*2
    END DO
    call my_concatv_inplace(data,revcnt,displ,parai%cp_inter_grp)
    ! ==--------------------------------------------------------------==
    RETURN
    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_grp_redist_array_r2_z
  ! ==================================================================
  SUBROUTINE cp_grp_redist_array_r2_d(data,ld,n)
    ! ==--------------------------------------------------------------==
    INTEGER :: ld,n
    REAL*8 :: data(ld,*)
    INTEGER :: n_states_proc(3,parai%cp_nogrp),revcnt(parai%cp_nogrp),displ(parai%cp_nogrp),i
    ! ==--------------------------------------------------------------==
    IF(parai%CP_NOGRP.EQ.1) RETURN 
    DO i = 1,parai%cp_nogrp
       CALL part_1d_get_blk_bounds(n,i-1,parai%cp_nogrp,n_states_proc(1,i),n_states_proc(2,i))
       n_states_proc(3,i)=n_states_proc(2,i)- n_states_proc(1,i)+1
    END DO
    DO i = 1,parai%cp_nogrp
       revcnt(i)=n_states_proc(3,i)*ld
    END DO
    displ=0
    DO i= 2,parai%cp_nogrp
       displ(i)=displ(i-1)+n_states_proc(3,i-1)*ld
    END DO
    CALL my_concatv_inplace(data,revcnt,displ,parai%cp_inter_grp)
    ! ==--------------------------------------------------------------==
    RETURN
    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_grp_redist_array_r2_d
  ! ==--------------------------------------------------------------==
  SUBROUTINE cp_grp_split_states(nstate,nstate_grp,first_state,last_state,norbpe,norbpe_cp,nst12,nst12_cp,nstate_cp,nst12_grp)
    IMPLICIT NONE
    INTEGER,INTENT(IN)                     :: nstate
    INTEGER,INTENT(OUT),optional           :: nstate_grp,first_state,last_state,norbpe_cp,norbpe,&
                                              nst12(2,0:parai%nproc-1),& ! nstate splitted among allgrp
                                              nst12_cp(2,0:parai%nproc-1),& ! local chunck from nst12 for members of cp_inter_grp
                                              nst12_grp(2,0:parai%nproc-1),& ! like nst12 but nstate_grp splitted, eg. min = first_state, max = last_state
                                              nstate_cp(3,0:parai%cp_nogrp-1)
    REAL(real_8)                           :: xstates,xsnow,xsaim
    INTEGER                                :: i,first,last,sum,state12(2,0:parai%nproc-1),state12_cp(2,0:parai%nproc-1),state_grps(3,0:parai%cp_nogrp-1),isub
    CHARACTER(*), PARAMETER                :: procedureN = 'cp_grp_split_states'
    CALL tiset(procedureN,isub)

    sum=nstate
    first=1
    last=nstate
    IF (parai%cp_nogrp .GT. 1) THEN
       CALL part_1d_get_blk_bounds(nstate,parai%cp_inter_me,parai%cp_nogrp,first,last)
       sum = last - first +1
    END IF

    CALL dist_size(sum,parai%nproc,state12)
    DO I=parai%NPROC,1,-1
       state12(1,i-1)=state12(1,i-1)+first-1
       state12(2,i-1)=state12(2,i-1)+first-1
    END DO

    IF (PRESENT(nst12_grp)) nst12_grp=state12
    IF (PRESENT (nstate_grp)) nstate_grp=sum
    IF (PRESENT (first_state)) first_state=first
    IF (PRESENT (last_state)) last_state=last

    CALL dist_size(nstate,parai%nproc,state12)
    IF (PRESENT(nst12)) nst12=state12
    IF (PRESENT(norbpe)) norbpe=state12(2,parai%mepos)-state12(1,parai%mepos)+1

    DO i=0,parai%nproc-1
       CALL part_1d_get_blk_bounds(state12(2,i)-state12(1,i)+1,parai%cp_inter_me,parai%cp_nogrp,first,last)
       state12_cp(1,i)=state12(1,i)+first-1
       state12_cp(2,i)=state12(1,i)+last-1
    END DO

    IF (PRESENT(nst12_cp)) nst12_cp=state12_cp
   
    IF (PRESENT(norbpe_cp)) norbpe_cp=state12_cp(2,parai%mepos)-state12_cp(1,parai%mepos)+1

    IF (PRESENT(nstate_cp)) THEN
       nstate_cp(1,:)=1
       nstate_cp(2,:)=nstate
       nstate_cp(3,:)=nstate
       IF (parai%cp_nogrp .GT. 1) THEN
          DO i = 0,parai%cp_nogrp-1
             CALL part_1d_get_blk_bounds(nstate,i,parai%cp_nogrp, &
                  nstate_cp(1,i),nstate_cp(2,i))
             nstate_cp(3,i)=nstate_cp(2,i)- nstate_cp(1,i)+1
          END DO
       END IF
    END IF
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE cp_grp_split_states
  ! ==--------------------------------------------------------------==
  SUBROUTINE cp_grp_split_atoms(na_split)
    INTEGER,INTENT(OUT)                    :: na_split(2,ions1%nsp,parai%cp_nogrp)
    INTEGER                                :: is,ia,iv,num(maxsys%nhxs,parai%cp_nogrp+1),&
                                              extra(parai%cp_nogrp), grp, aim, &
                                              distributed(parai%cp_nogrp), max_distributed, &
                                              counter, i, residual, isub, ierr
    LOGICAL,SAVE                           :: first=.true.
    INTEGER,ALLOCATABLE,SAVE               :: na_work(:,:,:)
    CHARACTER(*), PARAMETER                :: procedureN = 'cp_grp_split_atoms'
    ! ==------------------------------------------------------------==
    CALL tiset(procedureN,isub)

    IF (FIRST) THEN
       ALLOCATE(na_work(2,ions1%nsp,parai%cp_nogrp), stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate na_work',& 
            __LINE__,__FILE__)

       !determine the number of atoms with the same nlps_com%ngh(is)
       num=0
       DO iv=1,maxsys%nhxs
          DO is=1,ions1%nsp
             IF(nlps_com%ngh(is).EQ.iv) num(iv,1)=num(iv,1)+ions0%na(is)
          END DO
       END DO

       !divide the sum of atoms with the same nlps_com%ngh(is) into chunks for each group
       DO iv=1,maxsys%nhxs
          DO grp=1,parai%cp_nogrp
             num(iv,grp+1)=FLOOR(DBLE(num(iv,1))/DBLE(parai%cp_nogrp))
          END DO
       END DO

       !check for missing atoms
       distributed=0
       extra=0
       max_distributed=0
       DO iv=1,maxsys%nhxs
          counter=0
          DO grp=1,parai%cp_nogrp
             counter=counter+num(iv,grp+1)
          END DO
          IF (counter.LT.num(iv,1) ) THEN
             residual=num(iv,1)-counter
             !distribute es evenly as possible
             DO i=1,2
                max_distributed=minval(extra,1)
                DO grp=1,parai%cp_nogrp
                   !This group already got an extra atom
                   IF (extra(grp).GT.max_distributed) CYCLE
                   !Some atoms need to be distributed
                   IF (residual.GT.0) THEN
                      extra(grp)=extra(grp)+1
                      residual=residual-1               
                      num(iv,grp+1)=num(iv,grp+1)+1
                   END IF
                END DO
             END DO
          END IF
       END DO

       !sanity check
       DO iv=1,maxsys%nhxs
          counter=0
          DO grp=1,parai%cp_nogrp
             counter=counter+num(iv,grp+1)
          END DO
          IF (counter.NE.num(iv,1) ) CALL stopgm(procedureN, 'something went wront',& 
            __LINE__,__FILE__)
       END DO

       na_work(1,:,:)=0
       na_work(2,:,:)=-1

       !build up na_work
       DO iv=1,maxsys%nhxs
          DO grp=1,parai%cp_nogrp
             counter=0
             DO is=1,ions1%nsp
                IF (nlps_com%ngh(is).EQ.iv) THEN
                   aim=num(iv,grp+1)
                   max_distributed=0
                   IF (grp.GT.1) THEN
                      DO i=1,grp -1
                         max_distributed=max_distributed+num(iv,i+1)
                      END DO
                   END IF
                   DO ia=1,ions0%na(is)
                      IF(max_distributed.LE.counter.AND.aim.GE.(counter - max_distributed +1)) THEN
                         IF(na_work(1,is,grp).EQ.0) na_work(1,is,grp)=ia
                         na_work(2,is,grp)=ia
                         counter=counter+1
                      ELSE
                         counter=counter+1
                      END IF
                   END DO
                END IF
             END DO
          END DO
       END DO
       first=.FALSE.
    END IF
    
    na_split=na_work
    
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE cp_grp_split_atoms
  ! ==--------------------------------------------------------------==
  SUBROUTINE cp_grp_redist_dfnl_fnl(redist_fnl,redist_dfnl,nstate,ikind)
    LOGICAL,INTENT(IN)                     :: redist_fnl, redist_dfnl
    INTEGER,INTENT(IN)                     :: nstate, ikind
    CHARACTER(*), PARAMETER                :: procedureN = 'cp_grp_redist_dfnl_fnl'
    INTEGER                                :: na_grp(2,ions1%nsp,parai%cp_nogrp),&
                                              worksum(parai%cp_nogrp), ld_temp, il_temp(4),&
                                              is, ia, isa, isa0, ia0, ia_start, ia_end, &
                                              ia_sum, i, k, im, iv, ierr, grp, sendcnt, isub,&
                                              ii
    REAL(real_8), ALLOCATABLE              :: temp(:,:,:,:)
    ! ==------------------------------------------------------------==
    IF (parai%cp_nogrp.EQ.1) RETURN
    CALL tiset(procedureN,isub)
    !get cp_grp atomwise (d)fnl distribution
    CALL cp_grp_split_atoms(na_grp)
    
    worksum=0
    !count data sent to other groups
    DO grp=1,parai%cp_nogrp
       DO is=1,ions1%nsp
          worksum(grp)=worksum(grp)+(na_grp(2,is,grp)-na_grp(1,is,grp)+1)*nlps_com%ngh(is)
       END DO
    END DO   
    IF (redist_fnl) THEN       
       il_temp(1)=imagp
       il_temp(2)=MAXVAL(worksum)
       il_temp(3)=nstate
       il_temp(4)=parai%cp_nogrp
       ALLOCATE(temp(il_temp(1),il_temp(2),il_temp(3),il_temp(4)), stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate temp',& 
            __LINE__,__FILE__)
       !pack local fnl chunk
       !$omp parallel do private(i,isa0,ia0,is,ia_sum,iv,ia_start,ia_end,isa,ia,im)
       DO i=1,nstate
          isa0=0
          ia0=0
          DO is=1,ions1%nsp
             ia_sum=na_grp(2,is,parai%cp_inter_me+1)-na_grp(1,is,parai%cp_inter_me+1)+1
             IF(ia_sum.GT.0)THEN
                DO iv=1,nlps_com%ngh(is)
                   ia_start=ia0+1
                   ia_end=ia0+ia_sum                   
                   isa=isa0+na_grp(1,is,parai%cp_inter_me+1)-1
                   DO ia=ia_start,ia_end
                      isa=isa+1
                      DO im=1,imagp
                         temp(im,ia,i,parai%cp_inter_me+1)=fnl(im,isa,iv,i,ikind)
                      END DO
                   END DO
                   ia0=ia0+ia_sum
                END DO
             END IF
             isa0=isa0+ions0%na(is)
          END DO
       END DO
       sendcnt=il_temp(1)*il_temp(2)*il_temp(3)
       !exchange data
       CALL my_concat_inplace(temp,sendcnt,parai%cp_inter_grp)
       !unpack back to local fnl
       DO grp=1,parai%cp_nogrp
          IF (grp.EQ.parai%cp_inter_me+1) CYCLE
          !$omp parallel do private(i,isa0,ia0,is,ia_sum,iv,ia_start,ia_end,isa,ia,im)
          DO i=1,nstate
             isa0=0
             ia0=0
             DO is=1,ions1%nsp
                ia_sum=na_grp(2,is,grp)-na_grp(1,is,grp)+1
                IF(ia_sum.GT.0)THEN
                   DO iv=1,nlps_com%ngh(is)
                      ia_start=ia0+1
                      ia_end=ia0+ia_sum                   
                      isa=isa0+na_grp(1,is,grp)-1
                      DO ia=ia_start,ia_end
                         isa=isa+1
                         DO im=1,imagp
                            fnl(im,isa,iv,i,ikind)=temp(im,ia,i,grp)
                         END DO
                      END DO
                      ia0=ia0+ia_sum
                   END DO
                END IF
                isa0=isa0+ions0%na(is)
             END DO
          END DO
       END DO
       DEALLOCATE(temp, stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate temp',& 
            __LINE__,__FILE__)
    END IF

    IF (redist_dfnl) THEN       
       il_temp(1)=imagp
       il_temp(2)=MAXVAL(worksum)*3
       il_temp(3)=norbpe
       il_temp(4)=parai%cp_nogrp
       ALLOCATE(temp(il_temp(1),il_temp(2),il_temp(3),il_temp(4)), stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate temp',& 
            __LINE__,__FILE__)
       !pack local dfnl chunk
       !$omp parallel do private(i,ii,k,isa0,ia0,is,ia_sum,iv,ia_start,ia_end,isa,ia,im)
       DO i=parap%nst12(1,parai%mepos),parap%nst12(2,parai%mepos)
          ii=i-parap%nst12(1,parai%mepos)+1
          DO k=1,3
             isa0=0
             ia0=0
             DO is=1,ions1%nsp
                ia_sum=na_grp(2,is,parai%cp_inter_me+1)-na_grp(1,is,parai%cp_inter_me+1)+1
                IF(ia_sum.GT.0)THEN
                   DO iv=1,nlps_com%ngh(is)
                      ia_start=ia0+1
                      ia_end=ia0+ia_sum                   
                      isa=isa0+na_grp(1,is,parai%cp_inter_me+1)-1
                      DO ia=ia_start,ia_end
                         isa=isa+1
                         DO im=1,imagp
                            temp(im,ia,ii,parai%cp_inter_me+1)=dfnl(im,isa,iv,k,ii,ikind)
                         END DO
                      END DO
                      ia0=ia0+ia_sum
                   END DO
                END IF
                isa0=isa0+ions0%na(is)
             END DO
          END DO
       END DO
       sendcnt=il_temp(1)*il_temp(2)*il_temp(3)
       !exchange data
       CALL my_concat_inplace(temp,sendcnt,parai%cp_inter_grp)
       !unpack to local dfnl
       DO grp=1,parai%cp_nogrp
          IF (grp.EQ.parai%cp_inter_me+1) CYCLE
          !$omp parallel do private(i,ii,k,isa0,ia0,is,ia_sum,iv,ia_start,ia_end,isa,ia,im)
          DO i=parap%nst12(1,parai%mepos),parap%nst12(2,parai%mepos)
             ii=i-parap%nst12(1,parai%mepos)+1
             DO k=1,3               
                isa0=0
                ia0=0
                DO is=1,ions1%nsp
                   ia_sum=na_grp(2,is,grp)-na_grp(1,is,grp)+1
                   IF(ia_sum.GT.0)THEN
                      DO iv=1,nlps_com%ngh(is)
                         ia_start=ia0+1
                         ia_end=ia0+ia_sum                   
                         isa=isa0+na_grp(1,is,grp)-1
                         DO ia=ia_start,ia_end
                            isa=isa+1
                            DO im=1,imagp
                               dfnl(im,isa,iv,k,ii,ikind)=temp(im,ia,ii,grp)
                            END DO
                         END DO
                         ia0=ia0+ia_sum
                      END DO
                   END IF
                   isa0=isa0+ions0%na(is)
                END DO
             END DO
          END DO
       END DO
       DEALLOCATE(temp, stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate temp',& 
            __LINE__,__FILE__)
    END IF

    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE cp_grp_redist_dfnl_fnl
  ! ==--------------------------------------------------------------==
  
  ! ==================================================================
END MODULE cp_grp_utils
! ==================================================================

SUBROUTINE cp_grp_zero_g(c0,n,m,first_g,last_g)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
  INTEGER                                    :: n, m
  COMPLEX(real_8)                            :: c0(n,m)
  INTEGER                                    :: first_g, last_g

  COMPLEX(real_8), PARAMETER :: zzero = (0.0_real_8,0.0_real_8)

  INTEGER                                    :: i, j

  IF (first_g.EQ.1.AND.last_g.EQ.n) RETURN
  !$omp parallel do &
  !$omp            private(i) &
  !$omp            shared(first_g,last_g,N)
  DO j = 1,m
     DO i = 1,first_g-1
        c0(i,j) = zzero
     END DO
     DO i = last_g+1,n
        c0(i,j) = zzero
     END DO
  END DO
  !$omp end parallel do
  ! ==--------------------------------------------------------------==
  RETURN
  ! ==--------------------------------------------------------------==
END SUBROUTINE cp_grp_zero_g
! ==================================================================

SUBROUTINE cp_grp_zero_states(c0,n,m,ibeg,iend)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE parac , ONLY:parai
  IMPLICIT NONE
  INTEGER                                    :: n, m
  COMPLEX(real_8)                            :: c0(n,m)
  INTEGER                                    :: ibeg(0:parai%cp_nproc-1), &
                                                iend(0:parai%cp_nproc-1)

  COMPLEX(real_8), PARAMETER :: zzero = (0.0_real_8,0.0_real_8)

  INTEGER                                    :: nn

  IF (ibeg(parai%cp_me).GT.1) THEN
     nn = n * ( ibeg(parai%cp_me) - 1 )
     CALL zcopy(nn,zzero,0,c0(1,1),1)
  ENDIF
  IF (iend(parai%cp_me).LT.m) THEN
     nn = n * ( m - iend(parai%cp_me) - 1 )
     CALL zcopy(nn,zzero,0,c0(1,iend(parai%cp_me)+1),1)
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
  ! ==--------------------------------------------------------------==
END SUBROUTINE cp_grp_zero_states
! ==================================================================
SUBROUTINE cp_grp_copy_wfn_to_local(c0,ld,C0_l,LD_l,ibeg,M_l,n)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
  INTEGER                                    :: ld
  COMPLEX(real_8)                            :: c0(ld,*)
  INTEGER                                    :: LD_l
  COMPLEX(real_8)                            :: C0_l(LD_l,*)
  INTEGER                                    :: ibeg, M_l, n

  INTEGER                                    :: i

  IF (LD_l.GT.0.AND.ld.GT.0) THEN
     ! copy the right C0 block into the local array
     DO i=1,n
        CALL dcopy(2*M_l,c0(ibeg,i),1,C0_l(1,i),1)
     END DO
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
  ! ==--------------------------------------------------------------==
END SUBROUTINE cp_grp_copy_wfn_to_local
! ==================================================================
SUBROUTINE cp_grp_copy_local_to_wfn(C0_l,LD_l,c0,ld,ibeg,M_l,n)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
  INTEGER                                    :: LD_l
  COMPLEX(real_8)                            :: C0_l(LD_l,*)
  INTEGER                                    :: ld
  COMPLEX(real_8)                            :: c0(ld,*)
  INTEGER                                    :: ibeg, M_l, n

  INTEGER                                    :: i

  IF (LD_l.GT.0.AND.ld.GT.0) THEN
     ! copy back to the right C0 block
     DO i=1,n
        CALL dcopy(2*M_l,C0_l(1,i),1,c0(ibeg,i),1)
     END DO
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
  ! ==--------------------------------------------------------------==
END SUBROUTINE cp_grp_copy_local_to_wfn
! ==================================================================
SUBROUTINE cp_grp_add_local_to_wfn(C0_l,LD_l,c0,ld,ibeg,M_l,n)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
  INTEGER                                    :: LD_l
  COMPLEX(real_8)                            :: C0_l(LD_l,*)
  INTEGER                                    :: ld
  COMPLEX(real_8)                            :: c0(ld,*)
  INTEGER                                    :: ibeg, M_l, n

  INTEGER                                    :: i, j

  IF (LD_l.GT.0.AND.ld.GT.0) THEN
     ! copy back to the right C0 block
     !$omp parallel do private(i,j) shared(N,M_l,ibeg)
     DO i=1,n
        DO j=1,M_l
           c0(ibeg+j-1,i) = c0(ibeg+j-1,i) + C0_l(j,i)
        END DO
     END DO
     !$omp end parallel do
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
  ! ==--------------------------------------------------------------==
END SUBROUTINE cp_grp_add_local_to_wfn

