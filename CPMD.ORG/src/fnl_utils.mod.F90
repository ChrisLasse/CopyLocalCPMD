MODULE fnl_utils
  USE cp_grp_utils,                    ONLY: cp_grp_split_atoms  
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE nlps,                            ONLY: nlps_com,&
                                             imagp
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE sfac,                            ONLY: dfnl,fnl,&
                                             eigr
  USE system,                          ONLY: iatpt,&
                                             parap
  USE timer,                           ONLY: tihalt,&
                                             tiset
  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: unpack_fnl
  PUBLIC :: unpack_dfnl
  PUBLIC :: split_atoms_btw_buffers
  
CONTAINS
  ! ==================================================================
  SUBROUTINE unpack_dfnl(packed,unpacked_k,unpacked,na_grp,kpnt,ikind)
    ! ==--------------------------------------------------------------==
    ! == unpack packed_dfnl, mainly to get rid of pointer attributes  ==
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: na_grp(:,:,:),ikind
    LOGICAL,INTENT(IN)                       :: kpnt    
    REAL(real_8),INTENT(IN)                  :: packed(:,:)
    REAL(real_8),INTENT(INOUT)               :: unpacked_k(:,:,:,:,:,:),unpacked(:,:,:,:)
    INTEGER                                  :: i,ii,offset,isa0,is,iv,ia,isa,k,im

    IF(kpnt) THEN
       !$omp parallel do private(i,ii,offset,k,isa0,is,iv,ia,isa,im) proc_bind(close)
       DO i=parap%nst12(1,parai%mepos),parap%nst12(2,parai%mepos)
          ii=i-parap%nst12(1,parai%mepos)+1
          offset=0
          isa0=0
          DO is=1,ions1%nsp
             DO iv=1,nlps_com%ngh(is)
                DO k=1,3
                   DO ia=na_grp(1,is,parai%cp_inter_me +1),&
                        na_grp(2,is,parai%cp_inter_me +1)
                      isa=isa0+ia
                      DO im=1,imagp
                         offset=offset+1
                         unpacked_k(im,isa,iv,k,ii,ikind)=packed(offset,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
             isa0=isa0+ions0%na(is)
          ENDDO
       ENDDO
    ELSE
       !$omp parallel do private(i,ii,offset,k,isa0,is,iv,ia,isa) proc_bind(close)
       DO i=parap%nst12(1,parai%mepos),parap%nst12(2,parai%mepos)
          ii=i-parap%nst12(1,parai%mepos)+1
          offset=0
          isa0=0
          DO is=1,ions1%nsp
             DO iv=1,nlps_com%ngh(is)
                DO k=1,3
                   DO ia=na_grp(1,is,parai%cp_inter_me +1),&
                        na_grp(2,is,parai%cp_inter_me +1)
                      isa=isa0+ia
                      offset=offset+1
                      unpacked(isa,iv,k,ii)=packed(offset,i)
                   ENDDO
                ENDDO
             ENDDO
             isa0=isa0+ions0%na(is)
          ENDDO
       ENDDO
    END IF
  END SUBROUTINE unpack_dfnl
  ! ==================================================================
  SUBROUTINE unpack_fnl(packed,unpacked_k,unpacked,na_grp,nst_start,nst_end,kpnt,ikind)
    ! ==--------------------------------------------------------------==
    ! == unpack packed_dfnl, mainly to get rid of pointer attributes  ==
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: na_grp(:,:,:),nst_start,nst_end,ikind
    LOGICAL,INTENT(IN)                       :: kpnt
    REAL(real_8),INTENT(IN)                  :: packed(:,:)
    REAL(real_8),INTENT(INOUT)               :: unpacked_k(:,:,:,:,:),unpacked(:,:,:)
    INTEGER                                  :: i,offset,isa0,is,iv,ia,isa,k,im
    
    IF(kpnt) THEN
       !$omp parallel do private(i,offset,isa0,is,iv,ia,isa,im) proc_bind(close)
       DO i=nst_start,nst_end
          isa0=0
          offset=0
          DO is=1,ions1%nsp
             DO iv=1,nlps_com%ngh(is)
                DO ia=na_grp(1,is,parai%cp_inter_me +1),&
                     na_grp(2,is,parai%cp_inter_me +1)
                   isa=isa0+ia
                   DO im=1,imagp
                      offset=offset+1
                      unpacked_k(im,isa,iv,i,ikind)=packed(offset,i)
                   ENDDO
                ENDDO
             ENDDO
             isa0=isa0+ions0%na(is)
          ENDDO
       ENDDO
    ELSE
       !$omp parallel do private(i,offset,isa0,is,iv,ia,isa) proc_bind(close)
       DO i=nst_start,nst_end
          isa0=0
          offset=0
          DO is=1,ions1%nsp
             DO iv=1,nlps_com%ngh(is)
                DO ia=na_grp(1,is,parai%cp_inter_me +1),&
                     na_grp(2,is,parai%cp_inter_me +1)
                   isa=isa0+ia
                   offset=offset+1
                   unpacked(isa,iv,i)=packed(offset,i)
                ENDDO
             ENDDO
             isa0=isa0+ions0%na(is)
          ENDDO
       ENDDO
    END IF

  END SUBROUTINE unpack_fnl
  ! ==================================================================
  SUBROUTINE  split_atoms_btw_buffers(buffer_count,fractions,na_split,map_work,only_uspp)
    ! ==--------------------------------------------------------------==
    ! == Splits the atoms into buffer_count chunks according to       ==
    ! == fractions, sometimes only uspp have to be considered, in     ==
    ! == this case, specify uspp = .true.                             ==
    ! ==--------------------------------------------------------------==

    INTEGER, INTENT(IN)                      :: buffer_count
    REAL(real_8),INTENT(IN)                  :: fractions(buffer_count)
    INTEGER, INTENT(OUT)                     :: na_split(2,ions1%nsp,buffer_count,parai%cp_nogrp)
    INTEGER, INTENT(OUT), OPTIONAL           :: map_work(ions1%nat+1,buffer_count,parai%cp_nogrp)
    LOGICAL, INTENT(IN), OPTIONAL            :: only_uspp
    INTEGER                                  :: tot_work(parai%cp_nogrp), &
                                                na_work(2,ions1%nsp,parai%cp_nogrp), &
                                                div(buffer_count,parai%cp_nogrp),&
                                                grp, is, iv, ia, isa, isa0, buff, aim, count, &
                                                map(ions1%nat+1,buffer_count,parai%cp_nogrp)
    LOGICAL                                  :: uspp
    !if cp_grps are active, fnl and dfnl are distributed atomwise, so get the na(is) for the groups:
    CALL cp_grp_split_atoms(na_work)

    IF (PRESENT(only_uspp)) THEN
       uspp=only_uspp
    ELSE
       uspp=.FALSE.
    END IF

    tot_work=0
    DO grp=1,parai%cp_nogrp
       DO is=1,ions1%nsp
          IF (pslo_com%tvan(is).OR..NOT.uspp) THEN
             tot_work(grp)=tot_work(grp)+(na_work(2,is,grp)-na_work(1,is,grp)+1)*nlps_com%ngh(is)
          END IF
       END DO
    END DO

    div=0
    DO grp=1,parai%cp_nogrp
       DO buff=1,buffer_count
          div(buff,grp)=CEILING(REAL(tot_work(grp),real_8)*fractions(buff))
       END DO
    END DO   
    
    map=0
    DO grp=1,parai%cp_nogrp
       count=0
       buff=1
       aim=div(buff,grp)
       isa0=0
       DO is=1,ions1%nsp
          IF (pslo_com%tvan(is).OR..NOT.uspp) THEN
             DO ia=na_work(1,is,grp),na_work(2,is,grp)
                isa=isa0+ia
                !fail-save last buffer takes any remaining atom
                IF (buff.EQ.buffer_count) THEN
                   map(isa,buff,grp)=1
                   count=count+nlps_com%ngh(is)
                   !accumulate atoms in this buffer
                ELSE IF(count+nlps_com%ngh(is).LE.aim) THEN
                   map(isa,buff,grp)=1
                   count=count+nlps_com%ngh(is)
                   !the current atom would not fit into this buffer, add it to the next buffer
                ELSE IF( count+nlps_com%ngh(is).GT.aim) THEN
                   buff=buff+1
                   aim=aim+div(buff,grp)
                   map(isa,buff,grp)=1
                   count=count+nlps_com%ngh(is)
                ENDIF
             ENDDO
          END IF
          !isa0/isa is always the 'real' isa
          isa0=isa0+ions0%na(is)
       ENDDO
    END DO
    
    na_split(1,:,:,:)=0
    na_split(2,:,:,:)=-1

    DO grp=1,parai%cp_nogrp
       DO buff=1,buffer_count
          isa0=0
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                isa=isa0+ia
                !check if this atom is included in this buffer
                IF (map(isa,buff,grp).EQ.1) then
                   !Check if we found the first atom of this species
                   IF (na_split(1,is,buff,grp).EQ.0) na_split(1,is,buff,grp)=ia
                   !Set the last atom of this species
                   na_split(2,is,buff,grp)=ia
                ENDIF
             ENDDO
             isa0=isa0+ions0%na(is)
          ENDDO
       ENDDO
    END DO

    IF (PRESENT(map_work)) THEN
       DO grp=1,parai%cp_nogrp
          DO buff=1,buffer_count
             count=0   
             DO isa=1,ions1%nat
                IF (map(isa,buff,grp).EQ.1) THEN
                   is=iatpt(2,isa)
                   count=count+nlps_com%ngh(is)
                END IF
             END DO
             map(ions1%nat+1,buff,grp)=count
          END DO
       END DO
       map_work=map
    END IF

  END SUBROUTINE split_atoms_btw_buffers

END MODULE fnl_utils
