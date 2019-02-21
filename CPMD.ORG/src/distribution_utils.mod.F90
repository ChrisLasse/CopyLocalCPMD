#include "cpmd_global.h"

MODULE distribution_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: dist_size

CONTAINS

    SUBROUTINE dist_size(size,num,dist_array,nblock,nbmax,nblocal,iloc,fw)
    ! ==--------------------------------------------------------------==
    ! == Distributes size into num equally parts
    INTEGER,INTENT(IN)                       :: size,num
    INTEGER,INTENT(OUT)                      :: dist_array(2,0:num-1)
    INTEGER,INTENT(IN),OPTIONAL              :: nblock, fw, iloc
    INTEGER,INTENT(OUT),OPTIONAL             :: nbmax, nblocal
    INTEGER                                  :: ip, n
    REAL(real_8)                             :: xaim, xsize
    CHARACTER(*), PARAMETER                  :: procedureN = 'distribution'
! ==--------------------------------------------------------------==

    
    IF (PRESENT(nblock) .AND. nblock*num .GE. size) THEN
       xsize=REAL(nblock,kind=real_8)
    ELSE
       xsize=REAL(size,kind=real_8)
    END IF

    dist_array(1,:)=0
    dist_array(2,:)=-1

    IF (PRESENT(FW).AND.FW.EQ.1) THEN
       DO ip=1,num      
          xaim = ip * xsize/num
          dist_array(2,ip-1)=NINT(xaim)
          IF (ip.EQ.1) THEN
             dist_array(1,ip-1)=1
          !failsave
          ELSE
             dist_array(1,ip-1)=dist_array(2,ip-2)+1
          END IF
          IF (NINT(xaim).GT.size) dist_array(2,ip-1)=size
          IF (ip.EQ.num) dist_array(2,ip-1)=size
          IF (dist_array(1,ip-1).GT.size) dist_array(1,ip-1)=size+1
       ENDDO
    ELSE
       DO ip=num,1,-1
          xaim = xsize-(num-ip+1)*xsize/num+1
          dist_array(1,ip-1)=NINT(xaim)
          IF (ip.EQ.num) THEN
             dist_array(2,ip-1)=size
          !failsave
          ELSE
             dist_array(2,ip-1)=dist_array(1,ip)-1
          END IF
          IF (ip.EQ.1) THEN
             dist_array(1,ip-1)=1
          ENDIF
       ENDDO
    END IF

    IF (PRESENT(nbmax)) THEN
       nbmax=0
       DO ip=0,num-1
          nbmax=MAX(nbmax,dist_array(2,ip)-dist_array(2,ip)+1)
       END DO
    END IF
    
    IF (PRESENT(nblocal)) THEN
       IF (PRESENT(iloc)) THEN
          IF (iloc.GE.num) CALL stopgm(procedureN,'Programming error',&
               __LINE__,__FILE__)         
          nblocal=dist_array(2,iloc)-dist_array(1,iloc)+1
       ELSE
          CALL stopgm(procedureN,'Programming error',&
               __LINE__,__FILE__)
       END IF
    END IF

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dist_size
  

END MODULE distribution_utils
