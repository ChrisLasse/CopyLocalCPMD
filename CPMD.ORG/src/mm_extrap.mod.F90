MODULE mm_extrap
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! storage for wavefunction extrapolation in BOMD.
  ! scold for Vanderbilt
  COMPLEX(real_8), ALLOCATABLE, SAVE :: cold(:,:,:,:), scold(:,:,:,:)

  INTEGER, SAVE :: numcold,nnow

END MODULE mm_extrap
