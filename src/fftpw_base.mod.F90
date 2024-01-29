!=----------------------------------------------------------------------=!
MODULE fftpw_base
!=----------------------------------------------------------------------=!

  USE fftpw_stick_base,      ONLY: sticks_map
  USE fftpw_types,           ONLY: PW_fft_type_descriptor

  IMPLICIT NONE

  TYPE( PW_fft_type_descriptor ) :: dfft
  TYPE( PW_fft_type_descriptor ) :: dfftp

  TYPE (sticks_map) :: smap


  SAVE
  PRIVATE

  PUBLIC :: dfft
  PUBLIC :: dfftp
  PUBLIC :: smap

END MODULE fftpw_base
