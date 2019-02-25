C hello emacs, this is -*- fortran -*-
C These are additional unit names and numbers
#ifdef WRITE_FORCES
      DATA UNAME(JNTRJF)     /  'TRJF' /
      DATA IUNUM(JNTRJF)     /   IOTRJF /
#endif
#ifdef READ_COORD
      DATA UNAME(JNRDCO)     /  'RDCO' /
      DATA IUNUM(JNRDCO)     /   IORDCO /
#endif
