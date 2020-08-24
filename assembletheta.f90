!            **************************
!            **  SUBROUTINE ASSEMBLETheta  **
!            **************************
!******************************************************************************
!     *--------------------------------------
!     * Ce sous programme assemble la contribution de la
!     * la matrice elementaire locale Aloc du triangle jt

!     * Fonctionnement du programme
!     *----------------------------
!     * Utilise les modules longr, imprime, parmmmage
!******************************************************************************
SUBROUTINE ASSEMBLETHETA( A)
  !--------
  ! Modules
  !--------
  USE longr
  USE imprime
  USE parmmage

  IMPLICIT NONE

  !--------------------------
  ! Declaration des arguments
  !--------------------------
  TYPE(MatCreux)       :: A
  !----------------------------------
  ! Declaration des variables locales
  !----------------------------------
  CHARACTER(len=6)      :: oldprf
  INTEGER               :: i, j, k,  jt 
  REAL(kind=long), DIMENSION(3)     :: x, y


  !-------------------
  ! Debut du programme
  !-------------------
  oldprf = prefix
  prefix = 'ASSEMT'

  !------
  ! Corps
  !------
Do jt = 1, Nbt 
  CALL Ajout (jt, jt, Theta*AireK(jt), A )
END Do
  !------------
  ! Impressions
  !------------
!!$  IF (iprint >= 2) THEN
!!$     CALL prvari(uprint, ' Numero triangle = ', jt) 
!!$     WRITE (uprint,*) ( (matloc(iloc,jloc), jloc=1,3),iloc=1,3)
!!$  END IF

  !-----------------
  ! Fin du programme
  !-----------------
  prefix = oldprf

  RETURN
END SUBROUTINE assembletheta








