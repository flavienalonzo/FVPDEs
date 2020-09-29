!            **************************
!            **  SUBROUTINEUbord CLIMITE**
!            **************************
!*****************************************************************************
!     *--------------------------------------
!     * Ce sous programme donne les valeurs des conditions 
!     * aux limites dirichlet aux sommets du bord ou au milieu des bords
!******************************************************************************
SUBROUTINE UBORD(chxpb, chxsch, temps)
  !--------
  ! Modules
  !--------
  USE longr
  USE imprime
  USE parmmage
  USE fsourcemod
  use fsourcebreast

  IMPLICIT NONE

  !--------------------------
  ! Declaration des arguments
  !--------------------------
  Integer, intent(in) :: chxpb, chxsch
  REAL(kind=long), intent(in)      :: temps
  !--------------------------------
  ! Declaration des variables locales
  !----------------------------------
  CHARACTER(len=6)      :: oldprf
  integer               :: i

  !-------------------
  ! Debut du programme
  !-------------------
  oldprf = prefix
  prefix = 'UBORD '

  !------
  ! Corps
  !------

  select case  (Chxsch)
  case(VF4)

  case(P1Sommets)
     !! il faut GbU et GbV !!! a corriger later
     DO i = NsInt+1, Nbs
        Select case ( NtypS(i) ) 
        case (dirichlet)
           Gb(i) = gbord( temps, coordS(1,i), coordS(2,i),chxpb  )
        case(neumann)
           !! On fait rien 
        End Select
     END DO
  case(P1milieux,P1milieuxmonotone)
     DO i = Nsegint+1, Nseg
        Gb(i) = gbord( temps, coordMSeg(1,i), coordMSeg(2,i),chxpb  )
     END DO
  end select

  !print*,'max Gb dans ubord = ', maxval(Gb)
  !print*,'min Gb dans ubord = ', minval(Gb)

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
END SUBROUTINE UBORD
