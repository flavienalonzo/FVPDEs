!            **************************
!            **  FUNCTION AMATLOC  **
!            **************************
!******************************************************************************
!     * ECN
!     *--------------------------------------
!     * Ce sous programme calcule sur un triangle (jt) donne
!     * la matrice elementaire locale pour l'aproximation P1Sommets ou P1Milieux
!     * MATLOC fourni une matrice 3X3
!******************************************************************************
FUNCTION amatloc(jt)
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
  INTEGER, INTENT(in)                 :: jt       ! numero du triangle
  REAL(kind=long), DIMENSION(3,3)     :: amatloc
  !----------------------------------
  ! Declaration des variables locales
  !----------------------------------
  CHARACTER(len=6)                  :: oldprf
  INTEGER                           :: i, j, k, iloc, jloc
  REAL(kind=long), DIMENSION(3)     :: x, y
  REAL(kind=long), DIMENSION(2,3)   :: DF  
  LOGICAL :: trouver
  !-------------------
  ! Debut du programme
  !-------------------
  oldprf = prefix
  prefix = 'AMATLO'

  !------
  ! Corps
  !------
  i = NuSoK(1,jt) ; j = NuSoK(2,jt) ; k = NuSoK(3,jt)   !! i,j,k numero globale
  !! 1,2,3 numero locale

  x(1) = coordS(1,i) ; y(1) = coordS(2,i)
  x(2) = coordS(1,j) ; y(2) = coordS(2,j)
  x(3) = coordS(1,k) ; y(3) = coordS(2,k)
  
  Select case(ChoixSchema)

  case(P1sommets)

     ! calcul des gradients des fonctions de base
     ! DF(1:2,1) : les composantes de grad phi_1
     DF(1,1) = y(3)-y(2)    ; DF(2,1) = - (x(3)-x(2))
     DF(1,2) = -(y(3)-y(1)) ; DF(2,2) = x(3)-x(1)
     DF(1,3) = y(2)-y(1)    ; DF(2,3) = -( x(2)-x(1) )

     DO jloc = 1, 3
        DO iloc = 1, 3
           amatloc(iloc,jloc) = ( (SxxK(jt)*DF(1,iloc)+SxyK(jt)*DF(2,iloc) )*DF(1,jloc)&
                + ( SxyK(jt)*DF(1,iloc)+SyyK(jt)*DF(2,iloc) )*DF(2,jloc))/(2.D0 * AireK(jt))
        END DO
     END DO

  case(P1milieux,P1milieuxmonotone)
     ! calcul des gradients des fonctions de base
     ! DF(1:2,1) : les composantes de grad phi_1

     DF(1,1) = y(3)-y(2) ; DF(2,1) = - (x(3)-x(2))
     DF(1,2) = -(y(3)-y(1)) ; DF(2,2) = x(3)-x(1)
     DF(1,3) = y(2)-y(1) ; DF(2,3) = -( x(2)-x(1) )


     DO jloc = 1, 3
        DO iloc = 1, 3
           amatloc(iloc,jloc) = ( (SxxK(jt)*DF(1,iloc)+SxyK(jt)*DF(2,iloc) )*DF(1,jloc)&
                + (SxyK(jt)*DF(1,iloc)+SyyK(jt)*DF(2,iloc) )*DF(2,jloc))/(AireK(jt))
        END DO
     END DO


  case default
     !stop'pb de choixschema dans AMATLOC'

  end select



  !------------
  ! Impressions
  !------------
  IF (iprint >= 5) THEN
     !CALL prvari(uprint, ' Amatloc Numero triangle = ', jt) 
     !DO iloc=1,3
     !   WRITE (uprint,110) (amatloc(iloc,jloc), jloc=1,3)
     !END DO
     trouver = .false.
     DO iloc=1,3
        DO jloc=1,3
           If (jloc==iloc .and. amatloc(iloc,jloc)<0.D0 ) then 
              trouver = .true.
           ENDIF
           If (jloc.ne.iloc .and. amatloc(iloc,jloc)>0.D0 ) then 
              trouver = .true.
           ENDIF
        END DO
     END DO
     IF (trouver ) then 
        CALL prvari(uprint, ' Amatloc Numero triangle = ', jt) 
        DO iloc=1,3
           WRITE (uprint,110) (amatloc(iloc,jloc), jloc=1,3)
        END DO
     END If
  END IF
110 FORMAT (3(E16.9, 2X))
  !-----------------
  ! Fin du programme
  !-----------------
  prefix = oldprf


END FUNCTION amatloc
