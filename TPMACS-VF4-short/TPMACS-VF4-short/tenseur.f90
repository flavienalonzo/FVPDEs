!            **************************
!            **  FUNCTION tenseurAMATLOC  **
!            **************************
!******************************************************************************
!     * Ecole Centrale de Nantes
! calcul sur chaque triangle le tenseur de diffusion 
! SxxK, SxyK, SyyK
!
!******************************************************************************
SUBROUTINE tenseur(p)
  !--------
  ! Modules
  !--------
  USE longr
  USE imprime
  USE parmmage

  IMPLICIT NONE

  !----------------------------------
  ! Declaration des variables locales
  !----------------------------------
  CHARACTER(len=6)                  :: oldprf
  Integer, intent(in)               :: p
  INTEGER                           :: i, j, k, jt
  REAL(kind=long), DIMENSION(2)     :: z
  REAL(kind=long), DIMENSION(3)     :: x, y
  REAL(kind=long) :: sally, a, kscalaire, merci1, merci2, merci4

  !-------------------
  ! Debut du programme
  !-------------------
  oldprf = prefix
  prefix = 'TENSEUR'

  !------
  ! Corps
  !------
  Allocate(SxxK(Nbt),SxyK(Nbt),SyyK(Nbt))
  !! 
  !!  On peut ajouter la matrice M pour V
  !!  Allocate(MxxK(Nbt), MxyK(Nbt), MyyK(Nbt))

  DO jt = 1, Nbt

     i = NuSoK(1,jt) ; j = NuSoK(2,jt) ; k = NuSoK(3,jt)  !! i,j,k num. globale
     !! 1,2,3 numero locale

     x(1) = coordS(1,i) ; y(1) = coordS(2,i)
     x(2) = coordS(1,j) ; y(2) = coordS(2,j)
     x(3) = coordS(1,k) ; y(3) = coordS(2,k)

     z(1) = (x(1)+x(2)+x(3))/ 3.D0
     z(2) = (y(1)+y(2)+y(3))/ 3.D0

     Select case (p)

     case(1) 
        !-------
        ! S = Id
        !-------
        SxxK(jt)= 1.D0
        SyyK(jt) = 1.D0
        SxyK(jt)= 0.D0 
        !write(*,*)SxxK(jt)

     case(2)
        !-------------------------
        !anisotropie exponentielle
        !-------------------------
        SxxK(jt)= kscalaire(z(1),z(2),choixkscalaire)
        SyyK(jt) = SxxK(jt)
        SxyK(jt)= 0.D0

     case(3)
        !--------------------
        ! anisotropie moderee
        !--------------------
        SxxK(jt) = deltax 
        SxyK(jt) = deltaxy
        SyyK(jt) = deltay

     case(4)
        !-----------------------
        ! anisotropie heterogene
        !-----------------------
        sally= 1.D0 / ( z(1)**2 + z(2)**2 )
        SxxK(jt) = sally*( delta*z(1)**2 + z(2)**2 )
        SxyK(jt) = -(1.D0 - delta)*sally* z(1) * z(2)
        SyyK(jt) = sally*( z(1)**2 + delta*z(2)**2)

     case(5)
        !-----------------------
        ! anisotropie discontinu
        !-----------------------
        If( z(2) > 0.5D0) then
           SxxK(jt) = 1.D0
           SxyK(jt) = 0.D0
           SyyK(jt) = 1.D0
        Else
           SxxK(jt) = 1.5D0
           SxyK(jt) = 1.D0
           SyyK(jt) = 1.D0
        End if

     case(6)
        !---------------
        !oblique barrier
        !---------------
        if ( z(2) - 0.2D0 *(z(1)-0.5D0) - 0.475D0 < 0) then
           SxxK(jt) = 1.D0
           SxyK(jt) = 0.D0
           SyyK(jt) = 1.D0
        else
           if ( z(2) - 0.2D0 *(z(1)-0.5D0) - 0.475D0 > 0 .AND. &
                & z(2) - 0.2D0 *(z(1)-0.5D0) - 0.525D0 < 0) then

              SxxK(jt) = 1.D-2
              SxyK(jt) = 0.D0
              SyyK(jt) = 1.D-2
           else
              if(  z(2) - 0.2D0 *(z(1)-0.5D0) - 0.525D0 > 0) then
                 SxxK(jt) = 1.D0
                 SxyK(jt) = 0.D0
                 SyyK(jt) = 1.D0

              end if
           end if
        end if

     case(7)
        !--------------------
        ! anisotropie moderee
        !--------------------
        SxxK(jt) = z(1) 
        SxyK(jt) = 0.D0
        SyyK(jt) = z(2)

     case(8)
         !------------------------
         ! Isotropie par morceaux
         !------------------------
      if (ntypt(jt)==100) then
         SxxK(jt) = deltax 
         SxyK(jt) = 0.D0
         SyyK(jt) = deltay
      else 
         SxxK(jt) = 5*deltax 
         SxyK(jt) = 0.D0
         SyyK(jt) = 5*deltay
      end if
     
        
     end select
  end DO

  !------------
  ! Impressions
  !------------


  !-----------------
  ! Fin du programme
  !-----------------
  prefix = oldprf


END SUBROUTINE tenseur
