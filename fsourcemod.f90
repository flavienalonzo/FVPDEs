Module fsourcemod

Contains

 !======================================================
!=====================================================
  FUNCTION gbord(x,y,m)
    ! Calcul la solution exacte de l'equation
    !--------
    ! Modules
    !--------
    USE longr

    IMPLICIT NONE

    !--------------------------
    ! Declaration des arguments
    !--------------------------
    REAL(kind=long), INTENT(in)     :: x, y
    integer :: m
    REAL(kind=long)                 :: gbord
    !----------------------------------
    ! Declaration des variables locales
    !----------------------------------
    CHARACTER(len=6)      :: oldprf
    REAL(kind=long)                 :: pi

    !-------------------
    ! Debut du programme
    !-------------------
    oldprf = prefix
    prefix = 'Gbord'

    !------
    ! Corps
    !------
    pi = 4.D0*ATAN(1.D0)

    select case(m)
    case(1,7)
       gbord = 1.
    case(2,8)
       gbord  = x+y
    case(3) 
       gbord  = x*x - y*y
    case(4)
       gbord = cos(5.*pi*(x+y))
    case(5)
       gbord = x * (1-x)*y*(1-y)
    case(6)
       gbord= sin(pi * x)*sin(pi *y)
    case default 
       print*, ' pb gbord'
       stop
    end select
    !-----------------
    ! Fin du programme
    !-----------------
    prefix = oldprf

    RETURN
  END FUNCTION gbord



  

FUNCTION vitesse(x,y, m) result(Veloc)
    !--------
    ! Modules
    !--------
    USE longr

    IMPLICIT NONE

    !--------------------------
    ! Declaration des arguments
    !--------------------------
    REAL(kind=long), INTENT(in)     :: x , y
    Integer :: m
    REAL(kind=long),  DIMENSION(2)     :: Veloc
    !----------------------------------
    ! Declaration des variables locales
    !----------------------------------
    CHARACTER(len=6)                :: oldprf
    REAL(kind=long)                 :: pi
    REAL(kind=long) ::lea,leaa
    !-------------------
    ! Debut du programme
    !-------------------
    oldprf = prefix
    prefix = 'VITESS'


    !------
    ! Corps
    !------
    pi=4.D0*ATAN(1.D0)
    Select case (m)
    case(1:6)
       Veloc(1) = 0.D0  ;  Veloc(1) = 0.D0
    case (7,8)
      Veloc(1)= 1.D0 ;  Veloc(2)=1.D0
    end Select
    !----------------
    ! Fin du programme
    !-----------------
    prefix = oldprf

    RETURN
  END FUNCTION vitesse

  

FUNCTION fsource(x,y, m)
    !--------
    ! Modules
    !--------
    USE longr

    IMPLICIT NONE

    !--------------------------
    ! Declaration des arguments
    !--------------------------
    REAL(kind=long), INTENT(in)     :: x , y
    Integer :: m
    REAL(kind=long)                 :: fsource
    !----------------------------------
    ! Declaration des variables locales
    !----------------------------------
    CHARACTER(len=6)                :: oldprf
    REAL(kind=long)                 :: pi
    REAL(kind=long) ::lea,leaa
    !-------------------
    ! Debut du programme
    !-------------------
    oldprf = prefix
    prefix = 'FSOUR'


    !------
    ! Corps
    !------
    pi=4.D0*ATAN(1.D0)
    Select case (m)
    case(1,7)
       fsource =   theta
    case (2)
       fsource  = theta *(x+y)
    case (3)
       fsource = theta * (x*x - y*y)
    case(4)
       fsource =   50. *pi**2 * cos(5.*pi*(x+y)) + theta*cos(5.*pi*(x+y))
    case(5)
       fsource = 2. * (x - x**2 + y - y**2) + theta*x*(1-x)*y*(1-y)
    case(6)
       fsource = (1+delta) * (pi**2) *sin(pi * x)*sin(pi *y) + theta * sin(pi * x)*sin(pi *y)
    case(8) ! ici V=1, et u=x+y
       fsource=2+ theta *(x+y)
    end Select
    !----------------
    ! Fin du programme
    !-----------------
    prefix = oldprf

    RETURN
  END FUNCTION fsource


END Module fsourcemod
