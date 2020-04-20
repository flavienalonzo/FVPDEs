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
    case(1:6,10)
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
    case(10) 
      if (0.45<=x.and.x<=0.55.and.0.65<=y.and.y<=0.75) then
         fsource = 0.7
      else if (0.45<=x.and.x<=0.55.and.0.15<=y.and.y<=0.25) then
         fsource = 0.7
      else if (0.7<=x.and.x<=0.8.and.0.45<=y.and.y<=0.55) then
         fsource = 0.7
      else if (0.2<=x.and.x<=0.3.and.0.45<=y.and.y<=0.55) then
         fsource = 0.7
      else
         fsource = 0
      end if
      !fsource = 5*exp(-10*norm2((/x,y/)-(/0.5,0.75/))) + 5*exp(-10*norm2((/x,y/)-(/0.5,0.25/))) &
      !& + 5*exp(-10*norm2((/x,y/)-(/0.25,0.5/))) + 5*exp(-10*norm2((/x,y/)-(/0.75,0.5/))) 
    end Select
    !----------------
    ! Fin du programme
    !-----------------
    prefix = oldprf

    RETURN
  END FUNCTION fsource

  function cond_ini(x,y) result(z)
   USE longr
      implicit none
      REAL(kind=long), INTENT(IN) :: x,y
      REAL(kind=long) :: z 
      if (0.45<=x.and.x<=0.55.and.0.45<=y.and.y<=0.55) then 
         z = 0.5 !z = 0.5*exp(-30*norm2((/x,y/)-(/0.5,0.5/))) 
      else 
         z = 0
      end if
  end function cond_ini



END Module fsourcemod
