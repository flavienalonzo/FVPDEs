Module fsourcebreast
  implicit none 
Contains

  FUNCTION fsourceu(t,x,y, m)
    !--------
    ! Modules
    !--------
    USE longr

    IMPLICIT NONE

    !--------------------------
    ! Declaration des arguments
    !--------------------------
    REAL(kind=long), INTENT(in)     :: t, x , y
    Integer :: m
    REAL(kind=long)                 :: fsourceu, omega
    !----------------------------------
    ! Declaration des variables locales
    !----------------------------------
    CHARACTER(len=6)                :: oldprf

    !-------------------
    ! Debut du programme
    !-------------------
    oldprf = prefix
    prefix = 'FSOUR'

    !------
    ! Corps
    !------
    Select case (m)
    case(0) 
       fsourceu = 0.D0
    case (1)
       fsourceu = alpha *(x+y)
    case (2)
       fsourceu = alpha * (x*x - y*y)
    case(3)
       omega    = 20.D0*atan(1.)
       fsourceu =  2* omega**2 *cos(omega*(x+y))
    case(4)
       omega    = 20.D0*atan(1.)
       fsourceu = 2.D0*(omega*t +omega**2) *cos(omega*(t*t+x+y))
    case(5) 
       fsourceu =   1+beta*(t+x+y)
    case(66)!! Cas degenere avec a(u) =u(1-u)
       fsourceu = -2.D0*(1.D0-2.D0*(x+y))
    case(77)!! Cas degenere avec a(u) =u et A(u) = u**2/2, u =  x+y
       fsourceu = -2.D0
    case(99)
       fsourceu = 0.D0
    case default
       stop 'fsourceu'
    end Select
    !-----------------
    ! Fin du programme
    !-----------------
    prefix = oldprf

    RETURN
  END FUNCTION fsourceu
  !============================
  !============================
  FUNCTION fsourcev(t,x,y, m)
    !--------
    ! Modules
    !--------
    USE longr

    IMPLICIT NONE

    !--------------------------
    ! Declaration des arguments
    !--------------------------
    REAL(kind=long), INTENT(in)     :: t,x , y
    Integer :: m
    REAL(kind=long)                 :: fsourcev, omega
    !----------------------------------
    ! Declaration des variables locales
    !----------------------------------
    CHARACTER(len=6)                :: oldprf

    !-------------------
    ! Debut du programme
    !-------------------
    oldprf = prefix
    prefix = 'FSOUR'

    !------
    ! Corps
    !------
    Select case (m)
    case(0)
       fsourcev = 0.D0
    case (1)
       fsourcev  = alpha *(x+y)
    case (2)
       fsourcev = alpha * (x*x - y*y)
    case(3)
       omega = 20.D0*atan(1.)
       fsourcev =  2.D0* omega**2 *cos(omega*(x+y)) + (beta-alpha)*cos(omega*(x+y))
    case(4)
       omega = 20.D0*atan(1.)
       fsourcev = 2.D0*(omega*t +omega**2) *cos(omega*(t*t+x+y)) + (beta-alpha)*cos(omega*(t*t+x+y))
    case(5)
       fsourcev =   1+beta*(t+x+y)
    case(66,77,99)
       fsourcev =   0.D0
    case default
       stop 'fsourcev'
    end Select
    !-----------------
    ! Fin du programme
    !-----------------
    prefix = oldprf

    RETURN
  END FUNCTION fsourcev
  !======================================================
  !=====================================================
  FUNCTION gbord(t,x,y,m)
    !--------
    ! Modules
    !--------
    USE longr

    IMPLICIT NONE

    !--------------------------
    ! Declaration des arguments
    !--------------------------
    REAL(kind=long), INTENT(in)     :: t, x, y
    integer :: m
    REAL(kind=long)                 :: gbord, omega
    !----------------------------------
    ! Declaration des variables locales
    !----------------------------------
    CHARACTER(len=6)      :: oldprf

    !-------------------
    ! Debut du programme
    !-------------------
    oldprf = prefix
    prefix = 'Gbord'

    !------
    ! Corps
    !------
    select case(m)
    case(0)
       gbord = 1.D0
    case(1)
       gbord  = x+y
    case(2) 
       gbord  = x*x - y*y
    case(3)
       omega = 20.D0*atan(1.)
       gbord = cos(omega*(x+y))
    case(4)
       omega = 20.D0*atan(1.)
       gbord = cos(omega*(t*t+x+y))
    case(5)
       gbord = t+ x+y
    case(66)! cas degenere avec a(u) =u(1-u)
       gbord = x+y
    case(77)! cas degenere avec a(u) =u
       gbord = x+y
    case(99)
       !!       gbord = -1000000.
       gbord = 0.D0
    case default
       stop 'gbord'
    end select
    !-----------------
    ! Fin du programme
    !-----------------
    prefix = oldprf

    RETURN
  END FUNCTION gbord


  !=======================
  function Adegen(y)
    !=======================
    USE longr
    REAL(kind=long), INTENT(in)     :: y
    Real(kind=long) ::  Adegen


    select case ( ChoixAdeg )
    case(1)   !!test 1
       Adegen = CoefdiffuAdeg * y
    case(2)   !!test 2
       Adegen = CoefdiffuAdeg * y * y/2.D0 
    case(3)   !! test 3
       Adegen = CoefdiffuAdeg * (y * y/2.D0 - y * y * y/3.D0 )
    case(4)   !! test 4
       IF ( y <= ubar ) THEN
          Adegen=CoefdiffuAdeg* y * (1.D0 + ( (gamma-1.D0)/(gamma+1.D0) ) * ( (y/ubar)**gamma) )
       ELSE
          Adegen=0.D0
       END IF
    case default   
       !stop'Adegen'
    end select


  end function Adegen


  !=======================
  function DerivAdegen( y )
    !=======================
    USE longr
    REAL( kind=long) :: y, DerivAdegen


    select case ( ChoixAdeg )
    case(1) !!test 1
       DerivAdegen = CoefdiffuAdeg 
    case(2) !!test 2
       DerivAdegen = CoefdiffuAdeg * y
    case(3) !! test 3
       DerivAdegen = CoefdiffuAdeg * y *( 1.D0 - y )  
    case(4) !! test 4
       IF ( y < ubar ) THEN
          DerivAdegen = CoefdiffuAdeg * ( 1.D0 + (gamma-1.D0) * ( (y/ubar)**gamma ) )        
       ELSE 
          DerivAdegen = 0.D0
       END IF
    case default  
       !stop'DerivAdegen'
    end select

  end function DerivAdegen


  !=======================
  function DiffAdegen( y )
    !=======================
    USE longr
    REAL(kind=long), INTENT(in)     :: y
    REAL( kind=long)                :: DiffAdegen


    select case ( ChoixAdeg )
    case(1) !!test 1
       DiffAdegen = 0.D0
    case(2) !!test 2
       DiffAdegen = CoefdiffuAdeg
    case(3) !! test 3
       DiffAdegen = CoefdiffuAdeg  * ( 1.D0 -2.D0 * y )  
    case(4) !! test 4
       IF ( y < ubar ) THEN
          DiffAdegen = CoefdiffuAdeg *   (gamma-1.D0) * ( (y**(gamma-1.D0))/(ubar**gamma ))        
       ELSE 
          DiffAdegen = 0.D0
       END IF
    case default  
       !stop'DiffAdegen'
    end select

  end function DiffAdegen

  !==================
  FUNCTION AKL(x,y,z)
    !=================
    !--------
    ! Modules
    !--------
    USE longr

    IMPLICIT NONE

    !--------------------------
    ! Declaration des arguments
    !--------------------------
    REAL(kind=long), INTENT(in)  :: x, y, z    
    REAL(kind = long)            :: AKL

    !------
    ! Corps
    !------
    SELECT CASE(choixdegenere)
    case(0) !! Pb non degenere
       IF ( x >= 0.D0 ) THEN
          AKL = CoefdiffuAdeg*max( Eta(y), Eta(z) )
       ELSE
          AKL = CoefdiffuAdeg*min( Eta(y), Eta(z) )
       ENDIF
    case(1) !! pb degenere
       Select Case(choixAdeg)
       case(1)
          AKL = CoefdiffuAdeg
       case(2)
          IF ( x >= 0.D0 ) THEN
             if (y <= z) then
                AKL = DerivAdegen(z)
             else 
                AKL = DerivAdegen(y)
             endif
          ELSE
             if (y <= z) then
                AKL = DerivAdegen(y)
             else 
                AKL = DerivAdegen(z)
             endif
          ENDIF
       case(3)
          IF ( x >= 0.D0 ) THEN
             if ( ( y <= z .AND. z <= (1.D0/2.D0) ) .OR. ( y >= z .AND. z >= (1.D0/2.D0) )) then
                AKL = DerivAdegen(z)
             elseif ( ( z <= y .AND. y <= (1.D0/2.D0) ) .OR. ( z >= y .AND. y >= (1.D0/2.D0) )) then
                AKL = DerivAdegen(y)
             else
                AKL = DerivAdegen(1.D0/2.D0)
             endif
          ELSE
             if ( ( y <= z .AND. z <= (1.D0/2.D0) ) .OR. ( y >= z .AND. z >= (1.D0/2.D0) )) then
                AKL = DerivAdegen(y)
             elseif ( ( z <= y .AND. y <= (1.D0/2.D0) ) .OR. ( z >= y .AND. y >= (1.D0/2.D0) )) then
                AKL = DerivAdegen(z)
             else
                AKL = min(DerivAdegen(y),DerivAdegen(z))
             endif
          ENDIF
       END Select
    END SELECT

  END FUNCTION AKL
  
  !=======================
  FUNCTION DerivAKL(x,y,z)
    !=====================
    !--------
    ! Modules
    !--------
    USE longr

    IMPLICIT NONE
    !--------------------------
    ! Declaration des arguments
    !--------------------------
    REAL(kind=long), INTENT(in)      :: x, y, z
    REAL(kind=long), DIMENSION(2)    :: DerivAKL   

    !------
    ! Corps
    !------
    DerivAKL=0
    SELECT CASE(choixdegenere)
    case(0) !! Pb non degenere
       IF ( x >= 0.D0 ) THEN
          if ( z <= y .AND. y < 1.D0 .AND. y > 0.D0 ) then
             DerivAKL(1) = CoefdiffuAdeg
          elseif ( y <= z .AND. z < 1.D0 .AND. z > 0.D0 ) then
             DerivAKL(2) = CoefdiffuAdeg
          endif
       ELSE
          if ( y <= z .AND. y < 1.D0 .AND. y > 0.D0 ) then
             DerivAKL(1) = CoefdiffuAdeg
          elseif ( z <= y .AND. z < 1.D0 .AND. z > 0.D0 ) then
             DerivAKL(2) = CoefdiffuAdeg
          endif
       ENDIF
    case(1) !! pb degenere
       Select Case(choixAdeg)
       case(1)
          DerivAKL = 0.D0
       case(2)
          IF ( x >= 0.D0 ) THEN
             if (y <= z) then
                DerivAKL(2) = DiffAdegen(z)
             else 
                DerivAKL(1) = DiffAdegen(y)
             endif
          ELSE
             if (y <= z) then
                DerivAKL(1) = DiffAdegen(y)
             else 
                DerivAKL(2) = DiffAdegen(z)
             endif
          ENDIF
       case(3)
          IF ( x >= 0.D0 ) THEN
             if ( ( y <= z .AND. z <= (1.D0/2.D0) ) .OR. ( y >= z .AND. z >= (1.D0/2.D0) )) then
                DerivAKL(2) = DiffAdegen(z)
             elseif ( ( z <= y .AND. y <= (1.D0/2.D0) ) .OR. ( z >= y .AND. y >= (1.D0/2.D0) )) then
                DerivAKL(1) = DiffAdegen(y)
             endif
          ELSE
             if ( ( y <= z .AND. z <= (1.D0/2.D0) ) .OR. ( y >= z .AND. z >= (1.D0/2.D0) )) then
                DerivAKL(1) = DiffAdegen(y)
             elseif ( ( z <= y .AND. y <= (1.D0/2.D0) ) .OR. ( z >= y .AND. y >= (1.D0/2.D0) )) then
                DerivAKL(2) = DiffAdegen(z)
             else
                if (DerivAdegen(y) > DerivAdegen(z) ) then
                   DerivAKL(2) = DiffAdegen(z)
                elseif (DerivAdegen(y) < DerivAdegen(z) ) then
                   DerivAKL(1) = DiffAdegen(y)
                else
                   DerivAKL(1) = DiffAdegen(y)
                   DerivAKL(2) = DiffAdegen(z)
                end if
             endif
          ENDIF
       END Select
    END SELECT

  END FUNCTION DerivAKL

  !==============
  function Eta(y)
    !============
    USE longr
    REAL( kind=long) :: y, Eta
    
    if ( y < 0 ) then
       Eta = 0
    elseif ( y < 1) then
       Eta = y
    else
       Eta = 1
    endif
    
  end function Eta
  

  !=============
  function p(y)
    !===========
    USE longr
    REAL( kind=long) :: y, p

    SELECT CASE(choixdegenere)
    case(0) !! Pb non degenere
       if ( y <= 0.D0 ) then
          p = y/epsilon
       elseif ( y <= 1.D0 ) then
          p = log(1.D0+y/epsilon)
       else
          p = log(1.D0+1.D0/epsilon) + (1/(1.D0+epsilon))*(y-1)
       endif
    case(1) !! pb degenere
       p = y
    END SELECT

     end function p
  
  !===================
  function Derivp( y )
    !=================
    USE longr
    REAL( kind=long) :: y, Derivp
    
    SELECT CASE(choixdegenere)
    case(0) !! Pb non degenere
       Derivp = 1.D0/(Eta(y)+epsilon)
    case(1) !! pb degenere
       Derivp = 1.D0
    END SELECT
  end function Derivp
  
  
  !====================
  FUNCTION EtaKL(x,y,z)
    !==================
    
    !--------
    ! Modules
    !--------
    USE longr

    IMPLICIT NONE

    !--------------------------
    ! Declaration des arguments
    !--------------------------
    REAL(kind=long), INTENT(in)  :: x, y, z    
    REAL(kind = long)            :: EtaKL

    !------
    ! Corps
    !------
    
    IF ( x >= 0.D0 ) THEN
       EtaKL = max( Eta(y), Eta(z) )
    ELSE
       EtaKL = min( Eta(y), Eta(z) )
    ENDIF
        
  END FUNCTION EtaKL
  
  !=======================
  FUNCTION DerivEtaKL(x,y,z)
    !=====================
    !--------
    ! Modules
    !--------
    USE longr
    
    IMPLICIT NONE
    !--------------------------
    ! Declaration des arguments
    !--------------------------
    REAL(kind=long), INTENT(in)      :: x, y, z
    REAL(kind=long), DIMENSION(2)    :: DerivEtaKL   
    
    !------
    ! Corps
    !------
    
    DerivEtaKL=0
    IF ( x >= 0.D0 ) THEN
       if ( z <= y .AND. y < 1.D0 .AND. y > 0.D0 ) then
          DerivEtaKL(1) = 1.D0
       elseif ( y <= z .AND. z < 1.D0 .AND. z > 0.D0 ) then
          DerivEtaKL(2) = 1.D0
       endif
    ELSE
       if ( y <= z .AND. y < 1.D0 .AND. y > 0.D0 ) then
          DerivEtaKL(1) = 1.D0
       elseif ( z <= y .AND. z < 1.D0 .AND. z > 0.D0 ) then
          DerivEtaKL(2) = 1.D0
       endif
    ENDIF
    
  END FUNCTION DerivEtaKL

  !=======================
  function chi(y)
    !=======================
    USE longr
    REAL( kind=long) :: y, chi


    select case(ChoixChi)
    case(1)
       chi = CoefTranspChi 
    case(2)
       chi = CoefTranspChi *y
    case(3)
       chi= CoefTranspChi *y*(1.D0-y) 
    case(4)
       chi= CoefTranspChi *y*(1.D0-y)**2
    case(5)
       IF ( y <= ubar ) THEN
          chi= CoefTranspChi *y*( 1.D0-(y/ubar)**gamma)
       ELSE
          chi= 0.D0
       END IF
       case(6)
          chi= CoefTranspChi *( y*(1.D0-y) )**2
    case default
       !stop'pb chi'
    end select

  end function chi



  !=======================
  function ChiCroit(y)
    !=======================
    USE longr
    implicit none
    REAL( kind=long) :: y, z, ChiCroit, umax
    !integer :: choix



    select case(ChoixChi)
    case(1)
       ChiCroit = CoefTranspChi 
    case(2)
       ChiCroit = CoefTranspChi *y
    case(3)
       umax= 0.5D0
       z=min(y,umax)
       ChiCroit=  CoefTranspChi * z*(1.D0-z)
    case(4)
       umax= 1.D0/3.D0
       z=min(y,umax)
       ChiCroit= CoefTranspChi * z*(1.D0-z)**2
    case(5)
       umax= (ubar)/( (gamma+1.D0)**(1.D0/gamma) )
       z=min(y,umax)
       IF ( z <= ubar ) THEN
          ChiCroit= CoefTranspChi * z*(1.D0-(z/ubar)**gamma)
       ELSE
          ChiCroit= 0.D0
       END IF
    case(6)
       umax= 0.5D0
       z=min(y,umax)
       ChiCroit=  CoefTranspChi * (z*(1.D0-z))**2
    case default
       !stop'pb ChiCroit'
    end select

  end function ChiCroit

  !=======================
  function DerivChiCroit(y)
    !=======================
    USE longr
    REAL( kind=long) :: y, DerivChiCroit, umax
    !integer :: choix
    !

    select case(ChoixChi)
    case(1)
       DerivchiCroit = 0.D0 
    case(2)
       DerivchiCroit = CoefTranspChi
    case(3)
       umax= 0.5D0
       if (y <=umax) then
          DerivchiCroit =  CoefTranspChi * (1-2.D0*y)
       else
          DerivchiCroit = 0.D0
       end if
    case(4)
       umax= 1.D0/3.D0
       if (y < umax) then
          DerivchiCroit = CoefTranspChi * ( (1.D0-y)**2 -2.D0*y*(1.D0-y) )
       else
          DerivchiCroit = 0.D0
       end if
    case(5)
       umax= (ubar)/( (gamma+1.D0)**(1.D0/gamma) )
       IF (y < umax .AND. y< ubar) THEN !! Although we know that y<=ubar
          DerivchiCroit = CoefTranspChi * (1.D0 - ( gamma+1.D0 )*( (y/ubar)**gamma) )
       ELSE
          DerivchiCroit = 0.D0
       END IF
    case(6)
       umax= 0.5D0
       if (y <=umax) then
          DerivchiCroit =  CoefTranspChi * 2.D0*y*(1.D0-y)*(1-2.D0*y)
       else
          DerivchiCroit = 0.D0
       end if
    case default
       !stop'pb DerivChiCroit'
    end select

  end function DerivChiCroit

  !=======================
  function ChiDeCroit(y)
    !=======================
    USE longr
    REAL( kind=long) :: y, z, ChiDecroit, umax
    !

    select case(ChoixChi)
    case(1)
       ChiDeCroit = 0.D0 
    case(2)
       ChiDecroit = 0.D0
    case(3)
       umax= 0.5D0
       z=max(y,umax)
       ChiDeCroit=  CoefTranspChi * ( z*(1.D0-z) - umax*(1.D0-umax) )
    case(4)
       umax= 1.D0/3.D0
       z=max(y,umax)
       ChiDecroit= CoefTranspChi * ( z*(1.D0-z)**2 -umax*(1.D0-umax)**2 ) 
    case(5)
       umax= (ubar)/( (gamma+1.D0)**(1.D0/gamma) )
       z=max(y,umax)
       IF (z <= ubar ) THEN
          ChiDeCroit= CoefTranspChi * ( z*(1.D0-(z/ubar)**gamma) )
       ELSE
          ChiDeCroit= 0.D0
       END IF
       IF ( umax <= ubar )  THEN
          ChiDeCroit= ChiDeCroit  - CoefTranspChi*umax*( 1.D0-(umax/ubar)**gamma ) 
       END IF
    case(6)
       umax= 0.5D0
       z=max(y,umax)
       ChiDeCroit=  CoefTranspChi * ( (z*(1.D0-z))**2 - (umax*(1.D0-umax))**2 )
    case default
      ! stop'pb ChiCroit'
    end select

  end function ChiDeCroit

  !=======================
  function DerivChiDecroit(y)
    !=======================
    USE longr
    REAL( kind=long) :: y, DerivChiDecroit, umax
    !

    select case(ChoixChi)
    case(1)
       DerivChiDecroit = 0.D0 
    case(2)
       DerivChiDecroit = 0.D0
    case(3)
       umax= 0.5D0
       if (y > umax) then
          DerivChiDecroit =  CoefTranspChi * (1.D0-2.D0*y)
       else
          DerivChiDecroit = 0.D0
       end if
    case(4)
       umax= 1.D0/3.D0
       if (y > umax) then
          DerivChiDecroit = CoefTranspChi * ( (1.D0-y)**2 -2.D0*y*(1.D0-y) )
       else
          DerivChiDecroit = 0.D0
       end if
    case(5)
       umax= (ubar)/( (gamma+1.D0)**(1.D0/gamma) )
       IF (y > umax .AND. y < ubar) THEN
          DerivChiDeCroit = CoefTranspChi* ( 1.D0-(gamma+1.D0)*( (y/ubar)**gamma ) )
       ELSE
          DerivChiDeCroit = 0.D0
       END IF
    case(6)
       umax= 0.5D0
       if (y > umax) then
          DerivChiDecroit =  CoefTranspChi * 2.D0*y*(1.D0-y)*(1.D0-2.D0*y)
       else
          DerivChiDecroit = 0.D0
       end if
    case default
       !stop'pb DerivChiDecroit'
    end select

  end function DerivChiDecroit


  !=======================
  function Derivchi(y)
    !=======================
    USE longr
    REAL( kind=long) :: y,Derivchi


    select case(ChoixChi)
    case(1)
       Derivchi = 0.D0
    case(2)
       Derivchi = CoefTranspChi *1.D0
    case(3)
       Derivchi = CoefTranspChi * ( 1.D0 - 2.D0*y )
    case(4)
       Derivchi = CoefTranspChi *( 1.D0 - y ) * ( 1.D0 - 3.D0*y )
    case(5)
       IF ( y < ubar ) THEN
          Derivchi= CoefTranspChi *( 1.D0-( gamma+ 1.D0 )*( (y/ubar)**gamma ) )
       ELSE
          Derivchi= 0.D0
       END IF
    case(6)
       Derivchi = CoefTranspChi * 2.D0 *y*(1.D0-y)*( 1.D0 - 2.D0*y )
    case default
      ! stop'pb derivchi'
    end select


  end function Derivchi


  !=======================
  function muKL(y)
    !=======================
    USE longr
    REAL( kind=long) :: y, muKL


    select case(ChoixChi)
    case(4)
       muKL= CoefTranspChi *sqrt(y)*(1.D0-y)/CoefdiffuAdeg
   !case(5)
   !   IF ( y <= ubar ) THEN
   !      muKL= CoefTranspChi *y*( 1.D0-(y/ubar)**gamma)
   !   ELSE
   !      muKL= 0.D0
   !   END IF
    case(6)
       muKL= CoefTranspChi *y*(1.D0-y)/CoefdiffuAdeg
    case default
       !stop'pb muKL'
    end select
    
  end function muKL

  !=======================
  function MuKLCroit(y)
    !=======================
    USE longr
    implicit none
    REAL( kind=long) :: y, z, muKLCroit, umax
    !integer :: choix
    
    
    
    select case(ChoixChi)
    case(4)
       umax= 1.D0/3.D0
       z=min(y,umax)
       muKLCroit= CoefTranspChi * z*(1.D0-z)/CoefdiffuAdeg
   !case(5)
   !   umax= (ubar)/( (gamma+1.D0)**(1.D0/gamma) )
   !   z=min(y,umax)
   !   IF ( z <= ubar ) THEN
   !      muKLCroit= CoefTranspChi * z*(1.D0-(z/ubar)**gamma)
   !   ELSE
   !      muKLCroit= 0.D0
   !   END IF
    case(6)
       umax= 0.5D0
       z=min(y,umax)
       muKLCroit=  CoefTranspChi *z*(1.D0-z)/CoefdiffuAdeg
    case default
       !stop'pb muKLCroit'
    end select

  end function MuKLCroit

 !=======================
  function DerivMuKLCroit(y)
    !=======================
    USE longr
    REAL( kind=long) :: y, DerivMuKLCroit, umax
   ! integer :: choix
    !

    select case(ChoixChi)
    case(4)
       umax= 1.D0/3.D0
       if (y < umax) then
          DerivMuKLCroit = CoefTranspChi * ( (1.D0-y)**2 -2.D0*y*(1.D0-y) )/CoefdiffuAdeg
       else
          DerivMuKLCroit = 0.D0
       end if
    case(5)
       umax= (ubar)/( (gamma+1.D0)**(1.D0/gamma) )
       IF (y < umax .AND. y< ubar) THEN !! Although we know that y<=ubar
          DerivMuKLCroit = CoefTranspChi * (1.D0 - ( gamma+1.D0 )*( (y/ubar)**gamma) )/CoefdiffuAdeg
       ELSE
          DerivMuKLCroit = 0.D0
       END IF
    case(6)
       umax= 0.5D0
       if (y <=umax) then
          DerivMuKLCroit =  CoefTranspChi * (1-2.D0*y)/CoefdiffuAdeg
       else
          DerivMuKLCroit = 0.D0
       end if
    case default
      ! stop'pb DerivChiCroit'
    end select

  end function DerivMuKLCroit


 !=======================
  function MuKLDeCroit(y)
    !=======================
    USE longr
    REAL( kind=long)    :: y, z, MuKLDeCroit, umax
    !

    select case(ChoixChi)
    case(4)
       umax= 1.D0/3.D0
       z=max(y,umax)
       MuKLDeCroit= CoefTranspChi * ( z*(1.D0-z)**2 -umax*(1.D0-umax)**2 ) /CoefdiffuAdeg
    case(5)
       umax= (ubar)/( (gamma+1.D0)**(1.D0/gamma) )
       z=max(y,umax)
       IF (z <= ubar ) THEN
          MuKLDeCroit= CoefTranspChi * ( z*(1.D0-(z/ubar)**gamma) )/CoefdiffuAdeg
       ELSE
          MuKLDeCroit= 0.D0
       END IF
       IF ( umax <= ubar )  THEN
          MuKLDeCroit= MuKLDeCroit/CoefdiffuAdeg  - CoefTranspChi*umax*( 1.D0-(umax/ubar)**gamma )/CoefdiffuAdeg 
       END IF
    case(6)
       umax= 0.5D0
       z=max(y,umax)
       MuKLDeCroit=  CoefTranspChi * ( z*(1.D0-z) - umax*(1.D0-umax) )/CoefdiffuAdeg
    case default
       !stop'pb ChiCroit'
    end select

  end function MuKLDeCroit

  !=======================
  function DerivMuKLDecroit(y)
    !=======================
    USE longr
    REAL( kind=long) :: y, DerivMuKLDecroit, umax
    !

    select case(ChoixChi)
    case(4)
       umax= 1.D0/3.D0
       if (y > umax) then
          DerivMuKLDecroit = CoefTranspChi * ( (1.D0-y)**2 -2.D0*y*(1.D0-y) )/CoefdiffuAdeg
       else
          DerivMuKLDecroit = 0.D0
       end if
    case(5)
       umax= (ubar)/( (gamma+1.D0)**(1.D0/gamma) )
       IF (y > umax .AND. y < ubar) THEN
          DerivMuKLDeCroit = CoefTranspChi* ( 1.D0-(gamma+1.D0)*( (y/ubar)**gamma ) )/CoefdiffuAdeg
       ELSE
          DerivMuKLDeCroit = 0.D0
       END IF
    case(6)
       umax= 0.5D0
       if (y > umax) then
          DerivMuKLDecroit =  CoefTranspChi * (1.D0-2.D0*y)/CoefdiffuAdeg
       else
          DerivMuKLDecroit = 0.D0
       end if
    case default
       !stop'pb DerivMuKLDecroit'
    end select

  end function DerivMuKLDecroit

 
  function chemotherapy(y)
   use longr 
   implicit none 
   real(kind=long), intent(in) :: y 
   real(kind=long) :: chemotherapy 
   integer :: i 
   chemotherapy = 0.D0
   do i=1,Nchemo
      if (y>=chemo_time(1,i).and.chemo_time(2,i)> y) then
         chemotherapy = dose_chemo
      end if
   end do
  end function chemotherapy 


  function radiotherapy(y)
   use longr
   implicit none 
   real(kind=long), intent(in) :: y 
   real(kind=long) :: radiotherapy 
   integer :: i
   radiotherapy = 0.D0
   do i=1,Nradio
      if (y>=radio_time(1,i).and.radio_time(2,i)> y) then
         radiotherapy = 1-EXP(-(radio_alpha*dose_radio+radio_beta*dose_radio**2))
      end if
   end do
end function radiotherapy

END Module fsourcebreast
