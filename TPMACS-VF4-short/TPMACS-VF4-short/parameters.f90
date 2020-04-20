module parameters
    use longr 
    use imprime 
    use parmmage 

    implicit none 
 

    contains

    function f_vol(x,equation) result(y)
        implicit none 
        real (kind = long), dimension(2), intent(in):: x
        real (kind = long) :: y
        integer, intent(in) :: equation
        select case(equation)
        case(1)
            y = 0
        case(2) 
            y = 5*exp(-norm2(x-(/0.,2.5/))) + 5*exp(-norm2(x-(/0.,-2.5/))) &
            & + 5*exp(-norm2(x-(/2.5,0./))) + 5*exp(-norm2(x-(/-2.5,0./))) 
        case default
            print*,'Problème d équation'
        end select
    end function f_vol


    function diffusion(x,e,entity) result(y)
        real (kind = long), intent(in) :: x,e
        character(len=6), intent(in) :: entity
        real (kind = long) :: y 
    select case(entity)
    case('Endoth')
        y = 0
    case('u_norm')
        y = Diff_u*( (1-e)**2*x**3/3 - (1-e)*x**4/2 + x**5/5)
    end select
    end function diffusion 

    function diffprime(x,e,entity) result(y)
        real (kind = long), intent(in) :: x,e
        character(len=6), intent(in) :: entity
        real (kind = long) :: y 
        select case(entity)
        case('Endoth')
            y = 0
        case('u_norm')
            y = Diff_u*x**2*(1-x-e)**2
        end select
    end function diffprime 

    function chemo(x,e,entity) result(y) 
        real (kind = long), intent(in) :: x,e
        character(len=6), intent(in) :: entity
        real (kind = long) :: y 
        select case(entity)
        case('Endoth')
            y = 0
        case('u_norm')
            y = chi*x*(1-x-e)**2
        end select
    end function chemo 

    function chemoprime(x,e,entity) result(y) 
        real (kind = long), intent(in) :: x,e
        character(len=6), intent(in) :: entity
        real (kind = long) :: y 
        select case(entity)
        case('Endoth')
            y = 0
        case('u_norm')
            y = chi*( (1-x-e)**2 - 2*x*(1-x-e) )
        end select
    end function chemoprime 

    function reaction(x,e,entity) result(y)
        real (kind = long), intent(in) :: x,e
        character(len=6), intent(in) :: entity
        real (kind = long) :: y 
        select case(entity)
        case('Endoth')
            y = 0
        case('u_norm')
            y = rate*x*(1-x-e)**2
        end select
    end function reaction 

    function reactionprime(x,e,entity) result(y)
        real (kind = long), intent(in) :: x,e
        character(len=6), intent(in) :: entity
        real (kind = long) :: y 
        select case(entity)
        case('Endoth')
            y = 0
        case('u_norm')
            y =rate*( (1-x-e)**2 - 2*x*(1-x-e) )
        end select
    end function reactionprime 


end module parameters 