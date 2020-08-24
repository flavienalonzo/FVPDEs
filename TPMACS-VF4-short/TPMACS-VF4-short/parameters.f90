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

    function num2string(num) result(str)
        implicit none
        integer, intent(in) :: num
        character(:), allocatable :: str
        character(range(num) + 2) :: tmp
        write (tmp,'(i0)') num
        str = trim(adjustl(tmp))
    end function num2string

    function diffusion(x,entity) result(y)
        real (kind = long), dimension(n_enty), intent(in) :: x
        character(len=6), intent(in) :: entity
        real (kind = long) :: y 
        select case(entity)
            case('Nutrim')
                y = Coef_diffusion
            case('Endoth')
                y = 0.
            case('u_norm')
                y = Diff_u*x(index_norm)*(1-x(index_norm))
            case('VasEGF')
                y = 0.
        end select
    end function diffusion 

    function diffprime(x,entity,derivative) result(y)
        real (kind = long), dimension(n_enty), intent(in) :: x
        character(len=6), intent(in) :: entity,derivative
        real (kind = long) :: y 
        select case(entity)
        case('Nutrim')
            select case(derivative)
            case default
                y = 0
            end select
        case('Endoth')
            select case(derivative)
            case default
                y = 0
            end select
        case('u_norm')
            select case(derivative)
            case('u_norm')
                y = Diff_u*(1-2*x(index_norm))
            case default
                y = 0
            end select
            
        case('VasEGF')
            select case(derivative)
                case default
                    y = 0
                end select
            end select
    end function diffprime 

    function chemo(x,entity) result(y) 
        real (kind = long), dimension(n_enty), intent(in) :: x
        character(len=6), intent(in) :: entity
        real (kind = long) :: y
        select case(entity)
        case('Nutrim')
            y = 0
        case('Endoth')
            y = 0.
        case('u_norm')
            y = chi*x(index_norm)*(1-x(index_norm))
        case('VasEGF')
            y = 0
        end select
    end function chemo 

    function chemoprime(x,entity,derivative) result(y) 
        real (kind = long), dimension(n_enty), intent(in) :: x
        character(len=6), intent(in) :: entity, derivative
        real (kind = long) :: y 
        select case(entity)
        case('Nutrim')
            select case(derivative)
            case default
                y = 0
            end select
        case('Endoth')
            select case(derivative)
            case default
                y = 0
            end select
        case('u_norm')
            select case(derivative)
            case('u_norm')
                y = chi*(1-2*x(index_norm))
            case default
                y = 0
            end select
            
        case('VasEGF')
            select case(derivative)
            case default
                y = 0
            end select
        end select
    end function chemoprime 

    function reaction(x,entity) result(y)
        real (kind = long), dimension(n_enty), intent(in) :: x
        character(len=6), intent(in) :: entity
        real (kind = long) :: y 
        select case(entity)
        case('Nutrim')
            y = 0
        case('Endoth')
            y = 0
        case('u_norm')
            if (x(index_nut)>seuil_hypo) then
                y = rat_pop*x(index_norm)
            else if(x(index_nut)<seuil_necro) then
                y = -x(index_norm)
            else 
                y = 0
            end if
        case('VasEGF')
            y = 0
        end select
    end function reaction 

    function reactionprime(x,entity,derivative) result(y)
        real (kind = long), dimension(n_enty), intent(in) :: x
        character(len=6), intent(in) :: entity, derivative
        real (kind = long) :: y 
        select case(entity)
        case('Nutrim')
            select case(derivative)
            case default 
                y = 0
            end select
        case('Endoth')
            select case(derivative)
            case default
                y = 0
            end select
        case('u_norm')
            select case(derivative)
            case('u_norm')
                if (x(index_nut)>seuil_hypo) then
                    y = rat_pop
                else if(x(index_nut)<seuil_necro) then
                    y = -1
                else 
                    y = 0
                end if
            case default 
                y = 0
            end select
            
        case('VasEGF')
            select case(derivative)
            case default 
                y = 0
            end select
        end select
    end function reactionprime 


end module parameters 