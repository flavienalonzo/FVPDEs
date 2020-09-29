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
                y = Diff_endo
            case('u_norm')
                y = Diff_u
            case('VasEGF')
                y = VEGF_dif
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
            if (x(index_endo)>satur_endo) then 
                y = chemo_endo*x(index_endo)**2*(1-x(index_endo)-x(index_norm))**2
            else 
                y = 0
            end if
        case('u_norm')
            if (x(index_nut)<=seuil_hypo) then
                y = chi_u*x(index_norm)**2*(1-x(index_endo)-x(index_norm))**2
            else 
                y = 0
            end if
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
            case('Endoth')
                if (x(index_endo)>satur_endo) then 
                    y = chemo_endo*(4*x(index_endo)**3-6*(1-x(index_norm))*x(index_endo)**2&
                    & + x(index_endo)**2*(1-x(index_norm))**2)
                else
                    y = 0
                end if
            case('u_norm')
                if (x(index_endo)>satur_endo) then 
                    y = chemo_endo*(-2*x(index_endo)**3+2*x(index_endo)**2*(1-x(index_norm)))
                else
                    y = 0
                end if
            case default
                y = 0
            end select
        case('u_norm')
            select case(derivative)
            case('u_norm')
                if (x(index_nut)<=seuil_hypo) then
                    y = chi_u*(4*x(index_norm)**3-6*(1-x(index_endo))*x(index_norm)**2&
                    & + x(index_norm)**2*(1-x(index_endo))**2)
                else 
                    y = 0
                end if
            case('Endoth')
                if (x(index_nut)<=seuil_hypo) then
                    y = -chi_u*(-2*x(index_norm)**3+2*x(index_norm)**2*(1-x(index_endo)))
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
    end function chemoprime 

    function reaction(x,entity) result(y)
        real (kind = long), dimension(n_enty), intent(in) :: x
        character(len=6), intent(in) :: entity
        real (kind = long) :: y 
        select case(entity)
        case('Nutrim')
            y = rat_pop*x(index_endo)-Coef_cons*x(index_nut)*x(index_norm) - nut_degra*x(index_nut)
        case('Endoth')
                y = rate_endo*x(index_endo) - degr_endo*x(index_endo)
        case('u_norm')
            if (x(index_nut)>seuil_hypo) then
                y = rate*x(index_norm) - apop*x(index_norm)
            else 
                y = - apop*x(index_norm)
            end if
        case('VasEGF')
            if (x(index_nut)<=seuil_hypo.and.x(index_nut)>seuil_necro) then
                y = VEGF_prod*x(index_norm) - VEGF_degr*x(index_vegf)-VEGF_cons*x(index_vegf)*x(index_endo)
            else 
                y = 0
            end if
        end select
    end function reaction 

    function reactionprime(x,entity,derivative) result(y)
        real (kind = long), dimension(n_enty), intent(in) :: x
        character(len=6), intent(in) :: entity, derivative
        real (kind = long) :: y 
        select case(entity)
        case('Nutrim')
            select case(derivative)
            case('Nutrim')
                y = -(nut_degra + Coef_cons*x(index_norm))
            case('Endoth')
                y = rat_pop
            case('u_norm')
                y = -Coef_cons*x(index_nut)
            case default 
                y = 0
            end select
        case('Endoth')
            select case(derivative)
            case('Endoth')
                    y = rate_endo - degr_endo
            case default
                y = 0
            end select
        case('u_norm')
            select case(derivative)
            case('u_norm')
                if (x(index_nut)>seuil_hypo) then
                    y = rate - apop
                else 
                    y = -apop 
                end if
            case default 
                y = 0
            end select
            
        case('VasEGF')
            select case(derivative)
            case('VasEGF')
                if (x(index_nut)<=seuil_hypo.and.x(index_nut)>seuil_necro) then
                    y = -(VEGF_degr + VEGF_cons*x(index_endo))
                else 
                    y = 0
                end if
            case('u_norm')
                if (x(index_nut)<=seuil_hypo.and.x(index_nut)>seuil_necro) then
                    y = VEGF_prod
                else 
                    y = 0 
                end if
            case('Endoth')
                if (x(index_nut)<=seuil_hypo.and.x(index_nut)>seuil_necro) then
                    y = -VEGF_cons*x(index_vegf)
                else 
                    y = 0 
                end if
            case default 
                y = 0
            end select
        end select
    end function reactionprime 


end module parameters 