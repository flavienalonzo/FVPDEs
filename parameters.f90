module parameters

    implicit none 

    type matsparse
    integer, dimension(:), allocatable :: intI, intJ, PosDiag
    double precision, dimension(:), allocatable :: ValMat, Diag
    end type matsparse

    

    contains

    function norm_vec(x) result(y)
        implicit none 
        double precision, dimension (:,:), intent(in) :: x
        double precision, dimension(1:size(x,1)) :: y 
        integer :: i,n 
        n=size(x,1)
        do i=1,n 
            y(i) = sqrt(sum(x(i,:)**2))
        end do 
    end function norm_vec

    function f_vol(x,equation) result(y)
        implicit none 
        double precision, dimension(2), intent(in):: x
        double precision :: y
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

    function g_bord(x,equation) result(y)
        implicit none 
        double precision, dimension(2), intent(in) :: x
        double precision :: y
        integer, intent(in) :: equation
        select case(equation)
        case(1)
            y = sum(x)
        case(2) 
            y = 0
        case default
            print*,'Problème d équation'
        end select
    end function g_bord

    FUNCTION mat_vec(A, X) result(Y)
        !     *--------------------------------------
        !     * Ce sous programme calcule le produit 
        !     * d'une matrice par un vecteur 
        !     * 
        !     * La matrice A est stockée sous forme sparse
        !     * 
    
        IMPLICIT NONE
    
        TYPE(matsparse), intent(in) :: A
        double precision, DIMENSION(size(A%Diag)),intent(in) :: X
        double precision, DIMENSION(size(A%Diag)) :: Y

        integer :: k
        
        !------
        ! Corps
        !------
        Y = 0.D0
        do k = 1, SIZE(A%intI)
            Y( A%intI(k) ) =  Y(A%intI(k))  + A%ValMat(k) * X(A%intJ(k))
        end do

    END FUNCTION mat_vec
    
    FUNCTION transposee_mat_vec(A, X) result(Y)
        !     *--------------------------------------
        !     * Ce sous programme calcule le produit 
        !     * de la transposee d'une matrice par un vecteur 
        !     * 
        !     * La matrice A est stockée sous forme sparse
        !     * 
    
        IMPLICIT NONE
    
        TYPE(matsparse), intent(in) :: A
        double precision, DIMENSION( SIZE(A%Diag)),intent(in) :: X
        double precision, DIMENSION( SIZE(A%Diag))            :: Y

        integer :: k
    
        
        !------
        ! Corps
        !------
        Y = 0.D0
        do k = 1, SIZE(A%intI)
            Y( A%intJ(k) ) =  Y(A%intJ(k))  + A%ValMat(k) * X(A%intI(k))
        end do
    
    END FUNCTION transposee_mat_vec

    function diffusion(x) result(y)
        double precision, intent(in) :: x 
        double precision :: y, D 
        D = 0.05
        y = D*x 
    end function diffusion 

    function diffprime(x) result(y)
        double precision, intent(in) :: x 
        double precision :: y, D 
        D = 0.05
        y = D 
    end function diffprime 

    function chemo(x) result(y) 
        double precision, intent(in) :: x 
        double precision :: y,chi 
        chi = 0.005
        y = chi*x 
    end function chemo 

    function chemoprime(x) result(y) 
        double precision, intent(in) :: x 
        double precision :: y,chi 
        chi = 0.005
        y = chi 
    end function chemoprime 

    function production(x) result(y)
        double precision , intent(in) :: x 
        double precision :: y ,r 
        r = 0.
        y = r*x 
    end function production 

    function prodprime(x) result(y)
        double precision , intent(in) :: x 
        double precision :: y ,r 
        r = 0.
        y = r
    end function prodprime


end module parameters 