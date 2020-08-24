subroutine vide( Mat )

    USE longr
    USE imprime
    USE parmmage
    IMPLICIT NONE

    !--------------------------
    ! Declaration des arguments
    !--------------------------
    TYPE(MatCreux)       :: Mat 

    integer :: ii, ij 
    double precision :: coef 

    do ii = 1, size(Mat%IndPL)-1
        do ij = Mat%IndPL(ii),Mat%IndPL( ii + 1 ) - 1
            !coef = Mat%TMat(ij)
            !call ajout(ii,Mat%Indc(ij),-coef,Mat)
            Mat%TMat(ij)=0
        end do 
    end do

    Mat%F = 0.D0 
    Mat%Bg = 0.D0

end subroutine vide 