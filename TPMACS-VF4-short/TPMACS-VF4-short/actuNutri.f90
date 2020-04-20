subroutine actuNutri(N,U,C,E)

    USE longr 
    USE imprime
    USE parmmage
    use intbigradc
    implicit none

    type(MatCreux), intent(in) :: N 
    double precision, dimension(Nbt), intent(in) :: U
    double precision, dimension(Nbt), intent(out) :: C 
    double precision, dimension(Nbt) :: N0, E 

    N0=0.D0


    N%Bg = 0.5*(E*AireK*Coef_prod - Coef_cons*AireK*U + abs(E*AireK*Coef_prod - Coef_cons*AireK*U))
    C = bigradient(N, N%Bg,N0,5.d-10)



end subroutine actuNutri 