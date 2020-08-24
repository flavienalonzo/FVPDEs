subroutine actuNutri(N,C,Tab_U)

    USE longr 
    USE imprime
    USE parmmage
    use intbigradc
    implicit none

    type(MatCreux), intent(in) :: N 
    double precision, dimension(n_enty,Nbt), intent(in) :: Tab_U
    double precision, dimension(Nbt), intent(out) :: C 
    double precision, dimension(Nbt) :: N0

    N0=0.D0


    N%Bg = 0.5*(Tab_U(index_endo,:)*AireK*Coef_prod - Coef_cons*AireK*Tab_U(index_norm,:) &
    & + abs(Tab_U(index_endo,:)*AireK*Coef_prod - Coef_cons*AireK*Tab_U(index_norm,:) ))
    C = bigradient(N, N%Bg,N0,5.d-10)



end subroutine actuNutri 