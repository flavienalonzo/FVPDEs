subroutine actubigmatrix(Tab_A,A_big)
    USE longr 
    USE imprime
    USE parmmage
    use intbigradc
    implicit none
    type(MatCreux), dimension(n_enty,n_enty), INTENT(IN) :: Tab_A
    type(MatCreux), INTENT(INOUT) :: A_big
    integer :: k, h, ll, cc

    call vide(A_big)
    do k=1,n_enty
        do h=1,n_enty
            do ll=1,Nbt
                do cc=Tab_A(k,h)%IndPL(ll),Tab_A(k,h)%IndPL(ll+1)-1
                    call AJOUT(ll+(k-1)*Nbt,Tab_A(k,h)%Indc(cc)+(h-1)*Nbt,Tab_A(k,h)%TMat(cc),A_big)
                end do
            end do
        end do
        A_big%Bg(1+(k-1)*Nbt:k*Nbt) = Tab_A(k,k)%Bg
    end do
    


end subroutine actubigmatrix