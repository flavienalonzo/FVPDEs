subroutine bigmatrix(A,B)
    use longr
    use parmmage
    use imprime
    use intmatvec
    implicit none
    type(MatCreux), intent(out) :: A 
    type(MatCreux),  intent(in) :: B
    integer :: is, k , h , ll

    allocate(A%IndPL(n_enty*Nbt+1),A%Indc(n_enty**2*size(B%Indc)),&
    & A%TMat(n_enty**2*size(B%Indc)),A%F(n_enty*Nbt),A%Bg(n_enty*Nbt))
    A%IndPL(1) = 1
    do ll = 0,n_enty-1
        do is=1,Nbt
            k = B%IndPL(is+1)-B%IndPL(is)
            A%IndPL(ll*Nbt+is+1) = A%IndPL(ll*Nbt+is) + k*n_enty 
            do h=0,n_enty-1
                A%Indc(A%IndPL(ll*Nbt+is)+h*k:A%IndPL(ll*Nbt+is)+(h+1)*k) &
                & = h*Nbt+B%Indc(B%IndPL(is):B%IndPL(is+1)-1)
            end do
        end do
    end do

    A%TMat = 0D0
    A%F = 0D0
    A%Bg = 0D0

end subroutine bigmatrix