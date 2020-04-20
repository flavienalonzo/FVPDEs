subroutine MethNewton(A,U,U_p,Nutriment,Endothelial,entity,tol)
    USE longr
    USE parmmage
    USE imprime
    USE intmatvec
    USE algebrelineaire
    USE intbigradc
    use plotvtkmod
    Use fsourcemod

    implicit none
    real(kind=long), dimension(Nbt), INTENT(OUT) :: U
    real(kind=long), dimension(Nbt), INTENT(IN) :: U_p,Nutriment, Endothelial
    real(kind=long), INTENT(IN) :: tol
    type(MatCreux), INTENT(INOUT) :: A 
    character(len=6), intent(in) :: entity 

    real(kind=long), dimension(Nbt) :: dU,I0
    real(kind=long) :: seuil
    integer :: i

    !!!!! Corps
    U = U_p
    I0 = 1.D0
    call vide(A)
    A%F = 0.D0 
    A%Bg = U-U_p
    call assemblediffusionKS(A,U,Endothelial,entity) 
    call assemblechemoKS(A,U,Nutriment,Endothelial,entity)
    call assemblereactionKS(A,U,Endothelial,entity)
    do i=1,Nbt
        call AJOUT(i,i,1.D0,A)
    end do 
    dU = bigradient(A,-A%Bg,I0,tol) 
    U = U + dU 
    seuil = sqrt(dot_product(dU,dU))
    do while(seuil>tol)
        call vide(A)
        A%Bg = U-U_p
        call assemblediffusionKS(A,U,Endothelial,entity) 
        call assemblechemoKS(A,U,Nutriment,Endothelial,entity)
        call assemblereactionKS(A,U,Endothelial,entity)
        do i=1,Nbt
            call AJOUT(i,i,1.D0,A)
        end do
        dU = bigradient(A,-A%Bg,I0,tol) 
        U = U + dU 
        seuil = sqrt(dot_product(dU,dU))
        !print*,seuil
    end do
    !print*,A%TMat
    !print*,A%Bg
    !print*, seuil

end subroutine MethNewton