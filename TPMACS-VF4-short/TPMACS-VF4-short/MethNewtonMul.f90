subroutine MethNewtonMul(Tab_A,Tab_U,Tab_U_p,Tab_entity,tol)
    USE longr
    USE parmmage
    USE imprime
    USE intmatvec
    USE algebrelineaire
    USE intbigradc
    use plotvtkmod
    use fsourcemod

    implicit none
    type(MatCreux), dimension(n_enty), intent(inout) :: Tab_A
    real(kind=long), dimension(n_enty,Nbt), intent(inout) :: Tab_U, Tab_U_p
    character(len=6), dimension(n_enty), intent(in) :: Tab_entity
    real(kind=long), intent(in) :: tol

    real(kind=long), dimension(Nbt) :: dU,I0
    real(kind=long), dimension(n_enty) :: seuil
    integer :: i, k , index_endo, index_nut, index_norm

    Tab_U = Tab_U_p
    I0 = 1.D0
    index_endo=0
    index_norm=0
    index_nut=0
    do k=1,n_enty
        select case(Tab_entity(k))
        case('Endoth')
            index_endo = k
        case('Nutrim')
            index_nut = k
        case('u_norm')
            index_norm = k
        end select
    end do
    do k=1,n_enty
        if (k==index_nut) then 
            call actuNutri(Tab_A(k),Tab_U(index_nut,:),Tab_U(k,:),Tab_U(index_endo,:))
            seuil(k) = 0
        else
            call vide(Tab_A(k))
            Tab_A(k)%F = 0.D0 
            Tab_A(k)%Bg = Tab_U(k,:)-Tab_U_p(k,:)
            call assemblediffusionKS(Tab_A(k),Tab_U(k,:),Tab_U(index_endo,:),Tab_entity(k)) 
            call assemblechemoKS(Tab_A(k),Tab_U(k,:),Tab_U(index_nut,:),Tab_U(index_endo,:),Tab_entity(k))
            call assemblereactionKS(Tab_A(k),Tab_U(k,:),Tab_U(index_endo,:),Tab_entity(k))
            do i=1,Nbt
                call AJOUT(i,i,1.D0,Tab_A(k))
            end do 
            dU = bigradient(Tab_A(k),-Tab_A(k)%Bg,I0,tol) 
            Tab_U(k,:) = Tab_U(k,:) + dU 
            seuil(k) = sqrt(dot_product(dU,dU))
        end if
    end do
    do while(maxval(seuil)>tol)
        do k=1,n_enty
            if (k==index_nut) then
                call actuNutri(Tab_A(k),Tab_U(index_nut,:),Tab_U(k,:),Tab_U(index_endo,:))
                seuil(k) = 0
            else
                call vide(Tab_A(k))
                Tab_A(k)%F = 0.D0 
                Tab_A(k)%Bg = Tab_U(k,:)-Tab_U_p(k,:)
                call assemblediffusionKS(Tab_A(k),Tab_U(k,:),Tab_U(index_endo,:),Tab_entity(k)) 
                call assemblechemoKS(Tab_A(k),Tab_U(k,:),Tab_U(index_nut,:),Tab_U(index_endo,:),Tab_entity(k))
                call assemblereactionKS(Tab_A(k),Tab_U(k,:),Tab_U(index_endo,:),Tab_entity(k))
                do i=1,Nbt
                    call AJOUT(i,i,1.D0,Tab_A(k))
                end do 
                dU = bigradient(Tab_A(k),-Tab_A(k)%Bg,I0,tol) 
                Tab_U(k,:) = Tab_U(k,:) + dU 
                seuil(k) = sqrt(dot_product(dU,dU))
            end if
        end do
    end do
    






end subroutine MethNewtonMul