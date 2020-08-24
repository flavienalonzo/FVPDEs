subroutine MethNewtonMul(Tab_A,Tab_U,Tab_U_p,Tab_entity,A_big,Tab_equa,Tab_chemo,equa,tol,tps)
    USE longr
    USE parmmage
    USE imprime
    USE intmatvec
    USE algebrelineaire
    USE intbigradc
    use plotvtkmod
    use fsourcemod

    implicit none
    type(MatCreux), dimension(n_enty,n_enty), intent(inout) :: Tab_A
    real(kind=long), dimension(n_enty,Nbt), intent(inout) :: Tab_U, Tab_U_p
    character(len=6), dimension(n_enty), intent(in) :: Tab_entity, Tab_equa
    character(len=6), intent(in) :: equa
    type(MatCreux), intent(in) :: A_big
    real(kind=long), intent(in) :: tol,tps
    integer, dimension(n_enty), intent(in) :: Tab_chemo

    real(kind=long), dimension(n_enty*Nbt) :: B,dU,I0
    real(kind=long), dimension(n_enty) :: seuil
    integer :: i, k, h, n, n_max

    n_max = 15
    Tab_U = Tab_U_p
    I0 = 1.D0
    do k=1,n_enty
        do h=1,n_enty
            call vide(Tab_A(k,h))
        end do
        if (Tab_equa(k)=='instat') then 
            Tab_A(k,k)%Bg = Tab_U(k,:)-Tab_U_p(k,:)
            do i=1,Nbt
                call AJOUT(i,i,1.D0,Tab_A(k,k))
            end do 
        end if
        call assemblediffusionKS(Tab_A(k,:),Tab_U,Tab_entity(k),Tab_entity,Tab_equa,Tab_equa(k),k) 
        if (Tab_chemo(k)/=0) then
            call assemblechemoKS(Tab_A(k,:),Tab_U,Tab_entity(k),Tab_entity,Tab_equa,Tab_equa(k),k,Tab_chemo(k))
        end if
        call assemblereactionKS(Tab_A(k,:),Tab_U,Tab_entity(k),Tab_entity,Tab_equa,Tab_equa(k),k)
    end do
    !call actubigmatrix(Tab_A,A_big)
    !du = bigradient(A_big,-A_big%Bg,I0,tol)
    do k=1,n_enty
        B(1+(k-1)*Nbt:k*Nbt) = -Tab_A(k,k)%Bg  
    end do
        dU = BiCGSTAB(Tab_A,B,tol)!BGbloc(Tab_A,B,Tab_U,tol,tps)
    do k=1,n_enty
        Tab_U(k,:) = Tab_U(k,:) + dU(1+(k-1)*Nbt:k*Nbt) 
        seuil(k) = sqrt(dot_product(dU,dU))
    end do
    n=1
    do while(maxval(seuil)>tol.and.n<=n_max)
        do k=1,n_enty
            do h=1,n_enty
                call vide(Tab_A(k,h))
            end do
                if (Tab_equa(k)=='instat') then
                    Tab_A(k,k)%Bg = Tab_U(k,:)-Tab_U_p(k,:)
                    do i=1,Nbt
                        call AJOUT(i,i,1.D0,Tab_A(k,k))
                    end do 
                end if
                call assemblediffusionKS(Tab_A(k,:),Tab_U,Tab_entity(k),Tab_entity,Tab_equa,Tab_equa(k),k) 
                if (Tab_chemo(k)/=0) then
                    call assemblechemoKS(Tab_A(k,:),Tab_U,Tab_entity(k),Tab_entity,Tab_equa,Tab_equa(k),k,Tab_chemo(k))
                end if
                call assemblereactionKS(Tab_A(k,:),Tab_U,Tab_entity(k),Tab_entity,Tab_equa,Tab_equa(k),k)
        end do
        !call actubigmatrix(Tab_A,A_big)
        !dU = bigradient(A_big,-A_big%Bg,I0,tol)
        do k=1,n_enty
            B(1+(k-1)*Nbt:k*Nbt) = -Tab_A(k,k)%Bg  
        end do
        du = BiCGSTAB(Tab_A,B,tol)!BGbloc(Tab_A,B,Tab_U,tol,tps)
        do k=1,n_enty
            Tab_U(k,:) = Tab_U(k,:) + dU(1+(k-1)*Nbt:k*Nbt) 
            seuil(k) = sqrt(dot_product(dU(1+(k-1)*Nbt:k*Nbt) ,dU(1+(k-1)*Nbt:k*Nbt) ))
        end do
        n=n+1
        !print*,n,maxval(seuil)
    end do
    
end subroutine MethNewtonMul