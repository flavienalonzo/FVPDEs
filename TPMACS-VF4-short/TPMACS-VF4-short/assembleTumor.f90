subroutine assembleTumor(Tab_A,Tab_U,Tab_U_p)

    use longr
    use imprime
    use parmmage
    USE intmatvec
    USE algebrelineaire
    USE intbigradc
    use plotvtkmod
    use fsourcemod

    implicit none
    type(MatCreux), intent(inout) :: Tab_A
    real(kind=long), dimension(1:2,1:Nbt), intent(inout) :: Tab_U, Tab_U_p
    character(len=6), dimension(2) :: Tab_entity
    real(kind=long), dimension(1:Nbt) :: B,dU,I0
    real(kind=long) :: seuil, tol
    integer :: i, k, h, n, n_max
    tol = 5.D-12
    Tab_entity = (/'u_norm','Nutrim'/)
    n_max = 15
    Tab_U = Tab_U_p
    I0 = 1.D0
    call vide(Tab_A)
    Tab_A%Bg = Tab_U(index_norm,1:Nbt)-Tab_U_p(index_norm,1:Nbt)
    do i=1,Nbt
        call AJOUT(i,i,1.D0,Tab_A)
    end do 
    call assemblediffusionKS(Tab_A,Tab_U,'u_norm',Tab_entity,'instat','instat',1) 
    call assemblechemoKS(Tab_A,Tab_U,'u_norm',Tab_entity,'instat','instat',1,index_nut)
    call assemblereactionKS(Tab_A,Tab_U,'u_norm',Tab_entity,'instat','instat',1)
    B = -Tab_A%Bg  
    dU = bigradient(Tab_A,B,I0,tol)
    Tab_U(index_norm,1:Nbt) = Tab_U(index_norm,1:Nbt) + dU
    seuil = sqrt(dot_product(dU,dU))
    n=1
    do while(seuil>tol.and.n<=n_max)
        call vide(Tab_A)
        Tab_A%Bg = Tab_U(index_norm,1:Nbt)-Tab_U_p(index_norm,1:Nbt)
        do i=1,Nbt
            call AJOUT(i,i,1.D0,Tab_A)
        end do 
        call assemblediffusionKS(Tab_A,Tab_U,'u_norm',Tab_entity,'instat','instat',1) 
        call assemblechemoKS(Tab_A,Tab_U,'u_norm',Tab_entity,'instat','instat',1,index_nut)
        call assemblereactionKS(Tab_A,Tab_U,'u_norm',Tab_entity,'instat','instat',1)
        B = -Tab_A%Bg  
        dU = bigradient(Tab_A,B,I0,tol)
        Tab_U(index_norm,1:Nbt) = Tab_U(index_norm,1:Nbt) + dU
        seuil = sqrt(dot_product(dU,dU))
        n=n+1
        !print*,n,maxval(seuil)
    end do

end subroutine assembleTumor