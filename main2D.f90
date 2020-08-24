program main2D
    use load_module
    use plotvtkmod
    use parameters 
    implicit none 

    !!!Variables
    integer, dimension(:,:), allocatable :: ConnectivityList,Eall,Typesegment,TN 
    double precision, dimension(:,:), allocatable :: Points,Centre_tri,Dist_tri,Coef_trans,Volume,Long_segment
    integer :: k,n , Ncoef 
    double precision, dimension(:), allocatable :: b, U,U0,U_p,Nutriment,Fonc,dU
    type(matsparse) :: A, GradFonc

    !!!Constantes
    character(len=100),parameter :: Filename='Tri_2D_0.05_5_1_2_@(x)[-x(2),x(1)].txt'
    integer, parameter :: uniform=1,condition=2,lambda=5,equation=2
    double precision, parameter :: h=0.05,tol=5.d-10,Delta_t=0.01,T=0.1
    double precision :: seuil
    !!!Initialisation 

    !Chargement du maillage
    call load(Filename,Points,Centre_tri,Dist_tri,Coef_trans,Volume,Long_segment,ConnectivityList,Eall,Typesegment,TN)
    n=size(Centre_tri,1)
    Ncoef = count(TN>0)+n 
    allocate(U(1:n), U0(1:n),b(1:n),Nutriment(1:n),Fonc(1:n))
    
    !Calcul des matrices A et b pour le problème du laplacien
    call matrice_tri(equation,Points,ConnectivityList,&
&        Eall,Typesegment,TN,Centre_tri,Dist_tri,Coef_trans,Volume,Long_segment,A,b)

    !Résolution du système Au=b par la méthode des bigradients conjugués
    
    Nutriment = bigradient(A,b,b,tol)
    !print*,Nutriment(1:20)
    call plot_vtk (Nutriment,'Nutriment','Nutriment',Points,ConnectivityList,Centre_tri)

    deallocate(A%intI,A%intJ,A%PosDiag,A%ValMat,A%Diag)

    allocate(GradFonc%intI(1:Ncoef),GradFonc%intJ(1:Ncoef),GradFonc%PosDiag(1:n),&
    & GradFonc%ValMat(1:Ncoef),GradFonc%Diag(1:n))

    U0=0.5*exp(-norm_vec(Centre_tri))
    print*,Centre_tri(1,:)
    print*,U0(1)
    k=1
    U_p =U0
    call plot_vtk (U0,'Uinitial','U0',Points,ConnectivityList,Centre_tri)

    do while(k*Delta_t<=T)
    
        call foncNewton(Points,ConnectivityList,Eall,Typesegment,Delta_t,&
        & TN,Centre_tri,Dist_tri,Coef_trans,Volume,Long_segment,Fonc,GradFonc,U_p,U0,Nutriment)

        dU = bigradient(GradFonc,Fonc,U_p,tol)

        seuil = sqrt(dot_product(dU,dU))

        U = U_p + dU 
        U_p=U

        do while(seuil>tol)
            call foncNewton(Points,ConnectivityList,Eall,Typesegment,Delta_t,&
        & TN,Centre_tri,Dist_tri,Coef_trans,Volume,Long_segment,Fonc,GradFonc,U_p,U0,Nutriment)
            dU = bigradient(GradFonc,Fonc,U_p,tol)
            seuil = sqrt(dot_product(dU,dU))
            U = U_p + dU 
            U_p=U
        end do 
        k = k+1 
        U0 = U 
        U_p=U 
        print*,Delta_t*k , sum(U)
    end do       

    call plot_vtk (U,'Ucalcule','U',Points,ConnectivityList,Centre_tri)

    deallocate(U,U0,b)

end program main2D 


!gfortran plotvtkmod.f90 load_module.f90 parameters.f90 -o main2D main2D.f90
!.\main2D.exe

!.\SOFTWARE\*.f90