module load_module
    use parameters
    implicit none 
    contains

    subroutine load(Filename,Points,Centre_tri,Dist_tri,Coef_trans,Volume,Long_segment,ConnectivityList,Eall,Typesegment,TN)
        implicit none 

        character(len=100),intent(in) :: Filename
        integer :: n_row,n_col,k
        integer, dimension(:,:), allocatable :: B
        integer, dimension(:,:), intent(out), allocatable :: ConnectivityList,Eall,Typesegment,TN 
        double precision, dimension(:,:), allocatable :: A
        double precision, dimension(:,:), intent(out), allocatable :: Points,Centre_tri,Dist_tri,Coef_trans,Volume,Long_segment
        !Filename = 'Tri_2D_0.05_5_1_2_@(x)[-x(2),x(1)].txt'
        open(unit=10,file = Filename)

        !!! Lecture des coordonnées des sommets des triangles
        read(10,*) n_row,n_col
        allocate(A(1:n_row,1:n_col)) 
        allocate(Points(1:n_row,1:n_col))
        do k=1,n_row,1
            read(10,*) A(k,:)
        end do 
        Points = A 
        deallocate(A)
        !print*,Points(n_row,:)

        !!! Lecture de la liste des connectivités des triangles
        read(10,*) n_row,n_col
        allocate(B(1:n_row,1:n_col))
        allocate(ConnectivityList(1:n_row,1:n_col))
        do k=1,n_row,1
            read(10,*) B(k,:)
        end do 
        ConnectivityList = B  
        deallocate(B) 
        !print*,ConnectivityList(n_row,:)

        !!! Lecture de la liste des segments des triangles
        read(10,*) n_row,n_col
        allocate(B(1:n_row,1:n_col))
        allocate(Eall(1:n_row,1:n_col))
        do k=1,n_row,1
            read(10,*) B(k,:)
        end do 
        Eall = B  
        deallocate(B) 
        !print*,Eall(n_row,:)

        !!! Lecture de la liste donnant le type des segments
        read(10,*) n_row,n_col
        allocate(B(1:n_row,1:n_col))
        allocate(Typesegment(1:n_row,1:n_col))
        do k=1,n_row,1
            read(10,*) B(k,:)
        end do 
        Typesegment = B  
        deallocate(B) 
        !print*,Typesegment(n_row,:)

        !!! Lecture de la liste des triangles voisins d'un triangle
        read(10,*) n_row,n_col
        allocate(B(1:n_row,1:n_col))
        allocate(TN(1:n_row,1:n_col))
        do k=1,n_row,1
            read(10,*) B(k,:)
        end do 
        TN = B  
        deallocate(B) 
        !print*,TN(n_row,:)

        !!! Lecture des coordonnées des centres des triangles
        read(10,*) n_row,n_col
        allocate(A(1:n_row,1:n_col)) 
        allocate(Centre_tri(1:n_row,1:n_col))
        do k=1,n_row,1
            read(10,*) A(k,:)
        end do 
        Centre_tri = A 
        deallocate(A)
        !print*,Centre_tri(n_row,:)

        !!! Lecture des distances entre le centre des triangles et leur voisinage
        read(10,*) n_row,n_col
        allocate(A(1:n_row,1:n_col)) 
        allocate(Dist_tri(1:n_row,1:n_col))
        do k=1,n_row,1
            read(10,*) A(k,:)
        end do 
        Dist_tri = A 
        deallocate(A)
        !print*,Dist_tri(n_row,:)

        !!! Lecture des coefficients de transmissibilité
        read(10,*) n_row,n_col
        allocate(A(1:n_row,1:n_col)) 
        allocate(Coef_trans(1:n_row,1:n_col))
        do k=1,n_row,1
            read(10,*) A(k,:)
        end do 
        Coef_trans = A 
        deallocate(A)
        !print*,Coef_trans(n_row,:)

        !!! Lecture des surfaces des triangles
        read(10,*) n_row,n_col
        allocate(A(1:n_row,1:n_col)) 
        allocate(Volume(1:n_row,1:n_col))
        do k=1,n_row,1
            read(10,*) A(k,:)
        end do 
        Volume = A 
        deallocate(A)
        !print*,Volume(n_row,:)

        !!! Lecture des longueurs des sommets
        read(10,*) n_row,n_col
        allocate(A(1:n_row,1:n_col)) 
        allocate(Long_segment(1:n_row,1:n_col))
        do k=1,n_row,1
            read(10,*) A(k,:)
        end do 
        Long_segment = A 
        deallocate(A)
        !print*,Long_segment(n_row,:)

        close(10)

        

    end subroutine load


    subroutine matrice_tri(equation,Points,ConnectivityList,&
    & Eall,Typesegment,TN,Centre_tri,Dist_tri,Coef_trans,Volume,Long_segment,A,b)
        implicit none
        !Entrées/Sorties
        integer, intent(in) :: equation
        double precision, dimension(:,:), allocatable, intent(in) :: Points,Centre_tri,Dist_tri,Coef_trans,Volume,Long_segment
        integer, dimension(:,:), allocatable, intent(in) :: ConnectivityList,Eall,Typesegment,TN
        type(matsparse), intent(out) :: A
        double precision, dimension(1:size(Volume,1)), intent(out) :: b 
        !Variables et constante
        integer :: i,j,k,n,Ncoef,kdiag
        double precision :: d,alpha,beta

        !Initialisation des variables
        n=size(Centre_tri,1)
        !allocate(b(1:n))
        Ncoef = count(TN>0)+n 
        k=1
        d = 0.0001
        alpha = 0.01
        beta = 0.05
        allocate(A%intI(1:Ncoef),A%intJ(1:Ncoef),A%PosDiag(1:n),A%ValMat(1:Ncoef),A%Diag(1:n))

        !Insert entries one-by-one
        do i=1,size(TN,1)
            b(i) = alpha*Volume(i,1)*f_vol(Centre_tri(i,:),equation)
            A%intI(k) = i 
            A%intJ(k) = i 
            A%PosDiag(i) = i 
            A%ValMat(k) = 0
            A%Diag(i) = 0
            kdiag=k 
            k=k+1
            do j=1,size(TN,2)
                if (TN(i,j)==0) then 
                    select case(equation)
                    case(1)
                        A%ValMat(kdiag) = A%ValMat(kdiag) + Coef_trans(i,j)
                        A%Diag(i) = A%Diag(i) + Coef_trans(i,j)
                        b(i) = b(i) + Coef_trans(i,j)*g_bord(0.5*(Points(Eall(3*(i-1)+j,1),:)+Points(Eall(3*(i-1)+j,2),:)),equation)
                    case(2)

                    case default
                        print*,'Equation non prise-en compte'
                    end select
                else 

                    select case(equation)
                    case(1)
                        A%intI(k) = i 
                        A%intJ(k) = TN(i,j) 
                        A%ValMat(k) = -Coef_trans(i,j) 
                        A%ValMat(kdiag) = A%ValMat(kdiag) + Coef_trans(i,j) 
                        A%Diag(i) = A%Diag(i) + Coef_trans(i,j) 
                        k=k+1
                    case(2)
                        A%intI(k) = i 
                        A%intJ(k) = TN(i,j) 
                        A%ValMat(k) = -d*Coef_trans(i,j) 
                        A%ValMat(kdiag) = A%ValMat(kdiag) + d*Coef_trans(i,j) + beta*Volume(i,1)
                        A%Diag(i) = A%Diag(i) + d*Coef_trans(i,j) + beta*Volume(i,1)
                        k=k+1
                    case default
                        print*,'Equation non prise-en compte'
                    end select
                    
                end if
            end do
        end do
    end subroutine matrice_tri

    subroutine matrice_laplacien(equation,A,b)
            implicit none
            !Entrées/Sorties
            integer, intent(in) :: equation
            type(matsparse), intent(out) :: A
            double precision, dimension(1:Nbt), intent(out) :: b 
            !Variables et constante
            integer :: i,j,k,n,Ncoef,kdiag
            double precision :: d,alpha,beta
    
            !Initialisation des variables
            n=size(Centre_tri,1)
            !allocate(b(1:n))
            Ncoef = count(TN>0)+n 
            k=1
            d = 0.0001
            alpha = 0.01
            beta = 0.05
            allocate(A%intI(1:Ncoef),A%intJ(1:Ncoef),A%PosDiag(1:n),A%ValMat(1:Ncoef),A%Diag(1:n))
    
            !Insert entries one-by-one
            do i=1,size(TN,1)
                b(i) = alpha*Volume(i,1)*f_vol(Centre_tri(i,:),equation)
                A%intI(k) = i 
                A%intJ(k) = i 
                A%PosDiag(i) = i 
                A%ValMat(k) = 0
                A%Diag(i) = 0
                kdiag=k 
                k=k+1
                do j=1,size(TN,2)
                    if (TN(i,j)==0) then 
                        select case(equation)
                        case(1)
                            A%ValMat(kdiag) = A%ValMat(kdiag) + Coef_trans(i,j)
                            A%Diag(i) = A%Diag(i) + Coef_trans(i,j)
                            b(i) = b(i) + Coef_trans(i,j)*g_bord(0.5*(Points(Eall(3*(i-1)+j,1),:)+Points(Eall(3*(i-1)+j,2),:)),equation)
                        case(2)
    
                        case default
                            print*,'Equation non prise-en compte'
                        end select
                    else 
    
                        select case(equation)
                        case(1)
                            A%intI(k) = i 
                            A%intJ(k) = TN(i,j) 
                            A%ValMat(k) = -Coef_trans(i,j) 
                            A%ValMat(kdiag) = A%ValMat(kdiag) + Coef_trans(i,j) 
                            A%Diag(i) = A%Diag(i) + Coef_trans(i,j) 
                            k=k+1
                        case(2)
                            A%intI(k) = i 
                            A%intJ(k) = TN(i,j) 
                            A%ValMat(k) = -d*Coef_trans(i,j) 
                            A%ValMat(kdiag) = A%ValMat(kdiag) + d*Coef_trans(i,j) + beta*Volume(i,1)
                            A%Diag(i) = A%Diag(i) + d*Coef_trans(i,j) + beta*Volume(i,1)
                            k=k+1
                        case default
                            print*,'Equation non prise-en compte'
                        end select
                        
                    end if
                end do
            end do
        end subroutine matrice_laplacien

    function bigradient(A,b, X0,tolerance) result(X)
        !     *--------------------------------------------------
        !     * Ce sous programme resoud le systeme lineaire 
        !     *    A x = b 
        !     * par la methode des bigradients conjugues 
        !     * le resultat de l'operation x = b/A 
        !--------------------------------------------------------    
        ! Modules
        !--------
        implicit none
        !--------------------------
        ! Declaration des arguments
        !--------------------------
        TYPE(matsparse), INTENT(in)                           :: A
        double precision, DIMENSION( SIZE(A%Diag)),INTENT(in) :: b,X0
        double precision, DIMENSION( SIZE(A%Diag))            :: X
        double precision , intent(in) :: tolerance
        !----------------------------------
        ! Declaration des variables locales
        !----------------------------------
        double precision, dimension (size(A%Diag))   :: Rst, Prj, Q, Rstb, Prjb,Qb
        double precision  :: alfa, bta,seuil,prscal, prscal9
        integer          :: nbigr, nbitermax
        ! ----------------------------
        ! Signification des variables
        ! ---------------------------
        ! Rst : reste = B - A * X
        ! Prj : projection
        ! Q   := A * Prj
        ! Qb  := transposee (A) * Prjb
        ! Rstb:= B - transposee (A) * X
        ! prscal : designe un produit scalaire
        !-------------------
        ! Debut du programme
        !-------------------
        !
        Rst = 0. ;  Prj = 0. ; Q = 0.
        Rstb = 0. ; Prjb = 0. ; Qb = 0.
        nbigr = 1
        !
        X = X0
        Rst =  b - mat_vec(A,X) 
        Rstb = b - transposee_mat_vec(A,X) 
        Prj = Rst
        Prjb = Rstb
        prscal = DOT_PRODUCT(Rstb,Rst)
        seuil = sqrt(DOT_PRODUCT(Rst, Rst))
        nbitermax = size(X0) + 500
     
        do  while((seuil > tolerance) .and. nbigr <nbitermax)! nbitermax
            
            !print*,nbigr,seuil,sum(X)
            Q =  mat_vec(A,Prj)
            alfa = prscal/DOT_PRODUCT(Q,Prjb)
            X = X + alfa * Prj
            Rst = Rst - alfa * Q
            Qb = transposee_mat_vec(A,Prjb)!transposee_mat_vec(A,Prjb) 
            Rstb = Rstb - alfa * Qb
            prscal9 = DOT_PRODUCT(Rstb,Rst)
            bta = prscal9/prscal
            prscal = prscal9 
            Prj = Rst + bta * Prj
            Prjb = Rstb + bta * Prjb
            seuil = sqrt(DOT_PRODUCT(Rst, Rst))
            nbigr = nbigr + 1
        end do 
        print*,alfa,seuil,nbigr
     
    end function bigradient

    function bcgstab(A,b, X0,tolerance) result(X)
        !     *--------------------------------------------------
        !     * Ce sous programme resoud le systeme lineaire 
        !     *    A x = b 
        !     * par la methode des bigradients conjugues 
        !     * le resultat de l'operation x = b/A 
        !--------------------------------------------------------    
        ! Modules
        !--------
        implicit none
        !--------------------------
        ! Declaration des arguments
        !--------------------------
        TYPE(matsparse), INTENT(in)                           :: A
        double precision, DIMENSION( SIZE(A%Diag)),INTENT(in) :: b,X0
        double precision, DIMENSION( SIZE(A%Diag))            :: X
        double precision , intent(in) :: tolerance
        !----------------------------------
        ! Declaration des variables locales
        !----------------------------------
        double precision, dimension (size(A%Diag))   :: r0, r, v, p,h,t, s
        double precision  :: alpha, beta, rho,seuil,omega,rho_p
        integer          :: nbigr, nbitermax
        X=X0
        r0 = b - mat_vec(A,X) 
        r = r0
        rho = 1; alpha = 1; omega = 1;
        v = 0. ; p = 0.;
        nbigr = 1; nbitermax=10
        seuil = sqrt(dot_product(r,r))
        rho_p = rho
        do while((seuil>tolerance).and.(nbigr<nbitermax))
            rho = dot_product(r0,r)
            beta = rho/rho_p * alpha/omega 
            p = r + beta*(p-omega*v)
            v = mat_vec(A,p)
            alpha=rho/dot_product(r0,v)
            h = x + alpha*p 
            s = r - alpha*v 
            t = mat_vec(A,s) 
            omega = dot_product(t,s)/dot_product(t,t)
            X = h + omega * s 
            r = s - omega * t 
            seuil = sqrt(dot_product(r,r))
            rho_p=rho 
            nbigr = nbigr + 1
        end do
        print*,alpha,seuil,nbigr


    end function bcgstab

    function gradopt(A,b, X0,tolerance) result(X)
        TYPE(matsparse), INTENT(in)                           :: A
        double precision, DIMENSION( SIZE(A%Diag)),INTENT(in) :: b,X0
        double precision, DIMENSION( SIZE(A%Diag))            :: X
        double precision , intent(in) :: tolerance

        double precision, dimension(size(A%Diag)) :: omega , Prod
        double precision :: seuil,rho
        integer :: iter

        X=X0
        omega = mat_vec(A,X)-b 
        seuil = sqrt(dot_product(omega,omega))
        iter = 1

        do while((seuil>tolerance).and.(iter<2000))
            Prod = mat_vec(A,omega)
            rho = dot_product(omega,omega)/dot_product(Prod,omega)
            X = X - rho*omega 
            iter = iter + 1
            omega = mat_vec(A,X)-b 
            seuil = sqrt(dot_product(omega,omega))
        end do
        print*,rho,seuil,iter
    end function gradopt

    subroutine foncNewton(Points,ConnectivityList,Eall,Typesegment,delta_t,&
        & TN,Centre_tri,Dist_tri,Coef_trans,Volume,Long_segment,Fonc,GradFonc,x,U,Nutriment)
        implicit none
        !Entrées/Sorties
        double precision, dimension(:,:), allocatable, intent(in) :: Points,Centre_tri,Dist_tri,Coef_trans,Volume,Long_segment
        integer, dimension(:,:), allocatable, intent(in) :: ConnectivityList,Eall,Typesegment,TN
        double precision, dimension(:), allocatable, intent(in) :: x,U,Nutriment 
        type(matsparse), intent(inout) :: GradFonc
        double precision, dimension(1:size(Volume,1)), intent(out) :: Fonc 
        double precision, intent(in) :: delta_t

        integer :: n,Ncoef,k, i, j ,kdiag
        !Initialisation des variables
        n=size(Centre_tri,1)
        Ncoef = count(TN>0)+n 
        k=1
        !allocate(GradFonc%intI(1:Ncoef),GradFonc%intJ(1:Ncoef),GradFonc%PosDiag(1:n),&
        !& GradFonc%ValMat(1:Ncoef),GradFonc%Diag(1:n))
        do i=1,size(TN,1)
            Fonc(i) = x(i)-delta_t*production(x(i)) - U(i)
            GradFonc%intI(k) = i 
            GradFonc%intJ(k) = i 
            GradFonc%PosDiag(i) = i 
            GradFonc%ValMat(k) = 1 - Delta_t*prodprime(x(i))
            GradFonc%Diag(i) = 1 - Delta_t*prodprime(x(i))
            kdiag=k 
            k=k+1
            do j=1,size(TN,2)
                if (TN(i,j)/=0) then 
                    GradFonc%intI(k) = i 
                    GradFonc%intJ(k) = TN(i,j) 

                    GradFonc%ValMat(k) = -Delta_t/Volume(i,1)*Coef_trans(i,j)*diffprime(x(TN(i,j))) & 
                    & - Delta_t/Volume(i,1)*Coef_trans(i,j)*chemo(x(TN(i,j)))*&
                    & 0.5*(-abs(Nutriment(TN(i,j))-Nutriment(i))+Nutriment(TN(i,j))-Nutriment(i))

                    GradFonc%ValMat(kdiag) = GradFonc%ValMat(kdiag) + Delta_t/Volume(i,1)*Coef_trans(i,j)*diffprime(x(i)) &
                    & + Delta_t/Volume(i,1)*Coef_trans(i,j)*chemoprime(x(i))*&
                    & 0.5*(abs(Nutriment(TN(i,j))-Nutriment(i))+Nutriment(TN(i,j))-Nutriment(i))

                    GradFonc%Diag(i) = GradFonc%Diag(i) + Delta_t/Volume(i,1)*Coef_trans(i,j)*diffprime(x(i))&
                    & + Delta_t/Volume(i,1)*Coef_trans(i,j)*chemoprime(x(i))*&
                    & 0.5*(abs(Nutriment(TN(i,j))-Nutriment(i))+Nutriment(TN(i,j))-Nutriment(i))

                    Fonc(i) = Fonc(i) + Delta_t/Volume(i,1)*Coef_trans(i,j)*(diffusion(x(i))-diffusion(x(TN(i,j))))&
                    & + Delta_t/Volume(i,1)*Coef_trans(i,j)*chemo(x(i))*&
                    & 0.5*(abs(Nutriment(TN(i,j))-Nutriment(i))+Nutriment(TN(i,j))-Nutriment(i)) &
                    & + Delta_t/Volume(i,1)*Coef_trans(i,j)*chemo(x(TN(i,j)))*&
                    & 0.5*(-abs(Nutriment(TN(i,j))-Nutriment(i))+Nutriment(TN(i,j))-Nutriment(i))

                    k=k+1
                end if
            end do
        end do
    end subroutine foncNewton

end module load_module