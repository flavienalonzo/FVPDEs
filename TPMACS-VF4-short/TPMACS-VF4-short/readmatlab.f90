subroutine readmatlab
    use parmmage 
    use longr
    use imprime
    implicit none 

    integer :: n_row,n_col,k
    integer, dimension(:,:), allocatable :: B
    integer, dimension(:), allocatable :: D
    double precision, dimension(:,:), allocatable :: A 
    double precision, dimension(:), allocatable :: C

    open(unit=10,file = nom_mesh, status = 'old')

    !!! Lecture des coordonnées des sommets des triangles
    read(10,*) n_row,n_col
    allocate(A(1:n_row,1:n_col)) 
    allocate(Coords(1:n_col,1:n_row))
    do k=1,n_row,1
        read(10,*) A(k,:)
    end do 
    Coords = transpose(A) !Points dans Matlab
    Nbs = n_row 
    deallocate(A)

    !!! Lecture de la liste des connectivités des triangles
    read(10,*) n_row,n_col
    allocate(B(1:n_row,1:n_col))
    allocate(NuSoK(1:n_col,1:n_row))
    do k=1,n_row,1
        read(10,*) B(k,:)
    end do 
    NuSoK = transpose(B)  !ConnectivityList dans matlab 
    deallocate(B) 

    !!! Lecture de la liste des segments des triangles
    read(10,*) n_row,n_col
    allocate(B(1:n_row,1:n_col))
    allocate(NuMSeg(1:n_col,1:n_row))
    do k=1,n_row,1
        read(10,*) B(k,:)
    end do 
    NuMSeg = transpose(B)  
    Nseg = n_row
    deallocate(B) 

    !!! Lecture de la liste donnant le type des segments
    read(10,*) n_row,n_col
    allocate(D(1:n_row))
    allocate(ntypseg(1:n_row))
    do k=1,n_row,1
        read(10,*) D(k)
    end do 
    ntypseg = D 
    deallocate(D) 

    !!! Lecture de la liste des triangles voisins d'un triangle
    read(10,*) n_row,n_col
    allocate(B(1:n_row,1:n_col))
    allocate(NuVoisK(1:n_col,1:n_row))
    do k=1,n_row,1
        read(10,*) B(k,:)
    end do 
    NuVoisK = transpose(B)  
    Nbt = n_row 
    deallocate(B) 

    !!! Lecture des coordonnées des centres des triangles
    read(10,*) n_row,n_col
    allocate(A(1:n_row,1:n_col)) 
    allocate(CoordK(1:n_col,1:n_row))
    do k=1,n_row,1
        read(10,*) A(k,:)
    end do 
    CoordK = transpose(A) 
    deallocate(A)

    !!! Lecture des distances entre le centre des triangles et leur voisinage
    read(10,*) n_row,n_col
    allocate(A(1:n_row,1:n_col)) 
    allocate(dKL(1:n_row*n_col))
    do k=1,n_row,1
        read(10,*) A(k,:)
    end do 
    dKL = reshape(A,(/n_row*n_col/)) 
    deallocate(A)

    !!! Lecture des coefficients de transmissibilité
    read(10,*) n_row,n_col
    allocate(A(1:n_row,1:n_col)) 
    allocate(TauKL(1:n_row*n_col))
    do k=1,n_row,1
        read(10,*) A(k,:)
    end do 
    TauKL = reshape(A,(/n_row*n_col/)) 
    deallocate(A)

    !!! Lecture des surfaces des triangles
    read(10,*) n_row,n_col
    allocate(C(1:n_row)) 
    allocate(AireK(1:n_row))
    do k=1,n_row,1
        read(10,*) C(k)
    end do 
    AireK = C 
    deallocate(C)

    !!! Lecture des longueurs des sommets
    read(10,*) n_row,n_col
    allocate(A(1:n_row,1:n_col)) 
    !allocate(Long_segment(1:n_row,1:n_col))
    do k=1,n_row,1
        read(10,*) A(k,:)
    end do 
    !Long_segment = A 
    deallocate(A)

    close(10)

end subroutine readmatlab 