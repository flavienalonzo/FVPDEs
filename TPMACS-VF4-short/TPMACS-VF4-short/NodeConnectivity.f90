subroutine NodeConnectivity
    ! Ce programme établie la liste des triangles connectés à chaque sommets
    use longr 
    use imprime 
    use parmmage 
    implicit none 

    integer :: kt , it, ks
    integer, dimension(1:Nbs) :: kTrS

    allocate(NbrTrS(1:Nbs),AireS(1:Nbs))
    NbrTrS = 0
    AireS = 0.D0

    do kt=1,Nbt
        do it=1,3
            NbrTrS(NuSoK(it,kt)) = NbrTrS(NuSoK(it,kt)) + 1
        end do
    end do

    NmaxT = MAXVAL(NbrTrS)

    allocate(NuTrS(Nbs,NmaxT))
    kTrS = 1
    do kt=1,Nbt
        do it=1,3
            NuTrS(NuSoK(it,kt),kTrS(NuSoK(it,kt))) = kt
            AireS(NuSoK(it,kt)) = AireS(NuSoK(it,kt)) + AireK(kt)
            kTrS(NuSoK(it,kt)) = kTrS(NuSoK(it,kt)) + 1
        end do
    end do

end subroutine NodeConnectivity