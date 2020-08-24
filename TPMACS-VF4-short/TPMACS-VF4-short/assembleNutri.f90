subroutine assembleNutri(N,E,U)
    USE longr
    USE parmmage
    USE imprime
    Use fsourcemod
    implicit none
    type(MatCreux), intent(inout) :: N
    real(kind=long), dimension(Nbt), intent(in) :: E, U

    integer :: i,iseg

    N%F = 0.D0    
    
    DO iseg = 1, Nseg
        VitesseSeg(1:2,iseg) = Vitesse(CoordS(1,iseg),CoordS(2,iseg),choixpb)
    END DO
    do i=1,Nbt
        N%Bg(i)= E(i)*Coef_prod*AireK(i)
    end do
    call assembleVF4( N )  
    call assembleVitesse(N) 
    call assembletheta(N)   
   

end subroutine assembleNutri