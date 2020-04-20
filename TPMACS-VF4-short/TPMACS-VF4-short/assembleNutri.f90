subroutine assembleNutri(N,E)
    USE longr
    USE parmmage
    USE imprime
    Use fsourcemod
    implicit none
    type(MatCreux), intent(inout) :: N
    real(kind=long), dimension(Nbt), intent(in) :: E 

    integer :: i,iseg

    N%F = 0.D0    
    
    DO iseg = 1, Nseg
        VitesseSeg(1:2,iseg) = Vitesse(CoordS(1,iseg),CoordS(2,iseg),choixpb)
    END DO

    N%Bg= E*AireK*Coef_prod
    call assembleVF4( N )  
    call assembleVitesse( N ) 
    call assembletheta(N)   
   

end subroutine assembleNutri