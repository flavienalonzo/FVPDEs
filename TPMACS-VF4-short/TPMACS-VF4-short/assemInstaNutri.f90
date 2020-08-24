subroutine assemInstaNutri(N,E,U,Nut)
    USE longr
    USE parmmage
    USE imprime
    Use fsourcemod
    implicit none
    real(kind=long), dimension(Nbt), intent(inout) :: N
    real(kind=long), dimension(Nbt), intent(in) :: E, U, Nut

    integer :: i,iseg, ii, is, js, ik, jL

    N = Nut 

    do iseg=1,Nseg

        ii = (NtypSeg(iseg))
        select case (ii) 
        case (0) 
        is = NuSeg(1,iseg); js=NuSeg(2,iseg)
        ik= NumTVoisSeg(1,iseg); jL =NumTVoisSeg(2,iseg)

        !! contribution dans la ligne ik
        !CALL Ajout (ik, ik, delta/AireK(ik)*Coef_diffusion*TauKL(iseg),  N )
        N(ik) = N(ik) - delta/AireK(ik)*Coef_diffusion*TauKL(iseg)*Nut(ik)
        !CALL Ajout (ik, jL, -delta/AireK(ik)*Coef_diffusion*TauKL(iseg), N )
        N(ik) = N(ik) + delta/AireK(ik)*Coef_diffusion*TauKL(iseg)*Nut(jL)
        !! contribution dans la ligne jL

        !CALL Ajout (jL, jL, delta/AireK(jL)*Coef_diffusion*TauKL(iseg), N )
        N(jL) = N(jL) - delta/AireK(jL)*Coef_diffusion*TauKL(iseg)*Nut(jL)
        !CALL Ajout (jL, ik, -delta/AireK(jL)*Coef_diffusion*TauKL(iseg), N )
        N(jL) = N(jL) + delta/AireK(jL)*Coef_diffusion*TauKL(iseg)*Nut(ik)
        case default
        end select 

    end do
    do i=1,Nbt
        !CALL Ajout (i,i, delta*Theta, N )
        N(i) = N(i) - delta*theta*Nut(i)
        N(i) = N(i) + delta*(U(i)+E(i))*Coef_prod
    end do

end subroutine assemInstaNutri