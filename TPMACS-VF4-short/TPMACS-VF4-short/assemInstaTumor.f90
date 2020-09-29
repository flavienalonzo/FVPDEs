subroutine assemInstaTumor(Tum,E,U,Nut)
    USE longr
    USE parmmage
    USE imprime
    Use fsourcemod
    use parameters

    implicit none
    real(kind=long), dimension(Nbt), intent(inout) :: Tum
    real(kind=long), dimension(Nbt), intent(in) :: E, U, Nut
    real(kind=long) :: approx, cplusi, cmoinsi, cplusj, cmoinsj
    integer :: i,iseg, ii, is, js, ik, jL

    Tum = U

    do iseg=1,Nseg
        ii = (NtypSeg(iseg))
        select case (ii) 
        case (0) 
        is = NuSeg(1,iseg); js=NuSeg(2,iseg)
        ik= NumTVoisSeg(1,iseg); jL =NumTVoisSeg(2,iseg)

        !! contribution dans la ligne ik
        !CALL Ajout (ik, ik, delta/AireK(ik)*Coef_diffusion*TauKL(iseg),  N )
        Tum(ik) = Tum(ik) - delta/AireK(ik)*Diff_u*TauKL(iseg)*U(ik)
        !CALL Ajout (ik, jL, -delta/AireK(ik)*Coef_diffusion*TauKL(iseg), N )
        Tum(ik) = Tum(ik) + delta/AireK(ik)*Diff_u*TauKL(iseg)*U(jL)
        !! contribution dans la ligne jL

        !CALL Ajout (jL, jL, delta/AireK(jL)*Coef_diffusion*TauKL(iseg), N )
        Tum(jL) = Tum(jL) - delta/AireK(jL)*Diff_u*TauKL(iseg)*U(jL)
        !CALL Ajout (jL, ik, -delta/AireK(jL)*Coef_diffusion*TauKL(iseg), N )
        Tum(jL) = Tum(jL) + delta/AireK(jL)*Diff_u*TauKL(iseg)*U(ik)

        call approxGrad(iseg,Nut,approx)

        cplusi = 0.5*(Nut(jL)-Nut(ik)+ abs(Nut(jL)-Nut(ik)))
        cmoinsi = 0.5*(Nut(jL)-Nut(ik) - abs(Nut(jL)-Nut(ik)))
        cplusj = 0.5*(Nut(ik)-Nut(jL) + abs(Nut(ik) - Nut(jL)))
        cmoinsj = 0.5*(Nut(ik)-Nut(jL) - abs(Nut(ik) - Nut(jL)))

        Tum(ik) = Tum(ik) - delta/AireK(ik)*TauKL(iseg)*chi_u*(U(ik)*cplusi&
        &+U(jL)*cmoinsi)/(sqrt(1+satur_nutri*approx))

        Tum(jL) = Tum(jL) - delta/AireK(jL)*TauKL(iseg)*chi_u*(U(jL)*cplusj&
        &+U(ik)*cmoinsj)/(sqrt(1+satur_nutri*approx))

        case default
        end select 

    end do
    do i=1,Nbt 
        Tum(i) = Tum(i) +delta*reaction((/U(i),Nut(i)/),'u_norm')
    end do


end subroutine assemInstaTumor