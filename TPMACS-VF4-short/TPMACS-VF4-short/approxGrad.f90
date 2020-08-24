subroutine approxGrad(iseg,quant,approx)

! Calcul une approximation de |nabla quant|^2 sur le milieu du segment iseg
use longr
use imprime
use parmmage 

implicit none 

real (kind=long), intent(out) :: approx
integer, intent(in) :: iseg
real (kind=long), dimension(1:Nbt), intent(in) :: quant 
integer :: S1, S2, T1, T2, k
real (kind=long) :: cK, cL, c1, c2, sigma12

! Les sommets définissant le segment
S1 = NuSeg(1,iseg)
S2 = NuSeg(2,iseg)
! Les triangles communs au segment 
T1 = NumTVoisSeg(1,iseg)
T2 = NumTVoisSeg(2,iseg)
! Valeurs aux points K et L et approximations de celles en 1 et 2
if (NTypSeg(iseg)==0) then 
    cK = quant(T1)
    cL = quant(T2)
else ! Si on est sur le bord la composante du gradient selon la normale est nulle
    cK = 0
    cL = 0 
end if
c1 = 0.D0
c2 = 0.D0 
do k=1,NbrTrS(S1)
    c1 = c1 + AireK(NuTrS(S1,k))*quant(NuTrS(S1,k))
end do
c1 = c1/AireS(S1)
do k=1,NbrTrS(S2)
    c2 = c2 + AireK(NuTrS(S2,k))*quant(NuTrS(S2,k))
end do
c2 = c2/AireS(S2)
!longueur du segment reliant c1 à c2
sigma12 = TauKL(iseg)*dKL(iseg)
!Calcul de l'approximation
approx =  ((c1-c2)/sigma12)**2 + ((cK-cL)/dKL(iseg))**2
end subroutine approxGrad