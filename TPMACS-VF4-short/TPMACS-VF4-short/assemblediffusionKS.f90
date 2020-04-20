subroutine assemblediffusionKS (A,U,E,entity)
    USE longr
  USE imprime
  USE parmmage
  USE fsourcemod
  use parameters 

  IMPLICIT NONE

  !--------------------------
  ! Declaration des arguments
  !--------------------------
  TYPE(MatCreux)       :: A
  real (kind = long), dimension(Nbt) :: U , E
  !----------------------------------
  ! Declaration des variables locales
  !----------------------------------
  CHARACTER(len=6)      :: oldprf, entity
  INTEGER               :: is, js, ik, jL, iseg, ii, jv,  NcoefMat
  REAL(kind=long)       :: coef, x1, y1 

  !-------------------
  ! Debut du programme
  !-------------------
  oldprf = prefix
  prefix = 'ASSEMB'

  !------
  ! Corps
  !------

  !  Segment int�reurs
  DO iseg = 1, Nseg
     ii = (NtypSeg(iseg)) 

     Select case (ii) 

     case (0)   !! segment � l'nterieur  
        is = NuSeg(1,iseg); js=NuSeg(2,iseg)
        ik= NumTVoisSeg(1,iseg); jL =NumTVoisSeg(2,iseg)

        !! contribution dans la ligne ik
        CALL Ajout (ik, ik, delta/AireK(ik)*TauKL(iseg)*diffprime(U(ik),E(ik),entity),  A )

        CALL Ajout (ik, jL, -delta/AireK(ik)*TauKL(iseg)*diffprime(U(jL),E(jL),entity), A )

        A%Bg(ik) = A%Bg(ik) + delta/AireK(ik)*TauKL(iseg)*(diffusion(U(ik),E(ik),entity)-diffusion(U(jL),E(jL),entity))

        !! contribution dans la ligne jL

        CALL Ajout (jL, jL, delta/AireK(jL)*TauKL(iseg)*diffprime(U(jL),E(jL),entity), A )
        CALL Ajout (jL, ik, -delta/AireK(jL)*TauKL(iseg)*diffprime(U(ik),E(ik),entity), A )

        A%Bg(jL) = A%Bg(jL) + delta/AireK(jL)*TauKL(iseg)*(diffusion(U(jL),E(jL),entity)-diffusion(U(ik),E(ik),entity))

     case(dirichlet)
        stop 'Pas de condition de Dirichlet prevue'
     case(neumann)
        ! conditions aux limites homogenes le flux est nul
        ! sinon il faut changer A%Bg
     case default 
        stop 'pb assemble'
     end select
  end DO


  !------------
  ! Impressions
  !------------
  NcoefMat = size(A%Tmat)
  !!WRITE(*,*) ( A%Tmat(is) , is=1, NcoefMat )

  !-----------------
  ! Fin du programme
  !-----------------
  prefix = oldprf

  RETURN



end subroutine assemblediffusionKS