subroutine assemblediffusionKS (A,U)
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
  double precision, dimension(:) :: U 
  !----------------------------------
  ! Declaration des variables locales
  !----------------------------------
  CHARACTER(len=6)      :: oldprf
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
        CALL Ajout (ik, ik, delta/AireK(ik)*TauKL(iseg)*diffusionprime(U(ik)),  A )

        CALL Ajout (ik, jL, -delta/AireK(ik)*TauKL(iseg)*diffusionprime(U(jL)), A )

        A%Bg(ik) = A%Bg(ik) + delta/AireK(ik)*TauKL(iseg)*(diffusion(U(ik))-diffusion(U(jL)))

        !! contribution dans la ligne jL

        CALL Ajout (jL, jL, delta/AireK(ik)*TauKL(iseg)*diffusionprime(U(jL)), A )
        CALL Ajout (jL, ik, -delta/AireK(ik)*TauKL(iseg)*diffusionprime(U(ik)), A )

        A%Bg(jL) = A%Bg(jL) + delta/AireK(ik)*TauKL(iseg)*(diffusion(U(jL))-diffusion(U(ik)))

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
  IF (iprint >=5) THEN
     !-------------------------------
     ! Calcul de la digonale de A
     !-------------------------------
     Allocate(A%Diag(1:size(A%IndPL)-1));Allocate(A%PosiDiag(1:size(A%IndPL)-1))
     DO is = 1, size(A%IndPL)-1
        DO jv = A%IndPL(is), A%IndPL(is+1)-1 
           If (A%Indc(jv) == is) then 
              A%diag(is) = A%Tmat(jv)
              A%PosiDiag(is) = jv
              exit
           endif
         END DO
     END DO
     !
     !
     CALL prvari(uprint,'NcoefMat = ', NcoefMat )
     WRITE(uprint,*) 'A%TMAT de A est '
     WRITE(uprint,*) ( A%Tmat(is) , is=1, NcoefMat )
     WRITE(uprint,*) 'A%IndC de A est '
     WRITE(uprint,*) ( A%Indc(is) , is=1, NcoefMat )
     WRITE(uprint,*) 'A%IndPL de A est '
     WRITE(uprint,*) ( A%IndPL(is), is=1, Nbt+1 )
     WRITE(uprint,*) 'A%DIAG la diagonale de A est '
     WRITE(uprint,*) ( A%Diag(is), is=1, Nbt )
     WRITE(uprint,*) 'A%PosiDiag Position de la diagonale de A est '
     WRITE(uprint,*) ( A%PosiDiag(is), is=1, Nbt )
  ENDIF

  !-----------------
  ! Fin du programme
  !-----------------
  prefix = oldprf

  RETURN



end subroutine assemblediffusionKS