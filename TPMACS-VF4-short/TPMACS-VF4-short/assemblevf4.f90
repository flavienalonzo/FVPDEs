!            **************************
!            **  SUBROUTINE ASSEMBLEVF4  **
!            **************************
!******************************************************************************
!     *--------------------------------------
!     * Ce sous programme assemble la contribution de diffusion
!     * Assemblage  par segement l'operateur 
!      - d (d^2u/dx^2 + d^2/dy^2)
! en utilisant la méthode des volumes finis VF4
!     *----------------------------
!     * Fonctionnement du programme
!     *----------------------------
!     * Utilise les modules longr, imprime, parmmmage
!******************************************************************************
SUBROUTINE ASSEMBLEVF4( A )
  ! 
  !
  !--------
  ! Modules
  !--------
  USE longr
  USE imprime
  USE parmmage
  USE fsourcemod

  IMPLICIT NONE

  !--------------------------
  ! Declaration des arguments
  !--------------------------
  TYPE(MatCreux)       :: A
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

  !  Segment intéreurs
  DO iseg = 1, Nseg
     ii = (NtypSeg(iseg)) 

     Select case (ii) 

     case (0)   !! segment à l'nterieur  
        is = NuSeg(1,iseg); js=NuSeg(2,iseg)
        ik= NumTVoisSeg(1,iseg); jL =NumTVoisSeg(2,iseg)

        !! contribution dans la ligne ik
        CALL Ajout (ik, ik, Coef_diffusion*TauKL(iseg),  A )

        CALL Ajout (ik, jL, -Coef_diffusion*TauKL(iseg), A )

        !! contribution dans la ligne jL

        CALL Ajout (jL, jL, Coef_diffusion*TauKL(iseg), A )
        CALL Ajout (jL, ik, -Coef_diffusion*TauKL(iseg), A )


     case(dirichlet)

        is = NuSeg(1,iseg); js=NuSeg(2,iseg)
        ik= NumTVoisSeg(1,iseg) 

        ! contribution dans la ligne ik

        CALL Ajout (ik, ik, Coef_diffusion*TauKL(iseg),  A )

        !
        ! valeur de la donnée de la solution sur le bord  au point x1, y1
        !
        x1= (CoordS(1,is)+CoordS(1,js))/2.D0 ; y1= (CoordS(2,is)+CoordS(2,js))/2.D0
        !
        ! Ajout des conditions aux limlites 
        A%Bg (ik) = A%Bg(ik) + Coef_diffusion*TauKL(iseg)*gbord(x1,y1,choixpb)
        !
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
END SUBROUTINE ASSEMBLEVF4








