!            **************************
!            **  SUBROUTINE ASSEMBLEVF4  **
!            **************************
!******************************************************************************
!     *--------------------------------------
!     * Ce sous programme assemble la contribution de diffusion
!     * Assemblage  par segement l'operateur 
!      - d (d^2u/dx^2 + d^2/dy^2)
! en utilisant la m�thode des volumes finis VF4
!     *----------------------------
!     * Fonctionnement du programme
!     *----------------------------
!     * Utilise les modules longr, imprime, parmmmage
!******************************************************************************
SUBROUTINE ASSEMBLEVITESSE( A )
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
  ! Declaration des variables local  !-----------------------------
  CHARACTER(len=6)      :: oldprf
  INTEGER               :: is, js, ik, jL, iseg, ii, jv,  NcoefMat
  REAL(kind=long)       :: coef, x1, y1 ,nx,ny,sigmaKL,nnorm,inner

  !-------------------
  ! Debut du programme
  !-------------------
  oldprf = prefix
  prefix = 'ASSEMV'

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

        !! calcul de la normale ext�rieure n = x(Jl)-x(ik)
        nx= CoordK(1, jL)- CoordK(1, ik); ny= CoordK(2, jL)- CoordK(2, ik);
        nnorm= sqrt(nx*nx+ny*ny); nx=nx/nnorm; ny =ny/nnorm

        sigmaKL= sqrt( (CoordS(1,is)-CoordS(1,js))**2+ (CoordS(2,is)-CoordS(2,js))**2 )
        inner = sigmaKL* (VitesseSeg(1,iseg)*nx+VitesseSeg(2,nseg)*ny)

        if (inner >=0) then 

           !! contribution dans le coefficient  (ik,ik)
           CALL Ajout (ik, ik, inner,  A )

           !! contribution dans le coefficient  (jL,ik)
           CALL Ajout (jL, ik, -inner, A )

        else
           !! contribution dans le corfficient  (ik,jL)
           CALL Ajout (ik, jL, inner, A )

           !! contribution dans le coefficient  (jL,iL)
           CALL Ajout (jL, jL, -inner, A )
        end  if

     case(dirichlet)

        is = NuSeg(1,iseg); js=NuSeg(2,iseg)
        ik= NumTVoisSeg(1,iseg) 

        !
        ! valeur de la donn�e de la solution sur le bord  au point x1, y1
        !
        x1= (CoordS(1,is)+CoordS(1,js))/2.D0 ; y1= (CoordS(2,is)+CoordS(2,js))/2.D0
        !
        ! calcul exacte de la projection orthogonale du centre ik sur le segment sur le bord
        !

         
        
        !! calcul de la normale ext�rieure n = x(Jl)-x(ik)
        nx= x1- CoordK(1, ik); ny= y1- CoordK(2, ik);
        nnorm= sqrt(nx*nx+ny*ny); nx=nx/nnorm; ny =ny/nnorm

        sigmaKL= sqrt( (CoordS(1,is)-CoordS(1,js))**2+ (CoordS(2,is)-CoordS(2,js))**2)
        inner = sigmaKL* (VitesseSeg(1,iseg)*nx+VitesseSeg(2,nseg)*ny)
        !
        if (inner >=0) then 

           !! contribution dans le coefficient  (ik,ik)
           CALL Ajout (ik, ik, inner,  A )
        else
           ! Ajout des conditions aux limlites 
           A%Bg (ik) = A%Bg(ik) - inner*gbord2(x1,y1,choixpb)
           !
        end if
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
END SUBROUTINE ASSEMBLEVITESSE
