subroutine assemblechemoKS (A,U,C,E,entity)
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
  real (kind = long), dimension(Nbt) :: U, C , E
  !----------------------------------
  ! Declaration des variables locales
  !----------------------------------
  CHARACTER(len=6)      :: oldprf, entity
  INTEGER               :: is, js, ik, jL, iseg, ii, jv,  NcoefMat
  REAL(kind=long)       :: coef, x1, y1, cplusi, cmoinsi, cplusj, cmoinsj

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
        cplusi = 0.5*(C(jL)-C(ik)+ abs(C(jL)-C(ik)))
        cmoinsi = 0.5*(C(jL)-C(ik) - abs(C(jL)-C(ik)))
        cplusj = 0.5*(C(ik)-C(jL) + abs(C(ik) - C(jL)))
        cmoinsj = 0.5*(C(ik)-C(jL) - abs(C(ik) - C(jL)))
        !! contribution dans la ligne ik
        CALL Ajout (ik, ik, delta/AireK(ik)*TauKL(iseg)*chemoprime(U(ik),E(ik),entity)*cplusi,  A )

        CALL Ajout (ik, jL, delta/AireK(ik)*TauKL(iseg)*chemoprime(U(jL),E(jL),entity)*cmoinsi, A )

        A%Bg(ik) = A%Bg(ik) + delta/AireK(ik)*TauKL(iseg)*(chemo(U(ik),E(ik),entity)*cplusi+chemo(U(jL),E(jL),entity)*cmoinsi)

        !! contribution dans la ligne jL

        CALL Ajout (jL, jL, delta/AireK(jL)*TauKL(iseg)*chemoprime(U(jL),E(jL),entity)*cplusj, A )
        CALL Ajout (jL, ik, delta/AireK(jL)*TauKL(iseg)*chemoprime(U(ik),E(ik),entity)*cmoinsj, A )

        A%Bg(jL) = A%Bg(jL) + delta/AireK(jL)*TauKL(iseg)*(chemo(U(jL),E(jL),entity)*cplusj+chemo(U(ik),E(ik),entity)*cmoinsj)

     case(dirichlet)
        stop 'Pas de condition de Dirichlet prevue'
     case(neumann)
        ! conditions aux limites homogenes le flux est nul
        ! sinon il faut changer A%Bg
     case default 
        stop 'pb assemble'
     end select
  end DO

  RETURN
end subroutine assemblechemoKS