subroutine assemblechemoKS (A,Tab_U,entity,Tab_entity,Tab_equa,equa,k_enty,k_chemo)
   USE longr
   USE imprime
   USE parmmage
   USE fsourcemod
   use parameters 

  IMPLICIT NONE

  !--------------------------
  ! Declaration des arguments
  !--------------------------
  TYPE(MatCreux), dimension(n_enty)      :: A
  real (kind = long), dimension(n_enty,Nbt) :: Tab_U
  character(len=6), dimension(n_enty) :: Tab_entity, Tab_equa
  character(len=6), intent(in) :: equa
  !----------------------------------
  ! Declaration des variables locales
  !----------------------------------
  CHARACTER(len=6)      :: oldprf, entity
  INTEGER               :: is, js, ik, jL, iseg, ii, jv,  NcoefMat, k_enty, k, k_chemo
  REAL(kind=long)       :: coef, x1, y1, cplusi, cmoinsi, cplusj, cmoinsj, approx
  real(kind=long), dimension(Nbt) :: C

  !-------------------
  ! Debut du programme
  !-------------------
  oldprf = prefix
  prefix = 'ASSEMB'

  !------
  ! Corps
  !------

   C = Tab_U(k_chemo,:)
  !  Segment int�reurs
  DO iseg = 1, Nseg
     ii = (NtypSeg(iseg)) 
      call approxGrad(iseg,C,approx)
     Select case (ii) 

     case (0)   !! segment � l'nterieur  
        is = NuSeg(1,iseg); js=NuSeg(2,iseg)
        ik= NumTVoisSeg(1,iseg); jL =NumTVoisSeg(2,iseg)
        cplusi = 0.5*(C(jL)-C(ik)+ abs(C(jL)-C(ik)))
        cmoinsi = 0.5*(C(jL)-C(ik) - abs(C(jL)-C(ik)))
        cplusj = 0.5*(C(ik)-C(jL) + abs(C(ik) - C(jL)))
        cmoinsj = 0.5*(C(ik)-C(jL) - abs(C(ik) - C(jL)))
      if (equa=='instat') then 
        !! contribution dans la ligne ik
        do k=1,n_enty
         if (Tab_equa(k)==equa) then
         CALL Ajout (ik, ik, delta/AireK(ik)*TauKL(iseg)*chemoprime(Tab_U(:,ik),entity,Tab_entity(k))*cplusi,  A(k) )

         CALL Ajout (ik, jL, delta/AireK(ik)*TauKL(iseg)*chemoprime(Tab_U(:,jL),entity,Tab_entity(k))*cmoinsi, A(k) )
         end if
        end do
        if (cplusi>0) then
         CALL Ajout (ik, ik, delta/AireK(ik)*TauKL(iseg)*chemo(Tab_U(:,ik),entity),  A(k_chemo) )
        else 
         CALL Ajout (ik, jL, delta/AireK(ik)*TauKL(iseg)*chemo(Tab_U(:,jL),entity), A(k_chemo) )
        end if
        A(k_enty)%Bg(ik) = A(k_enty)%Bg(ik) + delta/AireK(ik)*TauKL(iseg)*(chemo(Tab_U(:,ik),entity)*cplusi&
        &+chemo(Tab_U(:,jL),entity)*cmoinsi)

        !! contribution dans la ligne jL
        do k=1,n_enty
         if (Tab_equa(k)==equa) then
         CALL Ajout (jL, jL, delta/AireK(jL)*TauKL(iseg)*chemoprime(Tab_U(:,jL),entity,Tab_entity(k))*cplusj, A(k) )
         CALL Ajout (jL, ik, delta/AireK(jL)*TauKL(iseg)*chemoprime(Tab_U(:,ik),entity,Tab_entity(k))*cmoinsj, A(k) )
         end if
        end do

        if (cplusj>0) then
         CALL Ajout (jL, jL, delta/AireK(jL)*TauKL(iseg)*chemo(Tab_U(:,jL),entity), A(k_chemo) )
        else 
         CALL Ajout (jL, ik, delta/AireK(jL)*TauKL(iseg)*chemo(Tab_U(:,ik),entity), A(k_chemo) )
        end if

        A(k_enty)%Bg(jL) = A(k_enty)%Bg(jL) + delta/AireK(jL)*TauKL(iseg)*(chemo(Tab_U(:,jL),entity)*cplusj&
        &+chemo(Tab_U(:,ik),entity)*cmoinsj)

      else 
         !! contribution dans la ligne ik
        do k=1,n_enty
         if (Tab_equa(k)==equa) then 
         CALL Ajout (ik, ik, TauKL(iseg)*chemoprime(Tab_U(:,ik),entity,Tab_entity(k))*cplusi,  A(k) )

         CALL Ajout (ik, jL, TauKL(iseg)*chemoprime(Tab_U(:,jL),entity,Tab_entity(k))*cmoinsi, A(k) )
         end if
        end do
        if (cplusi>0) then
         CALL Ajout (ik, ik, TauKL(iseg)*chemo(Tab_U(:,ik),entity), A(k_chemo) )
        else
         CALL Ajout (ik, jL, TauKL(iseg)*chemo(Tab_U(:,jL),entity), A(k_chemo) )
        end if
        A(k_enty)%Bg(ik) = A(k_enty)%Bg(ik) + TauKL(iseg)*(chemo(Tab_U(:,ik),entity)*cplusi&
        &+chemo(Tab_U(:,jL),entity)*cmoinsi)

        !! contribution dans la ligne jL
        do k=1,n_enty
         if (Tab_equa(k)==equa) then
         CALL Ajout (jL, jL, TauKL(iseg)*chemoprime(Tab_U(:,jL),entity,Tab_entity(k))*cplusj, A(k) )
         CALL Ajout (jL, ik, TauKL(iseg)*chemoprime(Tab_U(:,ik),entity,Tab_entity(k))*cmoinsj, A(k) )
         end if
        end do
        if (cplusj>0) then
         CALL Ajout (jL, jL, TauKL(iseg)*chemo(Tab_U(:,jL),entity), A(k_chemo) )
        else
         CALL Ajout (jL, ik, TauKL(iseg)*chemo(Tab_U(:,ik),entity), A(k_chemo) )
        end if
        A(k_enty)%Bg(jL) = A(k_enty)%Bg(jL) + TauKL(iseg)*(chemo(Tab_U(:,jL),entity)*cplusj&
        &+chemo(Tab_U(:,ik),entity)*cmoinsj)
      end if
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