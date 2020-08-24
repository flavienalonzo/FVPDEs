subroutine assemblediffusionKS (A,Tab_U,entity,Tab_entity,Tab_equa,equa,k_enty)
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
  INTEGER               :: is, js, ik, jL, iseg, ii, jv,  NcoefMat, k_enty, k
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
        if (equa=='instat') then 
            do k=1,n_enty
               if (Tab_equa(k)==equa) then
               CALL Ajout (ik, ik, delta/AireK(ik)*TauKL(iseg)*&
               &diffprime(Tab_U(:,ik),entity,Tab_entity(k))*AireK(ik)/(AireK(ik)+AireK(jL))*Tab_U(k_enty,ik),  A(k) )
               CALL Ajout (ik, ik, -delta/AireK(ik)*TauKL(iseg)*&
               &diffprime(Tab_U(:,ik),entity,Tab_entity(k))*AireK(ik)/(AireK(ik)+AireK(jL))*Tab_U(k_enty,jL),  A(k) )
               CALL Ajout (ik, ik, delta/AireK(ik)*TauKL(iseg)*&
               &diffprime(Tab_U(:,ik),entity,Tab_entity(k))*AireK(ik)/(AireK(ik)+AireK(jL)),  A(k) )
               end if
            end do
            CALL Ajout (ik, ik, delta/AireK(ik)*TauKL(iseg)*diffusion(Tab_U(:,jL),entity)&
            &*AireK(ik)/(AireK(ik)+AireK(jL)),  A(k_enty) )

            do k=1,n_enty
               if (Tab_equa(k)==equa) then
                  CALL Ajout (ik, jL, delta/AireK(ik)*TauKL(iseg)*&
                  &diffprime(Tab_U(:,jL),entity,Tab_entity(k))*AireK(jL)/(AireK(ik)+AireK(jL))*Tab_U(k_enty,ik),A(k))
                  CALL Ajout (ik, jL, -delta/AireK(ik)*TauKL(iseg)*&
                  &diffprime(Tab_U(:,jL),entity,Tab_entity(k))*AireK(jL)/(AireK(ik)+AireK(jL))*Tab_U(k_enty,jL),A(k))
               end if
            end do
            CALL Ajout (ik, jL, -delta/AireK(ik)*TauKL(iseg)*diffusion(Tab_U(:,ik),entity)&
            &*AireK(ik)/(AireK(ik)+AireK(jL)), A(k_enty) )
            CALL Ajout (ik, jL, -delta/AireK(ik)*TauKL(iseg)*diffusion(Tab_U(:,jL),entity)&
            &*AireK(jL)/(AireK(ik)+AireK(jL)), A(k_enty) )

            A(k_enty)%Bg(ik) = A(k_enty)%Bg(ik) + delta/AireK(ik)*TauKL(iseg)/(AireK(ik)+AireK(jL))*(AireK(ik)*&
            &diffusion(Tab_U(:,ik),entity)+AireK(jL)*diffusion(Tab_U(:,jL),entity))*(Tab_U(k_enty,ik)-Tab_U(k_enty,jL))

            !! contribution dans la ligne jL

            do k=1,n_enty
               if (Tab_equa(k)==equa) then
               CALL Ajout (jL, jL, delta/AireK(jL)*TauKL(iseg)*&
               &diffprime(Tab_U(:,jL),entity,Tab_entity(k))*AireK(jL)/(AireK(ik)+AireK(jL))*Tab_U(k_enty,jL), A(k) )
               CALL Ajout (jL, jL, -delta/AireK(jL)*TauKL(iseg)*&
               &diffprime(Tab_U(:,jL),entity,Tab_entity(k))*AireK(jL)/(AireK(ik)+AireK(jL))*Tab_U(k_enty,ik), A(k) )
               CALL Ajout (jL, jL, delta/AireK(jL)*TauKL(iseg)*diffprime(Tab_U(:,jL),entity,Tab_entity(k))*&
               &AireK(jL)/(AireK(ik)+AireK(jL)), A(k) )
               end if
            end do
            CALL Ajout (jL, jL, delta/AireK(jL)*TauKL(iseg)*diffusion(Tab_U(:,ik),entity)*AireK(ik)/&
            &(AireK(ik)+AireK(jL)), A(k_enty) )

            do k=1,n_enty
               if (Tab_equa(k)==equa) then
               CALL Ajout (jL, ik, delta/AireK(jL)*TauKL(iseg)*&
               &diffprime(Tab_U(:,ik),entity,Tab_entity(k))*AireK(ik)/(AireK(ik)+AireK(jL))*Tab_U(k_enty,jL), A(k) )
               CALL Ajout (jL, ik, -delta/AireK(jL)*TauKL(iseg)*&
               &diffprime(Tab_U(:,ik),entity,Tab_entity(k))*AireK(ik)/(AireK(ik)+AireK(jL))*Tab_U(k_enty,ik), A(k) )
               end if
            end do
            CALL Ajout (jL, ik, -delta/AireK(jL)*TauKL(iseg)*diffusion(Tab_U(:,jL),entity)*AireK(jL)/&
            &(AireK(ik)+AireK(jL)), A(k_enty) )
            CALL Ajout (jL, ik, -delta/AireK(jL)*TauKL(iseg)*diffusion(Tab_U(:,ik),entity)*AireK(ik)/&
            &(AireK(ik)+AireK(jL)), A(k_enty) )

            A(k_enty)%Bg(jL) = A(k_enty)%Bg(jL) + delta/AireK(jL)*TauKL(iseg)/(AireK(ik)+AireK(jL))&
            &*(AireK(jL)*diffusion(Tab_U(:,jL),entity)+AireK(ik)*diffusion(Tab_U(:,ik),entity))*(Tab_U(k_enty,jL)-Tab_U(k_enty,ik))
      else 
         call ASSEMBLEVF4(A(k_enty))
!         do k=1,n_enty
!            if (Tab_equa(k)==equa) then
!            CALL Ajout (ik, ik, TauKL(iseg)*diffprime(Tab_U(:,ik),entity,Tab_entity(k))*AireK(ik)/(AireK(ik)+AireK(jL))&
!            &*Tab_U(k_enty,ik),  A(k) )
!            CALL Ajout (ik, ik, -TauKL(iseg)*diffprime(Tab_U(:,ik),entity,Tab_entity(k))*AireK(ik)/(AireK(ik)+AireK(jL))&
!            &*Tab_U(k_enty,jL),  A(k) )
!            CALL Ajout (ik, ik, TauKL(iseg)*diffprime(Tab_U(:,ik),entity,Tab_entity(k))*AireK(ik)/(AireK(ik)+AireK(jL)),A(k))
!            end if
!         end do
!         CALL Ajout (ik, ik, TauKL(iseg)*diffusion(Tab_U(:,jL),entity)*AireK(jL)/(AireK(ik)+AireK(jL)),  A(k_enty) )

!         do k=1,n_enty
!            if (Tab_equa(k)==equa) then
!          CALL Ajout (ik, jL, TauKL(iseg)*diffprime(Tab_U(:,jL),entity,Tab_entity(k))*AireK(jL)/(AireK(ik)+AireK(jL))*&
!          &Tab_U(k_enty,ik), A(k) )
!          CALL Ajout (ik, jL, -TauKL(iseg)*diffprime(Tab_U(:,jL),entity,Tab_entity(k))*AireK(jL)/(AireK(ik)+AireK(jL))*&
!          &Tab_U(k_enty,jL), A(k) )
!            end if
!         end do
!         CALL Ajout (ik, jL, -TauKL(iseg)*diffusion(Tab_U(:,ik),entity)*AireK(ik)/(AireK(ik)+AireK(jL)), A(k_enty) )
!         CALL Ajout (ik, jL, -TauKL(iseg)*diffusion(Tab_U(:,jL),entity)*AireK(jL)/(AireK(ik)+AireK(jL)), A(k_enty) )

!         A(k_enty)%Bg(ik) = A(k_enty)%Bg(ik) + TauKL(iseg)/(AireK(ik)+AireK(jL))*(AireK(ik)*diffusion(Tab_U(:,ik),entity)&
!         &+AireK(jL)*diffusion(Tab_U(:,jL),entity))*(Tab_U(k_enty,ik)-Tab_U(k_enty,jL))

         !! contribution dans la ligne jL

!         do k=1,n_enty
!            if (Tab_equa(k)==equa) then
!            CALL Ajout (jL, jL, TauKL(iseg)*diffprime(Tab_U(:,jL),entity,Tab_entity(k))*AireK(jL)/(AireK(ik)+AireK(jL))&
!            &*Tab_U(k_enty,jL), A(k) )
!            CALL Ajout (jL, jL, -TauKL(iseg)*diffprime(Tab_U(:,jL),entity,Tab_entity(k))*AireK(jL)/(AireK(ik)+AireK(jL))&
!            &*Tab_U(k_enty,ik), A(k) )
!            CALL Ajout (jL, jL, TauKL(iseg)*diffprime(Tab_U(:,jL),entity,Tab_entity(k))*AireK(jL)/(AireK(ik)+AireK(jL)), A(k) )
!            end if
!         end do
!         CALL Ajout (jL, jL, TauKL(iseg)*diffusion(Tab_U(:,ik),entity)*AireK(ik)/(AireK(ik)+AireK(jL)), A(k_enty) )

!         do k=1,n_enty
!            if (Tab_equa(k)==equa) then
!            CALL Ajout (jL, ik, TauKL(iseg)*diffprime(Tab_U(:,ik),entity,Tab_entity(k))*AireK(ik)/(AireK(ik)+AireK(jL))&
!            &*Tab_U(k_enty,jL), A(k) )
!            CALL Ajout (jL, ik, -TauKL(iseg)*diffprime(Tab_U(:,ik),entity,Tab_entity(k))*AireK(ik)/(AireK(ik)+AireK(jL))&
!            &*Tab_U(k_enty,ik), A(k) )
!            end if
!         end do
!         CALL Ajout (jL, ik, -TauKL(iseg)*diffusion(Tab_U(:,jL),entity)*AireK(jL)/(AireK(ik)+AireK(jL)), A(k_enty) )
!         CALL Ajout (jL, ik, -TauKL(iseg)*diffusion(Tab_U(:,ik),entity)*AireK(ik)/(AireK(ik)+AireK(jL)), A(k_enty) )

!         A(k_enty)%Bg(jL) = A(k_enty)%Bg(jL) + TauKL(iseg)/(AireK(ik)+AireK(jL))*(AireK(jL)*diffusion(Tab_U(:,jL),entity)+&
!         &AireK(ik)*diffusion(Tab_U(:,ik),entity))*(Tab_U(k_enty,jL)-Tab_U(k_enty,ik))
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


  !------------
  ! Impressions
  !------------
  !!WRITE(*,*) ( A%Tmat(is) , is=1, NcoefMat )

  !-----------------
  ! Fin du programme
  !-----------------
  prefix = oldprf

  RETURN



end subroutine assemblediffusionKS