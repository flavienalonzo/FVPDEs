SUBROUTINE assemblereactionKS( A, Tab_U,entity,Tab_entity,Tab_equa,equa,k_enty)
    !--------
    ! Modules
    !--------
    USE longr
    USE imprime
    USE parmmage
    use parameters 
  
    IMPLICIT NONE
  
    !--------------------------
    ! Declaration des arguments
    !--------------------------
    TYPE(MatCreux) , dimension(n_enty)      :: A
    real (kind = long), dimension(n_enty,Nbt) :: Tab_U
    character(len=6), dimension(n_enty) :: Tab_entity, Tab_equa
    character(len=6), intent(in) :: equa
    !----------------------------------
    ! Declaration des variables locales
    !----------------------------------
    CHARACTER(len=6)      :: oldprf,entity
    INTEGER               :: i, j, k,  jt , k_enty
    REAL(kind=long), DIMENSION(3)     :: x, y
  
  
    !-------------------
    ! Debut du programme
    !-------------------
    oldprf = prefix
    prefix = 'ASSEMT'
  
    !------
    ! Corps
    !------
    if (equa=='instat') then 
      Do jt = 1, Nbt 
        do k=1,n_enty
          if (Tab_equa(k)==equa) then
          CALL Ajout (jt, jt, - delta*reactionprime(Tab_U(:,jt),entity,Tab_entity(k)), A(k) )
          end if 
        end do
        A(k_enty)%Bg(jt) = A(k_enty)%Bg(jt) - delta*reaction(Tab_U(:,jt),entity)
      END Do
    else 
      Do jt = 1, Nbt 
        do k=1,n_enty
          if (Tab_equa(k)==equa) then
          CALL Ajout (jt, jt, - AireK(jt)*reactionprime(Tab_U(:,jt),entity,Tab_entity(k)), A(k) )
          end if
        end do
        A(k_enty)%Bg(jt) = A(k_enty)%Bg(jt) - AireK(jt)*reaction(Tab_U(:,jt),entity)
      END Do
    end if
   
    !-----------------
    ! Fin du programme
    !-----------------
    prefix = oldprf
  
    RETURN
  end subroutine assemblereactionKS