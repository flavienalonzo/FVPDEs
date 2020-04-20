SUBROUTINE assemblereactionKS( A, U )
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
    TYPE(MatCreux)       :: A
    double precision, dimension(:) :: U 
    !----------------------------------
    ! Declaration des variables locales
    !----------------------------------
    CHARACTER(len=6)      :: oldprf
    INTEGER               :: i, j, k,  jt 
    REAL(kind=long), DIMENSION(3)     :: x, y
  
  
    !-------------------
    ! Debut du programme
    !-------------------
    oldprf = prefix
    prefix = 'ASSEMT'
  
    !------
    ! Corps
    !------
  Do jt = 1, Nbt 
    CALL Ajout (jt, jt, - delta*reactionprime(U(jt)), A )
    A%Bg(jt) = A%Bg(jt) - delta*reaction(U(jt))
  END Do
    !------------
    ! Impressions
    !------------
  !!$  IF (iprint >= 2) THEN
  !!$     CALL prvari(uprint, ' Numero triangle = ', jt) 
  !!$     WRITE (uprint,*) ( (matloc(iloc,jloc), jloc=1,3),iloc=1,3)
  !!$  END IF
  
    !-----------------
    ! Fin du programme
    !-----------------
    prefix = oldprf
  
    RETURN
  end subroutine assemblereactionKS