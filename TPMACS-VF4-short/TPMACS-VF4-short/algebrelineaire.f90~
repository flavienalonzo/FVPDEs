MODULE algebrelineaire
  IMPLICIT NONE

CONTAINS

  FUNCTION gradconj(b , A, tol)
    !     *--------------------------------------------------
    !     * Ce sous programme resoud le systeme lineaire 
    !     *    A x = b 
    !     * par la methode des gradients conjugues 
    !     * le resultat de l'operation x = b/A 
    !--------------------------------------------------------
    !--------------------------
    ! Declaration des arguments
    !--------------------------
    USE longr 
    USE imprime
    USE parmmage
    USE intmatvec
    IMPLICIT NONE
    TYPE(MatCreux), INTENT(in)                                :: A
    REAL(kind=long), DIMENSION( SIZE(A%IndPL) -1 ),INTENT(in) :: b
    REAL(kind=long), DIMENSION( SIZE(A%IndPL) -1 )            :: gradconj
    !----------------------------------
    ! Declaration des variables locales
    !----------------------------------
    CHARACTER(len=6)      :: oldprf
    REAL(kind=long),  DIMENSION( SIZE(A%IndPL) -1 ) ::  R, P, Q, X
    REAL(kind=long)                       :: alpha, beta, seuil, seuil9,tol
    INTEGER                                :: ngrad
    !-------------------
    ! Debut du programme
    !-------------------
    oldprf = prefix
    prefix = 'GRADCG '
    !
    R = 0.D0 ; P = 0.D0 ; Q = 0.D0
    ngrad = 0
    !--------------------------------------------------------------------
    ! Initialisation de X  ?? attention il faut donner X0 quand il le faut
    !-------------------------------------------------------------------
    X = 0.D0       
    R = b - A*X
    P = R

    seuil = DOT_PRODUCT(R, R)
    DO
       !WRITE(*,*)'GRAD  NITER = ',ngrad,' SEUIL = ', SQRT(seuil) 
       ngrad = ngrad + 1                       ! nombre d'iteration
       Q =  A* P 

       IF (SQRT(seuil) <= tol) EXIT

       alpha = seuil/DOT_PRODUCT(Q,P)
       X = X + alpha * P
       R = R - alpha * Q
       seuil9 = DOT_PRODUCT(R,R)
       beta = seuil9/seuil
       seuil = seuil9 
       P = R + beta * P

    END DO

    gradconj = X

    WRITE(*,*)'GRAD  NITER = ',ngrad,' SEUIL = ', SQRT(seuil) 
    !------------
    ! Impressions
    !------------
    !   if (iprint >= 3) then
    !      call prvarr(uprint, 'seuil',seuil)
    !   end if

    !-----------------
    ! Fin du programme
    !-----------------
    prefix = oldprf
    !
    RETURN
  END FUNCTION  gradconj

!!!!!

END MODULE algebrelineaire




