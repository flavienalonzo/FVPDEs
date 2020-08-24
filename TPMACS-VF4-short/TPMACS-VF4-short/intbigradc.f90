MODULE intbigradc
  USE longr 
  USE imprime
  USE parmmage
  USE intmatvec
  IMPLICIT NONE

 ! INTERFACE OPERATOR (/)
 !    MODULE PROCEDURE gradconj
 ! END INTERFACE
CONTAINS
  FUNCTION gradconj(b , A)
    !     *--------------------------------------------------
    !     * Ce sous programme resoud le systeme lineaire 
    !     *    A x = b 
    !     * par la methode des gradients conjugues 
    !     * le resultat de l'operation x = b/A 
    !--------------------------------------------------------
    !--------------------------
    ! Declaration des arguments
    !--------------------------
    TYPE(MatCreux), INTENT(in)                                :: A
    REAL(kind=long), DIMENSION( SIZE(A%IndPL) -1 ),INTENT(in) :: b
    REAL(kind=long), DIMENSION( SIZE(A%IndPL) -1 )            :: gradconj
    !----------------------------------
    ! Declaration des variables locales
    !----------------------------------
    CHARACTER(len=6)      :: oldprf
    REAL(kind=long), DIMENSION( SIZE(A%IndPL) -1 ) ::  R, P, Q, X
    REAL (kind=long)                       :: alfa, bta, seuil, seuil9
    INTEGER                                :: ngrad, maxitergrad
    !-------------------
    ! Debut du programme
    !-------------------
    oldprf = prefix
    prefix = 'GRADCG'
    !
    R = 0.D0 ; P = 0.D0 ; Q = 0.D0
    ngrad = 0
    !---------------------------
    ! Initialisation de X  ?? attention il faut donner X0 quand il le faut
    !---------------------------
    X = 0.D0       
    R = b - A*X
    P = R

    seuil = DOT_PRODUCT(R, R)
    Maxitergrad = 2000
    DO ngrad = 1, Maxitergrad
                            ! nombre d'iteration
       Q =  A* P 

       IF (SQRT(seuil) < 1.D-10) EXIT

       alfa = seuil/DOT_PRODUCT(Q,P)
       X = X + alfa * P
       R = R - alfa * Q
       seuil9 = DOT_PRODUCT(R,R)
       bta = seuil9/seuil
       seuil = seuil9 
       P = R + bta * P

    END DO
    
     gradconj = X

    WRITE(6,*)'GRAD  NITER = ',ngrad,' SEUIL = ', SQRT(seuil) 
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
!
!====================================
!
!====================================
!
!
  function bigradient(A,b, X0,tolerence) result(X)
    !     *--------------------------------------------------
    !     * Ce sous programme resoud le systeme lineaire 
    !     *    A x = b 
    !     * par la methode des bigradients conjugues 
    !     * le resultat de l'operation x = b/A 
    !--------------------------------------------------------    
    ! Modules
    !--------
 
    implicit none

    !--------------------------
    ! Declaration des arguments
    !--------------------------
    TYPE(MatCreux), INTENT(in)                                :: A
    REAL(kind=long), DIMENSION( SIZE(A%IndPL) -1 ),INTENT(in) :: b
    REAL(kind=long), DIMENSION( SIZE(A%IndPL) -1 )            :: X , X0

     real(kind=long), intent(in) :: tolerence

    !----------------------------------
    ! Declaration des variables locales
    !----------------------------------
    character(len=6) :: oldprf
    real(kind=long), dimension (size(b))   :: Rst, Prj, Q, Rstb, Prjb,Qb
    real(kind=long)  :: alfa, bta,seuil,prscal, prscal9
    integer          :: nbigr, nbitermax
    ! ----------------------------
    ! Signification des variables
    ! ---------------------------
    ! Rst : reste = B - A * X
    ! Prj : projection
    ! Q   := A * Prj
    ! Qb  := transposee (A) * Prjb
    ! Rstb:= B - transposee (A) * X
    ! prscal : designe un produit scalaire
    !-------------------
    ! Debut du programme
    !-------------------
    oldprf = prefix
    prefix = 'BIGRAD'
    !
    Rst = 0. ;  Prj = 0. ; Q = 0.
    Rstb = 0. ; Prjb = 0. ; Qb = 0.
    nbigr = 0
    !
    X= X0
    Rst = B - A*X
    Rstb = B - transposee_mat_vec(A,X)
    Prj = Rst
    Prjb = Rstb
    prscal = DOT_PRODUCT(Rstb,Rst)

    nbitermax = size(X0) + 500

    DO  nbigr = 1, nbitermax
       seuil = sqrt(DOT_PRODUCT(Rst, Rst))

       if (seuil < tolerence) exit

       Q =  A*Prj
       alfa = prscal/DOT_PRODUCT(Q,Prjb)
       X = X + alfa * Prj
       Rst = Rst - alfa * Q
       Qb = transposee_mat_vec(A,Prjb) 
       Rstb = Rstb - alfa * Qb
       prscal9 = DOT_PRODUCT(Rstb,Rst)
       bta = prscal9/prscal
       prscal = prscal9 
       Prj = Rst + bta * Prj
       Prjb = Rstb + bta * Prjb
       !================
       ! Impression
       !================
!!$       write(6,*)'nbigr =', nbigr
!!$       write(6,*)'prscal=',prscal
!!$       write(6,*)'alfa=',alfa
!!$       write(6,*)'bta=',bta
!!$       write(6,*)'seuil=',seuil , 'tolerence = ', tolerence   
    END DO
    !write(6,*)'nbigr =', nbigr, 'seuil bigrad=',seuil
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

  end function bigradient






!====================================
!
function BGbloc(Tab_A, B, Tab_U,tolerence,tps) result(Tab_X)
    !     *--------------------------------------------------
    !     * Ce sous programme resoud le systeme lineaire 
    !     *    A x = b 
    !     * par la methode des bigradients conjugues 
    !     * La matrice A est formé par bloc
    !       ! A11 ...  A1n !
    !    A= ! ...      ... !
    !       ! An1 ...  Ann !
    !    b=(b1, ..., bn)
    !--------------------------------------------------------    
    ! Modules
    !--------
 
    implicit none

    !--------------------------
    ! Declaration des arguments
    !--------------------------
    TYPE(MatCreux), dimension(n_enty,n_enty), INTENT(in)      :: Tab_A
    real(kind=long), dimension(n_enty*Nbt), intent(in) :: B
    real(kind=long), dimension(n_enty,Nbt), intent(in) :: Tab_U
    REAL(kind=long), DIMENSION( n_enty*Nbt )  :: Tab_X,Tab_X0
    real(kind=long), intent(in) :: tolerence, tps

    !----------------------------------
    ! Declaration des variables locales
    !----------------------------------
    character(len=6) :: oldprf
    real(kind=long), dimension (n_enty*Nbt)   :: Rst,Prj,Q,Rstb,Prjb,Qb
    real(kind=long)  :: alfa, bta,seuil,prscal, prscal9, condi
    integer          :: nbigr, nbitermax, k, h
    ! ----------------------------
    ! Signification des variables
    ! ---------------------------
    ! Rst : reste = B - A * X
    ! Prj : projection
    ! Q   := A * Prj
    ! Qb  := transposee (A) * Prjb
    ! Rstb:= B - transposee (A) * X
    ! prscal : designe un produit scalaire
    !-------------------
    ! Debut du programme
    !-------------------
    oldprf = prefix
    prefix = 'BIGRAD'
    ! 
    Rst = 0.D0  ;  Prj  = 0.D0 ; Q  = 0.D0
    Rstb = 0.D0 ;  Prjb = 0.D0 ; Qb = 0.D0
    condi = 1.D0
    if (tps<2.D-1.and.tps>3.05D-1) then 
      Tab_X0 = 1.D0
    else 
      Tab_X0 = reshape(Tab_U,(/n_enty*Nbt/))
    end if
    nbigr = 0
    !
    Tab_X = Tab_X0
    do k=1,n_enty
      Rst(1+(k-1)*Nbt:k*Nbt) = B(1+(k-1)*Nbt:k*Nbt)
      Rstb(1+(k-1)*Nbt:k*Nbt) = B(1+(k-1)*Nbt:k*Nbt)
      do h=1,n_enty
        Rst(1+(k-1)*Nbt:k*Nbt) = Rst(1+(k-1)*Nbt:k*Nbt) - Tab_A(k,h)*Tab_X(1+(h-1)*Nbt:h*Nbt)
        Rstb(1+(k-1)*Nbt:k*Nbt) = Rstb(1+(k-1)*Nbt:k*Nbt) -&
        & transposee_mat_vec(Tab_A(k,h),Tab_X(1+(h-1)*Nbt:h*Nbt))
      end do
    end do
    !Rst1  = b1 - A11*X1-A12*X2 ;  Rst2 = b2 - A21*X1-A22*X2 
    !Rstb1 = b1 - transposee_mat_vec(A11,X1)-transposee_mat_vec(A21,X2)
    !Rstb2 = b2 - transposee_mat_vec(A12,X1)-transposee_mat_vec(A22,X2)
    Prj = Rst*condi; Prjb = Rstb*condi
    prscal = DOT_PRODUCT(Rstb,Rst)*condi
    !Prj1   = Rst1  ;  Prj2 = Rst2   
    !Prjb1  = Rstb1 ;  Prjb2= Rstb2 
    !prscal = DOT_PRODUCT(Rstb1,Rst1) + DOT_PRODUCT(Rstb2,Rst2)

    nbitermax = 200/n_enty*size(Tab_X0) + 500

    DO  nbigr = 1, nbitermax
      seuil = sqrt( DOT_PRODUCT(Rst,Rst) )
       !seuil = sqrt( DOT_PRODUCT(Rst1, Rst1) + DOT_PRODUCT(Rst2, Rst2) )
      if (seuil < tolerence) exit
      do k=1,n_enty
        Q(1+(k-1)*Nbt:k*Nbt)=0.D0
        do h=1,n_enty
          Q(1+(k-1)*Nbt:k*Nbt) = Q(1+(k-1)*Nbt:k*Nbt) + Tab_A(k,h)*Prj(1+(h-1)*Nbt:h*Nbt)
        end do
      end do
       !Q1 =  A11*Prj1 + A12*Prj2 
       !Q2 =  A21*Prj1 + A22*Prj2 
      alfa = prscal/(DOT_PRODUCT(Q,Prjb))
      ! alfa = prscal/(DOT_PRODUCT(Q1,Prjb1)+DOT_PRODUCT(Q2,Prjb2))
      Tab_X = Tab_X + alfa*Prj
      ! X1 = X1 + alfa * Prj1 
      ! X2 = X2 + alfa * Prj2
      Rst = Rst - alfa*Q
      ! Rst1 = Rst1 - alfa * Q1
      ! Rst2 = Rst2 - alfa * Q2
      do k=1,n_enty
        Qb(1+(k-1)*Nbt:k*Nbt) = 0.D0
        do h=1,n_enty
          Qb(1+(k-1)*Nbt:k*Nbt) = Qb(1+(k-1)*Nbt:k*Nbt) + transposee_mat_vec(Tab_A(k,h),Prjb(1+(h-1)*Nbt:h*Nbt))
        end do
      end do 
      ! Qb1 = transposee_mat_vec(A11,Prjb1) +  transposee_mat_vec(A21,Prjb2) 
      ! Qb2 = transposee_mat_vec(A12,Prjb1) +  transposee_mat_vec(A22,Prjb2) 
      Rstb = Rstb - alfa*Qb
      ! Rstb1 = Rstb1 - alfa * Qb1
      ! Rstb2 = Rstb2 - alfa * Qb2
      prscal9 = DOT_PRODUCT(Rstb,Rst)*condi
      ! prscal9 = DOT_PRODUCT(Rstb1,Rst1) + DOT_PRODUCT(Rstb2,Rst2)
      bta = prscal9/prscal
      prscal = prscal9 
      Prj = condi*Rst + bta * Prj 
      ! Prj1 = Rst1 + bta * Prj1
      ! Prj2 = Rst2 + bta * Prj2
      Prjb = condi*Rstb + bta * Prjb 
      ! Prjb1 = Rstb1 + bta * Prjb1
      ! Prjb2 = Rstb2 + bta * Prjb2
       !================
       ! Impression
       !================
  !!$       write(6,*)'nbigr =', nbigr
  !!$       write(6,*)'prscal=',prscal
  !!$       write(6,*)'alfa=',alfa
  !!$       write(6,*)'bta=',bta
  !!$       write(6,*)'seuil=',seuil , 'tolerence = ', tolerence   
  END DO
    !write(6,*)'nbigr =', nbigr, 'seuil bigrad=',seuil
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

end function BGbloc

function BiCGSTAB(Tab_A,B,tolerence) result(Tab_X)
implicit none
type(MatCreux), dimension(n_enty,n_enty), intent(in) :: Tab_A
real(kind=long), dimension(n_enty*Nbt), intent(in) :: B
real(kind=long), intent(in) :: tolerence
real(kind=long), dimension(n_enty*Nbt) :: Tab_X

real(kind=long) :: rho, alpha, omega , beta
real(kind=long), dimension(n_enty*Nbt) :: r0,r, v, p, s, t, g
integer :: n,k,h,n_max

n= 0 ; n_max =200*Nbt + 500
rho = 1.D0 ; alpha = 1.D0 ; omega = 1.D0
v = 0.D0 ; p = 0.D0
Tab_X = 1.D0
r0 = 1.D0
r = B
!Initialisation :
do k=1,n_enty
  do h=1,n_enty
    r(1+(k-1)*Nbt : k*Nbt) = r(1+(k-1)*Nbt : k*Nbt) - Tab_A(k,h)*Tab_X(1+(h-1)*Nbt:h*Nbt)
  end do
end do
do while(n<=n_max)
  beta = (DOT_PRODUCT(r,r0)/rho)*(alpha/omega)
  rho = DOT_PRODUCT(r,r0)
  p = r + beta*(p-omega*v)
  v = 0.D0
  do k=1,n_enty
    do h=1,n_enty
      v(1+(k-1)*Nbt:k*Nbt) = v(1+(k-1)*Nbt:k*Nbt) + Tab_A(k,h)*p(1+(h-1)*Nbt:h*Nbt)
    end do
  end do
  alpha = rho/DOT_PRODUCT(r0,v)
  g = Tab_X + alpha*p
  if (DOT_PRODUCT(g,g)<tolerence/1.D5) then
    Tab_X = g
    exit
  end if 
  s = r -alpha*v
  t = 0.D0
  do k=1,n_enty
    do h=1,n_enty
      t(1+(k-1)*Nbt:k*Nbt) = t(1+(k-1)*Nbt:k*Nbt)+Tab_A(k,h)*s(1+(h-1)*Nbt:h*Nbt)
    end do
  end do
  omega =DOT_PRODUCT(t,s)/DOT_PRODUCT(t,t)
  Tab_X = g + omega*s
  if (DOT_PRODUCT(omega*s,omega*s)<tolerence) then
    !print*,n
    exit
  end if
  r = s-omega*t
  n=n+1
  !print*,n,DOT_PRODUCT(omega*s,omega*s)
end do

end function BiCGSTAB




END MODULE intbigradc
