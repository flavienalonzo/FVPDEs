SUBROUTINE Newtonconditioninitiale(A,Vold, V,U,Ue, ndim, choixf,  temps)
    !--------
    ! Modules
    !--------
    USE longr
    USE imprime
    USE parmmage
    USE intmatvec
    USE intbigradc
    USE fsourcemod
    use fsourcebreast
  
    IMPLICIT NONE
  
    !--------------------------
    ! Declaration des arguments
    !--------------------------
    
    TYPE(MatCreux)                                   :: A
    Integer, intent(in)                              :: ndim, choixf
    REAL(kind=long), intent(in)                      :: temps
    REAL(kind=long), dimension (ndim), intent(in)    :: Vold, U, Ue
    REAL(kind=long), dimension (ndim), intent(out )  :: V
  
    !----------------------------------
    ! Declaration des variables locales
    !----------------------------------
  
    CHARACTER(len=6)      :: oldprf
    INTEGER               :: jt, ik, j,k,jL, iseg, ii, kiter, kitermax, is , js ,i
    REAL(kind=long)       :: CoefAjout, x1, y1 
    REAL(kind = long)     :: Adegenbord
    REAL(kind = long)     :: dVKLmoins, dVKLplus
    REAL(kind = long)     :: dVKL, dVLK, Ubord, Vbord,temp
    REAL(kind = long), dimension(size(V))            :: Xk,Yk,X0, FF
    REAL(kind = long), dimension(:), ALLOCATABLE     :: T
    REAL(kind = long), dimension(3,3)                :: Aloc
    REAL(kind = long), dimension(2)                  :: DerEtaKL
    
    !--------------------------------------------------------
    !--------------------------------------------------------
    INTERFACE
       FUNCTION amatloc(jt)
         USE longr
         IMPLICIT NONE
         INTEGER, INTENT(in)               :: jt       ! numero du triangle
         REAL(kind=long), DIMENSION(3,3)   :: amatloc
       END FUNCTION amatloc
    END INTERFACE
    !-------------------
    ! Debut du programme
    !-------------------
    oldprf = prefix
    prefix = 'ASSP1S'
  
    !------
    ! Corps
    !------
  
    kitermax = 100
    Xk = Vold
  
    DO kiter = 1, kitermax
  
       A%TMAT = 0.D0
       
       !!===============================
       !! calcul de FF(Xk) et FF'(Xk)=AV
       !!===============================
  
       !--------------------
       !--1-- terme en temps
       !--------------------
  
       DO i = 1, NsInt
          FF(i) = - AireDsommet(i)*(Coef_prod*Ue(i) - nut_degra*Xk(i) &
          & -Coef_cons*U(i)*Xk(i) )
          !! derivee de FF par rapport a Xk
          call ajout(i,i, AireDsommet(i)*(Coef_cons*U(i) + nut_degra), A ) 
       END DO
       !-----------------------------------------
       ! --2-- terme laplacien  -div(D(x) grad v)
       !-----------------------------------------
  
       DO jt =1,Nbt
          i = NuSoK(1,jt) ; j = NuSoK(2,jt) ; k = NuSoK(3,jt)
  
          Aloc = - Amatloc(jt)*CoefDiffV !! We multiply by CoefDiffV
          
          !-----------------------------
          ! contribution dans la ligne i
          !-----------------------------
          
          IF ( NtypS(i) /= 1) THEN
             IF (NtypS(j) /= 1 ) THEN 
                temp = EtaKL(Aloc(1,2), Xk(i), Xk(j))
                FF(i) = FF(i) + Aloc(1,2)*temp*(p(Xk(i)) - p(Xk(j)))
                DerEtaKL = DerivEtaKL(Aloc(1,2), Xk(i), Xk(j))
                !! derivee de FF(i) par rapport a Xk(i)
                CALL Ajout (i, i, Aloc(1,2)*( DerEtaKL(1)*(p(Xk(i)) - p(Xk(j))) + temp*Derivp(Xk(i)) ), A )
                !! derivee de FF(i) par rapport a Xk(j)
                CALL Ajout (i, j, Aloc(1,2)*( DerEtaKL(2)*(p(Xk(i)) - p(Xk(j))) - temp*Derivp(Xk(j)) ), A )
             ELSE
                Select case ( NtypS(j) ) 
                case (dirichlet+11)
                   temp = EtaKL(Aloc(1,2), Xk(i), Gb(j))
                   FF(i) = FF(i) + Aloc(1,2)*temp*(p(Xk(i)) - p(Gb(j)))
                   DerEtaKL=DerivEtaKL(Aloc(1,2), Xk(i), Gb(j))
                   CALL Ajout (i, i, Aloc(1,2)*( DerEtaKL(1)*(p(Xk(i)) - p(Gb(j))) + temp*Derivp(Xk(i))) , A )
                case(1)
                   !! On fait rien 
                End Select
             ENDIF  ! fin (i,j) 
             !
             IF (NtypS(k) /= 1 ) THEN
                temp = EtaKL(Aloc(1,3), Xk(i), Xk(k))
                FF(i) = FF(i) + Aloc(1,3)*temp*(p(Xk(i)) - p(Xk(k)))
                DerEtaKL = DerivEtaKL(Aloc(1,3), Xk(i), Xk(k))
                !! derivee de FF(i) par rapport a Xk(i)
                CALL Ajout (i, i, Aloc(1,3)*( DerEtaKL(1)*(p(Xk(i)) - p(Xk(k))) + temp*Derivp(Xk(i)) ), A )
                !! derivee de FF(i) par rapport a Xk(k)
                CALL Ajout (i, k, Aloc(1,3)*( DerEtaKL(2)*(p(Xk(i)) - p(Xk(k))) - temp*Derivp(Xk(k)) ), A )
             ELSE
                Select case ( Ntyps(k) ) 
                case (dirichlet+11)
                   temp = EtaKL(Aloc(1,3), Xk(i), Gb(k))
                   FF(i) = FF(i) + Aloc(1,3)*temp*(p(Xk(i)) - p(Gb(k)))
                   DerEtaKL = DerivEtaKL(Aloc(1,3), Xk(i), Gb(k))
                   CALL Ajout (i, i, Aloc(1,3)*( DerEtaKL(1)*(p(Xk(i)) - p(Gb(k))) + temp*Derivp(Xk(i)) ), A )
                case(1)
                   !! on fait rien 
                End Select
                
             ENDIF ! fin (i,k) interieur
          END IF  ! fin i 
          
          ! contribution dans la ligne j 
          !-----------------------------
          
          IF ( NtypS(j) /= 1) THEN
             IF (NtypS(i) /= 1 ) THEN
                temp = EtaKL(Aloc(2,1), Xk(j), Xk(i)) 
                FF(j) = FF(j) + Aloc(2,1)*temp*(p(Xk(j)) - p(Xk(i)))
                DerEtaKL = DerivEtaKL(Aloc(2,1), Xk(j), Xk(i))
                !! derivee de FF(j) par rapport a Xk(j)
                CALL Ajout (j, j, Aloc(2,1)*( DerEtaKL(1)*(p(Xk(j)) - p(Xk(i))) + temp*Derivp(Xk(j)) ), A )
                !! derivee de FF(j) par rapport a Xk(i)
                CALL Ajout (j, i, Aloc(2,1)*( DerEtaKL(2)*(p(Xk(j)) - p(Xk(i))) - temp*Derivp(Xk(i)) ), A )
             ELSE
                Select case ( Ntyps(i) ) 
                case (dirichlet+11)
                   temp = EtaKL(Aloc(2,1), Xk(j), Gb(i))
                   FF(j) = FF(j) + Aloc(2,1)*temp*(p(Xk(j)) - p(Gb(i)))
                   DerEtaKL = DerivEtaKL(Aloc(2,1), Xk(j), Gb(i))
                   !! derivee de FF(j) par rapport a Xk(j)
                   CALL Ajout (j, j, Aloc(2,1)*( DerEtaKL(1)*(p(Xk(j)) - p(Gb(i))) + temp*Derivp(Xk(j)) ), A )
                case(1)
                   !! on fait rien 
                End Select
             ENDIF  ! fin (j,i) 
             !
             IF (NtypS(k) /= 1 ) THEN
                temp = EtaKL(Aloc(2,3), Xk(j), Xk(k))
                FF(j) = FF(j) + Aloc(2,3)*temp*(p(Xk(j)) - p(Xk(k)))
                DerEtaKL = DerivEtaKL(Aloc(2,3), Xk(j), Xk(k))
                !! derivee de FF(j) par rapport a Xk(j)
                CALL Ajout (j, j, Aloc(2,3)*( DerEtaKL(1)*(p(Xk(j)) - p(Xk(k))) + temp*Derivp(Xk(j)) ), A )
                !! derivee de FF(j) par rapport a Xk(k)
                CALL Ajout (j, k, Aloc(2,3)*( DerEtaKL(2)*(p(Xk(j)) - p(Xk(k))) - temp*Derivp(Xk(k)) ), A )
             ELSE
                Select case ( Ntyps(k) ) 
                case (dirichlet+11)
                   temp = EtaKL(Aloc(2,3), Xk(j), Gb(k))
                   FF(j) = FF(j) + Aloc(2,3)*temp*(p(Xk(j)) - p(Gb(k)))
                   DerEtaKL = DerivEtaKL(Aloc(2,3), Xk(j), Gb(k))
                   !! derivee de FF(j) par rapport a Xk(j)
                   CALL Ajout (j, j, Aloc(2,3)*( DerEtaKL(1)*(p(Xk(j)) - p(Gb(k))) + temp*Derivp(Xk(j)) ), A )
                case(1)
                   !! on fait rien 
                End Select
             ENDIF ! fin (j,k) 
          END IF  ! fin j
          
          ! contribution dans la ligne k
          !-----------------------------
          
          IF ( NtypS(k) /= 1) THEN
             IF (NtypS(i) /= 1 ) THEN
                temp = EtaKL(Aloc(3,1), Xk(k), Xk(i))
                FF(k) = FF(k) + Aloc(3,1)*temp*(p(Xk(k)) - p(Xk(i)))
                DerEtaKL = DerivEtaKL(Aloc(3,1), Xk(k), Xk(i))
                !! derivee de FF(k) par rapport a Xk(k)
                CALL Ajout (k, k, Aloc(3,1)*( DerEtaKL(1)*(p(Xk(k)) - p(Xk(i))) + temp*Derivp(Xk(k)) ), A )
                !! derivee de FF(k) par rapport a Xk(i)
                CALL Ajout (k, i, Aloc(3,1)*( DerEtaKL(2)*(p(Xk(k)) - p(Xk(i))) - temp*Derivp(Xk(i)) ), A )
             ELSE
                Select case ( Ntyps(i) ) 
                case (dirichlet+11)
                   temp = EtaKL(Aloc(3,1), Xk(k), Gb(i))
                   FF(k) = FF(k) + Aloc(3,1)*temp*(p(Xk(k)) - p(Gb(i)))
                   DerEtaKL = DerivEtaKL(Aloc(3,1), Xk(k), Gb(i))
                   !! derivee de FF(k) par rapport a Xk(k)
                   CALL Ajout (k, k, Aloc(3,1)*( DerEtaKL(1)*(p(Xk(k)) - p(Gb(i))) + temp*Derivp(Xk(k)) ), A )
                case(1)
                   !! on fait rien 
                End Select
             ENDIF  ! fin (k,i) 
             !
             IF (NtypS(j) /= 1 ) THEN
                temp = EtaKL(Aloc(3,2), Xk(k), Xk(j))
                FF(k) = FF(k) + Aloc(3,2)*temp*(p(Xk(k)) - p(Xk(j)))
                DerEtaKL = DerivEtaKL(Aloc(3,2), Xk(k), Xk(j))
                !! derivee de FF(k) par rapport a Xk(k)
                CALL Ajout (k, k, Aloc(3,2)*( DerEtaKL(1)*(p(Xk(k)) - p(Xk(j))) + temp*Derivp(Xk(k)) ), A )
                !! derivee de FF(k) par rapport a Xk(j)
                CALL Ajout (k, j, Aloc(3,2)*( DerEtaKL(2)*(p(Xk(k)) - p(Xk(j))) - temp*Derivp(Xk(j)) ), A )
             ELSE
                Select case ( Ntyps(j) ) 
                case (dirichlet+11)
                   temp = EtaKL(Aloc(3,2), Xk(k), Gb(j))
                   FF(k) = FF(k) + Aloc(3,2)*temp*(p(Xk(k)) - p(Gb(j)))
                   DerEtaKL = DerivEtaKL(Aloc(3,2), Xk(k), Gb(j))
                   !! derivee de FF(k) par rapport a Xk(k)
                   CALL Ajout (k, k, Aloc(3,2)*( DerEtaKL(1)*(p(Xk(k)) - p(Gb(j))) + temp*Derivp(Xk(k)) ), A )
                case(1)
                   !! il fait après dans le cas de neumann car un point spécifique
                End Select
             ENDIF ! fin (k,j) 
          END IF  !! fin k 
       END DO
       !!
       !!============================
       !!=============================
       !!print*,'debut sys lin'
       !Yk = -FF/A
       X0= 0.D0 
       Yk = bigradient(A,-FF, X0,Tolerencegradient)
       !write(*,*)'Yk=',Yk
       !write(*,*)'FF(1)=',FF(1)
       !print*,'kiter = ', kiter,'erreur NEWTONv sqrt(sum(Yk*Yk))',sqrt(sum(Yk*Yk))
       
       If (sqrt(dot_product(Yk,Yk)) <TolerenceNewton) exit
       Xk=Xk+Yk
       !print*,'max Xk', maxval(Xk), 'min Xk', minval(Xk)
       !print*,'max FF', maxval(FF), 'min FF', minval(FF)
       
       !! impression de A etde F
       !     DO i = 1, size(A%IndPL)-1
       !           print*,'iseg = ', i, 'sum line = ', SUM(A%Tmat(A%IndPL(i): A%IndPL(i+1)-1))
       !    END DO
       
    end Do
  
    V=Xk+Yk
  
    !print*,'fin newton'
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
  
  
  END SUBROUTINE Newtonconditioninitiale