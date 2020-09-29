SUBROUTINE NewtonEndP1sommets(A, F, Uold, U, S, AA, ndim, choixf, temps)
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
    REAL(kind=long), dimension (ndim), intent(in)    :: F, Uold, S, AA
    REAL(kind=long), dimension (ndim), intent(out )  :: U
  
    !----------------------------------
    ! Declaration des variables locales
    !----------------------------------
    CHARACTER(len=6)      :: oldprf
    INTEGER               :: jt, ik, j, k, jL, iseg, ii, kiter, kitermax, is , js ,i
    REAL(kind=long)       :: CoefAjout, x1, y1 
    REAL(kind = long)     :: Adegenbord
    REAL(kind = long)     :: dVKLmoins, dVKLplus
    REAL(kind = long)     :: dVKL, dVLK, Ubord, Fbord, temp
    REAL(kind = long), dimension(size(F))            :: Xk, Yk, X0, FF
    REAL(kind = long), dimension(:), ALLOCATABLE     :: T
    REAL(kind = long), dimension(3,3)                :: Aloc
    REAL(kind = long), dimension(2)                  :: DerAKL
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
    Xk = Uold
    !print*,'max Xk initial',maxval(Xk)
    DO kiter = 1, kitermax
  
       A%TMAT = 0.D0
  
       !!=============================
       !! calcul de FF(Xk) et FF'(Xk)=AU
       !!=============================
       !
       !--------------------
       !--1-- terme en temps
       !--------------------
       !

       DO i = 1, NsInt
          FF(i) = AireDsommet(i)*( Xk(i)-Uold(i) )/dt -AireDsommet(i)*( rate_endo*Xk(i)*(1-Xk(i)) - degr_endo*Xk(i) )
          !! derivee de FF par rapport a Xk
          call ajout(i,i, AireDSommet(i)/dt, A )
          call ajout(i,i, - AireDsommet(i)*rate_endo*(1-2*Xk(i)), A )
          call ajout(i,i, AireDSommet(i)*degr_endo,A)
       END DO
       
 
       !--------------------------------------------
       ! --2-- terme laplacien  -div(S(x) grad A(u))
       !--------------------------------------------
      ! print*,'max FF apres terme en temps',maxval(FF)
       DO jt =1,Nbt
          i = NuSoK(1,jt) ; j = NuSoK(2,jt) ; k = NuSoK(3,jt)
  
          Aloc = - amatloc(jt) 
  
          ! contribution dans la ligne i
          !-----------------------------
         
          IF ( Ntyps(i) /= 1) THEN
             IF (Ntyps(j) /= 1 ) THEN 
                temp = AKL(Aloc(1,2), Xk(i), Xk(j))
                FF(i) = FF(i) + Aloc(1,2)*temp*(p(Xk(i)) - p(Xk(j)))
                DerAKL = DerivAKL(Aloc(1,2), Xk(i), Xk(j))
                !! derivee de FF(i) par rapport a Xk(i)
                CALL Ajout (i, i, Aloc(1,2)*( DerAKL(1)*(p(Xk(i)) - p(Xk(j))) + temp*Derivp(Xk(i)) ), A )
                !! derivee de FF(i) par rapport a Xk(j)
                CALL Ajout (i, j, Aloc(1,2)*( DerAKL(2)*(p(Xk(i)) - p(Xk(j))) - temp*Derivp(Xk(j)) ), A )
             ELSE
                Select case ( NtypS(j) ) 
                case (dirichlet)
                   temp = AKL(Aloc(1,2), Xk(i), Gb(j))
                   FF(i) = FF(i) + Aloc(1,2)*temp*(p(Xk(i)) - p(Gb(j)))
                   DerAKL=DerivAKL(Aloc(1,2), Xk(i), Gb(j))
                   CALL Ajout (i, i, Aloc(1,2)*( DerAKL(1)*(p(Xk(i)) - p(Gb(j))) + temp*Derivp(p(Xk(i)))) , A )
                case(neumann)
                   !! On fait rien 
                End Select
             ENDIF  ! fin (i,j) 
             !
             IF (Ntyps(k) /= 1 ) THEN
                temp = AKL(Aloc(1,3), Xk(i), Xk(k))
                FF(i) = FF(i) + Aloc(1,3)*temp*(p(Xk(i)) - p(Xk(k)))
                DerAKL = DerivAKL(Aloc(1,3), Xk(i), Xk(k))
                !! derivee de FF(i) par rapport a Xk(i)
                CALL Ajout (i, i, Aloc(1,3)*( DerAKL(1)*(p(Xk(i)) - p(Xk(k))) + temp*Derivp(Xk(i)) ), A )
                !! derivee de FF(i) par rapport a Xk(k)
                CALL Ajout (i, k, Aloc(1,3)*( DerAKL(2)*(p(Xk(i)) - p(Xk(k))) - temp*Derivp(Xk(k)) ), A )
             ELSE
                Select case ( Ntyps(k) ) 
                case (dirichlet+11)
                   temp = AKL(Aloc(1,3), Xk(i), Gb(k))
                   FF(i) = FF(i) + Aloc(1,3)*temp*(p(Xk(i)) - p(Gb(k)))
                   DerAKL = DerivAKL(Aloc(1,3), Xk(i), Gb(k))
                   CALL Ajout (i, i, Aloc(1,3)*( DerAKL(1)*(p(Xk(i)) - p(Gb(k))) + temp*Derivp(Xk(i)) ), A )
                case(1)
                   !! on fait rien 
                End Select
  
             ENDIF ! fin (i,k) interieur
          END IF  ! fin i 
          !print*,'max FF apres i',maxval(FF)
          ! contribution dans la ligne j 
          !-----------------------------
          
          IF ( Ntyps(j) /= 1) THEN
             IF (Ntyps(i) /= 1 ) THEN
                temp = AKL(Aloc(2,1), Xk(j), Xk(i)) 
                FF(j) = FF(j) + Aloc(2,1)*temp*(p(Xk(j)) - p(Xk(i)))
                DerAKL = DerivAKL(Aloc(2,1), Xk(j), Xk(i))
                !! derivee de FF(j) par rapport a Xk(j)
                CALL Ajout (j, j, Aloc(2,1)*( DerAKL(1)*(p(Xk(j)) - p(Xk(i))) + temp*Derivp(Xk(j)) ), A )
                !! derivee de FF(j) par rapport a Xk(i)
                CALL Ajout (j, i, Aloc(2,1)*( DerAKL(2)*(p(Xk(j)) - p(Xk(i))) - temp*Derivp(Xk(i)) ), A )
             ELSE
                Select case ( Ntyps(i) ) 
                case (dirichlet+11)
                   temp = AKL(Aloc(2,1), Xk(j), Gb(i))
                   FF(j) = FF(j) + Aloc(2,1)*temp*(p(Xk(j)) - p(Gb(i)))
                   DerAKL = DerivAKL(Aloc(2,1), Xk(j), Gb(i))
                   !! derivee de FF(j) par rapport a Xk(j)
                   CALL Ajout (j, j, Aloc(2,1)*( DerAKL(1)*(p(Xk(j)) - p(Gb(i))) + temp*Derivp(Xk(j)) ), A )
                case(1)
                   !! on fait rien 
                End Select
             ENDIF  ! fin (j,i) 
             !
             IF (Ntyps(k) /= 1 ) THEN
                temp = AKL(Aloc(2,3), Xk(j), Xk(k))
                FF(j) = FF(j) + Aloc(2,3)*temp*(p(Xk(j)) - p(Xk(k)))
                DerAKL = DerivAKL(Aloc(2,3), Xk(j), Xk(k))
                !! derivee de FF(j) par rapport a Xk(j)
                CALL Ajout (j, j, Aloc(2,3)*( DerAKL(1)*(p(Xk(j)) - p(Xk(k))) + temp*Derivp(Xk(j)) ), A )
                !! derivee de FF(j) par rapport a Xk(k)
                CALL Ajout (j, k, Aloc(2,3)*( DerAKL(2)*(p(Xk(j)) - p(Xk(k))) - temp*Derivp(Xk(k)) ), A )
             ELSE
                Select case ( Ntyps(k) ) 
                case (dirichlet+11)
                   temp = AKL(Aloc(2,3), Xk(j), Gb(k))
                   FF(j) = FF(j) + Aloc(2,3)*temp*(p(Xk(j)) - p(Gb(k)))
                   DerAKL = DerivAKL(Aloc(2,3), Xk(j), Gb(k))
                   !! derivee de FF(j) par rapport a Xk(j)
                   CALL Ajout (j, j, Aloc(2,3)*( DerAKL(1)*(p(Xk(j)) - p(Gb(k))) + temp*Derivp(Xk(j)) ), A )
                case(1)
                   !! on fait rien 
                End Select
             ENDIF ! fin (j,k) 
          END IF  ! fin j
          !print*,'max FF apres j',maxval(FF)
          ! contribution dans la ligne k
          !-----------------------------
          
          IF ( Ntyps(k) /= 1) THEN
             IF (Ntyps(i) /= 1 ) THEN
                temp = AKL(Aloc(3,1), Xk(k), Xk(i))
                FF(k) = FF(k) + Aloc(3,1)*temp*(p(Xk(k)) - p(Xk(i)))
                DerAKL = DerivAKL(Aloc(3,1), Xk(k), Xk(i))
                !! derivee de FF(k) par rapport a Xk(k)
                CALL Ajout (k, k, Aloc(3,1)*( DerAKL(1)*(p(Xk(k)) - p(Xk(i))) + temp*Derivp(Xk(k)) ), A )
                !! derivee de FF(k) par rapport a Xk(i)
                CALL Ajout (k, i, Aloc(3,1)*( DerAKL(2)*(p(Xk(k)) - p(Xk(i))) - temp*Derivp(Xk(i)) ), A )
             ELSE
                Select case ( Ntyps(i) ) 
                case (dirichlet+11)
                   temp = AKL(Aloc(3,1), Xk(k), Gb(i))
                   FF(k) = FF(k) + Aloc(3,1)*temp*(p(Xk(k)) - p(Gb(i)))
                   DerAKL = DerivAKL(Aloc(3,1), Xk(k), Gb(i))
                   !! derivee de FF(k) par rapport a Xk(k)
                   CALL Ajout (k, k, Aloc(3,1)*( DerAKL(1)*(p(Xk(k)) - p(Gb(i))) + temp*Derivp(Xk(k)) ), A )
                case(1)
                   !! on fait rien 
                End Select
             ENDIF  ! fin (k,i) 
             !
             IF (Ntyps(j) /= 1 ) THEN
                temp = AKL(Aloc(3,2), Xk(k), Xk(j))
                FF(k) = FF(k) + Aloc(3,2)*temp*(p(Xk(k)) - p(Xk(j)))
                DerAKL = DerivAKL(Aloc(3,2), Xk(k), Xk(j))
                !! derivee de FF(k) par rapport a Xk(k)
                CALL Ajout (k, k, Aloc(3,2)*( DerAKL(1)*(p(Xk(k)) - p(Xk(j))) + temp*Derivp(Xk(k)) ), A )
                !! derivee de FF(k) par rapport a Xk(j)
                CALL Ajout (k, j, Aloc(3,2)*( DerAKL(2)*(p(Xk(k)) - p(Xk(j))) - temp*Derivp(Xk(j)) ), A )
             ELSE
                Select case ( Ntyps(j) ) 
                case (dirichlet+11)
                   temp = AKL(Aloc(3,2), Xk(k), Gb(j))
                   FF(k) = FF(k) + Aloc(3,2)*temp*(p(Xk(k)) - p(Gb(j)))
                   DerAKL = DerivAKL(Aloc(3,2), Xk(k), Gb(j))
                   !! derivee de FF(k) par rapport a Xk(k)
                   CALL Ajout (k, k, Aloc(3,2)*( DerAKL(1)*(p(Xk(k)) - p(Gb(j))) + temp*Derivp(Xk(k)) ), A )
                case(1)
                   !! il fait après dans le cas de dirichlet car un point spécifique
                End Select
             ENDIF ! fin (k,j) 
          END IF  !! fin k 
       END DO
       !print*,'max FF apres k',maxval(FF)
       !-------------------------------
       ! ajout des termes de convection
       !-------------------------------
      
       !----------------------------
       !contribution dans la ligne i
       !----------------------------
  
       DO jt = 1,Nbt
        
          Aloc =  - amatloc(jt)  
          
          i = NuSoK(1,jt) ; j = NuSoK(2,jt) ; k = NuSoK(3,jt)
                  
          IF ( Ntyps(i) /= 1) THEN       
             IF (Ntyps(j) /= 1 ) THEN 
                dVKL = Aloc(1,2)* (F(j)-F(i))
                dVKLplus= max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0) 
  
                temp = AKL(Aloc(1,2), Xk(i), Xk(j))
                DerAKL = DerivAKL(Aloc(1,2), Xk(i), Xk(j))

                FF(i) = FF(i) + temp*( dVKLplus*( MuKLCroit(Xk(i)) + MuKLDecroit(Xk(j)) ) + &
                     & dVKLmoins*( MuKLCroit(Xk(j)) + MuKLDecroit(Xk(i))) )
                !! derivee de FF(i) par rapport a Xk(i)
                CoefAjout = DerAKL(1)*( dVKLplus*( MuKLCroit(Xk(i)) + MuKLDecroit(Xk(j)) ) + &
                     & dVKLmoins*( MuKLCroit(Xk(j)) + MuKLDecroit(Xk(i))) )+ &
                     & temp*(dVKLplus*DerivMuKLCroit(Xk(i)) + dVKLmoins*DerivMuKLDecroit(Xk(i)))
                call ajout(i,i, CoefAjout, A)
                !! derivee de FF(i) par rapport a Xk(j)
                CoefAjout = DerAKL(2)*( dVKLplus*( MuKLCroit(Xk(i)) + MuKLDecroit(Xk(j)) ) + &
                     & dVKLmoins*( MuKLCroit(Xk(j)) + MuKLDecroit(Xk(i))) )+&
                     & temp*( dVKLplus*DerivMuKLDecroit(Xk(j)) + dVKLmoins*DerivMuKLCroit(Xk(j)) )
                call ajout(i,j, CoefAjout, A)
             ElSE
                Select case ( Ntyps(j) ) 
                case (dirichlet+11) 
                   Ubord= Gb(j) ; Fbord = Gb(j)
                   dVKL = Aloc(1,2)* ( Fbord- F(i)) 
                   dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)   
                   FF(i)=FF(i)+ dVKLplus*( MuKLCroit(Xk(i)) + MuKLDecroit(Ubord) ) +&
                        & dVKLmoins*( MuKLCroit(Ubord) + MuKLDecroit(Xk(i)))
                   !! derivee de FF(i) par rapport a Xk(i)
                   CoefAjout = dVKLplus*DerivMuKLCroit(Xk(i)) + dVKLmoins*DerivMuKLDecroit(Xk(i))
                   call ajout(i,i, CoefAjout, A)
                case(1)
                   !! on fait rien  
                End Select
             END IF
             !!
             IF (Ntyps(k) /= 1 ) THEN 
                dVKL = Aloc(1,3)* (F(k)-F(i))
                dVKLplus= max(dVKL, 0.D0) ; dVKLmoins=min(dVKL, 0.D0)
                !
                temp = AKL(Aloc(1,3), Xk(i), Xk(k))
                DerAKL = DerivAKL(Aloc(1,3), Xk(i), Xk(k))
  
                FF(i) = FF(i) + temp*( dVKLplus*( MuKLCroit(Xk(i)) + MuKLDecroit(Xk(k)) ) + &
                     & dVKLmoins*( MuKLCroit(Xk(k)) + MuKLDecroit(Xk(i))) )
                !! derivee de FF(i) par rapport a Xk(i)
                CoefAjout = DerAKL(1)*( dVKLplus*( MuKLCroit(Xk(i)) + MuKLDecroit(Xk(k)) ) + &
                     & dVKLmoins*( MuKLCroit(Xk(k)) + MuKLDecroit(Xk(i))) )+&
                     & temp*(dVKLplus*DerivMuKLCroit(Xk(i)) + dVKLmoins*DerivMuKLDecroit(Xk(i)))
                call ajout(i,i, CoefAjout, A)
                !! derivee de FF(i) par rapport a Xk(k)
                CoefAjout = DerAKL(2)*( dVKLplus*( MuKLCroit(Xk(i)) + MuKLDecroit(Xk(k)) ) + &
                     & dVKLmoins*( MuKLCroit(Xk(k))+MuKLDecroit(Xk(i))) )+&
                     & temp*( dVKLplus*DerivMuKLDecroit(Xk(k)) + dVKLmoins*DerivMuKLCroit(Xk(k)) )
                call ajout(i,k, CoefAjout, A)
             Else
                Select case ( Ntyps(k) ) 
                case (dirichlet+11)   
                   Ubord = Gb(k) ; Fbord = Gb(k)
                   dVKL = Aloc(1,3)* (Fbord-F(i))  
                   dVKLplus= max(dVKL, 0.D0) ; dVKLmoins=min(dVKL, 0.D0)
                   FF(i) = FF(i)+ dVKLplus*( MuKLCroit(Xk(i) )+MuKLDecroit(Ubord) ) +&
                        & dVKLmoins*( MuKLCroit(Ubord )+MuKLDecroit(Xk(i)) )
                   !! derivee de FF(i) par rapport a Xk(i)
                   CoefAjout = dVKLplus*DerivMuKLCroit(Xk(i)) + dVKLmoins*DerivMuKLDecroit(Xk(i))
                   call ajout(i,i, CoefAjout, A)
                case(1)
                   !! rien
                End Select
             END IF
          END IF
   ! fin i 
          !print*,'max FF apres i',maxval(FF)
          !------------------------------
          !  Contribution dans la ligne j
          !------------------------------
  
          IF ( Ntyps(j) /= 1) THEN 
             IF (Ntyps(i) /= 1 ) THEN 
                dVKL = Aloc(2,1)* (F(i)-F(j))
                dVKLplus= max(dVKL, 0.D0) ; dVKLmoins=min(dVKL, 0.D0)
                
                temp = AKL(Aloc(2,1), Xk(j), Xk(i))
                DerAKL = DerivAKL(Aloc(2,1), Xk(j), Xk(i))
  
                FF(j) = FF(j) + temp*( dVKLplus*( MuKLCroit(Xk(j)) + MuKLDecroit(Xk(i)) ) + &
                     & dVKLmoins*( MuKLCroit(Xk(i)) + MuKLDecroit(Xk(j))) )
                !! derivee de FF(j) par rapport a Xk(j)
                CoefAjout = DerAKL(1)*( dVKLplus*( MuKLCroit(Xk(j)) + MuKLDecroit(Xk(i)) ) + &
                     & dVKLmoins*( MuKLCroit(Xk(i)) + MuKLDecroit(Xk(j))) ) + &
                     & temp*(dVKLplus*DerivMuKLCroit(Xk(j)) + dVKLmoins*DerivMuKLDecroit(Xk(j)))
                call ajout(j,j, CoefAjout, A)
                !! derivee de FF(j) par rapport a Xk(i)
                CoefAjout = DerAKL(2)*( dVKLplus*( MuKLCroit(Xk(j)) + MuKLDecroit(Xk(i)) ) + &
                     & dVKLmoins*( MuKLCroit(Xk(i)) + MuKLDecroit(Xk(j))) )+&
                     & temp*( dVKLplus*DerivMuKLDecroit(Xk(i)) + dVKLmoins*DerivMuKLCroit(Xk(i)) )
                call ajout(j,i, CoefAjout, A)
             Else
                Select case ( Ntyps(i) ) 
                Case (dirichlet+11)   
                   Ubord = Gb(i) ; Fbord = Gb(i)
                   dVKL = Aloc(2,1)* (Fbord-F(j))
                   !  
                   dVKLplus= max(dVKL, 0.D0) ; dVKLmoins=min(dVKL, 0.D0) 
                   FF(j) =FF(j)+ dVKLplus*( MuKLCroit(Xk(j)) + MuKLDecroit(Ubord) ) +&
                        & dVKLmoins*( MuKLCroit(Ubord) + MuKLDecroit(Xk(j)) )
                   !! derivee de FF(j) par rapport a Xk(j)
                   CoefAjout = dVKLplus*DerivMuKLCroit(Xk(j)) + dVKLmoins*DerivMuKLDecroit(Xk(j))
                   call ajout(j,j, CoefAjout, A)
                Case(1)
                   !! rien
                End Select
  
             END IF
             IF (Ntyps(k) /= 1 ) THEN 
                dVKL = Aloc(2,3)* (F(k)-F(j)) 
                dVKLplus= max(dVKL, 0.D0) ; dVKLmoins=min(dVKL, 0.D0)
  
                temp = AKL(Aloc(2,3), Xk(j), Xk(k))
                DerAKL = DerivAKL(Aloc(2,3), Xk(j), Xk(k))
  
                FF(j) = FF(j) + temp*( dVKLplus*( MuKLCroit(Xk(j)) + MuKLDecroit(Xk(k)) ) + &
                     & dVKLmoins*( MuKLCroit(Xk(k)) + MuKLDecroit(Xk(j))) )
                !! derivee de FF(j) par rapport a Xk(j)
                CoefAjout = DerAKL(1)*( dVKLplus*( MuKLCroit(Xk(j)) + MuKLDecroit(Xk(k)) ) + &
                     & dVKLmoins*( MuKLCroit(Xk(k))+MuKLDecroit(Xk(j))) )+&
                     & temp*(dVKLplus*DerivMuKLCroit(Xk(j)) + dVKLmoins*DerivMuKLDecroit(Xk(j)))
                call ajout(j,j, CoefAjout, A)
                !! derivee de FF(j) par rapport a Xk(k)
                CoefAjout = DerAKL(2)*( dVKLplus*( MuKLCroit(Xk(j)) + MuKLDecroit(Xk(k)) ) + &
                     & dVKLmoins*( MuKLCroit(Xk(k)) + MuKLDecroit(Xk(j))) )+&
                     & temp*( dVKLplus*DerivMuKLDecroit(Xk(k)) + dVKLmoins*DerivMuKLCroit(Xk(k)) )
                call ajout(j,k, CoefAjout, A)
             Else
                Select case ( Ntyps(k) ) 
                case (dirichlet+11)   
                   Ubord = Gb(k) ;  Fbord = Gb(k)
                   dVKL = Aloc(2,3)* (Fbord-F(j))  
                   dVKLplus= max(dVKL, 0.D0) ; dVKLmoins=min(dVKL, 0.D0) 
                   FF(j) =FF(j)+ dVKLplus*( MuKLCroit(Xk(j)) + MuKLDecroit(Ubord) ) +&
                        & dVKLmoins*( MuKLCroit(Ubord) + MuKLDecroit(Xk(j)) )
                   !! derivee de FF(j) par rapport a Xk(j)
                   CoefAjout = dVKLplus*DerivMuKLCroit(Xk(j)) + dVKLmoins*DerivMuKLDecroit(Xk(j))
                   call ajout(j,j, CoefAjout, A)
                case(1)
                   !! rien
                End Select
  
             END IF
          END IF  ! fin j
         ! print*,'max FF apres j',maxval(FF)
  
          !-----------------------------
          ! Contribution dans la ligne k
          !-----------------------------
  
          IF ( Ntyps(k) /= 1) THEN 
             IF (Ntyps(i) /= 1 ) THEN 
                dVKL = Aloc(3,1)* (F(i)-F(k)) 
                dVKLplus= max(dVKL, 0.D0) ; dVKLmoins=min(dVKL, 0.D0)
  
                temp = AKL(Aloc(3,1), Xk(k), Xk(i))
                DerAKL = DerivAKL(Aloc(3,1), Xk(k), Xk(i))
  
                FF(k) = FF(k) + temp*( dVKLplus*( MuKLCroit(Xk(k)) + MuKLDecroit(Xk(i)) ) + &
                     & dVKLmoins*( MuKLCroit(Xk(i)) + MuKLDecroit(Xk(k))) )
                !! derivee de FF(k) par rapport a Xk(k)
                CoefAjout = DerAKL(1)*( dVKLplus*( MuKLCroit(Xk(k)) + MuKLDecroit(Xk(i)) ) + &
                     & dVKLmoins*( MuKLCroit(Xk(i)) + MuKLDecroit(Xk(k))) )+&
                     & temp*(dVKLplus*DerivMuKLCroit(Xk(k)) + dVKLmoins*DerivMuKLDecroit(Xk(k)))
                call ajout(k,k, CoefAjout, A)
                !! derivee de FF(k) par rapport a Xk(i)
                CoefAjout = DerAKL(2)*( dVKLplus*( MuKLCroit(Xk(k)) + MuKLDecroit(Xk(i)) ) + &
                     & dVKLmoins*( MuKLCroit(Xk(i)) + MuKLDecroit(Xk(k))) )+&
                     & temp*( dVKLplus*DerivMuKLDecroit(Xk(i)) + dVKLmoins*DerivMuKLCroit(Xk(i)) )
                call ajout(k,i, CoefAjout, A)
             ELSE 
                Select case ( Ntyps(i) ) 
                case (dirichlet+11)   
                   Ubord = Gb(i) ; Fbord = Gb(i)
                   dVKL = Aloc(3,1)* (Fbord-F(k))  
                   dVKLplus= max(dVKL, 0.D0) ; dVKLmoins=min(dVKL, 0.D0)
                   FF(k) =FF(k)+ dVKLplus*( MuKLCroit(Xk(k)) + MuKLDecroit(Ubord) ) +&
                        & dVKLmoins*( MuKLCroit(Ubord) + MuKLDecroit(Xk(k)) )
                   !! derivee de FF(k) par rapport a Xk(k)
                   CoefAjout = dVKLplus*DerivMuKLCroit(Xk(k)) + dVKLmoins*DerivMuKLDecroit(Xk(k))
                   call ajout(k,k, CoefAjout, A)
                case(1)
                   !!  rien
                End Select
             END IF
             !print*,'max FF apres ki',maxval(FF)
             !
             IF (Ntyps(j) /= 1 ) THEN 
                dVKL = Aloc(3,2)* (F(j)-F(k))
                dVKLplus= max(dVKL, 0.D0) ; dVKLmoins=min(dVKL, 0.D0)
  
                temp = AKL(Aloc(3,2), Xk(k), Xk(j))
                DerAKL = DerivAKL(Aloc(3,2), Xk(k), Xk(j))
  
                FF(k) = FF(k) + temp*( dVKLplus*( MuKLCroit(Xk(k)) + MuKLDecroit(Xk(j)) ) + &
                     & dVKLmoins*( MuKLCroit(Xk(j)) + MuKLDecroit(Xk(k))) )
                !! derivee de FF(k) par rapport a Xk(k)
                CoefAjout = DerAKL(1)*( dVKLplus*( MuKLCroit(Xk(k)) + MuKLDecroit(Xk(j)) ) + &
                     & dVKLmoins*( MuKLCroit(Xk(j)) + MuKLDecroit(Xk(k))) )+&
                     & temp*(dVKLplus*DerivMuKLCroit(Xk(k)) + dVKLmoins*DerivMuKLDecroit(Xk(k)))
                call ajout(k,k, CoefAjout, A)
                !! derivee de FF(k) par rapport a Xk(j)
                CoefAjout = DerAKL(2)*( dVKLplus*( MuKLCroit(Xk(k)) + MuKLDecroit(Xk(j)) ) + &
                     & dVKLmoins*( MuKLCroit(Xk(j))+MuKLDecroit(Xk(k))) )+&
                     & temp*( dVKLplus*DerivMuKLDecroit(Xk(j)) + dVKLmoins*DerivMuKLCroit(Xk(j)) )
                call ajout(k,j, CoefAjout, A)
             Else
                Select case ( Ntyps(j) ) 
                case (dirichlet+11)   
                   Ubord = Gb(j);  Fbord =  Gb(j)
                   dVKL = Aloc(3,2)* (Fbord-F(k))  
                   dVKLplus= max(dVKL, 0.D0) ; dVKLmoins=min(dVKL, 0.D0) 
                   FF(k) =FF(k)+ dVKLplus*( MuKLCroit(Xk(k) ) + MuKLDecroit(Ubord ) ) +&
                        & dVKLmoins*( MuKLCroit(Ubord ) + MuKLDecroit(Xk(k)) )
                   !! derivee de FF(k) par rapport a Xk(k)
                   CoefAjout = dVKLplus*DerivMuKLCroit(Xk(k)) + dVKLmoins*DerivMuKLDecroit(Xk(k))
                   call ajout(k,k, CoefAjout, A)
                case(1)
                   !!!!!!print*,'neumann',neumann
                   !!!!!!print*,'max FF neumann',maxval(FF)
                   !!  rien
                End Select
             END IF
            ! print*,'max FF apres kj',maxval(FF)
             !
          END IF  ! fin k
         ! print*,'max FF apres k',maxval(FF)
       END DO   !! fin boucle sur jt
       !!
       !!============================
       !!=============================
       !!print*,'debut sys lin'
       !Yk = -FF/A
       !print*,'max FF avant bigradient', maxval(FF), 'min FF', minval(FF)
       X0= 0.D0 
       Yk = bigradient(A,-FF, X0,Tolerencegradient)
       !write(*,*)'Yk=',Yk
       !write(*,*)'F(1)=',F(1)
       !print*,'kiter = ', kiter,'erreur NEWTON sqrt(sum(Yk*Yk))',sqrt(sum(Yk*Yk))
  
       If (sqrt(dot_product(Yk,Yk)) <TolerenceNewton) exit
       Xk=Xk+Yk
       !print*,'max Xk', maxval(Xk), 'min Xk', minval(Xk)
       !print*,'max FF', maxval(FF), 'min FF', minval(FF)
  
       !! impression de A etde F
       !     DO i = 1, size(A%IndPL)-1
       !           print*,'iseg = ', i, 'sum line = ', SUM(A%Tmat(A%IndPL(i): A%IndPL(i+1)-1))
       !    END DO
  
    end Do
  
    U=Xk+Yk
  
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
  
  
  END SUBROUTINE NewtonEndP1sommets