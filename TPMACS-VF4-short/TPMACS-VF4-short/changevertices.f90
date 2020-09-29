subroutine changevertices

    USE longr
    USE imprime
    USE parmmage
    IMPLICIT NONE 
    CHARACTER(len=6)      :: oldprf
    INTEGER               :: n, k, ks,Kvois,Lvois,is1, is2, ismin, ismax
    INTEGER               :: isegmin, isegmax, jt, iseg, next,typeseg
    INTEGER               :: m, ii, ik, jL, is, js , ls, i , j , kt, ncentrevernoi
    INTEGER               :: Nbstriang, nusommet1, nusommet2, nextint,nextext
  
  
    REAL(kind=long) , DIMENSION(:,:), ALLOCATABLE ::  CoordSbis
    Integer , DIMENSION(:,:), ALLOCATABLE         ::  NuSoKbis, NumTVoisSegbis
    INTEGER, DIMENSION( : ), ALLOCATABLE          :: NombresommetstrouveK
    INTEGER, DIMENSION( :,: ), ALLOCATABLE        ::  NuSegbis
    INTEGER, DIMENSION(:), ALLOCATABLE            :: PermutSegment,PermutSommet,ntypsbis,NTypSegbis
  
  
    real(kind=long), dimension(2)        :: xs1, xs2, xL,xK, S1,S2,S3,P,v0,v1,v2
    real(kind=long), dimension(3)        :: x, y 
    real(kind=long), dimension(17)       :: ligne
  
  
    real(kind=long)        :: x1,y1, ps, Denom, ubc, vbc
    real(kind=long)        :: as1,as2,as3,angle1,angle2,angle3,pi
    Logical                :: Zvernoi, trouve
    CHARACTER(len=14)      :: chaine
  
    !-------------------
    ! Debut du programme
    !-------------------
    oldprf = prefix
    prefix = 'SEGMT '

    !--------------------------------------------------------------
  ! calcul de NsInt le nombre de sommets ou la solution est calcul�e
  !-----------------------------------------------------------------
  NsInt = 0
  DO i= 1, Nbs
     IF( Ntyps(i) /= 1) NsInt = NsInt +1 
  END DO
  Nbord = Nbs - NsInt
  write(*,*)'le nombre de sommets interieurs est:',NsInt
  write(*,*)'le nombre de sommets exterieurs est:',Nbord
  write(*,*)'le nombre de sommets total est:',Nbs

  CALL prvari(uprint,'nombre des sommets interieurs  (NsInt) ', NsInt )
  CALL prvari(uprint,'nombre des sommets sur le bord (Nbord) ', Nbord )
  CALL prvari(uprint,'nombre des sommets total (Nbs) ', Nbs )
  !
  !! Calcule de Nsegint et Nsegbord
  !
  Nsegint =0
  Do iseg =1,Nseg
     if(NtypSeg(iseg) == 0) Then
        Nsegint =Nsegint+1
     end if
  end do

  Nsegbord= Nseg-Nsegint
  write(*,*)'le nombre de segments interieurs est:',Nsegint
  write(*,*)'le nombre de segments exterieurs est:',Nsegbord
  write(*,*)'le nombre de segments total est:',Nseg

  CALL prvari(uprint,'le nombre de segments interieurs (Nsegint) ', Nsegint)
  CALL prvari(uprint,'le nombre de segments exterieur (Nsegbord) ', Nsegbord)
  CALL prvari(uprint,'le nombre de segments total             (Nseg) ', Nseg)
  !!---------------------------------------------------
  !! renumerotation des sommets int�rieurs de 1 � NsInt
  !! --------------------------------------------------
  Allocate(PermutSommet(Nbs))
  Nextint = 1 ; Nextext = 1
  DO is = 1, Nbs
     IF( Ntyps(is) /= 1) then 
        PermutSommet(is)=Nextint
        Nextint=Nextint+1
     else 
        PermutSommet(is)=Nsint+Nextext
        Nextext=Nextext+1
     end IF
  Enddo
  !!write(*,*)'PermutSommet=',PermutSommet

  !!----------------------------------------------------------------------
  !! permuter les coordonn�es des sommets et type des sommets COORDS NTYPS
  !!----------------------------------------------------------------------

  Allocate(CoordSbis(2, Nbs),ntypSbis(Nbs))
  CoordSbis = CoordS
  NtypSbis  = NtypS
  DO is = 1, Nbs
     CoordS(1,PermutSommet(is))=  CoordSbis(1,is)
     CoordS(2,PermutSommet(is))=  CoordSbis(2,is)
     NtypS(PermutSommet(is)) = NtypSbis(is)
  EndDO
  !write(*,*)'CoordSbis avant permutation=',CoordSbis
  !write(*,*)'CoordS apres permutation=',CoordS

  !write(*,*)'ntypsbis avant permutation =', NtypSbis
  !write(*,*)'ntyps apres permutation =', NtypS
  DEALLOCATE(CoordSbis,Ntypsbis)

  !! ------------------------------------------------------------------------------
  !! Mise � jour des numeros des sommets des triangles: NuSoK
  !! -----------------------------------------------------------------------------
  Allocate(NuSoKbis(1:3,Nbt))
  NuSoKbis = NuSoK
  DO kt = 1, Nbt
     !! is = NuSoKbis(1,kt); js = NusoKbis(2,kt); ks = NusoKbis(3,kt)
     Do m = 1, typegeomaille
        NuSoK (m,kt) = permutsommet(NuSoKbis(m,kt))
     END Do
  END DO
  !write(*,*)'NuSoKbis avant permutation=', NuSoKbis
  !write(*,*)'NuSoK apres permutation=', NuSoK
  Deallocate(NuSoKbis)

  !!--------------------------------------------------
  !! Mise � jour des numeros des sommets des segments : NuSeg
  !! --------------------------------------------------
  Allocate(NuSegbis(2,Nseg))
  NuSegbis = NuSeg
  Do iseg = 1,Nseg
     !! is = NuSegbis(1,iseg); js = NuSegbis(2,iseg)
     NuSeg (1,iseg) = permutsommet(NuSegbis(1,iseg))
     NuSeg (2,iseg) = permutsommet(NuSegbis(2,iseg))
  End do

  !write(*,*)'NuSegbis avant permutation en Sommets=',NuSegbis
  Deallocate(NuSegbis)
  !write(*,*)'NuSeg apres permutation en Sommets =',NuSeg


  !! Les autres tableaux : 
  !! NuVoisK, Ntypseg, NumTVoisSeg, NombVoisSeg, CoordK : ces tableaux sont inchang�s

  !!============================================================================
  !! NUMEROTATION DES SEGMENTS A L'INTERIEUR
  !!============================================================================
  !! A ce stade nous avons renumerot� les inconnues � l'int�rieur 
  !! en gardant la nouvelle num�rotation
  !! 
  !! Maintenant on va renumeroter les segments 
  !! a l'interieur et mettre � jour les tableau
  !! 

  !!--------------------------------------------------
  !! renumerotation des segmets interieurs de 1 a  Nsegint
  !! --------------------------------------------------
  Allocate (PermutSegment(1:Nseg))
  Nextint = 1 ; Nextext = 1
  DO iseg= 1, Nseg
     IF( NtypSeg(iseg) == 0) then
        PermutSegment(iseg)=Nextint
        Nextint=Nextint+1
     else
        PermutSegment(iseg)=Nextext+Nsegint
        Nextext=Nextext+1
     end IF
  Enddo
  !write(*,*)'PermutSegment=',PermutSegment

  !!---------------------------------------------------------------------------------
  !! Permuter les types des segment et numeros des sommets 
  !!---------------------------------------------------------------------------------

  Allocate(NtypSegbis(Nseg), NuSegbis(2,Nseg))
  NtypSegbis  = NtypSeg 
  NuSegbis = NuSeg
  DO iseg = 1, Nseg
     Nuseg(1, PermutSegment(iseg)) =  Nusegbis(1,iseg)
     Nuseg(2, PermutSegment(iseg)) =  Nusegbis(2,iseg)
     NtypSeg(PermutSegment(iseg)) = NtypSegbis(iseg)
  END DO
  !write(*,*)'NuSegbis avant permutation en Segments=', Nusegbis
  !write(*,*)'NuSeg apres permutation en segments=', NuSeg

  !write(*,*)'NtypSegbis avant  permutation =', NtypSegbis
  !write(*,*)'NtypSeg apres permutation =', NtypSeg

  DEALLOCATE(NtypSegbis)

  !------------------------------------
  ! Remettre a jour nombVoisSeg
  !-----------------------------------
  DO iseg=1, Nseg
     If (NTypSeg(iseg) == 0)  Then 
        NombVoisSeg(iseg)= 2
     Else 
        NombVoisSeg(iseg)= 1
     endif
  ENDDO

  !!-------------------------------------
  !! Remettre a jour NumTVoisSeg(2,Nseg)
  !!-------------------------------------

  Allocate(NumTVoisSegbis(2,Nseg))
  NumTVoisSegbis  = NumTVoisSeg
  DO iseg = 1, Nseg
     NumTVoisSeg(1:2,PermutSegment(iseg)) = NumTVoisSegbis(1:2,iseg)
  EndDO
  !write(*,*)'NumTVoisSegbis avant permutation =', NumTVoisSegbis
  !write(*,*)'NumTVoisSeg apres permutation =', NumTVoisSeg
  DEALLOCATE(NumTVoisSegbis)



  !----------------------------------------------------------------------------------
  !! Calcul des coordonn�es aux milieux des segments
  !----------------------------------------------------------------------------------
  ALLOCATE(coordMSeg(1:2,1:Nseg))
  DO iseg = 1,Nseg
     is = NuSeg(1,iseg) ; js = NuSeg(2,iseg)
     coordMSeg(1,iseg)=(coordS(1,is)+coordS(1,js))/2.
     coordMSeg(2,iseg)=(coordS(2,is)+coordS(2,js))/2.
  End do


!!$  !-----------------------------------
!!$  ! Construction de NuMSeg(1:3,1:nbt)
!!$  !-----------------------------------
!!$  ALLOCATE(NuMSeg(1:3,1:Nbt))
!!$  Do k = 1, Nbt
!!$     i = NuSoK(1,k);j = NuSoK(2,k);l = NuSoK(3,k);
!!$     Do ii = 1 , Nseg
!!$        If(NuSeg(1,ii) == i .AND. NuSeg(2,ii) == j)Then
!!$           NuMSeg(1,k)= ii
!!$        End if
!!$        If(NuSeg(1,ii) == j .AND. NuSeg(2,ii) == i)Then
!!$           NuMSeg(1,k)= ii
!!$        End if
!!$        If(NuSeg(1,ii) == j .AND. NuSeg(2,ii) == l)Then
!!$           NuMSeg(2,k)= ii
!!$        End if
!!$        If(NuSeg(1,ii) == l .AND. NuSeg(2,ii) == j)Then
!!$           NuMSeg(2,k)= ii
!!$        End if
!!$        If(Nuseg(1,ii) == i .AND. NuSeg(2,ii) == l)Then
!!$           NuMSeg(3,k)= ii
!!$        End if
!!$        If(Nuseg(1,ii) == l .AND. NuSeg(2,ii) == i)Then
!!$           NuMSeg(3,k)= ii
!!$        End if
!!$     End Do
!!$  End do
!!$
!!$  write(*,*)'NuMSeg des segments=',NuMSeg

  !-----------------------------------------
  ! Construction de NuMSeg(1:3,1:nbt) pour que  (maz)
  !-----------------------------------------
  ALLOCATE(NuMSeg(1:3,1:Nbt))
  Do kt = 1, Nbt
     is = NuSoK(1,kt); js = NuSoK(2,kt); ls = NuSoK(3,kt)
     !  NuMSeg(1, jt) : le numero du segment en face du sommet is
     !  NuMSeg(2, jt) : le numero du segment en face du sommet js
     !  NuMSeg(3, jt) : le numero du segment en face du sommet ls
     Do ii = 1 , Nseg
        If (NuSeg(1,ii) == is .AND. NuSeg(2,ii) == js) Then
           NuMSeg(3,kt)= ii
        End if
        If (NuSeg(1,ii) == js .AND. NuSeg(2,ii) == is) Then
           NuMSeg(3,kt)= ii
        End if
        If (NuSeg(1,ii) == js .AND. NuSeg(2,ii) == ls) Then
           NuMSeg(1,kt)= ii
        End if
        If (NuSeg(1,ii) == ls .AND. NuSeg(2,ii) == js) Then
           NuMSeg(1,kt)= ii
        End if
        If (NuSeg(1,ii) == is .AND. NuSeg(2,ii) == ls) Then
           NuMSeg(2,kt)= ii
        End if
        If (Nuseg(1,ii) == ls .AND. NuSeg(2,ii) == is)Then
           NuMSeg(2,kt)= ii
        End if
     End Do
  End do

  !write(*,*)'NuMSeg des segments=',NuMSeg



  If ( Typegeomaille == 3 )  then  !!! maillage triangulaire
     !============================
     !=============================
     !! ------------------------------------------------------------------------
     !!  calcul AireK des sommets pour le maillage triangulaire
     !!--------------------------------------------------------------------------
     Do jt = 1, Nbt
        i = NuSoK(1,jt) ; j = NuSoK(2,jt) ; k = NuSoK(3,jt) !! i,j,k numero globale

        x(1) = coordS(1,i) ; y(1) = coordS(2,i)
        x(2) = coordS(1,j) ; y(2) = coordS(2,j)
        x(3) = coordS(1,k) ; y(3) = coordS(2,k)
        !!
        AireK(jt) = ABS( (x(2)-x(1))* (y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)) ) /2.
     END Do
     CALL prvarr(uprint, 'Max aireK = ', maxval(AireK))
     CALL prvarr(uprint, 'Min aireK = ', minval(AireK))
     !write(*,*)'AireK=',AireK

     !------------------------------
     ! Calcul des aires des diamands
     !------------------------------
     ALLOCATE(AireD(Nseg))
     Do iseg = 1, Nseg
        Kvois = NumTVoisSeg(1,iseg)
        js = Nuseg(1,iseg) ; ks = Nuseg(2,iseg)

        x(1) = coordK(1,Kvois) ; y(1) = coordK(2,Kvois)
        x(2) = coordS(1,js) ; y(2) = coordS(2,js)
        x(3) = coordS(1,ks) ; y(3) = coordS(2,ks)
        !!
        AireD(iseg) = ABS( (x(2)-x(1))* (y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)) ) /2.D0
        If (NombvoisSeg(iseg)==2) then

           Lvois = NumTVoisSeg(2,iseg)
           js = NuSeg(1,iseg) ; ks = NuSeg(2,iseg)

           x(1) = coordK(1,Lvois) ; y(1) = coordK(2,Lvois)
           x(2) = coordS(1,js) ; y(2) = coordS(2,js)
           x(3) = coordS(1,ks) ; y(3) = coordS(2,ks)
           !!
           AireD(iseg) = AireD(iseg)+ &
                & ABS( (x(2)-x(1))* (y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)) )/2.D0
        end If
     End Do
     CALL prvarr(uprint, 'Max aireD = ', maxval(AireD))
     CALL prvarr(uprint, 'Min aireD = ', minval(AireD))
     !write(*,*)'AireK=',AireK
     !!write(*,*)'AireDiamant=',AireD

  end If

  ! -------------------------------------------------------------------
  !!  calcul AireDSommet associe a chaque sommet maillage de Donald
  !!----------------------------------------------------------------
  Allocate (AireDSommet(NbS))
  AireDSommet=0.D0
  Do jt = 1, Nbt
     i = NuSoK(1,jt) ; j = NuSoK(2,jt) ; k = NuSoK(3,jt) 

     AireDSommet(i) = AireDSommet(i) + AireK(jt)/3.D0
     AireDSommet(j) = AireDSommet(j) + AireK(jt)/3.D0
     AireDSommet(k) = AireDSommet(k) + AireK(jt)/3.D0
  END Do
  CALL prvarr(uprint, 'Max AireDSommet = ', maxval(AireDSommet))
  CALL prvarr(uprint, 'Min AireDSommet = ', minval(AireDSommet))
  !write(*,*)'AireDSommet=',AireDSommet


  !=============
  ! un test pour voir si ca bien marche
  !==============
  Do iseg = 1, Nseg
     If(NTypSeg(iseg)==0 .and. NombVoisSeg(iseg)/= 2) then
        print*,'PROBLEME segment interieur iseg = ', iseg
        print*,'Nuseg', Nuseg(:,iseg)
        print*,'Ntypseg', Ntypseg(iseg)
        print*,'NumTVoisSeg',NumTVoisSeg(:, iseg)
        stop
     endif
     If(NTypSeg(iseg)/=0 .and. NombVoisSeg(iseg)/= 1) then
        print*,'PROBLEME segment bord = ',iseg
        stop
     endif
  END Do

  !! -----------------------
  !!  calcul dKL la distance entre deux mailles voisines,
  !! c a d en prenant un segment, les deux triangles voisins
  !! SKL designe la longueur du segment
  !!------------------------
  Allocate (SKL(Nseg))
  !  Segment int�reurs
  Do iseg = 1, Nseg
     ii = (NtypSeg(iseg)) 

     Select case (ii) 

     case (0)   !! segment � l'nterieur  
        is = NuSeg(1,iseg); js=NuSeg(2,iseg)
        ik= NumTVoisSeg(1,iseg); jL =NumTVoisSeg(2,iseg)
        dKL(iseg) = sqrt((CoordK(1,ik)-coordK(1,jL))**2 +(CoordK(2,ik)-coordK(2,jL))**2) 
        SKL(iseg) = sqrt((CoordS(1,is)-coordS(1,js))**2 +(CoordS(2,is)-coordS(2,js))**2)

        TauKL(iseg) = SKL(iseg)/dKL(iseg)     

     case (dirichlet)
        is = NuSeg(1,iseg); js=NuSeg(2,iseg)
        ik= NumTVoisSeg(1,iseg) 
        SKL(iseg) = sqrt((CoordS(1,is)-coordS(1,js))**2 +(CoordS(2,is)-coordS(2,js))**2)
        ! milieu du segment 
        x1= (CoordS(1,is)+CoordS(1,js))/2.D0 ; y1= (CoordS(2,is)+CoordS(2,js))/2.D0 ; 
        dKL(iseg) = sqrt((CoordK(1,ik)-x1)**2 +(CoordK(2,ik)-y1)**2) 

        TauKL(iseg) = SKL(iseg)/dKL(iseg)       
     case(neumann)
        !stop ' Segments : condition de newmann reflechir'
        TauKL(iseg) = -100000
     case default
        print*,'type segment =', ii
        !stop'segments prob taukl' 
     End Select
  end Do



end subroutine changevertices