SUBROUTINE   segments
  !------------------------------------------------
  ! Nbs  : nombre des sommets total
  ! Nbt  : nombre des triangle
  ! Ns   : nombre des sommets interieur au domaine 
  !        nombre des sommets ou la solutions est calculee
  ! Nb   : nombre des sommets au bord
  ! Nbs = Ns + Nb 
  ! NuMSeg : numero arete triangle, numero milieu segment
  ! CoordMseg : coordonnee du popint milieu du segment
  !----------------------------------
  ! Declaration des variables locales
  !----------------------------------
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
  INTEGER, DIMENSION( : ), ALLOCATABLE          :: NombVoisSeg, NombresommetstrouveK
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


  If (maillageregulier == 1) then 

     Typegeomaille = 4                 ! type geometrie de la maille (rectangle)
     !
     ! Construction de la liste des segments 
     ! les numeros des triangles voisins
     ! calcul des aires des triangles
     !
     Nseg = NX*(Ny+1) + Ny*(Nx+1)
     ALLOCATE ( NuSeg(1:2,1:Nseg), NTypSeg(Nseg) )
     NTypSeg = -1000
     ALLOCATE(NumTVoisSeg(1:2, NSeg))  ! Numero des triangles voisins
     ALLOCATE(NombVoisSeg(NSeg))       ! Nombre de voisins au plus 2; 
     NumTVoisSeg = -10000
     !  sinon 1 (cest a dire segment sur le bord
     !
     ! Pour un triangle kt = i+(j-1)*Nx les sommets sont
     ! kt+j-1  : coin bas gauche
     ! kt + j  : coin bas droit
     ! kt+j+Nx : coin haut gauche
     ! kt+j+Nx+1 : coin haut droit
     !

     iseg = 0
     ! numéroter les segments à gauche et à droite 
     DO j = 1, NY
        i = 1 ;  kt = i+(j-1)*Nx ; iseg = iseg +1
        NuSeg(:,iseg ) = (/ kt+j-1, kt+j+Nx/)
        NtypSeg(iseg) = CLgauche ; NombVoisSeg(iseg) = 1 ; NumTVoisSeg(1, iseg) = kt
        DO  i = 2, NX
           iseg = iseg +1
           kt = i+(j-1)*Nx
           NuSeg(:,iseg ) = (/ kt+j-1, kt+j+Nx/)
           NtypSeg(iseg)=0 ;       NombVoisSeg(iseg) = 2 
           NumTVoisSeg(1, iseg) = kt; NumTVoisSeg(2, iseg) = kt-1 
        END DO
        i = NX ;  kt = i+(j-1)*Nx ; iseg = iseg +1
        NuSeg(:,iseg ) = (/ kt+j, kt+j+Nx+1/)
        NtypSeg(iseg) = CLdroite; NombVoisSeg(iseg) = 1
        NumTVoisSeg(1, iseg) = kt
     END DO
     !!     print*,'iseg = ', iseg
     ! numéroter les segments en bas et en haut
     DO i = 1, NX
        j = 1 ;  kt = i+(j-1)*Nx ; iseg = iseg +1
        NuSeg(:,iseg ) = (/ kt+j-1, kt+j/)
        NtypSeg(iseg) = CLbas
        NombVoisSeg(iseg) = 1 ; NumTVoisSeg(1, iseg) = kt

        DO  j = 2, NY
           iseg = iseg +1
           !!     print*,'iseg = ', iseg
           kt = i+(j-1)*Nx
           NuSeg(:,iseg ) = (/ kt+j-1, kt+j/)
           NtypSeg(iseg) = 0 ; NombVoisSeg(iseg) = 2 
           NumTVoisSeg(1, iseg) = kt ; NumTVoisSeg(2, iseg) = kt-NX
        END DO

        j = NY ;  kt = i+(j-1)*Nx ; iseg = iseg +1
        NuSeg(:,iseg ) = (/ kt+j+Nx, kt+j+Nx+1/)
        NtypSeg(iseg) = CLhaut ; NombVoisSeg(iseg) = 1
        NumTVoisSeg(1, iseg) = kt
     END DO


     !! ----------------------------------------------
     !!  calcul AireK pour le maillage rectangulaire
     !!------------------------------------------------
     Allocate (AireK(Nbt))
     dx = Lx/Nx ; dy = Ly/Ny
     AireK = dx*dy
     CALL prvarr(uprint, 'Max aaireK = ', maxval(AireK))
     CALL prvarr(uprint, 'Min aireK = ', minval(AireK))

  ELSE 
     Typegeomaille = 3 ! type geometrie de la maille (triangle)
     !! maillage triangulaire

     SELECT CASE (Maillage_Format)
     CASE (1) 
        !! maillage issus du logiciel triangle

        ! lecture de 
        ! liste des segments NuSeg et NtypSeg
        ! construction de la liste des triangles voisins NumTvoisSeg, NomVoisSeg
        ! calcul des aires des triangles
        ! ------------------------------------------------------------- 
        !-- lecture de la liste des segments et des types nom.1.edge
        ! -----------------------------------------------------------
        OPEN (uedge,     file = nom_edge   ,  status = 'old')

        READ (uedge,*) Nseg, typeseg

        CALL prvari(uprint,'nombre des segments   (Nseg) ', Nseg )

        ALLOCATE ( NuSeg(1:2,1:Nseg), NTypSeg(Nseg) )

        DO iseg = 1,Nseg
           READ(uedge,*) n, NuSeg(1,iseg), NuSeg(2,iseg), NTypSeg(iseg)
        END DO

        close(uedge)

        ! ==============
        ! ==============

        PRINT*,'ATTENTION LE TYPE DU BORD EST IMPOSE Neumann'
        WRITE(uprint,*)'ATTENTION LE TYPE DU BORD EST IMPOSE Neumann'
        DO iseg = 1,Nseg
           IF  (NTypSeg(iseg)/=0) NTypSeg(iseg) = Neumann
        END DO
        !===================================================================
        ! Construction  de la liste des deux triangles 
        !===================================================================

        ALLOCATE( NumTVoisSeg(1:2, NSeg) )  ! Numero des triangles voisins
        ALLOCATE( NombVoisSeg(NSeg) )       ! Nombre de voisins au plus 2; 
        !sinon 1 (cest a dire segment sur le bord


        NombVoisSeg = 0 ; NumTVoisSeg = -1
        DO jt = 1, Nbt
           DO k = 1,typegeomaille
              next = MOD(k,typegeomaille)+1
              is1 = NuSoK(k,   jt)
              is2 = NuSoK(next,jt)
              ismin = MIN(is1,is2)
              ismax = MAX(is1,is2)
              !!
              DO iseg = 1, Nseg
                 isegmin = MIN(NuSeg(1,iseg), NuSeg(2,iseg))
                 isegmax = MAX(NuSeg(1,iseg), NuSeg(2,iseg))

                 IF(ismin == isegmin.AND. ismax == isegmax) Then
                    NombVoisSeg(iseg)=NombVoisSeg(iseg)+1  !! Nombre de Voisins
                    NumTVoisSeg(NombVoisSeg(iseg), iseg)=jt            !! 
                    EXIT
                 END IF
              ENDDO
           END DO
        END DO

     CASE (2)
        !! format de type boyer
        ! ------------------------------------------------------------- 
        !-- lecture de la liste des segments et des types maillage.arete
        ! -----------------------------------------------------------
        OPEN (uedge,     file = nom_edge   ,  status = 'old')

        DO i = 1,17
           READ (uedge,*) 
        ENDDO

        READ (uedge,*) chaine, Nseg  !!, typeseg

        CALL prvari(uprint,'nombre des segments   (Nseg) ', Nseg )

        ALLOCATE (NuSeg(1:2,1:Nseg), NTypSeg(Nseg)) ! Numero segment et typesegment
        ALLOCATE(NumTVoisSeg(1:2, NSeg))  ! Numero des triangles voisins au segment
        ALLOCATE(NombVoisSeg(NSeg))       ! Nombre de voisins au plus 2;
        ! sinon 1 (cest a dire segment sur le bord

        DO iseg = 1,Nseg
           READ(uedge,*) (ligne(j),j=1,17)

           NuSeg(1,iseg) = int(ligne(1))
           NuSeg(2,iseg) = int(ligne(2))
           NumTVoisSeg(1,iseg) = int(ligne(3))
           NumTVoisSeg(2,iseg) = int(ligne(4))
           NTypSeg(iseg) = int(ligne(17))
        END DO

        close(uedge)

        !! affecter le nombre des voisins
        DO iseg=1, Nseg
           If (NTypSeg(iseg) == 0)  Then 
              NombVoisSeg(iseg)= 2
           Else 
              NombVoisSeg(iseg)= 1
           endif
        ENDDO
        !
        PRINT*,'ATTENTION LE TYPE DU BORD EST IMPOSE Neumann'
        WRITE(uprint,*)'ATTENTION LE TYPE DU BORD EST IMPOSE Neumann'
        DO iseg = 1,Nseg
           IF  (NTypSeg(iseg)/=0) NTypSeg(iseg) = Neumann
        END DO
        !! ============================================================
        !! calcul de Ntyps dans le cas de maillage boyer 
        !! =============================================================
        Allocate (Ntyps(1:Nbs) )  
        Ntyps = 0 
        DO iseg = 1,Nseg
           is = NuSeg(1,iseg) ; js = NuSeg(2,iseg)
           IF  (NTypSeg(iseg)/=0) then 
              Ntyps(is) = Neumann ; Ntyps(js) = Neumann
           Endif
        Enddo
        !!
        !! Construction de NuSoK pour calculer l'aire de K
        !!
        !! On parcout les segments et connaissant les numeros des triangles voisins
        !! on déduit que les deux sommets du segment sont deux sommets du triangle
        !! on arête dès que on trouve les trois sommets. 

        Allocate (NuSoK(1:3,Nbt))

        Allocate ( NombresommetstrouveK(1:Nbt) )
        NombresommetstrouveK = 0

        DO iseg = 1, Nseg
           kt = NumTVoisSeg(1,iseg)
           If (kt >0) THEN
              nusommet1= NuSeg(1,iseg) ;  nusommet2= NuSeg(2,iseg) 

              select case (NombresommetstrouveK(kt))
              case (0) 
                 NuSoK(1,kt) = nusommet1 ; NusoK(2,kt)= nusommet2
                 Nbstriang = 2
              case (2)
                 ! recherche si nusommet est dans la liste des sommets de K
                 If (nusommet1 == NuSoK(1,kt) .or. nusommet1 == NuSoK(2,kt) ) then 
                    trouve = .true. 
                 else
                    NuSoK(3,kt) = nusommet1
                    Nbstriang = 3 
                 end If
                 If (nusommet2 == NuSoK(1,kt) .or. nusommet2 == NuSoK(2,kt) ) then 
                    trouve = .true. 
                 else
                    NuSoK(3,kt) = nusommet2
                    Nbstriang = 3 
                 end If

              case (3)
                 !! rien a faire
                 Nbstriang = 3 
              case default 
                 stop'pb nu'              
              end select
              NombresommetstrouveK(kt)= Nbstriang
           ENDIF
           !!
           !! On regarde le deuxième voisin 
           !!
           kt = NumTVoisSeg(2,iseg)
           If (kt >0) THEN
              nusommet1= NuSeg(1,iseg) ;  nusommet2= NuSeg(2,iseg) 

              select case (NombresommetstrouveK(kt))
              case (0) 
                 NuSoK(1,kt) = nusommet1 ; NusoK(2,kt)= nusommet2
                 Nbstriang = 2
              case (2)
                 ! recherche si nusommet est dans la liste des sommets de K
                 If (nusommet1 == NuSoK(1,kt) .or. nusommet1 == NuSoK(2,kt) ) then 
                    trouve = .true. 
                 else
                    NuSoK(3,kt) = nusommet1
                    Nbstriang = 3 
                 end If
                 If (nusommet2 == NuSoK(1,kt) .or. nusommet2 == NuSoK(2,kt) ) then 
                    trouve = .true. 
                 else
                    NuSoK(3,kt) = nusommet2
                    Nbstriang = 3 
                 end If

              case (3)
                 !! rien a faire
                 Nbstriang = 3 
              case default 
                 stop'pb segment NuTvois'              
              end select
              NombresommetstrouveK(kt)= Nbstriang

           ENDIF
        END DO
        DEALLOCATE ( NombresommetstrouveK )
     case default
        print*,'segment : format maillage inconnu'

     end select
  End If

  !! resume : apres maillage 
  !! on dispose pour tout les type de maillage de : 
  !!               regulier  triangle boyer   amdba 
  !! NuSoK :            1       1          
  !! NuVoisK            1       1
  !! Nbs                1       1      1     1
  !! Nbt                1       1      1     1
  !! NombVoisSeg        
  !! CoordS             1       1      1     1
  !! ntyps              1       1      1
  !! CoordK             1       1      1     1         

  !! resume : segment 
  !! on dispose pour tout les type de maillage de : 
  !!               regulier  triangle boyer  amdba
  !! NuSoK :                          1      
  !! NuSeg :             1    1       1  
  !! NTypSeg :           1    1       1
  !! NumTVoisSeg         1    1       1
  !! NombVoisSeg         0    1       1
  !! aireK               1




  !--------------------------------------------------------------
  ! calcul de NsInt le nombre de sommets ou la solution est calculée
  !-----------------------------------------------------------------
  NsInt = 0
  DO i= 1, Nbs
     IF( Ntyps(i) == 0) NsInt = NsInt +1 
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
  !! renumerotation des sommets intérieurs de 1 à NsInt
  !! --------------------------------------------------
  Allocate(PermutSommet(Nbs))
  Nextint = 1 ; Nextext = 1
  DO is = 1, Nbs
     IF( Ntyps(is) == 0) then 
        PermutSommet(is)=Nextint
        Nextint=Nextint+1
     else 
        PermutSommet(is)=Nsint+Nextext
        Nextext=Nextext+1
     end IF
  Enddo
  !!write(*,*)'PermutSommet=',PermutSommet

  !!----------------------------------------------------------------------
  !! permuter les coordonnées des sommets et type des sommets COORDS NTYPS
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
  !! Mise à jour des numeros des sommets des triangles: NuSoK
  !! -----------------------------------------------------------------------------
  Allocate(NuSoKbis(typegeomaille,Nbt))
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
  !! Mise à jour des numeros des sommets des segments : NuSeg
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
  !! NuVoisK, Ntypseg, NumTVoisSeg, NombVoisSeg, CoordK : ces tableaux sont inchangés

  !!============================================================================
  !! NUMEROTATION DES SEGMENTS A L'INTERIEUR
  !!============================================================================
  !! A ce stade nous avons renumeroté les inconnues à l'intérieur 
  !! en gardant la nouvelle numérotation
  !! 
  !! Maintenant on va renumeroter les segments 
  !! a l'interieur et mettre à jour les tableau
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
  !! Calcul des coordonnées aux milieux des segments
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
     Allocate (AireK(Nbt))
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



  !----------------
  ! Impressions eventuelles
  !------------------------
  IF (iprint >=2) THEN 
     write(uprint,*)'  kt , (NuSoK(j,kt),j=1,2, 3, ...) '
     DO kt=1,nbt
        WRITE(uprint,113) kt , (NuSoK(j,kt),j=1,typegeomaille)
     ENDDO
  END IF
  !------------------------

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
  Allocate (dKL(Nseg), SKL(Nseg), tauKL(Nseg))
  !  Segment intéreurs
  Do iseg = 1, Nseg
     ii = (NtypSeg(iseg)) 

     Select case (ii) 

     case (0)   !! segment à l'nterieur  
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
        stop'segments prob taukl' 
     End Select
  end Do

  !===========================================
  ! Verification si le maillage est de Voroni
  !==========================================
  ! 1.  verifer si le segment reliant les centres des triangles est orthogonale
  !     au segment
  ! 
  Select case ( ChoixSchema)
  case(VF4,P1milieux)

     Zvernoi = .false.
     Do iseg = 1, Nseg
        If(NombVoisSeg(iseg) == 2 ) then 
           is1 = NuSeg(1, iseg) ; is2 = NuSeg(2,iseg)   ! numeros sommets
           XS1 = coordS(:,is1)  ; XS2 = coordS(:,is2)
           xL = CoordK(:,NumTVoisSeg(1,iseg)) ;  xK = CoordK(:,NumTVoisSeg(2,iseg))  
           ps = (xs2(1)-xs1(1))*(xL(2)-xk(2))-(xs2(2)-xs1(2))*(xL(1)-xk(1))
           If(abs(ps)< 1.d-9) then
              Zvernoi = .true.
              write(uprint,*)' pb vornoi : segment et distance entre les deux centres PAS ORTHOGONALE ', ' iseg = ',iseg
              write(uprint,*)' les coordonnees des sommets du segment ', ' iseg =', iseg
              write(uprint,*) ' sommet is1 = ',is1, ' xs1 = ', xs1
              write(uprint,*) ' sommet is2 = ',is2, ' xs2 = ', xs2
              write(uprint,*) ' triangle voisin 1 = ',NumTVoisSeg(1,iseg), ' xL = ', xL
              write(uprint,*) ' triangle voisin 2 = ',NumTVoisSeg(2,iseg), ' xk = ', xK
           endif
        end If
     end Do

  case(P1Sommets)
     Zvernoi = .false.
     Do iseg = 1, Nseg
        If(NombVoisSeg(iseg) == 2 ) then
           is1 = NuSeg(1, iseg) ; is2 = NuSeg(2,iseg)   ! numeros sommets
           XS1 = coordS(:,is1)  ; XS2 = coordS(:,is2)
           xL = CoordK(:,NumTVoisSeg(1,iseg)) ;  xK = CoordK(:,NumTVoisSeg(2,iseg))
           ps = (xs2(1)-xs1(1))*(xL(2)-xk(2))-(xs2(2)-xs1(2))*(xL(1)-xk(1))
           If(abs(ps)< 1.d-9) then
              Zvernoi = .true.
              write(uprint,*)' pb vornoi : segment et distance entre les deux centres PAS ORTHOGONALE ', ' iseg = ',iseg
              write(uprint,*)' les coordonnees des sommets du segment ', ' iseg =', iseg
              write(uprint,*) ' sommet is1 = ',is1, ' xs1 = ', xs1
              write(uprint,*) ' sommet is2 = ',is2, ' xs2 = ', xs2
              write(uprint,*) ' triangle voisin 1 = ',NumTVoisSeg(1,iseg), ' xL = ', xL
              write(uprint,*) ' triangle voisin 2 = ',NumTVoisSeg(2,iseg), ' xk = ', xK
           endif
        end If
     end Do
  case default

  end select

  !
  ! 2. verifier si la distance entre les centres est positive
  ! 
  Do iseg = 1, Nseg
     ii = (NtypSeg(iseg))

     !!print*,'NtypSeg(',iseg,')=',NtypSeg(iseg)

     Select case (ii) 

     case (0)   !! segment à l'nterieur  
        is = NuSeg(1,iseg); js=NuSeg(2,iseg)
        ik= NumTVoisSeg(1,iseg); jL =NumTVoisSeg(2,iseg)
        dKL = sqrt((CoordK(1,ik)-coordK(1,jL))**2+(CoordK(2,ik)-coordK(2,jL))**2) 
        If (abs(dKL(iseg))<1.d-6) then
           Zvernoi = .true.
           write(uprint,*)'SEGT Centre vernoi confondu dKL=0','iseg = ', iseg
        endif
!!!
     case (dirichlet) 
        !!
        is = NuSeg(1,iseg); js=NuSeg(2,iseg)
        ik= NumTVoisSeg(1,iseg) 

        ! milieu du segment 
        x1= (CoordS(1,is)+CoordS(1,js))/2.; y1= (CoordS(2,is)+CoordS(2,js))/2. 
        !dKL = sqrt((CoordK(1,ik)-x1)**2 +(CoordK(2,ik)-y1)**2) 

        If (abs(dKL(iseg))<1.d-6) then
           zvernoi= .true. 
           print*,'SEGT centre milieu segment sur le BORD ','iseg = ', iseg
        endif
     case (neumann)
     case default
        stop'SEGT pb segt  Vernoi'
     end select
  end Do
  !!
  !! verifier si le centre du triangle (centre de vernoi) est dans le triangle
  !!
  !! pour le maillage triangulaire
  If (maillageregulier == 0) then
     ncentrevernoi=0
     pi=4.D0*atan(1.D0)
     DO jt = 1,Nbt
        !! verifier si le point (xpoint(ii), ypoint(ii)) dans le triangle jt
        ! Compute vectors
        !! S1 = A=
        S1=coordS(:,NuSoK(1,jt)); S2=coordS(:,NuSoK(2,jt));S3 =coordS(:,NuSoK(3,jt)) 
        P=(/coordk(1,jt), coordk(2,jt)/)
        v0 = S1 - S3
        v1 = S2 - S3
        v2 = P - S3
        !
        Denom = v0(1)*v1(2)-v0(2)*v1(1)
        ubc = ( v1(2)*v2(1) -v1(1)*v2(2) )/Denom
        vbc = ( -v0(2)*v2(1)+v0(1)*v2(2) )/Denom
        ! Check if point is in triangle
        if (ubc < 0.D0 .or. vbc < 0.D0 .or. (ubc + vbc > 1.D0)) then 
           ncentrevernoi = ncentrevernoi+1
           write(uprint,*)'centre du triangle vernoi',jt,'pas dans le triangle'
           ! calcul des angles de ce triangle car il y a forcement un angle obtu
           ! formule de AL-Kashi
           ! un triangle de sommet A, B, C et les longueurs des cotes sont a, b, c
           ! le cote a est en face du sommet A
           ! A^=arccos((b^2+c^2-a^2)/2bc)
           as1 = sqrt(dot_product(S2-S3, S2-S3))
           as2 = sqrt(dot_product(S1-S3, S1-S3))
           as3 = sqrt(dot_product(S1-S2, S1-S2))
           angle1 = acos( (as2*as2+as3*as3-as1*as1)/(2.D0*as2*as3) )
           angle2 = acos( (as1*as1+as3*as3-as2*as2)/(2.D0*as1*as3) )
           angle3 = acos( (as2*as2+as1*as1-as3*as3)/(2.D0*as2*as1) )
           write(uprint,*)'angles : ', angle1, angle2,angle3, (angle1+angle2+angle3),(angle1+angle2+angle3)/2 
           write(uprint,*)'anglessdegre : ', angle1*180/pi, angle2*180/pi,angle3*180/pi, (angle1+angle2+angle3)*180/pi
        end if

     end DO
     write(uprint,*)'Le nombre de triangle non vernoi est ', ncentrevernoi
  endif







  !---------------------------------------------
  ! Fichier des segments a trace avec gnuplot
  !  plot 'SEGT.out' w l
  !---------------------------------------------
  !OPEN(unit = 2, file = 'SEGT.OUT',status='unknown')
  !DO iseg = 1, Nseg
  !   WRITE(2,300) CoordS(1, NuSeg(1, iseg)),  CoordS(2, NuSeg(1, iseg)), NtypSeg(iseg)
  !   WRITE(2,300) CoordS(1, NuSeg(2, iseg)),  CoordS(2, NuSeg(2, iseg)), NtypSeg(iseg)
  !END DO
  !CLOSE(2)
  !------------------------
!!$  ! Impressions eventuelles
!!$  !------------------------
  IF (iprint >=5) THEN 
     CALL prvari(uprint, 'Nseg = ', Nseg)
     Write(uprint,*)'iseg, NombVoisSeg(iseg), NumTVoisSeg(1,iseg), NumTVoisSeg(2,iseg), Ntyp'
     DO iseg = 1 , Nseg 
        WRITE(uprint,200) iseg, NombVoisSeg(iseg), NumTVoisSeg(1,iseg), NumTVoisSeg(2,iseg),NTypSeg(iseg)
     ENDDO
     !!
     CALL prvari(uprint, 'NBS = ', Nbs)
     CALL prvari(uprint, 'NsInt = ', NsInt)

     Write(uprint,*)'is, CoordS(1,is), CoordS(2,is), Ntyps(is)'
     DO is = 1 , Nbs
        WRITE(uprint,500) is, CoordS(1,is), CoordS(2,is), Ntyps(is)
     END DO
  ENDIF

  !!
  !! Verifier si le centre des triangles est dans le triangle
  !!
  If (Zvernoi) then 
     print*,' PB DE MAILLAGE VERNOI (segments)'
     !stop
  endif

300 FORMAT (2(E12.4,2x), I5)
113 FORMAT (5i8)
200 FORMAT (5(i8,2x))

500 FORMAT (I6, 2(E12.4, 2x), I5)
  !-----------------
  ! Fin du programme
  !----------------
  prefix = oldprf

END SUBROUTINE segments
