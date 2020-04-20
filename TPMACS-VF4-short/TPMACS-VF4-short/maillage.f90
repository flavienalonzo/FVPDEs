SUBROUTINE   maillage
  !------------------------------------------------
  ! nbs  : nombre des sommets total
  ! nbt  : nombre des triangle
  ! Ns   : nombre des sommets interieur au domaine 
  !        nombre des sommets ou la solutions est calculee
  ! Nb   : nombre des sommets au bord
  ! nbs = Ns + Nb 
  !----------------------------------
  ! Declaration des variables locales
  !----------------------------------
  USE longr
  USE imprime
  USE parmmage
  IMPLICIT NONE 
  CHARACTER(len=6)                                :: oldprf
  CHARACTER(len=10)                                :: chaine

  REAL (kind = long)                              :: errmach, r1, r2

  INTEGER :: i, j , kt, n, is, bidon,bidon1,bidon2,bidon3,Nbt2,typet, nums
  Integer ::  Dimespace,Typedom,Typsommet  
  !-------------------
  ! Debut du programme
  !-------------------
  oldprf = prefix
  prefix = 'MAILAG'
  IF (maillageregulier== 1) then 
     Typegeomaille = 4  ! type geometrie de la maille (rectangle)
     !! maillage rectangulaire
     ! affecter  Les numeros des sommets des triangles NuSoK
     ! affecter  Numeros des triangles voisins de chauqe triangles NuVois 
     ! affecter  les coordonnées de chaque sommet
     ! affecter  les centres de Voronoi de chaque triangle

     Nbt = Nx*Ny
     ! les rectnagles sont numeroté de i=1,Nx, j=1,NY 
     ! le numero dur rectangle est kt = i+(j-1)*Nx
     ! Les sommets sont numerotés de i=1,Nx+1; j=1, NY+1
     ! Les sommets d'un triangle kt sont :
     ! kt+j-1, kt+j, kt+(j-1)+Nx+1; kt+j+Nx+1 
     !
     ! Numeros des sommets des mailles


     ALLOCATE ( NuSoK(1:typegeomaille,1:Nbt) )
     DO i=1,Nx
        DO j=1,NY
           kt = i+(j-1)*Nx
           NuSoK(:,kt) = (/ kt+j-1, kt+j, kt+(j-1)+Nx+1, kt+j+Nx+1 /)
        END DO
     END DO


     ! ------------------------------------- 
     ! -- Liste des voisins de chaque maille
     ! -------------------------------------  
     ALLOCATE(NuVoisK(typegeomaille,nbt))
     DO i=1,Nx
        DO j=1,NY
           kt = i+(j-1)*Nx
           NuVoisK(:,kt) = (/ kt-1, kt+1, kt-Nx, kt+Nx /)
           ! Les mailles en dehors de la geometrie ne sont pas mis à -1
        END DO
     END DO

     ! ------------------------------------------------------------
     ! -- Les coordonnées de chaque sommet et le type du sommet
     ! ------------------------------------------------------------ 
     NbS = (Nx+1)*(Ny+1)
     ALLOCATE(coordS(2,nbs) , ntyps(1:nbs) )
     dx = LX/NX  ; dy = LY/NY
     DO i = 1,NX+1
        DO j = 1,NY+1
           nums = i + (j-1)*(NX+1)
           coordS(:,nums)= (/ (i-1)*dx, (j-1)*dy/)
        END DO
     END DO
     ntyps = 0
     !! conditions aux limites 
     !! coté bas et haut
     DO i = 1, Nx+1
        ntyps(i) = -1
        ntyps(i + NY*(Nx+1))=-2
     ENDDO
     !! cote gauche et droite 
     DO j = 2, NY
        ntyps(1+(j-1)*(Nx+1)) = -3
        ntyps(Nx+1+(j-1)*(Nx+1))=-4
     ENDDO

     ! ------------------------------------------------ 
     ! -- coordonnées du centroide de chaque maille  
     ! ------------------------------------------------ 

     ALLOCATE(coordK(2,nbt))
     DO i = 1,NX
        DO j = 1,NY
           kt = i + (j-1)*NX
           coordK(:,kt)= (/ (i-0.5D0)*dx, (j-0.5D0)*dy/)
        END DO
     END DO








  ELSE 
     Typegeomaille = 3 ! type geometrie de la maille (triangle)
     !! maillage triangulaire
     If (Maillage_Format == 1) then 
        !! maillage de format triangle 
        !! 
        ! ------------------------------------- 
        !-- lecture des numeros des sommets de chaque triangle nom.1.ele
        ! ------------------------------------- 
        OPEN (uele,     file = nom_ele   ,  status = 'old')

        READ (uele,*) nbt, bidon, typet

        CALL prvari(uprint,'nombre des traingles   (Nbt) ', Nbt )

        ALLOCATE ( NuSoK(1:3,1:Nbt) )

        DO kt = 1,Nbt
           READ(uele,*) n, ( NuSoK(j,kt),j=1,3 )
        END DO

        close(uele)

        ! ------------------------------------- 
        ! -- lecture de la liste des voisins de chaque triangle nom.1.neigh
        ! ------------------------------------- 
        OPEN (uneigh,     file = nom_neigh   ,  status = 'old')

        READ (uneigh,*) nbt2, bidon

        If (nbt2/=nbt) then 
           print*,'pb de fichier .neigh incompatible'
        end If
        CALL prvari(uprint,'nombre des traingles   (Nbt) ', Nbt )



        ALLOCATE(NuVoisK(3,nbt))

        DO kt = 1,Nbt
           READ(uneigh,*) n,(NuVoisK(j,kt),j=1,3)
        END DO
        !! attention quand est le numero est -1; c'est a dire pas de voisin

        close(uneigh)

        ! ---------------------------------------------------------------- 
        ! -- lecture des coordonnées de chaque sommet 
        ! Le fichier nom.1.node
        ! Dans le fichier on peut lire aussi 
        ! type du sommet, type du triangle, type du segment
        ! car un sommet appartient a tout cela a la fois CE N'EST PAS CLAIR
        ! -----------------------------------------------------------------
        ! 
        OPEN (unode,     file = nom_node   ,  status = 'old')

        READ (unode,*) Nbs, Dimespace,Typedom,Typsommet  

        CALL prvari(uprint,'nombre des sommets     (Nbs) ', Nbs )

        ALLOCATE(coordS(2,nbs) , ntyps(1:nbs) )
        Select Case(Typedom)
        case(0)
           DO is = 1,nbs
              READ(unode,*) n,(coordS(j,is),j=1,2), ntyps(is)
           END DO
        Case(1) 
           DO is = 1,nbs
              READ(unode,*) n,(coordS(j,is),j=1,2), r1, ntyps(is)
           END DO
        Case(2) 
           DO is = 1,nbs
              READ(unode,*) n,(coordS(j,is),j=1,2), r1, r2, ntyps(is)
           END DO
        end Select


        !! Il faut se mettre d'accord sur le numero de sommet

        close(unode)
        !====================
        ! ---------------------------------------------------------------- 
        ! -- lecture des coordonnées du centroide de chaque triangle  
        ! Le fichier nom.1.v.node    : c'est le fichier de Vernoui
        ! Dans le fichier on peut lire le type des triangle et  ....? 
        ! -----------------------------------------------------------------
        ! 
        OPEN (uvnode,     file = nom_v_node   ,  status = 'old')

        READ (uvnode,*) nbt2, bidon1,bidon2,bidon3  
        If (nbt2/=nbt) then 
           print*,'pb de fichier .v.node incompatible'
        end If

        ALLOCATE(coordK(2,nbt))

        DO kt = 1,nbt
           READ(uvnode,*) n,(coordK(j,kt),j=1,2)
        END DO

        !! on remet les coordonnees barycentriques ici 
        DO kt = 1,nbt
           CoordK(1,kt) = (CoordS(1,NuSoK(1,kt)) + CoordS(1,NuSoK(2,kt)) + CoordS(1,NuSoK(3,kt)) )/3.D0 
           CoordK(2,kt) = (CoordS(2,NuSoK(1,kt)) + CoordS(2,NuSoK(2,kt)) + CoordS(2,NuSoK(3,kt)) )/3.D0 
        END DO

        close(uvnode)
        !====================
        !------------------------
        ! Impressions eventuelles
        !------------------------
        IF (iprint >=2) THEN 
           write(uprint,*)' is ,coordS(1,is),coordS(2,is),ntyps(is)'
           DO is =1,nbs
              WRITE(uprint,112) is ,coordS(1,is),coordS(2,is),ntyps(is)
           ENDDO
           write(uprint,*)'  kt , (NuSoK(j,kt),j=1,3) '
           DO kt=1,nbt
              WRITE(uprint,113) kt , (NuSoK(j,kt),j=1,3)
           ENDDO
        ENDIF
     ENDIF

    !----- !!!!!!!!!!!!!!!!!!!!!
    !----- !!!!!!!!!!!!!!!!!!!!!
    !-----  !!!!!!!!!!!!!!!!!!!!!
     IF (Maillage_Format == 2) THEN
        ! fichier de type boyer
        ! ---------------------------------------------------------------- 
        ! -- lecture des coordonnées de chaque sommet 
        ! -----------------------------------------------------------------
        ! 
        OPEN (unode,     file = nom_node   ,  status = 'old')
        DO  kt=1,4
           READ(unode,*)
        ENDDO
        READ (unode,*)chaine, Nbs

        CALL prvari(uprint,'nombre des sommets     (Nbs) ', Nbs )

        ALLOCATE(coordS(2,Nbs) )

        DO is = 1,Nbs
           READ(unode,*) (coordS(j,is),j=1,2)
        END DO
        close(unode)
        ! ---------------------------------------------------------------- 
        ! -- lecture des coordonnées du centroide de chaque triangle  
        ! Le fichier nom.1.v.node    : c'est le fichier de Vernoui
        ! -----------------------------------------------------------------
        ! 
        OPEN (uvnode,     file = nom_v_node   ,  status = 'old')
        DO  kt=1,4
           READ(uvnode,*)
        ENDDO
        READ (uvnode,*)chaine, nbt

        !! pour ke schema p1 milieu, on a besoin du barycentre qui est different
        !! des points de voronoi
        ALLOCATE(coordK(2,nbt))

        DO kt = 1,nbt
           READ(uvnode,*) (coordK(j,kt),j=1,2)
        END DO


        close(uvnode)
        !! le type du sommet n'est pas utilisé pour l'instant
        !ALLOCATE(ntyps(1:nbs) )
        !ntyps = -9999
     END IF

  END IF
 IF (Maillage_Format == 3) THEN
    stop 'verifier ce format .amdba '
   OPEN (uma, file = nom_maillage ,  status = 'old')

  READ (uma,*) Nbs, Nbt
  CALL prvari(uprint,'nombre des sommets     (Nbs) ', Nbs )
  CALL prvari(uprint,'nombre des traingles   (Nbt) ', Nbt )

  ALLOCATE (coordS(2,1:Nbs) , ntyps(1:Nbs) )

  DO i = 1,Nbs
     READ(uma,*)n,coordS(1,i),coordS(2,i),ntyps(i)
  ENDDO

  ALLOCATE(NuSoK(3,nbt))

  DO i = 1,Nbt
     READ(uma,*) n,(NuSoK(j,i),j=1,3)
  END DO

  CLOSE(uma)

end if

  !------------------------
  ! Impressions eventuelles
  !------------------------
  IF (iprint >=2) THEN 
     write(uprint,*)' is ,coordS(1,is),coordS(2,is),ntyps(is)'
     DO is =1,nbs
        WRITE(uprint,112) is ,coordS(1,is),coordS(2,is)  !!,ntyps(is)
     ENDDO
     !write(uprint,*)'  kt , (NuSoK(j,kt),j=1,3) '
     !DO kt=1,nbt
     !   WRITE(uprint,113) kt , (NuSoK(j,kt),j=1,3)
     !ENDDO
  ENDIF


112 FORMAT (i5,2(E12.4,1x),i5)
113 FORMAT (5i8)
  !-----------------
  ! Fin du programme
  !----------------
  prefix = oldprf

END SUBROUTINE maillage



