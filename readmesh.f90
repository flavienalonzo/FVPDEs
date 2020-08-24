SUBROUTINE   readmesh
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
  CHARACTER(len=6)           :: oldprf
  CHARACTER(len=10)          :: chaine
  CHARACTER(len=len_buffer) :: buffer

  REAL (kind = long)        :: MaxaireK, MinaireK

  INTEGER ::  j, kt, n, is, iseg

  !-------------------
  ! Debut du programme
  !-------------------
  oldprf = prefix
  prefix = 'MESH'



  OPEN (umesh,     file = nom_mesh,  status = 'old')

  buffer=lireligne(umesh)
  READ(buffer, *, err = 10) Nbs
  CALL prvari (uprint, 'Nbs          : ', Nbs)

  buffer=lireligne(umesh)
  READ(buffer, *, err = 10) Nbt
  CALL prvari (uprint, 'Nbt          : ', Nbt)

  buffer=lireligne(umesh)
  READ(buffer, *, err = 10) Nseg
  CALL prvari (uprint, 'Nseg         : ', Nseg)

  buffer=lireligne(umesh)
  READ(buffer, *, err = 10) Nsint
  CALL prvari (uprint, 'Nsint        : ', Nsint)

  buffer=lireligne(umesh)
  READ(buffer, *, err = 10) Nbord
  CALL prvari (uprint, 'Nbord        : ', Nbord)


  buffer=lireligne(umesh)
  READ(buffer, *, err = 10) MaxaireK
  CALL prvarr (uprint, 'MaxaireK     : ', MaxaireK)

  buffer=lireligne(umesh)
  READ(buffer, *, err = 10) MinaireK
  CALL prvarr (uprint, 'MinaireK     : ', MinaireK)

  !! lecture des coordonnées sommets

  READ(umesh, *)

  ALLOCATE(coordS(2,Nbs) , ntyps(1:Nbs) )
  DO is = 1,nbs
     READ(umesh,*) n,(coordS(j,is),j=1,2), ntyps(is)
   !  write(*,*) n,(coordS(j,is),j=1,2), ntyps(is)
  END DO

  !! lecture des numeros des sommets pour chaque triangle, coord centre K, AireK

  READ(umesh, *)
  ALLOCATE ( NuSoK(1:3,1:Nbt), CoordK(1:2,Nbt), AireK(Nbt) )
  DO kt = 1,Nbt
     READ(umesh,*) n, ( NuSoK(j,kt),j=1,3 ), (CoordK(j,kt),j=1,2), AireK(kt)
    ! write(*,*) n, ( NuSoK(j,kt),j=1,3 ), (CoordK(j,kt),j=1,2), AireK(kt)
  END DO

  !! lecture des connectivités entre segments et voisins

  READ(umesh, *)
  !
  ! NuSeg : numero de deux sommets du segment iseg
  ! Nombre Voisin du segment iseg,
  ! 
  ALLOCATE ( NuSeg(1:2,1:Nseg),  NombVoisSeg(1:Nseg),  NumTVoisSeg(1:2,1:Nseg) )
  ALLOCATE ( NtypSeg(1:Nseg),  TauKL(1:Nseg), dKL(1:Nseg))

  DO iseg = 1,Nseg
     READ(umesh,*) n,  (NuSeg(j,iseg),j=1,2),  NombVoisSeg(iseg), (NumTVoisSeg(j,iseg),j=1,2),&
          NtypSeg(iseg),  TauKL(iseg), dKL(iseg) 
     !write(*,*)    n,  (NuSeg(j,iseg),j=1,2),  NombVoisSeg(iseg), (NumTVoisSeg(j,iseg),j=1,2),&
     !     NtypSeg(iseg),  TauKL(iseg), dKL(iseg)  
  END DO


  PRINT*,'ATTENTION LE TYPE DU BORD EST IMPOSE Dirichlet'
  WRITE(uprint,*)'ATTENTION LE TYPE DU BORD EST IMPOSE Dirichlet'
  DO iseg = 1,Nseg
     IF  (NTypSeg(iseg)/=0) NTypSeg(iseg) = Dirichlet
  END DO

  

  close(umesh)
  !====================
  !------------------------
  ! Impressions eventuelles
  !------------------------
  IF (iprint >=6) THEN 
     write(uprint,*)' is ,coordS(1,is),coordS(2,is),ntyps(is)'
     DO is =1,nbs
        WRITE(uprint,112)  is ,coordS(1,is),coordS(2,is),ntyps(is)
     ENDDO
     write(uprint,*)  'kt , (NuSoK(j,kt),j=1,3) '
     DO kt=1,nbt
        WRITE(uprint,113) kt , (NuSoK(j,kt),j=1,3)
     ENDDO
     
     write(uprint,*) 'iseg ,  (NuSeg(j,iseg),j=1,2)'
     !
     DO iseg = 1,Nseg
     write(uprint,*) iseg,  (NuSeg(j,iseg),j=1,2),  NombVoisSeg(iseg), (NumTVoisSeg(j,iseg),j=1,2),&
          NtypSeg(iseg),  TauKL(iseg), dKL(iseg)  
  END DO
  ENDIF



112 FORMAT (i5,2(E12.4,1x),i5)
113 FORMAT (5i8)
  !-----------------
  ! Fin du programme
  !----------------
  prefix = oldprf

  RETURN

10 PRINT*,"Erreur dans l'entree des parametres"
  
  STOP

  RETURN
END SUBROUTINE readmesh



