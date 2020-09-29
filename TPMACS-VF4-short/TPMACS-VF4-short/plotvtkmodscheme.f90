module plotvtkmodscheme

  Use longr
  Use imprime
  Use parmmage

  Implicit none

contains

  Subroutine plot_vtk_scheme (vec,chaine,nomchamps, typesh)

    Real(kind=long), dimension(:),intent(in)   :: vec
    Character(len=*),intent(in)                :: chaine,nomchamps
    Integer, intent(in)                        :: typesh

    Integer                    :: kt, Lt, is, js, ks, iseg, jseg, kseg, countk, i, im, jm, km
    Real(kind=long), dimension(:), allocatable :: WW
    Integer, dimension(:,:), allocatable       :: NuSoDiam

    !print*,"Creation du fichier d'impression"
    uplotvtk = 69
    FPLOTVTK = chaine//'.vtk'

    open (unit=uplotvtk,file=FPLOTVTK,status='replace')
    write(uplotvtk,'(A)') '# vtk DataFile Version 3.0'
    write(uplotvtk,'(A)') 'LAPLACIEN  2D'
    write(uplotvtk,'(A)') 'ASCII'
    write(uplotvtk,'(A)') 'DATASET UNSTRUCTURED_GRID'


    If (maillageregulier == 1)  then 
       ! maillage rectangulaire
       write(uplotvtk,*) 'POINTS',NBS,' float'
       do is = 1, NBS
          write (uplotvtk,*) CoordS(1,is), CoordS(2,is), 0.D0
       end do

       write(uplotvtk,*) 'CELLS ',Nbt, 5*Nbt

       do kt=1,Nbt
          write (uplotvtk,*) 4,NuSoK(1,kt)-1, NuSoK(2,kt)-1, NuSoK(3,kt)-1, NuSoK(4,kt)-1
       end do

       write(uplotvtk,*) 'CELL_TYPES ',Nbt
       DO kt=1, Nbt
          write(uplotvtk,*) 8
       END DO

       WRITE(uplotvtk,*) 'CELL_DATA',Nbt 

       WRITE(uplotvtk,*) 'SCALARS ',nomchamps,' float'
       WRITE(uplotvtk,*) 'LOOKUP_TABLE default'

       DO kt=1, Nbt
          write (uplotvtk,500) vec(kt)
       END DO

       close (uplotvtk)


    else 

       select case ( typesh)
       case ( VF4 )
          ! maillage triangulaire
          write(uplotvtk,*) 'POINTS',NBS,' float'
          do is = 1, NBS
             write (uplotvtk,*) CoordS(1,is), CoordS(2,is), 0.D0
          end do

          write(uplotvtk,*) 'CELLS ',Nbt, 4*Nbt

          do kt=1,Nbt
             write (uplotvtk,*) 3,NuSoK(1,kt)-1, NuSoK(2,kt)-1, NuSoK(3,kt)-1
          end do

          write(uplotvtk,*) 'CELL_TYPES ',Nbt
          DO kt=1, Nbt
             write(uplotvtk,*) 5
          END DO

          WRITE(uplotvtk,*) 'CELL_DATA',Nbt 

          WRITE(uplotvtk,*) 'SCALARS ',nomchamps,' float'
          WRITE(uplotvtk,*) 'LOOKUP_TABLE default'

          DO kt=1, Nbt
             write (uplotvtk,500) max(vec(kt), 1.0D-30)
          END DO

       case ( P1Sommets )
          select case (ChoixPlot)
          case ( 0 )
             !----------------------------------!
             ! maillage formé par les triangles ! 
             !----------------------------------!

             write(uplotvtk,*) 'POINTS',NBS,' float'

             do is = 1, Nbs
                write (uplotvtk,*) CoordS(1,is), CoordS(2,is), 0.D0
             end do

             write(uplotvtk,*) 'CELLS ',Nbt, 4*Nbt

             do kt=1,Nbt
                write (uplotvtk,*) 3,NuSoK(1,kt)-1, NuSoK(2,kt)-1, NuSoK(3,kt)-1
             end do

             write(uplotvtk,*) 'CELL_TYPES ',Nbt
             DO kt=1, Nbt
                write(uplotvtk,*) 5
             END DO

             WRITE(uplotvtk,*) 'CELL_DATA',Nbt 

             WRITE(uplotvtk,*) 'SCALARS ',nomchamps,' float'
             WRITE(uplotvtk,*) 'LOOKUP_TABLE default'

             !! la solution est calculee aux sommets, 
             !! on fait une interpolation pour donner la solution au centre
             ALLOCATE(WW(Nbt))
             WW = 0.D0
             DO kt = 1,Nbt
                is = NuSoK(1,kt); js = NuSoK(2,kt); ks = NuSoK(3,kt)
                countk = 0
                if (is<=Nsint) then
                   WW(kt) = WW(kt)+vec(is) ; countk=countk+1
                endif
                if (js<=Nsint) then
                   WW(kt) = WW(kt)+vec(js) ; countk=countk+1
                endif
                if (ks<=Nsint) then
                   WW(kt) = WW(kt)+ vec(ks) ; countk=countk+1
                endif
                If (countk == 0) then 
                   !write(*,*) ' triangle entier sur le bord : PLOTVTKMOD'
                   !write(*,*)'ww(',kt,')=',ww(kt), 'countk = ', countk 
                   ww(kt)= -10.
                else
                   WW(kt)=WW(kt)/countk
                   !!write(*,*)'ww(',kt,')=',ww(kt), 'countk = ', countk 
                end If
             END DO

             DO kt = 1, Nbt
                write (uplotvtk,500)  max(WW(kt), 1.D-30)   !!WW(kt)
             END DO
             deallocate(WW)
             close(uplotvtk)

          case ( 1 )
             !-----------------------------------------------!
             ! maillage formé par les sous-mailles de Donald ! 
             !-----------------------------------------------!

             write(uplotvtk,*) 'POINTS',Nbs+Nseg+Nbt,' float'

             ! 1. Construction de sous-mailles de Donald:
             !-------------------------------------------

             ! 1.1. Definir les coordonnes de sommets des sous-mailles: 

             DO is = 1, Nbs
                write (uplotvtk,*) CoordS(1,is), CoordS(2,is), 0.D0
             END DO

             DO iseg = 1, Nseg
                is = NuSeg(1,iseg); js=NuSeg(2,iseg)
                write (uplotvtk,*) (CoordS(1,is)+CoordS(1,js))/2.D0, (CoordS(2,is)+CoordS(2,js))/2.D0, 0.D0
             END DO

             DO kt = 1, Nbt
                write (uplotvtk,*) CoordK(1,kt), CoordK(2,kt), 0.D0
             END DO

             ! 1.2. Definir le numero de sommets des sous-mailles: NuSoDiam (1:4, 3*Nbt) 
             !--------------------------------------------------------------------------

             Allocate (NuSoDiam (1:4, 3*Nbt) )
             do kt = 1, Nbt
                is = NuSoK(1,kt);  js = NuSoK(2,kt);  ks = NuSoK(3,kt)
                im = NuMSeg(1,kt); jm = NuMSeg(2,kt); km = NuMSeg(3, kt)

                NuSoDiam (1, kt) = is
                NuSoDiam (2, kt) = Nbs+jm
                NuSoDiam (3, kt) = Nbs+Nseg+kt
                NuSoDiam (4, kt) = Nbs+km

                NuSoDiam (1, kt+Nbt) = js
                NuSoDiam (2, kt+Nbt) = Nbs+km
                NuSoDiam (3, kt+Nbt) = Nbs+Nseg+kt
                NuSoDiam (4, kt+Nbt) = Nbs+im

                NuSoDiam (1, kt+2*Nbt) = ks
                NuSoDiam (2, kt+2*Nbt) = Nbs+im
                NuSoDiam (3, kt+2*Nbt) = Nbs+Nseg+kt
                NuSoDiam (4, kt+2*Nbt) = Nbs+jm 
             end do

             ! 1.3. Joindre les sommets de chaque sous-maille:
             !-------------------------------------------------

             write(uplotvtk,*) 'CELLS ',3*Nbt, 15*Nbt
             DO kt = 1, Nbt          
                write (uplotvtk,*) 4, NuSoDiam(1,kt)-1, NuSoDiam(2,kt)-1, &
                     & NuSoDiam(3,kt)-1, NuSoDiam(4,kt)-1
             END DO
             DO kt = Nbt+1, 2*Nbt          
                write (uplotvtk,*) 4, NuSoDiam(1,kt)-1, NuSoDiam(2,kt)-1, &
                     & NuSoDiam(3,kt)-1, NuSoDiam(4,kt)-1
             END DO
             DO kt = 2*Nbt+1, 3*Nbt          
                write (uplotvtk,*) 4, NuSoDiam(1,kt)-1, NuSoDiam(2,kt)-1, &
                     & NuSoDiam(3,kt)-1, NuSoDiam(4,kt)-1
             END DO
             Deallocate (NuSoDiam)

             write(uplotvtk,*) 'CELL_TYPES ',3*Nbt

             DO i = 1, 3*Nbt
                write(uplotvtk,*) 9
             END DO

             WRITE(uplotvtk,*) 'CELL_DATA', 3*Nbt 
             WRITE(uplotvtk,*) 'SCALARS ',nomchamps,' float'
             WRITE(uplotvtk,*) 'LOOKUP_TABLE default'

             ! 2. Attribuer les valeurs aux sous-mailles:
             !-------------------------------------------

             Allocate (WW(3*Nbt))
             ! Pour les sous-mailles admettant des sommets a l''interieur
             DO kt = 1,Nbt
                is = NuSoK(1,kt); js = NuSoK(2,kt); ks = NuSoK(3,kt)
                If (is<=Nsint)   WW(kt) = vec(is)
                If (js<=Nsint)   WW(kt+Nbt) = vec(js)
                If (ks<=Nsint)   WW(kt+2*Nbt) = vec(ks)
             END DO
             ! Pour les sous-mailles admettant des sommets sur le bord
             DO kt = 1,Nbt
                countk = 0
                DO i = 1, 3 
                   If (NuSoK(i,kt) > NsInt ) countk = countk+1
                END DO
                is = NuSoK(1,kt); js = NuSoK(2,kt); ks = NuSoK(3,kt)
                select case (countk)
                case(1)
                   If (is > NsInt)    WW(kt) = (vec(js)+vec(ks))/2.D0
                   If (js > NsInt)    WW(kt+Nbt) = (vec(is)+vec(ks))/2.D0
                   If (ks > NsInt)    WW(kt+2*Nbt) = (vec(is)+vec(js))/2.D0
                case(2)
                   if (is <= NsInt) then
                      WW(kt+Nbt)    = vec(is)
                      WW(kt+2*Nbt)  = vec(is)
                   endif
                   if (js <= NsInt) then
                      WW(kt)        = vec(js)
                      WW(kt+2*Nbt)  = vec(js)
                   endif
                   if (ks <= NsInt) then
                      WW(kt)        = vec(ks)
                      WW(kt+Nbt)    = vec(ks)
                   endif
                case(3)
                WW(kt) = 0.D0
                WW(kt+Nbt) = 0.D0
                WW(kt+2*Nbt) = 0.D0
                   !!stop 'Triangle entier sur le bord, Probleme d''interpolation'
                end select
             END DO

             ! 3. Creation des fichiers pour visit:
             !-------------------------------------

             DO kt = 1, Nbt
                write (uplotvtk,500)  max(WW(kt), 1.D-30) 
             END DO
             DO kt = 1, Nbt
                write (uplotvtk,500)  max(WW(kt+Nbt), 1.D-30) 
             END DO
             DO kt = 1, Nbt
                write (uplotvtk,500)  max(WW(kt+2*Nbt), 1.D-30) 
             END DO
             deallocate(WW)
             close (uplotvtk)
          end select
       case(P1Milieux, P1Milieuxmonotone)

          ! maillage formé des diamand 

          write(uplotvtk,*) 'POINTS',Nbs+Nbt+Nsegbord,' float'
          !!ms           write(uplotvtk,*) 'POINTS',Nbs+Nbt,' float'

          ! on est obligé d'ajouter un sommet fictif sur le bord
          ! maillage diamand
          ! 
          ! Construction des diamaond : 
          ! NuSDiam (1:4, Nseg) ; 

          ! numerotation des sommets des diama,ds voir NuSoDiam
          do is = 1, Nbs
             write (uplotvtk,*) CoordS(1,is), CoordS(2,is), 0.D0
          end do
          do kt = 1, Nbt
             write (uplotvtk,*) CoordK(1,kt), CoordK(2,kt), 0.D0
          end do

          countk=0
          DO iseg = 1, Nseg
             Lt = NumTVoisSeg(2,iseg)
             !!write(*,*)'NumTVoisSeg(2,iseg)', NumTVoisSeg(2,iseg)
             if (Lt < 0) then
                !!write (uplotvtk,*) CoordMseg(1,iseg), CoordMseg(2,iseg), 0.D0
                is = NuSeg(1,iseg)
                write (uplotvtk,*) CoordS(1,is), CoordS(2,is), 0.D0
                countk = countk+1
             end if
          END DO
          write(*,*)'nb seg bord via plot est ',countk

          ! NuSDiam (1:4, Nseg) ; 
          ! on numerote d'abord les sommets du maillage,
          !  ensuite on translate les numeros des traingles de Nbs et 
          !  enfin on ajoute un sommet au demi diamand 
          Allocate (NuSoDiam (1:4, Nseg) )
          countk = 1 
          do iseg = 1, Nseg
             is = NuSeg(1,iseg); js=NuSeg(2,iseg)
             kt = NumTVoisSeg(1,iseg); Lt= NumTVoisSeg(2,iseg)
             if (Lt > 0) then
                NuSoDiam (1, iseg)= is
                NuSoDiam (2, iseg)= kt+Nbs
                NuSoDiam (3, iseg)= js
                NuSoDiam (4, iseg)= Lt+Nbs            
             else
                NuSoDiam (1, iseg)= is
                NuSoDiam (2, iseg)= kt+Nbs
                NuSoDiam (3, iseg)= js
                NuSoDiam (4, iseg)= Nbt+Nbs + countk
                countk = countk+1 
             end if
          end do

          write(uplotvtk,*) 'CELLS ',Nseg, 5*Nseg
          DO iseg = 1,Nseg          
             write (uplotvtk,*) 4, NuSoDiam(1,iseg)-1, NuSoDiam(2,iseg)-1, &
                  & NuSoDiam(3,iseg)-1, NuSoDiam(4,iseg)-1
          END DO
          Deallocate (NuSoDiam)

          write(uplotvtk,*) 'CELL_TYPES ',Nseg

          DO iseg=1, Nseg
             write(uplotvtk,*) 9
          END DO

          WRITE(uplotvtk,*) 'CELL_DATA',Nseg 
          WRITE(uplotvtk,*) 'SCALARS ',nomchamps,' float'
          WRITE(uplotvtk,*) 'LOOKUP_TABLE default'

          !! affectation des condions aux bords

          Allocate (WW(Nseg))

          Do iseg = 1, Nsegint
             WW(iseg) = vec(iseg)
          End do

          !! interpolation des conditions aux bords

          DO kt = 1, Nbt
             countk = 0
             DO i = 1, 3 
                If (NuMSeg(i,kt) > Nsegint ) countk = countk+1
             END DO
             iseg= NuMSeg(1,kt) ; jseg= NuMSeg(2,kt) ; kseg= NuMSeg(3,kt)
             select case (countk)
             case(1)     ! iseg sur le bord
                If (iseg > Nsegint)  WW(iseg) = (vec(jseg)+vec(kseg))/2.D0
                If (jseg >Nsegint)   WW(jseg) = (vec(iseg)+vec(kseg))/2.D0
                If (kseg >Nsegint)   WW(kseg) = (vec(iseg)+vec(jseg))/2.D0
             case(2)    ! iseg et jseg sur le bord
                if (iseg <=Nsegint) then 
                   WW(jseg)= vec(iseg) ; WW(kseg)= vec(iseg)
                end if
                if (jseg <=Nsegint) then 
                   WW(iseg)= vec(jseg) ; WW(kseg)= vec(jseg)
                end if
                if (kseg <=Nsegint) then 
                   WW(iseg)= vec(kseg) ; WW(jseg)= vec(kseg)
                end if
             case(3)
                WW(iseg) = 0.D0
                WW(jseg) = 0.D0
                WW(kseg) = 0.D0
                !! stop 'pb dans interpolation bord dans plot'
             end select
          END DO

          !         DO iseg = Nsegint+1, Nseg
          !            WW(iseg) = 0.D0 
          !         END DO

          DO iseg = 1, Nseg
             write (uplotvtk,500)  max(WW(iseg), 1.D-30)   !! max(max(WW(iseg),0.D0),1.D-30)
          END DO
          deallocate(WW)
          close (uplotvtk)
       case default 
          print*, 'pb shema dans CREATION PLOT'
       end select
    endif

   ! print*,"OK"
   ! print*," "

400 format (E10.5)
500 format (E30.20)

  end subroutine plot_vtk_scheme




end module plotvtkmodscheme
