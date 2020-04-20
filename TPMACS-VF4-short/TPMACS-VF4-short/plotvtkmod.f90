module plotvtkmod

use longr
use imprime
use parmmage

implicit none

contains

  subroutine plot_vtk (vec,chaine,nomchamps)

    real(kind=long), dimension(:),intent(in) :: vec
    character(len=*),intent(in) :: chaine,nomchamps

    Integer :: kt,Lt,is,js,ks,iseg, jseg, kseg, countk, i
    real(kind=long), dimension(:),allocatable :: WW
    Integer, dimension(:,:),allocatable :: NuSoDiam

    print*,"Creation du fichier d'impression"
    uplotvtk = 63
    FPLOTVTK = chaine//'.vtk'

    open (unit=uplotvtk,file=FPLOTVTK,status='replace')
    write(uplotvtk,'(A)') '# vtk DataFile Version 3.0'
    write(uplotvtk,'(A)') 'LAPLACIEN  2D'
    write(uplotvtk,'(A)') 'ASCII'
    write(uplotvtk,'(A)') 'DATASET UNSTRUCTURED_GRID'



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
       write (uplotvtk,500) vec(kt)
    END DO

    close(uplotvtk)


    print*,"OK"
    print*," "

400 format (E10.5)
500 format (E30.20)
600 format (4(I6,3X))

  end subroutine plot_vtk




end module plotvtkmod
