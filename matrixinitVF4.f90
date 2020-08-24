!            **************************
!            **  SUBROUTINE MATRIXINIT**
!            **************************
!******************************************************************************
!     * M. Saad ECN
!     *--------------------------------------
!     * Ce sous programme donne la structure generale de 
!     * de la matrice A%TMAT ainsi le remplissage de 
!     * A%IndPL  
!     * A%Indc
!****************************************************************************
SUBROUTINE matrixinitVF4(Mat)
  !--------
  ! Modules
  !--------
  USE longr
  USE parmmage
  USE imprime
  IMPLICIT NONE

  !--------------------------
  ! Declaration des arguments
  !--------------------------
  TYPE(MatCreux),           INTENT(OUT) :: Mat
  !----------------------------------
  ! Declaration des variables locales
  !----------------------------------
  CHARACTER(len=6)      :: oldprf
  INTEGER                             :: iseg, is, js, i1, jv, kv,jmin,mloc
  INTEGER                             :: NcoefMat
  INTEGER, DIMENSION(Nbt)              :: IndPL

  !-------------------
  ! Debut du programme
  !-------------------
  oldprf = prefix
  prefix = 'MATINIT'
  !------
  ! Corps
  !------
  NcoefMat = Nbt

  DO iseg = 1, Nseg
     IF(NtypSeg(iseg) == 0)   NcoefMat = NcoefMat + 2
  END DO
  ALLOCATE(Mat%IndPL(1:Nbt+1), Mat%Indc(NcoefMat), Mat%TMat(NcoefMat))
  ALLOCATE( Mat%F ( Nbt ), Mat%Bg ( Nbt ) )
   
  Mat%TMat = 0.0D0 ;  Mat%F = 0.D0 ;  Mat%Bg = 0.D0

  DO is = 1, Nbt+1
     Mat%IndPL(is) = is
  END DO

  DO iseg = 1, Nseg

     IF(NtypSeg(iseg)==0) then 
        is = NumTVoisSeg(1,iseg) ; js = NumTVoisSeg(2, iseg)
        Mat%IndPL(is+1:Nbt+1) = Mat%IndPL(is+1:Nbt+1) + 1
        Mat%IndPL(js+1:Nbt+1) = Mat%IndPL(js+1:Nbt+1) + 1
     end IF
  END DO

  IndPL(1:Nbt) = Mat%IndPL(1:Nbt)

  DO is = 1, Nbt
     Mat%Indc( IndPL(is) ) = is
     IndPL(is)           = IndPL(is) + 1
  ENDDO

  DO iseg = 1, Nseg
     IF (NtypSeg(iseg)==0) then 
        is = NumTVoisSeg(1,iseg) ;  js = NumTVoisSeg(2, iseg)
        Mat%Indc(IndPL(is)) = js
        Mat%Indc(IndPL(js)) = is
        IndPL(is)           = IndPL(is) + 1
        IndPL(js)           = IndPL(js) + 1
     end IF
  END DO

  DO is = 1, Nbt 
     DO jv = Mat%IndPL(is), Mat%IndPL(is+1)-1 
        jmin = jv
        do   kv = jv +1,  Mat%IndPL(is+1)-1
           IF (Mat%Indc(kv) < Mat%Indc(jmin)) jmin = kv
        END DO
        mloc =  Mat%Indc(jmin )
        Mat%Indc(jmin ) = Mat%Indc(jv ) 
        Mat%Indc(jv ) = mloc 
     END DO
  END DO
  !------------
  ! Impressions
  !------------
  IF (iprint >=5) THEN
     CALL prvari(uprint,'NcoefMat = ', NcoefMat )
     WRITE(uprint,*) 'matrixinit : MAT%TMAT est '
     WRITE(uprint,*) ( Mat%Tmat(is) , is=1, NcoefMat )
     WRITE(uprint,*) 'matrixinit : MAT%IndC est '
     WRITE(uprint,*) ( Mat%Indc(is) , is=1, NcoefMat )
     WRITE(uprint,*) 'matrixinit : MAT%IndPL est '
     WRITE(uprint,*) ( Mat%IndPL(is), is=1, Nbt+1 )
  ENDIF
  !-----------------
  ! Fin du programme
  !-----------------
  prefix = oldprf

END SUBROUTINE matrixinitVF4


  
