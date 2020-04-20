module parmmage
  !------------------------------------------------------------------
  ! On definit ici les parametres du maillage
  !------------------------------------------------------------------
  !  Nbs                 : nombres des sommets du maillages
  !  Nbt                 : nombres des triangles du maillage
  !  Nseg                : Nombre de SEGments
  !  coordS(1:2,1:nbs)   : les coordonnees des sommets du maillage
  !  coordK(1:2,1:nbt)   : les coordonnees des centroides du maillage
  !  NuSoK (1:3,1:nbt)   : les Numeros des SOmmets des Triangles
  ! NumTVoisSeg          : Numero des Triangles voisins du segment
  !------------------------------------------------------------------
  use longr
  implicit none

  REAL (kind = long), DIMENSION(:,:), POINTER :: CoordS, CoordK
  INTEGER,            DIMENSION(:,:), POINTER   :: NuSoK
  Integer , DIMENSION(:,:), ALLOCATABLE         :: NuMSeg
  INTEGER,            DIMENSION(:,:), POINTER   :: NuVoisK
  INTEGER,            DIMENSION(:)  , POINTER   :: ntyps,ntypt, ntypseg, NombVoisSeg
  INTEGER                                       :: Nbs,Nbt,Ns,Nb, Nseg,NsInt,Nbord,Nsegint,Nsegbord 
  INTEGER,            DIMENSION(:,:), POINTER :: NuSeg, NumTVoisSeg
  REAL (kind = long), DIMENSION(:),   ALLOCATABLE ::  TauKL, dKL, AireK



 !------------------------------------------------------------------- 
 ! Stockage du systeme lineaire a resoudre
 !-------------------------------------------------------------------
 ! IndPL      : INDice du Premier non coef non nul de la Ligne i
 ! Indc       : INDice de la colonne du coef (j)
 ! TMat       : les coefficients non nuls de la matrice stokes alors 
 !            : par ligne et puis ordre croissant par colonne
 !  Diag      : Extrait la diagonale de A
 ! PosiDiag   : dans la position de la diagonale dans TMat
 ! F          : le second membre physique 
 ! Bg         : le second membre du systeme a resoudre Tmat * x = bg
 !            : Bg = F + Gb   ;  Bg la donnee au bord physique  
 !-------------------------------------------------------------------
  TYPE MatCreux
     INTEGER, DIMENSION(:), POINTER :: IndPL, Indc, PosiDiag
     REAL (kind = long), DIMENSION(:), POINTER    :: TMat
     REAL (kind = long), DIMENSION(:), POINTER    :: Bg, F, Diag
  END TYPE MatCreux
end module parmmage


