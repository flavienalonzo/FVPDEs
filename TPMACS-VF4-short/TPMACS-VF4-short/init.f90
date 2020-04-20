SUBROUTINE  INIT
  !****************************************************************
  !     * Ce sous programme lit le fichier uread
  !     * le fichier d'entree
  !****************************************************************
  !--------
  ! Modules
  !--------
  USE longr
  USE imprime
  USE parmmage

  IMPLICIT NONE

  !--------------------------
  ! Declaration des arguments
  !--------------------------
  !----------------------------------
  ! Declaration des variables locales
  !----------------------------------
  CHARACTER(len=6)      :: oldprf
  INTEGER               :: i
  REAL(kind = long)    :: gbord
  CHARACTER(len=len_buffer) :: buffer
  !-------------------
  ! Debut du programme
  !-------------------
  oldprf = prefix
  prefix = 'INIT'

  uread   = 9   ! fichier de donnee pourra etre ajoute par la suite 
  uprint = 12   ! Unite de listage pour verifier les donnees
  umesh= 39   ! unite  fichier maillage

  
  ! Le fichier UPRINT sert à l'ecriture de certainses données pour valider le code
  open(unit=uprint, file='UPRINT',status='unknown')

  !---------------------------------------------------------------------
  !     Donnees pour le probleme  de laplace sur un domaine qcq
  !     par  la methode des volumes finis sur un maillage admissible.
  !     equation :
  !       - c div( grad u) +theta * u + div(Vu)= F    dans omega
  !             u = Gbord        sur le bord
  !


  !-----------------------------------------------
  ! maillage triangulaire : fichier contenant le maillage  et les coonectivités

  ! il y a six maillages dispo 'MAILLAGEGEOx', x=1,...,6
  
   nom_mesh =  'MAILLAGEGEO4'   
   !nom_mesh = 'Tri_2D_0.05_5_1_2_@(x)[-x(2),x(1)].txt' 

  !--------------------------------------------------------------------

   iprint = 6  ! niveau d'impression


   n_enty = 3
   !Pour le nutriment :
   Coef_diffusion = 0.0001D0
   theta = 0.05D0
   Coef_prod = 0.05D0

   !Pour les cellules tumorales :
    Diff_u = 0.05D0
   chi = 0.005D0
   rate = 0.1D0 
    Coef_cons = 0.05D0

   delta = 0.001
   Tf = 30.1

  
  !--------------------------------------------
  ! Choix du probleme à resoudre (ChoixPb)
  ! ChoixPb = 1, Uexacte =  1.
  ! ChoixPb = 2, Uexacte = x+y          
  ! ChoixPb = 3, Uexacte = x*x - y*y    
  ! ChoixPb = 4, Uexacte =  cos(5.*pi*(x+y)) 
  ! ChoixPb = 5, Uexacte =  x(1-x)y(1-y)
  ! Choixpb = 6, Uexacte = sin(pi*x)sin(pi*y)
  !----------------------------------------------
  !
  ! La fonction gbord contient la solution exacte 

  ChoixPb = 10        



  !-----------------
  ! Fin du programme
  !-----------------
  prefix = oldprf
  RETURN

10 PRINT*,"Erreur dans l'entree des parametres"

  STOP


  RETURN
END SUBROUTINE INIT








