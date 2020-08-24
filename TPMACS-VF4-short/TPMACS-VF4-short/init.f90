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
  
   nom_mesh =  'MAILLAGEGEO'   
   !nom_mesh = 'Tri_2D_0.05_5_1_2_@(x)[-x(2),x(1)].txt' 

  !--------------------------------------------------------------------

   iprint = 6  ! niveau d'impression


   n_enty = 2
   index_norm=1;index_nut=2;index_endo=3;index_vegf=4;
   !Pour le nutriment :
    Coef_diffusion = 1.D-3
    theta = 1.D-1
    Coef_prod = 5.D0
    seuil_hypo = 1.D-1
    seuil_necro = 1.D-2
    satur_nutri = 1.D1

   !Pour les cellules tumorales :
    Diff_u = 1.D-4
   chi = 1.D-2
   rate = 6.D-2 
   apop = 1.D-2
   rat_pop = rate - apop
    Coef_cons = 5.D-1
    satur_norm = 0.0D0

  !Pour les cellules endotheliales
    Diff_endo = 1.D-5
    chemo_endo = 1.D-1
    satur_endo = 4.D-1
    rate_endo = 0.D0

    !Pour le VEGF
  VEGF_prod = Coef_prod
  VEGF_dif = Coef_diffusion
  VEGF_cons = Coef_cons
  VEGF_degr = theta


  WhichPb = 29
  delta = 1.D-5
  Tf = 20.1

  
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








