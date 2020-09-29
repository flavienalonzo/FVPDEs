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
  uele      =  15  ! unite de lecture du maillage
  uneigh    =  16  ! unite de lecture du maillage
  unode     =  17  ! unite de lecture du maillage
  uplotvtk  =  20  ! fichier de sauvgarde

  
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
  
   nom_mesh =  'MAILLAGEBRAIN'   
   !nom_mesh = 'Tri_2D_0.05_5_1_2_@(x)[-x(2),x(1)].txt' 

  !--------------------------------------------------------------------

   iprint = 6  ! niveau d'impression


   n_enty = 4
   index_norm=1;index_nut=2;index_endo=3;index_vegf=4;
   !Pour le nutriment :
    Coef_diffusion = 1.D-3
    theta = 1.D-1
    Coef_prod = 5.D-2
    seuil_hypo = 8.D0
    seuil_necro = 6.D0
    satur_nutri = 1.D1
    nut_degra = 5.D-2

   !Pour les cellules tumorales :
    Diff_u = 1.D-4
   chi_u = 5.D-3
   rate = 1.D-1 
   apop = 5.D-2
   rat_pop = rate - apop
    Coef_cons = 5.D-1
    satur_norm = 0.D0
    choixkscalaireu = 1;
    deltau = 5.D-3; deltaxu = 3.5D-1; deltaxyu = 50.D-3; deltayu = 3.5D-1;

  !Pour les cellules endotheliales
    Diff_endo = 1.D-5
    chemo_endo = 5.D-1
    satur_endo = 0.D-1
    rate_endo = 5.D-2
    degr_endo = 1.D-2

    !Pour le VEGF
  VEGF_prod = Coef_prod
  VEGF_dif = Coef_diffusion
  VEGF_cons = Coef_cons
  VEGF_degr = nut_degra

  !Traitements
  !!Chirurgie 
  time_surg_d= 0.05 ; time_surg_f=0.06; seuil_surg=1.D-1 ;
  !!Chimiotherapie
  dose_chemo = 6.D-1; chemo_time(1,1:3) = (/1.0,2.5,4.0/); chemo_time(2,1:3) = (/1.4,2.9,4.4/)
  Nchemo = 3
  !!Radiotherapie
  Nradio = 6; radio_time(1,1:6) = (/1.0,1.4,2.5,2.9,4.0,4.4/); radio_time(2,1:6) = (/1.2,1.6,2.7,3.1,4.2,4.6/)
  dose_radio=60.D0 ; radio_beta=0.027D-1 ; radio_alpha =0.027D0 ;

  WhichPb = 33
  dt = 5.D-3
  Tf = 10.001
  TolerenceGradient = 1.D-13
  TolerenceNewton = 1.D-11
  ChoixPlot = 1 ; ChoixDegenere = 1 ; ChoixUpwind = 0 ; ChoixSchema = 2;
  choixanisu = 8 ; choixanis = 1; ChoixAdeg = 3 ; ChoixChi = 6;
  epsilon = 1.D-10
  
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

  ChoixPb = 99        



  !-----------------
  ! Fin du programme
  !-----------------
  prefix = oldprf
  RETURN

10 PRINT*,"Erreur dans l'entree des parametres"

  STOP


  RETURN
END SUBROUTINE INIT








