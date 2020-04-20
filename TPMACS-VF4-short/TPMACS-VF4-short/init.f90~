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
  
   nom_mesh =  'MAILLAGEGEO3'    

  !--------------------------------------------------------------------

   iprint = 6  ! niveau d'impression

   Coef_diffusion = 1.D0

   theta = 0.D0

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

  ChoixPb = 1        



  !-----------------
  ! Fin du programme
  !-----------------
  prefix = oldprf
  RETURN

10 PRINT*,"Erreur dans l'entree des parametres"

  STOP


  RETURN
END SUBROUTINE INIT








