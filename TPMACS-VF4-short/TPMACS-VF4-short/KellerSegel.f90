program KellerSegel

  USE longr
  USE parmmage
  USE imprime
  USE intmatvec
  USE algebrelineaire
  USE intbigradc
  use plotvtkmod
  Use fsourcemod
  use parameters
  IMPLICIT NONE
  !==========================
  ! Declaration des tableaux
  !==========================
  TYPE(MatCreux)       :: A,N,Mat_E,A_vegf, Big_A
  type(Matcreux), dimension(:,:), allocatable :: Tab_A
  !==================================
  ! Declaration des variables locales
  !==================================
  INTEGER                             :: jt, i, j, is, jv,iseg,  kiter ,k,h,kplot
  REAL(kind=long), DIMENSION(:), ALLOCATABLE  :: U,U0,dU, U_p , Uexacte, Nutriment, N0 , Endothelial, Cells_density, Vasegf&
  &, Nutact, Tum
   real (kind=long), dimension(:,:), allocatable :: Tab_U, Tab_U_p
   character(len=6) , dimension(:), allocatable :: Tab_entity,Tab_equa
   integer, dimension(:), allocatable :: Tab_chemo
  REAL(kind = long)               :: tol, seuil, approx
  character(:), ALLOCATABLE :: str, strPb
  logical :: usemethodmul

  !===================
  ! Debut du programme
  !===================
  prefix = 'LAPLAC  '
  !
  CALL init                    !* Initialisation
  print*,  'init ok '

  print*, 'Keller Segel'

  ! Lecture du maillage � partir d'un fichier : MAILLAGEGEOx, x=1...6

  CALL readmesh
  !CALL readmatlab

  print*,  'readmesh ok '

  call NodeConnectivity
  print*, 'NodeConnectivity ok', NmaxT

  usemethodmul=.false.



  ALLOCATE(Nutriment(Nbt), N0(Nbt),Nutact(1:Nbt))
  allocate(Endothelial(Nbt),Vasegf(Nbt))
  allocate(Cells_density(Nbt))
  ALLOCATE(VitesseSeg(1:2,1:Nseg))
   allocate(Tab_A(n_enty,n_enty))

   do k=1,n_enty
      do h=1,n_enty
         call matrixinitVF4 ( Tab_A(k,h) ) 
      end do
   end do
   call bigmatrix(Big_A,Tab_A(1,1))

  do i=1,Nbt 
      Endothelial(i) = fsource( coordK(1,i), coordK(2,i), choixpb ) 
  end do

  ! resolution du syst�me lin�aire  par la mathode du gradient conjugu�
  tol = 5.D-12
  N0(:)=1.D0

   allocate(U(1:Nbt),U0(1:Nbt),dU(1:Nbt),U_p(1:Nbt),Tum(1:Nbt))
   U0 = 0.D0 
   do i=1,Nbt
      U0(i) = cond_ini(CoordK(1,i),CoordK(2,i))
      !Nutriment(i) = 0.5*fsource(CoordK(1,i),CoordK(2,i), ChoixPb)
   end do 
   U = U0

   call matrixinitVF4( N )
   call assembleNutri(N,Endothelial,U0)
   Nutriment = bigradient(N, N%Bg,N0,tol)
   call vide( N )
   call assembleNutri(N,U0,Endothelial)
   Vasegf = bigradient(N, N%Bg, N0,tol)

   allocate(Tab_U(n_enty,Nbt),Tab_U_p(n_enty,Nbt),Tab_entity(n_enty),Tab_equa(n_enty),Tab_chemo(n_enty))
   Tab_U(1,:) = U0;              Tab_U_p(1,:) = U0;            Tab_entity(1)='u_norm'; Tab_equa(1)='instat';
   Tab_chemo(1) = index_nut;
   Tab_U(2,:) = Nutriment;       Tab_U_p(2,:) = Nutriment;     Tab_entity(2)='Nutrim'; Tab_equa(2)='instat';
   Tab_chemo(2) = 0;
   !Tab_U(3,:) = Endothelial;     Tab_U_p(3,:) = Endothelial;   Tab_entity(3)='Endoth'; Tab_equa(3)='instat';
   !Tab_chemo(3) = index_vegf;
   !Tab_U(4,:) = Vasegf;          Tab_U_p(4,:) = Vasegf;        Tab_entity(4)='VasEGF'; Tab_equa(4)='instat';
   !Tab_chemo(4) = 0;

   Cells_density = Tab_U(index_norm,:)! + Tab_U(index_endo,:) 

   strPb = num2string(WhichPb)

   call plot_vtk(U0,'TC_Pb'//strPb//'_0','Tumor_cells')
   call plot_vtk(Nutriment,'N_Pb'//strPb//'_0','Nutrients')
   call plot_vtk(Cells_density,'CD_Pb'//strPb//'_0','Cells_density')
   !call plot_vtk(Tab_U(index_vegf,:),'VEGF_Pb'//strPb//'_0','VEGF')
   !call plot_vtk(Endothelial,'E_Pb'//strPb//'_0','Endothelial_cells')
   U = U0
   
   if (usemethodmul.eqv..true.) then
   print*, 'matrix init ok' 

   call MethNewtonMul(Tab_A,Tab_U,Tab_U_p,Tab_entity,Big_A,Tab_equa,Tab_chemo,'instat',tol,0)
   Tab_U_p=Tab_U
 print*,'methnewton ok'
 i=1
 kplot = 1
 do while(i*delta<=Tf)
   call MethNewtonMul(Tab_A,Tab_U,Tab_U_p,Tab_entity,Big_A,Tab_equa,Tab_chemo,'instat',tol,i*delta)
   Cells_density = Tab_U(index_norm,:) + Tab_U(index_endo,:) 
   Tab_U_p = Tab_U
   i=i+1
   if (modulo(i,10)==1) then
   print*,i*delta
   end if
   if (modulo(i,200)==1) then 
      str = num2string(kplot)
      call plot_vtk(Tab_U(index_nut,:),'N_Pb'//strPb//'_'//str,'Nutrients')
      call plot_vtk(Tab_U(index_norm,:),'TC_Pb'//strPb//'_'//str,'Tumor_cells')
     call plot_vtk(Cells_density,'CD_Pb'//strPb//'_'//str,'Cells_density')
     call plot_vtk(Tab_U(index_vegf,:),'VEGF_Pb'//strPb//'_'//str,'VEGF')
     call plot_vtk(Tab_U(index_endo,:),'E_Pb'//strPb//'_'//str,'Endothelial_cells')
      print*,i*delta
      kplot = kplot + 1
   end if
 end do

else 
   i=1
   kplot = 1
   do while(i*delta<=Tf) 
   !Actualisation des nutriments
   call assemInstaNutri(Nutact,Endothelial,U,Nutriment)
   !Actualisation des cellules tumorales
   call assemInstaTumor(Tum,Endothelial,U,Nutriment) 
   Nutriment = Nutact
   U = Tum
   
   !Les plots
   !if (modulo(i,10)==0) then
   !   print*,i*delta
   !end if
   if (modulo(i,5000)==0) then 
         str = num2string(kplot)
         call plot_vtk(Nutriment,'N_Pb'//strPb//'_'//str,'Nutrients')
         call plot_vtk(U,'TC_Pb'//strPb//'_'//str,'Tumor_cells')
         !call plot_vtk(Cells_density,'CD_Pb'//strPb//'_'//str,'Cells_density')
         print*,i*delta
         kplot = kplot + 1
   end if
         i=i+1

   end do 

end if


100 FORMAT(10(E10.3,2x))
200 FORMAT(6(E14.6,2x))
  CLOSE (uprint)
  PRINT*,'fin du travail laplacien'
  print*,'CHOIX PB = ', choixpb


end program KellerSegel 