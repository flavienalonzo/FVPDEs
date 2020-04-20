program KellerSegel

  USE longr
  USE parmmage
  USE imprime
  USE intmatvec
  USE algebrelineaire
  USE intbigradc
  use plotvtkmod
  Use fsourcemod
  IMPLICIT NONE
  !==========================
  ! Declaration des tableaux
  !==========================
  TYPE(MatCreux)       :: A,N,Mat_E
  type(Matcreux), dimension(:), allocatable :: Tab_A
  !==================================
  ! Declaration des variables locales
  !==================================
  INTEGER                             :: jt, i, j, is, jv,iseg,  kiter
  REAL(kind=long), DIMENSION(:), ALLOCATABLE  :: U,U0,dU, U_p , Uexacte, Nutriment, N0 , Endothelial, Cells_density
   real (kind=long), dimension(:,:), allocatable :: Tab_U, Tab_U_p
   character(len=6) , dimension(:), allocatable :: Tab_entity
  REAL(kind = long)               :: tol, seuil

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

  ALLOCATE(Nutriment(Nbt), N0(Nbt))
  allocate(Endothelial(Nbt))
  allocate(Cells_density(Nbt))
  ALLOCATE(VitesseSeg(1:2,1:Nseg))

  CALL matrixinitVF4 ( N ) 
   call matrixinitVF4 ( A )
   call matrixinitVF4 ( Mat_E )

  do i=1,Nbt 
      Endothelial(i) = fsource( coordK(1,i), coordK(2,i), choixpb ) 
  end do

  call plot_vtk(Endothelial,'E_Pb20_0','Endothelial_cells')

  call assembleNutri(N,Endothelial)

  ! resolution du syst�me lin�aire  par la mathode du gradient conjugu�
  tol = 5.d-14
  N0(:)=1.D0
  Nutriment = bigradient(N, N%Bg,N0,tol)

   allocate(U(1:Nbt),U0(1:Nbt),dU(1:Nbt),U_p(1:Nbt))
   U0 = 0.D0 
   do i=1,Nbt
      U0(i) = cond_ini(CoordK(1,i),CoordK(2,i))
   end do 
   
   allocate(Tab_A(n_enty),Tab_U(n_enty,Nbt),Tab_U_p(n_enty,Nbt),Tab_entity(n_enty))
   Tab_A(1) = N;Tab_U(1,:) = Nutriment;Tab_U_p(1,:) = Nutriment;Tab_entity(1)='Nutrim'
   Tab_A(2) = Mat_E;Tab_U(2,:) = Endothelial;Tab_U_p(2,:) = Endothelial;Tab_entity(2)='Endoth'
   Tab_A(3) = A;Tab_U(3,:) = U0;Tab_U_p(3,:) = U0;Tab_entity(3)='u_norm'

   Cells_density = U0 + Endothelial

   call plot_vtk(U0,'TC_Pb20_0','Tumor_cells')
   call plot_vtk(Nutriment,'N_Pb20_0','Nutrients')
   call plot_vtk(Cells_density,'CD_Pb20_0','Cells_density')
   U = U0
   
   print*, 'matrix init ok' 

 !call MethNewton(A,U,U0,Nutriment,Endothelial,'u_norm',tol)
 !U_p = U
   call MethNewtonMul(Tab_A,Tab_U,Tab_U_p,Tab_entity,tol)
   Tab_U_p=Tab_U
 print*,'methnewton ok'
 i=1
 do while(i*delta<=Tf)
   !call actuNutri(N,U,Nutriment,Endothelial)
   !call MethNewton(A,U,U_p,Nutriment,Endothelial,'u_norm',tol)
   call MethNewtonMul(Tab_A,Tab_U,Tab_U_p,Tab_entity,tol)
   Cells_density = Tab_U(2,:) + Tab_U(3,:)
   Tab_U_p = Tab_U
   i=i+1
   if (i==500) then
     call plot_vtk(Tab_U(1,:),'N_Pb20_1','Nutrients')
     call plot_vtk(Tab_U(3,:),'TC_Pb20_1','Tumor_cells')
     call plot_vtk(Cells_density,'CD_Pb20_1','Cells_density')
     print*,i*delta
   else if(i==2000) then
      call plot_vtk(Tab_U(1,:),'N_Pb20_2','Nutrients')
      call plot_vtk(Tab_U(3,:),'TC_Pb20_2','Tumor_cells')
      call plot_vtk(Cells_density,'CD_Pb20_2','Cells_density')
      print*,i*delta
   else if(i==5000) then
      call plot_vtk(Tab_U(1,:),'N_Pb20_3','Nutrients')
      call plot_vtk(Tab_U(3,:),'TC_Pb20_3','Tumor_cells')
      call plot_vtk(Cells_density,'CD_Pb20_3','Cells_density')
      print*,i*delta
   else if(i==10000) then
      call plot_vtk(Tab_U(1,:),'N_Pb20_4','Nutrients')
      call plot_vtk(Tab_U(3,:),'TC_Pb20_4','Tumor_cells')
      call plot_vtk(Cells_density,'CD_Pb20_4','Cells_density')
      print*,i*delta
   else if(i==15000) then 
      call plot_vtk(Tab_U(1,:),'N_Pb20_5','Nutrients')
      call plot_vtk(Tab_U(3,:),'TC_Pb20_5','Tumor_cells')
      call plot_vtk(Cells_density,'CD_Pb20_5','Cells_density')
      print*,i*delta
   else if(i==20000) then
      call plot_vtk(Tab_U(1,:),'N_Pb20_6','Nutrients')
      call plot_vtk(Tab_U(3,:),'TC_Pb20_6','Tumor_cells')
      call plot_vtk(Cells_density,'CD_Pb20_6','Cells_density')
      print*,i*delta
   else if(i==25000) then
      call plot_vtk(Tab_U(1,:),'N_Pb20_7','Nutrients')
      call plot_vtk(Tab_U(3,:),'TC_Pb20_7','Tumor_cells')
      call plot_vtk(Cells_density,'CD_Pb20_7','Cells_density')
      print*,i*delta
   else if(i==30000) then
      call plot_vtk(Tab_U(1,:),'N_Pb20_8','Nutrients')
      call plot_vtk(Tab_U(3,:),'TC_Pb20_8','Tumor_cells')
      call plot_vtk(Cells_density,'CD_Pb20_8','Cells_density')
      print*,i*delta
   end if
   if (modulo(i,100)==1) then 
      print*,i*delta
   end if
 end do


100 FORMAT(10(E10.3,2x))
200 FORMAT(6(E14.6,2x))
  CLOSE (uprint)
  PRINT*,'fin du travail laplacien'
  print*,'CHOIX PB = ', choixpb


end program KellerSegel 