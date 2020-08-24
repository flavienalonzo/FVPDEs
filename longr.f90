module longr
  implicit none
  integer, parameter     :: long = 8
  character (len=6)      :: prefix
  integer                :: uread, uprint, iprint,umesh, uplotvtk
!  integer                :: uma,uele, uneigh, unode, Maillage_Format
  real (kind = long)     :: theta,delta, Coef_diffusion
  CHARACTER (len = 50)   :: nom_mesh
  !nom_maillage,nom_ele, nom_neigh, nom_node, nom_v_node, &
  !      & nom_edge, FPLOTVTK, dossiermaillage
  CHARACTER (len = 50)   :: FPLOTVTK
  integer :: ChoixPb         ! probleme à resoudre 

  integer, parameter     :: Neumann=2, Dirichlet=1

 ! integer, parameter     :: VF4=1, P1Sommets=2, P1milieux=3, P1milieuxmonotone=4

  INTEGER,PARAMETER      :: len_buffer=80

end module longr
