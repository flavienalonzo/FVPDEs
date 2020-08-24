module longr
  implicit none
  integer, parameter     :: long = 8
  character (len=6)      :: prefix
  integer                :: uread, uprint, iprint,umesh, uplotvtk , n_enty, index_endo, index_nut, index_norm, index_vegf, WhichPb
!  integer                :: uma,uele, uneigh, unode, Maillage_Format
  real (kind = long)     :: theta,delta, Coef_diffusion, chi , Coef_prod, rate, Diff_u, Tf, Coef_cons, seuil_hypo, seuil_necro
  real (kind = long)     :: apop, rat_pop, VEGF_prod,VEGF_dif, VEGF_cons, chemo_endo, VEGF_degr, satur_endo, Diff_endo, satur_norm,&
                            & rate_endo, satur_nutri
  CHARACTER (len = 50)   :: nom_mesh
  !nom_maillage,nom_ele, nom_neigh, nom_node, nom_v_node, &
  !      & nom_edge, FPLOTVTK, dossiermaillage
  CHARACTER (len = 50)   :: FPLOTVTK
  integer :: ChoixPb         ! probleme � resoudre 

  integer, parameter     :: Neumann=2, Dirichlet=1

 ! integer, parameter     :: VF4=1, P1Sommets=2, P1milieux=3, P1milieuxmonotone=4

  INTEGER,PARAMETER      :: len_buffer=80

end module longr
