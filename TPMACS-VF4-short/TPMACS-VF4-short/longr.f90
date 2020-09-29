module longr
  implicit none
  integer, parameter     :: long = 8
  character (len=6)      :: prefix
  integer                :: uread, uprint, iprint,umesh, uplotvtk , n_enty, index_endo, index_nut, index_norm, index_vegf, WhichPb
  integer                :: uma,uele, uneigh, unode, Maillage_Format, Nchemo, Nradio
  real (kind = long)     :: theta,dt, Coef_diffusion, chi_u , Coef_prod, rate, Diff_u, Tf, Coef_cons, seuil_hypo, seuil_necro
  real (kind = long)     :: apop, rat_pop, VEGF_prod,VEGF_dif, VEGF_cons, chemo_endo, VEGF_degr, satur_endo, Diff_endo, satur_norm,&
                            & rate_endo, satur_nutri, nut_degra, degr_endo, epsilon
  real (kind = long)     :: TolerenceGradient, TolerenceNewton
  real (kind = long)     :: delta, deltau, deltax, deltay, deltaxy, deltaxu, deltayu, deltaxyu, a33, b33, c33, d33
  real (kind = long)     :: mu, uc, ubar, gamma, us, vs, epsu, epsv, NBnu, NBmu, NBnv, NBmv
  real (kind = long)     :: alpha, beta, CoefDiffV, CoefDiffuAdeg, CoefTranspChi, dose_chemo, seuil_surg, time_surg_d, time_surg_f,&
                            & radio_alpha, radio_beta, dose_radio
  real (kind = long), dimension(2,3) :: chemo_time
  real (kind = long), dimension(2,6) :: radio_time
  CHARACTER (len = 50)   :: nom_mesh,nom_maillage,nom_ele, nom_neigh, nom_node, nom_v_node, &
        & nom_edge, FPLOTVTK, dossiermaillage
  integer :: ChoixPb         ! probleme � resoudre 

  integer, parameter     :: transparent=3, Neumann=2, Dirichlet=1
  integer, parameter     :: VF4=1, P1Sommets=2, P1milieux=3, P1milieuxmonotone=4
  integer, parameter     :: len_buffer=80, maillageregulier=0, typegeomaille=3
  integer                :: choixkscalaire, choixkscalaireu,ChoixPlot, ChoixDegenere, ChoixUpwind, &
  & choixanisu, choixanis, NbInc, NbTotal, ChoixChi, ChoixAdeg, ChoixSchema 



end module longr
