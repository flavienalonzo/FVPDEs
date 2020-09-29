!!!!! Informations sur le maillage

la lecture d'un maillage est réalisée dans le fichier readmesh.f90, le type des données récupérées est précisé dans parmmage.f90 :

On a des informations sur les sommets du maillage :
	- Nbs (integer) est le nombre de sommets dans le maillage. Pour chaque sommet, on a :
		- CoordS (real, dimension(1:2,Nbs)), les coordonnées dans le plan de chaque sommet, par convention CoordS(:,i) donne les coordonnées du i-ème sommet, la première ligne correspond aux abscisses et la seconde aux ordonnées.
		- ntyps (integer), le type du sommet vis-à-vis des conditions aux bords (0 : le sommet n'est pas sur le bord, 1 : appartient à une condition de type Dirichlet, 2 : appatient à une condition de type Neumann)

On a des informations sur les mailles (triangles) :
	- Nbt (integer) est le nombre de triangles formés dans le maillage, pour chaque triangle, on a :
		- NuSoK (integer, dimension (1:3,Nbt)) correspond à la liste des connectivités entre les triangles et les sommets du maillage, ainsi NuSoK(:,i) donne le numéro des 3 sommets qui définissent le i-ème triangle.
		- CoordK (real, dimension(1:2,Nbt)) donne les coordonnées des centres des triangles, ainsi CoordK(:,i) est le centre du i-ème triangle.
		- AireK (real, dimension(Nbt)) donne l'aire de chaque triangle, ainsi AireK(i) est l'aire du i-ème triangle.

On a des informations sur les segments :
	- Nseg (integer) est le nombre de segment présent dans le maillage. Pour chaque segment, on a :
		-NuSeg (integer, dimension(1:2,Nseg)), donne le numéro des sommets qui définissent un segment, ainsi NuSeg(:,i) contient les 2 sommets du i-ème segment, à priori on ne privilégie pas d'ordre pour le numéro des sommets.
		-NombVoisSeg (integer,dimension(Nseg)) donne le nombre de triangles que chaque sommet touche, ainsi NombVoisSeg(i) contient le nombre de triangles qui ont i comme sommet, 1 si le segment est sur le bord du maillage, 2 si le segment est à l'intérieur.
		-NumTVoisSeg (integer,dimension (1:2,Nseg)) donne les numéros des triangles que chaque sommet touche, ainsi NumTVoisSeg(:,i) contient les numéros des 2 triangles qui admettent i comme segment. Par convention s'il n'y a qu'un seul triangle son numéro est dans NumTVoisSeg(1,i) et on complète NumTVoisSeg(2,i) par un entier quelconque négatif.
		-NtypSeg (integer,dimension(Nseg)) donne le type de chaque segment vis-à-vis des conditions aux bords, ainsi NtypSeg(i) contient 0 si le i-ème segment n'est pas sur le bord, 1 si le i-ème segment appartient à une condition de type Dirichlet, 2 si le i-ème segment appartient à une condition de type Neumann.
		-TauKL (real, dimension(Nseg)) donne la valeur du coefficient de transmissibilité de chaque segment vis-à-vis des triangles qui touchent chaque segment. Ainsi TauKL(i) donne le coefficient de transmissibilité du triangle NumTVoisSeg(1,i) au triangle NumTVoisSeg(2,i) (TauKL(i) est aussi celui du triangle NumTVoisSeg(2,i) au triangle NumTVoisSeg(1,i)). Si le i-ème segment n'appartient qu'à un seul triangle on a une valeur ...
		-dKL (real, dimension(Nseg)) donne la distance entre les centres des triangles qui ont un segment commun. Ainsi dKL(i) est la distance entre les centres des triangles NumTVoisSeg(1,i) et NumTVoisSeg(2,i). Si un segment n'appartient qu'à un seul triangle, on a une valeur ...

En complément de ces valeurs, on peut définir :
	-Nsint : le nombre de sommets qui sont à l'intérieur du maillage, 
	-Nbord : le nombre de sommets sur le bord du maillage
		
!!! Conditions aux bords :

Les conditions sont données dans le fichier readmesh.f90, on peut changer à la main le type de condition aux bords en mettant directement 'Dirichlet' ou 'Neumann'. C'est information correspond bien à l'entier 1 ou 2 parce que ces valeurs ont été assignées dans longr.f90

!!! Stockage des matrices sparses :

Un type est défini pour stocker les informations d'une matrice sparse : type(MatCreux)
Si on note A une matrice sparse on pourra alors avoir :
	-A%IndPL (integer,dimension(Nbt+1)) 	permet de retrouver l'entrée i
	-A%Indc (integer,dimension(Ncoef)) 	permet de retrouver l'entrée j
	-A%TMat (real, dimension(Ncoef))	permet de retrouver le coefficient Aij
Pour toute ligne k de la matrice, [A%IndPL(k),A%IndPL(k+1)-1] correspond aux indices des tableaux A%Indc et A%TMat relatifs à la ligne k de la matrice, ainsi A%Indc([A%IndPL(k),A%IndPL(k+1)-1]) et A%TMat([A%IndPL(k),A%IndPL(k+1)-1]) sont respectivement les indices des colonnes des coefficients sur la k-ème ligne et la valeur des coefficients A(kj).
	-A%Diag (real, dimension(Nbt)) donne les coefficients de la diagonale de A, ainsi A%Diag(k) donne A(kk)
	-A%PosiDiag (integer,dimension(Nbt)) donne les indices dans %TMat des coefficients diagonaux de A, ainsi A%PosiDiag(k) donne l'indice du coefficient A(kk) dans A%TMat.
	-A%F (real,dimension(Nbt)) qui est un attribut propre à A donc par exemple on peut y associé les forces volumiques.
	-A%Bg (real,dimension(Nbt)) est aussi un atribut propre à A donc par exemble on peut y associé les termes de bords.

Le programme matrixinitVF4.f90 permet de créer une matrice sparse associée à un maillage en initialisant tous les coefficients à 0. La subroutine ajout.f90 permet ultérieurement d'ajouter à A un coefficient à une ligne i et colonne j.

Les programmes assemblevf4.f90, assembleVitesse.f90 et assembletheta.f90 permettent respectivement d'associées la discrétisation à A du terme de diffusion, du terme de convection et du terme de réaction.
	
