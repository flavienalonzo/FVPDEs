program mainKellerSegel

use longr
use parmmage

use plotvtkmod
use load_module
use parameters

type(MatCreux)       :: A

INTEGER                             :: jt, i, j, is, jv,iseg,  kiter
REAL(kind=long), DIMENSION(:), ALLOCATABLE  :: U,U0, Uexacte

REAL(kind = long)               :: tol, seuil

call readmatlab('Tri_2D_0.05_5_1_2_@(x)[-x(2),x(1)].txt')

!call matrixinitVF4(A)



end program mainKellerSegel