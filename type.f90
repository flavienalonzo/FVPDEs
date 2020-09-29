module type

implicit none 
type matsparse
integer, dimension(:), pointer :: intI, intJ, PosDiag
double precision, dimension(:), pointer :: ValMat, Diag
end type matsparse


end module type