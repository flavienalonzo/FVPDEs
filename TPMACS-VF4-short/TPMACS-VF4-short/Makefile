FORTRAN=gfortran
#COPTS= -c  -O -g -pedantic -fbounds-check 
COPTS= -c  -O3 -g 
#-c -checkbound -O -r8 -v
PROG=mainvf4
SUFFIXES=.f90.o
.SUFFIXES : .f90 .o
$(SUFFIXES):
	$(FORTRAN) $(COPTS) $*.f90

OBJS	      = longr.o\
		parmmage.o \
		imprime.o \
		intmatvec.o \
		algebrelineaire.o\
		intbigradc.o\
		plotvtkmod.o\
		fsourcemod.o\
		init.o \
		matrixinitVF4.o \
		assembletheta.o \
		assemblevf4.o \
		assembleVitesse.o \
		readmesh.o \
        ajout.o \
		parameters.o \
		vide.o \
		assemblediffusionKS.o \
		assemblereactionKS.o \
		assemblechemoKS.o \
		MethNewton.o \
		plotvtktime.o \
		actuNutri.o \
		assembleNutri.o \
		MethNewtonMul.o \
		KellerSegel.o 

$(PROG): $(OBJS)
	$(FORTRAN) $(OBJS) -o $(PROG)

clean:;	rm -f *.o *.mod objets/*



