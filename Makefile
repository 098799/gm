F90     = gfortran
#OPTIONS = -fcheck=all -Wall -Wextra -Warray-temporaries -Wrealloc-lhs-all -pedantic -std=f2008
OPTIONS = -O3
BLAS    = -llapack -lblas

NAME = gm
OBJ  = precision.o commontypes.o inputread.o diis.o integrals.o scf_driver.o f12.o main.o

%.o : %.f90
	$(F90) -c $(OPTIONS) $<

$(NAME) : $(OBJ)
	$(F90) $(OPTIONS) -o $@ $(OBJ) $(BLAS)

precision.o   :
commontypes.o : precision.o
inputread.o   : precision.o commontypes.o
diis.o        : precision.o commontypes.o
integrals.o   : precision.o commontypes.o
scf_driver.o  : precision.o commontypes.o diis.o integrals.o
f12.o         : precision.o commontypes.o
main.o        : precision.o commontypes.o inputread.o diis.o scf_driver.o f12.o

clean:
	rm *.o *.mod $(NAME)
