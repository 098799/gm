F90 = gfortran
OPTIONS = -fcheck=all -Wall -Wextra -Warray-temporaries -Wrealloc-lhs-all -pedantic -std=f2008

NAME = gm
OBJ  = precision.o commontypes.o inputread.o integrals.o scf.o f12.o main.o

%.o : %.f90
	$(F90) -c $(OPTIONS) $<

$(NAME) : $(OBJ)
	$(F90) $(OPTIONS) -o $@ $(OBJ)

precision.o   :
commontypes.o : precision.o
inputread.o   : precision.o commontypes.o
integrals.o   : precision.o commontypes.o
scf.o         : precision.o commontypes.o integrals.o
f12.o         : precision.o commontypes.o
main.o        : precision.o commontypes.o inputread.o

clean:
	rm *.o *.mod $(NAME)
