F90   = pgf90
F90OPTFLAGS= -O3 -mp
F90FLAGS =$(F90OPTFLAGS)

.SUFFIXES:
.SUFFIXES: .o .f .f90 .c 
#
.f90.o:; $(F90) -c $(F90FLAGS)  $(F90OPTFLAG) $<
.f.o:; $(F90) -c -loglist $(F90FLAGS)  $(F90OPTFLAG) $<
#
OBJS = \
Multigrid_1D_Neumann_BC.o Multigrid_1D_Neumann_BC_OMP.o Main.o \

TARGET = Run

all: $(TARGET)

$(TARGET): $(OBJS) 
	$(F90) $(F90FLAGS) -o $(TARGET) \
	$(OBJS) \
	$(F90FLAGS)

clean:
	rm -f *.o $(TARGET) *.mod *~ PI* *.log *.lst *.txt *.plt
