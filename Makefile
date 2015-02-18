FC      = gfortran
FFLAGS  = -Wall -Wextra -O3 -fimplicit-none -march=native -ffast-math
#FFLAGS  = -Wall -Wextra -O3 -fimplicit-none
#FFLAGS += -pedantic -fbounds-check -fmax-errors=1 -g
#FFLAGS += $(shell pkg-config --cflags plplotd-f95)
#LDFLAGS = $(shell pkg-config --libs plplotd-f95)
LDFLAGS = $(HOME)/iccp/XPS/lib/libxps.a -lX11 
LIBS    =

COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) 

TARGET = main.exe       # Name of final executable to produce
#OBJS = md_plot.o 
OBJS =
OBJS += simParam.o
OBJS += initial_conditions.o
OBJS += plot_part.o
OBJS += physical_quantities.o
OBJS += output_file.o
OBJS += simulation.o
OBJS += main.o # List of object dependencies

$(TARGET): $(OBJS)
	$(LINK) -o $@ $^ $(LIBS) $(LDFLAGS)

%.o:%.f95
	$(COMPILE) -c $<

.PHONY:clean
clean:
	rm -r *.mod *.dat *.exe *.o

