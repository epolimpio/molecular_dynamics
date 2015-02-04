FC      = gfortran
FFLAGS  = -Wall -Wextra -O3 -fimplicit-none -march=native 
#FFLAGS += -pedantic -fbounds-check -fmax-errors=1 -g
#FFLAGS += $(shell pkg-config --cflags plplotd-f95)
#LDFLAGS = $(shell pkg-config --libs plplotd-f95)
LDFLAGS = 
LIBS    =

COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) $(LDFLAGS)

TARGET = main.exe       # Name of final executable to produce
#OBJS = md_plot.o 
OBJS =
OBJS += initial_conditions.o
OBJS += md.o # List of object dependencies

$(TARGET): $(OBJS)
	$(LINK) -o $@ $^ $(LIBS)

%.o:%.f95
	$(COMPILE) -c $<
