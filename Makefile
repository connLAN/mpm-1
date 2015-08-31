CC = gcc
#MPICC = mpicc
MPICC = /opt/openmpi-1.8/bin/mpicc
INCLUDES = -I./include
S_OBJECTS = ./build/initmpm.o ./build/mpmfunctions.o ./build/commfunctions.o ./build/mpmio.o
P_OBJECTS = ./build/initmpm.o ./build/mpmfunctions.o ./build/commfunctions.o ./build/mpmio.o
#FLAGS = -Wall -g
FLAGS = -Wall

serial: $(S_OBJECTS)
	$(MPICC) $(INCLUDES) $(S_OBJECTS) serial_mpm.c -o ./build/serial_mpm

parallel: $(P_OBJECTS)
	$(MPICC) $(FLAGS) $(INCLUDES) $(P_OBJECTS) parallel_mpm.c -o ./build/parallel_mpm

./build/initmpm.o:
	$(MPICC) $(FLAGS) $(INCLUDES) -c initmpm.c -o ./build/initmpm.o

./build/mpmfunctions.o:
	$(MPICC) $(FLAGS) $(INCLUDES) -c mpmfunctions.c -o ./build/mpmfunctions.o

./build/commfunctions.o:
	$(MPICC) $(FLAGS) $(INCLUDES) -c commfunctions.c -o ./build/commfunctions.o

./build/mpmio.o:
	$(MPICC) $(FLAGS) $(INCLUDES) -c mpmio.c -o ./build/mpmio.o

clean:
	rm ./build/*

clean_debug:
	rm -rf ./build/parallel_mpm.dSYM
	rm ./build/parallel_mpm
	rm $(P_OBJECTS)

clean_all:
	find ./build/ -maxdepth 1 -name '*.vtk' -delete
	rm ./build/*
