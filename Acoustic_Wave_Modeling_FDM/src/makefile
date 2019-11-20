CC = gcc

LDFLAGS = -lm

OBJ1 = hete_aco_wave_2D.o seis_wave.o
OBJ2 = homo_aco_wave_1D.o seis_wave.o
OBJ3 = homo_aco_wave_2D.o seis_wave.o
OBJ4 = txt2sac_1D.o sacio.o

BIN1 = ../bin/hete_aco_wave_2D
BIN2 = ../bin/homo_aco_wave_1D
BIN3 = ../bin/homo_aco_wave_2D
BIN4 = ../bin/txt2sac_1D

all: $(BIN1) $(BIN2) $(BIN3) $(BIN4)
.PHONY: all

$(BIN1): $(OBJ1)
	$(CC) -o $(BIN1) $(OBJ1) $(LDFLAGS)
$(BIN2): $(OBJ2)
	$(CC) -o $(BIN2) $(OBJ2) $(LDFLAGS)
$(BIN3): $(OBJ3)
	$(CC) -o $(BIN3) $(OBJ3) $(LDFLAGS)
$(BIN4): $(OBJ4)
	$(CC) -o $(BIN4) $(OBJ4) $(LDFLAGS)

$(OBJ1): seis_wave.h
$(OBJ2): seis_wave.h
$(OBJ3): seis_wave.h
$(OBJ4): sacio.h

clean:
	rm *.o
