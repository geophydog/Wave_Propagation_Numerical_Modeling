CC = gcc

LDFLAGS = -lm

OBJ = hete_aco_wave_2D.o seis_wave.o

BIN = ../bin/hete_aco_wave_2D

all: $(BIN)
.PHONY: all

$(BIN): $(OBJ)
	$(CC) -o $(BIN) $(OBJ) $(LDFLAGS)

$(OBJ): seis_wave.h

clean:
	rm *.o
