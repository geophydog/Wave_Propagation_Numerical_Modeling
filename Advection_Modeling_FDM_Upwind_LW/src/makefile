CC = g++

LDMFLAG = -lm

OBJ = advection_forward.o
OBJ1 = ad1D_run.o $(OBJ)
OBJ2 = ad2D_run.o $(OBJ)

BIN1 = ../bin/ad1D_run
BIN2 = ../bin/ad2D_run

all: $(BIN1) $(BIN2)
.PHONY: all

$(BIN1): $(OBJ1)
	$(CC) -o $(BIN1) $(OBJ1) $(LDMFLAG)
$(BIN2): $(OBJ2)
	$(CC) -o $(BIN2) $(OBJ2) $(LDMFLAG) $(LDFFTFLAG)

$(OBJ1) $(BIN2): advection_forward.h

clean:
	rm *.o
