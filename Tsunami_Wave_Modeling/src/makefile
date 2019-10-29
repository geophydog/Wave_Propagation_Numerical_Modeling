CC = gcc

LDMFLAG = -lm
LDFFTFLAG = -lfftw3

OBJ = wave_forward.o
OBJ1 = tsunami_wave_FD22_2D_hete.o $(OBJ)
OBJ2 = tsunami_wave_FD22_2D_homo.o $(OBJ)
OBJ3 = tsunami_wave_FD24_2D_hete.o $(OBJ)
OBJ4 = tsunami_wave_FD24_2D_homo.o $(OBJ)
OBJ5 = tsunami_wave_PSM_2D_homo.o $(OBJ)

BIN1 = tsunami_wave_FD22_2D_hete
BIN2 = tsunami_wave_FD22_2D_homo
BIN3 = tsunami_wave_FD24_2D_hete
BIN4 = tsunami_wave_FD24_2D_homo
BIN5 = tsunami_wave_PSM_2D_homo

all: $(BIN1) $(BIN2) $(BIN3) $(BIN4) $(BIN5)
.PHONY: all

$(BIN1): $(OBJ1)
	$(CC) -o $(BIN1) $(OBJ1) $(LDMFLAG)
$(BIN2): $(OBJ2)
	$(CC) -o $(BIN2) $(OBJ2) $(LDMFLAG)
$(BIN3): $(OBJ3)
	$(CC) -o $(BIN3) $(OBJ3) $(LDMFLAG)
$(BIN4): $(OBJ4)
	$(CC) -o $(BIN4) $(OBJ4) $(LDMFLAG)
$(BIN5): $(OBJ5)
	$(CC) -o $(BIN5) $(OBJ5) $(LDMFLAG) $(LDFFTFLAG)

$(OBJ1) $(BIN2) $(BIN3) $(BIN4) $(BIN5): wave_forward.h

clean:
	rm *.o
