CC=mpicc
RUN=mpirun
FLAGS=-lm
SOURCE=$(wildcard *.c) 
OBJ=$(basename $(SOURCE))
NP=2
MACHINEFILE=machinefile
PARAMS=-np $(NP)
MFILE = -machinefile $(MACHINEFILE)


$(OBJ): $(SOURCE)
	$(CC) -o $(OBJ) $(SOURCE) $(FLAGS)

run: $(OBJ)
	$(RUN) $(PARAMS) $(OBJ)

run_remote: $(OBJ)
	$(RUN) $(MFILE) $(PARAMS) $(OBJ)

clean: 
	rm $(OBJ) 