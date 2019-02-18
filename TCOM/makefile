CC=mpicc
RUN=mpirun
FLAGS=-lm
SOURCE=$(wildcard *.c) 
OBJ=$(basename $(SOURCE))
NP=4
MACHINEFILE=machinefile
PARAMS=-np $(NP)
MFILE = -machinefile $(MACHINEFILE)


$(OBJ): $(SOURCE)
	$(CC) -o $(OBJ) $(SOURCE) $(FLAGS)

run: 
	$(RUN) $(PARAMS) $(OBJ)

run_remote:
	$(RUN) $(MFILE) $(PARAMS) $(OBJ)

clean: 
	rm $(OBJ) 