CC=mpicc
RUN=mpirun
SOURCE=$(wildcard *.c) 
OBJ=$(basename $(SOURCE))
NP=4
MACHINEFILE=machinefile
PARAMS=-np $(NP)
MFILE = -machinefile $(MACHINEFILE)
ARGS = 5


$(OBJ): $(SOURCE)
	$(CC) -o $(OBJ) $(SOURCE)

run: 
	$(RUN) $(PARAMS) $(OBJ) $(ARGS)

run_remote:
	$(RUN) $(MFILE) $(PARAMS) $(OBJ) $(ARGS)

clean: 
	rm $(OBJ) 