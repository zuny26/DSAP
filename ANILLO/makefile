CC=mpicc
RUN=mpirun
SOURCE=$(wildcard *.c) 
OBJ=$(basename $(SOURCE))
NP=4
MACHINEFILE=machinefile
PARAMS=-np $(NP)
MFILE = -machinefile $(MACHINEFILE)


$(OBJ): $(SOURCE)
	$(CC) -o $(OBJ) $(SOURCE)

run: 
	$(RUN) $(PARAMS) $(OBJ)

run_remote:
	$(RUN) $(MFILE) $(PARAMS) $(OBJ)

clean: 
	rm $(OBJ) 