CC=mpicc
RUN=mpirun
FLAGS= -lm
SOURCE=$(wildcard *.c) 
OBJ=$(basename $(SOURCE))
NP=4
MACHINEFILE=machinefile
PARAMS=-np $(NP)
MFILE = -machinefile $(MACHINEFILE)
ARGS = 


$(OBJ): $(SOURCE)
	$(CC) -o $(OBJ) $(SOURCE) $(FLAGS)

run: 
	$(RUN) $(PARAMS) $(OBJ) $(ARGS)

run_remote:
	$(RUN) $(MFILE) $(PARAMS) $(OBJ) $(ARGS)

clean: 
	rm $(OBJ) 