CC=mpicc
RUN=mpirun
FLAGS= -lm
SOURCE=$(wildcard *.c) 
OBJ=$(basename $(SOURCE))
NP=8
MACHINEFILE=machinefile
PARAMS=-np $(NP)
MFILE = -machinefile $(MACHINEFILE)
ARGS = 


$(OBJ): $(SOURCE)
	$(CC) -o $(OBJ) $(SOURCE) $(FLAGS)

run: $(OBJ)
	$(RUN) $(PARAMS) $(OBJ) $(ARGS)

run_remote: $(OBJ)
	$(RUN) $(MFILE) $(PARAMS) $(OBJ) $(ARGS)

clean: 
	rm $(OBJ) *.pgm