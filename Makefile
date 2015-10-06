CC=g++
CCOPT=-O3 -Wall

all : Delta
	echo `date` @ `uname -n` > VERSION

Delta : Delta.o Moose.o Organism.o Gamete.o
	$(CC) $(CCOPT) -o delta $^ -lgsl -lgslcblas -lm

Delta.o : Delta.cpp Moose.h
	$(CC) $(CCOPT) -c Delta.cpp

Moose.o : Moose.cpp Moose.h Organism.h
	$(CC) $(CCOPT) -c Moose.cpp

Organism.o : Organism.cpp Organism.h Gamete.h
	$(CC) $(CCOPT) -c Organism.cpp

Gamete.o : Gamete.cpp Gamete.h
	$(CC) $(CCOPT) -c Gamete.cpp

clean :
	rm delta *.o
