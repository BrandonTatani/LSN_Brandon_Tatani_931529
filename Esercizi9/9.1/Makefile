CC = g++
CFLAGS = -Wall -O3 --std=c++11

Genetic_Algorithm.exe : main.cpp Genetic_Algorithm.o random.o
	$(CC) random.o Genetic_Algorithm.o main.cpp -o Genetic_Algorithm.exe
Genetic_Algorithm.o : Genetic_Algorithm.cpp
	$(CC) -c Genetic_Algorithm.cpp -o Genetic_Algorithm.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o Monte_Carlo_ISING_1D.exe
output:
	rm -rf *.dat
	rm -rf *.config
config:
	rm -rf config.final
