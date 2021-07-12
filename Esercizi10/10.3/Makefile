CC = g++
CFLAGS = -Wall -O3 --std=c++11

Genetic_Algorithm.exe : main.cpp Genetic_Algorithm.o random.o
	$(CC) random.o Genetic_Algorithm.o main.cpp -o Genetic_Algorithm.exe
Genetic_Algorithm.o : Genetic_Algorithm.cpp
	$(CC) -c Genetic_Algorithm.cpp -o Genetic_Algorithm.o $(CFLAGS)

exec:
	rm *.dat
	rm percorso.config
	./Genetic_Algorithm.exe	

	
	
MPICXX:	main2.cpp Genetic_Algorithm.o random.o
	mpicxx Genetic_Algorithm.o random.o main2.cpp

MPIEXEC:
	rm *.dat
	rm percorso.config
	mpiexec -np 2 a.out	

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

output:
	rm -rf *.dat
	rm -rf *.config

