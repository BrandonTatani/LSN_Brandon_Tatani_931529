CC = g++
CFLAGS = -Wall -O3 --std=c++11

Simulated_Annealing.exe : main.cpp Simulated_Annealing.o random.o
	$(CC) random.o Simulated_Annealing.o main.cpp -o Simulated_Annealing.exe
Simulated_Annealing.o : Simulated_Annealing.cpp
	$(CC) -c Simulated_Annealing.cpp -o Simulated_Annealing.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

output:
	rm -rf *.dat
	rm -rf *.config
