CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o
	$(CC) main.o -o main.exe
main.o : MolDyn_NVE.cpp
	$(CC) -c MolDyn_NVE.cpp -o main.o $(CFLAGS)
	
start:
	cp config.fcc config.0
	./main.exe 0
	
restart:
	cp config.final config.0
	cp old.final old.0
	./main.exe 1
	
clean:
	rm *.dat
