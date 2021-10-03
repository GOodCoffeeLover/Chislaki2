all: P1D

P1D : main.o P1D.o 
	g++ P1D.o main.o  -o P1D

main.o: main.cpp Progonka.hpp P1D.h
	g++ -c main.cpp 

P1D.o: P1D.cpp P1D.h Progonka.hpp
	g++ -c P1D.cpp 	

clean:
	rm P1D.o main.o

